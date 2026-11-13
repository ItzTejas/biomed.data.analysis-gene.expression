"""Survival statistics used by the audit, implemented with numpy only.

Everything here is cross-checked against lifelines in tests/test_survival.py.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


# ---------------------------------------------------------------- concordance

def comparable_pairs(time: np.ndarray, event: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Index pairs (i, j) where patient i is known to fail before patient j.

    Harrell's definition: i had an event and t_i < t_j, or t_i == t_j with i an
    event and j censored. Higher risk score for i means the pair is concordant.
    """
    time = np.asarray(time, float)
    event = np.asarray(event, bool)
    ti, tj = time[:, None], time[None, :]
    ei, ej = event[:, None], event[None, :]
    mask = ei & ((ti < tj) | ((ti == tj) & ~ej))
    i, j = np.nonzero(mask)
    return i.astype(np.int32), j.astype(np.int32)


def concordance_from_pairs(score: np.ndarray, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """C-index for one score vector (n,) or many (n, k). Ties count 0.5."""
    s = np.asarray(score, float)
    if len(i) == 0:
        return np.full(s.shape[1:] if s.ndim > 1 else (), np.nan)
    si, sj = s[i], s[j]
    return ((si > sj).sum(0) + 0.5 * (si == sj).sum(0)) / len(i)


def concordance(time, event, score) -> float:
    i, j = comparable_pairs(time, event)
    return float(concordance_from_pairs(score, i, j))


def bootstrap_concordance(time, event, score, n_boot: int = 200, seed: int = 0,
                          groups: np.ndarray | None = None) -> np.ndarray:
    """Bootstrap distribution of the C-index (patients resampled with replacement)."""
    rng = np.random.default_rng(seed)
    time, event, score = (np.asarray(a) for a in (time, event, score))
    n = len(time)
    out = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        out[b] = concordance(time[idx], event[idx], score[idx])
    return out


# ------------------------------------------------------------------ Cox model

@dataclass
class CoxResult:
    names: list[str]
    beta: np.ndarray
    se: np.ndarray
    loglik: float
    n: int
    events: int
    converged: bool

    @property
    def hr(self):
        return np.exp(self.beta)

    @property
    def z(self):
        return self.beta / self.se

    @property
    def p(self):
        return 2 * stats.norm.sf(np.abs(self.z))

    def ci(self, level: float = 0.95):
        q = stats.norm.ppf(0.5 + level / 2)
        return np.exp(self.beta - q * self.se), np.exp(self.beta + q * self.se)

    def row(self, name: str) -> dict:
        k = self.names.index(name)
        lo, hi = self.ci()
        return {"hr": float(self.hr[k]), "lo": float(lo[k]), "hi": float(hi[k]),
                "p": float(self.p[k]), "beta": float(self.beta[k]), "se": float(self.se[k])}


def cox_fit(time, event, X, names: list[str] | None = None, max_iter: int = 50,
            tol: float = 1e-9, ridge: float = 0.0) -> CoxResult:
    """Cox proportional hazards by Newton-Raphson, Breslow handling of ties."""
    time = np.asarray(time, float)
    event = np.asarray(event, bool)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    n, p = X.shape
    names = names or [f"x{k}" for k in range(p)]

    order = np.argsort(-time, kind="stable")          # descending time
    t, e, Xs = time[order], event[order], X[order]
    # for each position, the last index sharing its time (risk set = everyone up to it)
    last = np.searchsorted(-t, -t, side="right") - 1

    beta = np.zeros(p)
    loglik_old = -np.inf
    converged = False
    for _ in range(max_iter):
        eta = Xs @ beta
        eta -= eta.max()
        w = np.exp(eta)
        S0 = np.cumsum(w)[last]
        S1 = np.cumsum(w[:, None] * Xs, axis=0)[last]
        S2 = np.cumsum(w[:, None, None] * Xs[:, :, None] * Xs[:, None, :], axis=0)[last]
        ev = e
        loglik = float(np.sum(eta[ev] - np.log(S0[ev]))) - 0.5 * ridge * beta @ beta
        mean = S1[ev] / S0[ev, None]
        grad = (Xs[ev] - mean).sum(0) - ridge * beta
        hess = -(S2[ev] / S0[ev, None, None] - mean[:, :, None] * mean[:, None, :]).sum(0)
        hess -= ridge * np.eye(p)
        try:
            step = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            break
        # the Cox partial likelihood is concave, so plain Newton steps converge;
        # cap the step to keep near-separated data from blowing up
        step = np.clip(step, -5, 5)
        beta = beta - step
        if abs(loglik - loglik_old) < tol:
            converged = True
            break
        loglik_old = loglik
    try:
        cov = np.linalg.inv(-hess)
        se = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        se = np.full(p, np.nan)
    return CoxResult(names, beta, se, loglik, n, int(event.sum()), converged)


# -------------------------------------------------------------- meta-analysis

def random_effects_meta(log_hr: np.ndarray, se: np.ndarray) -> dict:
    """DerSimonian-Laird random-effects pooling of log hazard ratios."""
    y = np.asarray(log_hr, float)
    v = np.asarray(se, float) ** 2
    ok = np.isfinite(y) & np.isfinite(v) & (v > 0)
    y, v = y[ok], v[ok]
    k = len(y)
    if k == 0:
        return {"k": 0}
    w = 1 / v
    fixed = np.sum(w * y) / np.sum(w)
    q = float(np.sum(w * (y - fixed) ** 2))
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 else 0.0
    wr = 1 / (v + tau2)
    est = np.sum(wr * y) / np.sum(wr)
    se_est = np.sqrt(1 / np.sum(wr))
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 and k > 1 else 0.0
    return {"k": k, "hr": float(np.exp(est)), "lo": float(np.exp(est - 1.96 * se_est)),
            "hi": float(np.exp(est + 1.96 * se_est)),
            "p": float(2 * stats.norm.sf(abs(est / se_est))), "tau2": float(tau2), "i2": float(i2)}
