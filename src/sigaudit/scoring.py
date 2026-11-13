"""Turning a signature into a per-patient risk score, and the baseline scores
every signature is compared against."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .signature import Signature

# PAM50 proliferation index genes (Nielsen et al., Clin Cancer Res 2010), current symbols
PROLIFERATION_11 = ["BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80",
                    "MKI67", "PTTG1", "RRM2", "TYMS", "UBE2C"]
# anchor for the meta-PCNA-style score: PCNA plus well-established cell-cycle genes
PROLIFERATION_CORE = PROLIFERATION_11 + ["PCNA", "TOP2A", "MCM2", "MCM6", "CENPF", "BUB1"]


@dataclass
class Scored:
    score: pd.Series
    found: list[str]
    missing: list[str]

    @property
    def coverage(self) -> float:
        n = len(self.found) + len(self.missing)
        return len(self.found) / n if n else 0.0


def score_signature(sig: Signature, z: pd.DataFrame) -> Scored:
    """Weighted mean of per-gene z-scores: sum(w_g * z_g) / sum(|w_g|) over measured genes."""
    w = pd.Series(sig.weights, index=sig.genes, dtype=float)
    present = [g for g in sig.genes if g in z.index and z.loc[g].notna().all()]
    missing = [g for g in sig.genes if g not in present]
    if not present:
        return Scored(pd.Series(np.nan, index=z.columns), [], missing)
    wp = w[present]
    s = (z.loc[present].T @ wp) / wp.abs().sum()
    return Scored(s, present, missing)


def random_scores(sig: Signature, z: pd.DataFrame, k: int = 1000, seed: int = 0,
                  exclude: set[str] | None = None, chunk: int = 50) -> np.ndarray:
    """Scores of `k` random signatures the same size as `sig` (after gene matching).

    Each random signature reuses the real signature's weights, shuffled onto
    randomly drawn genes, so the null keeps the same size, sign mix and weight
    spread. Returns an (n_samples, k) array.
    """
    rng = np.random.default_rng(seed)
    exclude = exclude or set()
    zz = z.dropna(how="any")
    pool = np.array([g for g in zz.index if g not in exclude])
    Z = zz.loc[pool].to_numpy(np.float32)
    present = [g for g in sig.genes if g in zz.index]
    w = np.array([sig.weights[sig.genes.index(g)] for g in present] or sig.weights, np.float32)
    m = len(w)
    out = np.empty((Z.shape[1], k), np.float32)
    for start in range(0, k, chunk):
        kk = min(chunk, k - start)
        idx = np.stack([rng.choice(len(pool), m, replace=False) for _ in range(kk)])
        ww = np.stack([rng.permutation(w) for _ in range(kk)])
        out[:, start:start + kk] = np.einsum("kmn,km->nk", Z[idx], ww) / np.abs(w).sum()
    return out


def proliferation_score(z: pd.DataFrame, top_frac: float = 0.01) -> tuple[pd.Series, list[str]]:
    """Meta-PCNA-style proliferation score: mean z of the ~1% of genes most
    correlated with a core proliferation anchor.

    Venet et al. (2011) took the top 1% of genes correlated with PCNA in a
    normal-tissue compendium. Here the correlation is computed within the
    cohort itself (no outcome data are used), and the anchor is the mean of a
    17-gene proliferation core instead of the PCNA probe alone, because a
    single noisy probe (as on some arrays) gives a weak, misleading baseline.
    """
    zz = z.dropna(how="any")
    core = [g for g in PROLIFERATION_CORE if g in zz.index]
    if len(core) < 5:
        genes = [g for g in PROLIFERATION_11 if g in zz.index]
        return zz.loc[genes].mean(), genes
    x = zz.to_numpy(float)
    p = zz.loc[core].mean().to_numpy(float)
    xc = x - x.mean(1, keepdims=True)
    pc = p - p.mean()
    r = (xc @ pc) / (np.sqrt((xc ** 2).sum(1)) * np.sqrt((pc ** 2).sum()) + 1e-12)
    n_top = max(10, int(round(top_frac * len(r))))
    genes = list(zz.index[np.argsort(-r)[:n_top]])
    return zz.loc[genes].mean(), genes
