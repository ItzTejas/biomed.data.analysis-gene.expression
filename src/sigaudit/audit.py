"""The audit: the checks applied to one signature across one or more cohorts."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import survival as sv
from .cohort import Cohort
from .scoring import random_scores, proliferation_score, score_signature
from .signature import Signature

MIN_COVERAGE = 0.5        # below this fraction of genes measured, a cohort is skipped
MIN_SUBGROUP_EVENTS = 10  # subgroups with fewer events are reported but not judged
NULL_PAIRS = 300_000      # comparable pairs sampled when scoring the random null


@dataclass
class AuditSettings:
    n_random: int = 1000
    n_boot: int = 200
    horizon_years: float | None = 10.0
    seed: int = 0
    subgroups: tuple[str, ...] = ("ancestry", "race", "sex", "age_band")


@dataclass
class CohortAudit:
    cohort: str
    endpoint: str
    platform: str
    n: int
    events: int
    training: bool
    genes_found: list[str]
    genes_missing: list[str]
    c_index: float = np.nan
    c_lo: float = np.nan
    c_hi: float = np.nan
    hr: dict = field(default_factory=dict)            # univariate, per SD of score
    null_c: list[float] = field(default_factory=list)  # random-signature C-indices
    random_percentile: float = np.nan                  # of |C-0.5| among random |C-0.5|
    random_p: float = np.nan
    prolif_c: float = np.nan
    prolif_genes: list[str] = field(default_factory=list)
    hr_adj_prolif: dict = field(default_factory=dict)
    hr_adj_clinical: dict = field(default_factory=dict)
    adj_covariates: list[str] = field(default_factory=list)
    subgroups: list[dict] = field(default_factory=list)
    skipped: str = ""

    @property
    def coverage(self) -> float:
        n = len(self.genes_found) + len(self.genes_missing)
        return len(self.genes_found) / n if n else 0.0

    @property
    def beats_random(self) -> bool:
        return bool(self.random_p < 0.05) and self.c_index > 0.5

    @property
    def beats_proliferation(self) -> bool:
        return bool(self.hr_adj_prolif) and self.hr_adj_prolif["p"] < 0.05 \
            and self.hr_adj_prolif["hr"] > 1


@dataclass
class AuditResult:
    signature: Signature
    unmapped: list[str]
    cohorts: list[CohortAudit]
    meta: dict
    settings: AuditSettings

    @property
    def independent(self) -> list[CohortAudit]:
        return [c for c in self.cohorts if not c.training and not c.skipped]

    def verdict(self) -> dict:
        """Plain-language summary over the cohorts the signature was NOT built on."""
        ind = self.independent
        k = len(ind)
        if k == 0:
            return {"grade": "untested", "lines": ["No independent cohort could be scored."]}
        n_dir = sum(c.c_index > 0.5 for c in ind)
        n_rand = sum(c.beats_random for c in ind)
        n_prol = sum(c.beats_proliferation for c in ind)
        lines = [
            f"Points the right way (C > 0.5) in {n_dir} of {k} independent cohort(s).",
            f"Beats size-matched random gene sets in {n_rand} of {k}.",
            f"Still predicts after adjusting for proliferation in {n_prol} of {k}.",
        ]
        if self.meta.get("k", 0) >= 2:
            m = self.meta
            lines.append(f"Pooled HR per SD {m['hr']:.2f} ({m['lo']:.2f}-{m['hi']:.2f}), "
                         f"heterogeneity I² = {100 * m['i2']:.0f}%.")
        gaps = [s for c in ind for s in c.subgroups if s.get("flag")]
        for s in gaps:
            lines.append(f"{s['cohort']}: C-index {s['c']:.2f} in {s['attribute']}={s['group']} vs "
                         f"{s['ref_c']:.2f} in {s['ref']} (gap {s['gap']:+.2f}, 95% CI "
                         f"{s['gap_lo']:+.2f} to {s['gap_hi']:+.2f}).")
        if n_rand == k and n_prol >= max(1, k // 2 + k % 2):
            grade = "holds up"
        elif n_rand >= max(1, k // 2 + k % 2):
            grade = "mixed"
        else:
            grade = "does not hold up"
        return {"grade": grade, "lines": lines, "subgroup_gaps": len(gaps)}


def _age_band(age: pd.Series) -> pd.Series:
    return pd.cut(age, [0, 50, 70, 200], labels=["<50", "50-69", "70+"], right=False).astype(str)


def _subgroup_audit(cohort: Cohort, score: pd.Series, attrs, n_boot: int, seed: int) -> list[dict]:
    rows = []
    clin = cohort.clin.copy()
    if "age" in clin:
        clin["age_band"] = _age_band(clin["age"])
    rng = np.random.default_rng(seed)
    for attr in attrs:
        if attr not in clin or clin[attr].nunique(dropna=True) < 2:
            continue
        if attr == "sex" and clin["sex"].value_counts().min() < 0.05 * len(clin):
            continue  # effectively single-sex cancer
        groups = {}
        for g, sub in clin.groupby(attr):
            if g in ("nan", "None", "") or pd.isna(g):
                continue
            idx = sub.index
            t, e, s = sub["time"].to_numpy(), sub["event"].to_numpy(), score[idx].to_numpy()
            ev = int(e.sum())
            c = sv.concordance(t, e, s) if ev >= 2 else np.nan
            groups[g] = (t, e, s, ev, c, len(idx))
        if len(groups) < 2:
            continue
        ref = max(groups, key=lambda g: groups[g][5])  # largest group is the reference
        rt, re_, rs, rev, rc, rn = groups[ref]
        for g, (t, e, s, ev, c, n) in groups.items():
            row = {"cohort": cohort.name, "attribute": attr, "group": str(g), "n": n,
                   "events": ev, "c": float(c), "ref": str(ref), "ref_c": float(rc)}
            if g != ref and ev >= MIN_SUBGROUP_EVENTS and rev >= MIN_SUBGROUP_EVENTS:
                diffs = []
                for _ in range(n_boot):
                    a = rng.integers(0, n, n)
                    b = rng.integers(0, rn, rn)
                    diffs.append(sv.concordance(t[a], e[a], s[a]) - sv.concordance(rt[b], re_[b], rs[b]))
                diffs = np.array(diffs)
                lo, hi = np.nanpercentile(diffs, [2.5, 97.5])
                row.update(gap=float(c - rc), gap_lo=float(lo), gap_hi=float(hi),
                           flag=bool(hi < 0))  # reliably worse than the reference group
            else:
                row.update(gap=float(c - rc) if np.isfinite(c) else np.nan,
                           note="too few events to judge" if g != ref else "reference")
            rows.append(row)
    return rows


def audit_cohort(sig: Signature, cohort: Cohort, settings: AuditSettings) -> CohortAudit:
    cohort = cohort.truncate(settings.horizon_years)
    training = cohort.name in sig.training_cohorts
    z = cohort.zscores()
    sc = score_signature(sig, z)
    res = CohortAudit(cohort.name, cohort.endpoint, cohort.platform, cohort.n, cohort.events,
                      training, sc.found, sc.missing)
    if sc.coverage < MIN_COVERAGE:
        res.skipped = f"only {len(sc.found)}/{sig.size} genes measured"
        return res
    if cohort.events < 10:
        res.skipped = f"only {cohort.events} events"
        return res

    clin = cohort.clin
    t, e = clin["time"].to_numpy(), clin["event"].to_numpy().astype(bool)
    s = sc.score[clin.index].to_numpy()
    i, j = sv.comparable_pairs(t, e)
    res.c_index = float(sv.concordance_from_pairs(s, i, j))
    boots = sv.bootstrap_concordance(t, e, s, settings.n_boot, settings.seed)
    res.c_lo, res.c_hi = (float(x) for x in np.nanpercentile(boots, [2.5, 97.5]))

    sd = s.std(ddof=1) or 1.0
    res.hr = sv.cox_fit(t, e, s / sd, ["sig"]).row("sig")

    # random-gene null (Venet et al. 2011)
    null = random_scores(sig, z[clin.index], settings.n_random, settings.seed, set(sig.genes))
    rng = np.random.default_rng(settings.seed)
    if len(i) > NULL_PAIRS:
        pick = rng.choice(len(i), NULL_PAIRS, replace=False)
        ni, nj = i[pick], j[pick]
    else:
        ni, nj = i, j
    null_c = sv.concordance_from_pairs(null, ni, nj)
    res.null_c = [float(x) for x in null_c]
    mag, null_mag = abs(res.c_index - 0.5), np.abs(null_c - 0.5)
    res.random_percentile = float((null_mag < mag).mean() * 100)
    res.random_p = float((1 + (null_mag >= mag).sum()) / (1 + len(null_mag)))

    # proliferation
    prolif, pgenes = proliferation_score(z[clin.index])
    p = prolif[clin.index].to_numpy()
    res.prolif_genes = pgenes
    res.prolif_c = float(sv.concordance_from_pairs(p, i, j))
    X = np.column_stack([s / sd, (p - p.mean()) / (p.std(ddof=1) or 1)])
    res.hr_adj_prolif = sv.cox_fit(t, e, X, ["sig", "prolif"]).row("sig")

    # clinical adjustment where available
    cols, names = [s / sd], ["sig"]
    for c in ("age", "stage"):
        if c in clin and clin[c].notna().mean() > 0.8:
            v = clin[c].astype(float)
            v = v.fillna(v.median()).to_numpy()
            cols.append((v - v.mean()) / (v.std(ddof=1) or 1))
            names.append(c)
    if len(names) > 1:
        res.adj_covariates = names[1:]
        res.hr_adj_clinical = sv.cox_fit(t, e, np.column_stack(cols), names).row("sig")

    res.subgroups = _subgroup_audit(cohort, sc.score, settings.subgroups,
                                    settings.n_boot, settings.seed)
    return res


def run_audit(sig: Signature, cohorts: list[Cohort], settings: AuditSettings | None = None) -> AuditResult:
    settings = settings or AuditSettings()
    sig_n, unmapped = sig.normalized()
    results = [audit_cohort(sig_n, c, settings) for c in cohorts]
    ind = [r for r in results if not r.training and not r.skipped and r.hr]
    meta = sv.random_effects_meta(np.log([r.hr["hr"] for r in ind]), [r.hr["se"] for r in ind]) \
        if ind else {"k": 0}
    return AuditResult(sig_n, unmapped, results, meta, settings)
