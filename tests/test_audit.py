"""End-to-end audit on synthetic cohorts where the right answer is known."""
import numpy as np
import pandas as pd

from sigaudit import AuditSettings, Cohort, Signature, run_audit
from sigaudit.scoring import score_signature


def synthetic_cohort(name, n=400, n_genes=600, seed=0, signal=1.0, prolif_signal=0.0,
                     ancestry_noise=None):
    """Genes G0..G9 carry real signal; a PCNA-like module carries proliferation."""
    rng = np.random.default_rng(seed)
    genes = [f"G{k}" for k in range(n_genes)]
    # use real symbols so HGNC normalisation keeps them
    real = ["PCNA", "MKI67", "TOP2A", "CCNB1", "BIRC5", "TP53", "ESR1", "GATA3", "ERBB2", "EGFR",
            "KRT5", "FOXA1", "PGR", "BCL2", "MYC", "CDK1", "AURKA", "BUB1", "CCNE1", "RRM2"]
    from sigaudit.genes import symbol_map
    pool = [s for s in sorted(set(symbol_map().values())) if s not in real][: n_genes - len(real)]
    genes = real + pool
    prolif = rng.normal(size=n)
    risk = rng.normal(size=n)
    X = rng.normal(size=(len(genes), n))
    X[:5] += 2 * prolif                  # PCNA..BIRC5 track proliferation
    X[len(real):len(real) + 60] += 2 * prolif  # plus a wider proliferation module
    X[5:10] += risk                      # TP53..EGFR carry the true risk
    anc = rng.choice(["EUR", "AFR"], n, p=[0.7, 0.3])
    lp = signal * risk + prolif_signal * prolif
    if ancestry_noise is not None:       # signature works worse in AFR
        lp = np.where(anc == "AFR", ancestry_noise * rng.normal(size=n), lp)
    t = rng.exponential(np.exp(-lp))
    c = rng.exponential(np.median(t) * 2, n)
    clin = pd.DataFrame({"time": np.minimum(t, c), "event": (t <= c).astype(float),
                         "age": rng.normal(60, 10, n), "sex": "female", "ancestry": anc},
                        index=[f"P{k}" for k in range(n)])
    expr = pd.DataFrame(X, index=genes, columns=clin.index)
    return Cohort(name, expr, clin)


FAST = AuditSettings(n_random=200, n_boot=50, horizon_years=None)


def test_true_signature_beats_random():
    sig = Signature.from_genes(["TP53", "ESR1", "GATA3", "ERBB2", "EGFR"], "true")
    res = run_audit(sig, [synthetic_cohort("a", seed=1), synthetic_cohort("b", seed=2)], FAST)
    for c in res.cohorts:
        assert c.c_index > 0.6
        assert c.beats_random
        assert c.beats_proliferation
    assert res.verdict()["grade"] == "holds up"


def test_noise_signature_does_not_hold_up():
    sig = Signature.from_genes(["KRT5", "FOXA1", "PGR", "BCL2", "MYC"], "noise")
    res = run_audit(sig, [synthetic_cohort("a", seed=3), synthetic_cohort("b", seed=4)], FAST)
    assert res.verdict()["grade"] == "does not hold up"


def test_proliferation_proxy_is_caught():
    # outcome driven only by proliferation; signature = proliferation genes
    sig = Signature.from_genes(["MKI67", "TOP2A", "CCNB1", "BIRC5"], "prolif-proxy")
    # realistic gene count, so the meta-PCNA score averages ~20 module genes (as
    # it averages ~200 in real data) and is a sharper proliferation measure
    coh = [synthetic_cohort(x, n_genes=2000, seed=s, signal=0, prolif_signal=1.0)
           for x, s in (("a", 5), ("b", 6))]
    res = run_audit(sig, coh, FAST)
    for c in res.cohorts:
        assert c.beats_random            # it does predict outcome...
        assert not c.beats_proliferation  # ...but adds nothing beyond proliferation


def test_ancestry_gap_flagged():
    sig = Signature.from_genes(["TP53", "ESR1", "GATA3", "ERBB2", "EGFR"], "true")
    coh = synthetic_cohort("gap", n=800, seed=8, ancestry_noise=1.0)
    res = run_audit(sig, [coh], FAST)
    rows = [r for r in res.cohorts[0].subgroups if r["attribute"] == "ancestry" and r["group"] == "AFR"]
    assert rows and rows[0]["flag"]


def test_training_cohort_excluded_from_verdict():
    sig = Signature.from_genes(["TP53", "ESR1", "GATA3"], "t", training_cohorts=["a"])
    res = run_audit(sig, [synthetic_cohort("a", seed=1)], FAST)
    assert res.cohorts[0].training
    assert res.verdict()["grade"] == "untested"


def test_score_handles_missing_genes():
    coh = synthetic_cohort("a", n=50, seed=1)
    sc = score_signature(Signature.from_genes(["TP53", "NOT_A_GENE"]), coh.zscores())
    assert sc.found == ["TP53"] and sc.coverage == 0.5
