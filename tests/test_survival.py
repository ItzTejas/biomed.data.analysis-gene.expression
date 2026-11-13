import numpy as np
import pandas as pd
import pytest

from sigaudit import survival as sv

lifelines = pytest.importorskip("lifelines")


def _toy(n=300, seed=1, ties=False):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 2))
    t = rng.exponential(np.exp(-(0.7 * x[:, 0] - 0.3 * x[:, 1])))
    if ties:
        t = np.round(t, 1) + 0.1
    c = rng.exponential(2.0, n)
    time = np.minimum(t, c)
    event = (t <= c).astype(int)
    return time, event, x


@pytest.mark.parametrize("ties", [False, True])
def test_concordance_matches_lifelines(ties):
    from lifelines.utils import concordance_index
    time, event, x = _toy(ties=ties)
    ours = sv.concordance(time, event, x[:, 0])
    theirs = 1 - concordance_index(time, x[:, 0], event)  # lifelines: higher score = longer survival
    assert ours == pytest.approx(theirs, abs=1e-9)


def test_concordance_many_scores_vectorised():
    time, event, x = _toy()
    i, j = sv.comparable_pairs(time, event)
    both = sv.concordance_from_pairs(x, i, j)
    assert both[0] == pytest.approx(sv.concordance(time, event, x[:, 0]))
    assert both[1] == pytest.approx(sv.concordance(time, event, x[:, 1]))


def test_cox_matches_lifelines_breslow():
    from lifelines import CoxPHFitter
    time, event, x = _toy(ties=True)
    ours = sv.cox_fit(time, event, x, ["a", "b"])
    df = pd.DataFrame({"T": time, "E": event, "a": x[:, 0], "b": x[:, 1]})
    cph = CoxPHFitter(baseline_estimation_method="breslow")
    cph.fit(df, "T", "E")
    # lifelines uses Efron's tie correction, we use Breslow's; with heavy ties
    # Breslow shrinks estimates slightly toward zero, so allow ~10%
    np.testing.assert_allclose(ours.beta, cph.params_.values, rtol=0.1)
    np.testing.assert_allclose(ours.se, cph.standard_errors_.values, rtol=0.1)
    assert np.all(np.abs(ours.beta) <= np.abs(cph.params_.values) + 1e-9)
    assert ours.converged
    assert ours.beta[0] > 0.5 and ours.beta[1] < 0


def test_cox_no_ties_exact():
    from lifelines import CoxPHFitter
    time, event, x = _toy(seed=7)
    ours = sv.cox_fit(time, event, x[:, :1], ["a"])
    df = pd.DataFrame({"T": time, "E": event, "a": x[:, 0]})
    cph = CoxPHFitter().fit(df, "T", "E")
    assert ours.beta[0] == pytest.approx(cph.params_["a"], abs=1e-6)
    assert ours.se[0] == pytest.approx(cph.standard_errors_["a"], rel=1e-4)


def test_random_effects_meta_homogeneous():
    m = sv.random_effects_meta(np.log([1.5, 1.5, 1.5]), [0.1, 0.2, 0.15])
    assert m["hr"] == pytest.approx(1.5)
    assert m["tau2"] == 0 and m["i2"] == 0
