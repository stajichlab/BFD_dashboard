"""Fixture tests for stats_util (run: python3 -m pytest test_stats_util.py, or python3 test_stats_util.py)."""
import numpy as np
from scipy import stats

from stats_util import bh_qvalues, expected_coabsence


def test_bh_known_values_and_input_order():
    p = [0.01, 0.04, 0.03, 0.005]
    # sorted: .005 .01 .03 .04 -> x4/rank: .02 .02 .04 .04 -> monotone -> same; back in input order:
    assert np.allclose(bh_qvalues(p), [0.02, 0.04, 0.04, 0.02])


def test_bh_matches_scipy():
    rng = np.random.default_rng(3)
    p = rng.uniform(size=200) ** 3
    assert np.allclose(bh_qvalues(p), stats.false_discovery_control(p, method="bh"))


def test_bh_nan_and_cap():
    q = bh_qvalues([0.5, np.nan, 0.9])
    assert np.isnan(q[1]) and np.all(q[[0, 2]] <= 1.0) and np.isclose(q[2], 0.9)


def test_old_th4_formula_was_wrong():
    # the original code: sort p, scale, accumulate, assign in the unsorted row order
    p = np.array([0.01, 0.04, 0.03, 0.005])
    k = len(p)
    old = np.sort(p) * k / np.arange(1, k + 1)
    old = np.minimum.accumulate(old[::-1])[::-1]
    assert not np.allclose(old, bh_qvalues(p))


def test_bh_is_row_order_invariant():
    rng = np.random.default_rng(1)
    p = rng.uniform(size=50)
    perm = rng.permutation(50)
    assert np.allclose(bh_qvalues(p)[perm], bh_qvalues(p[perm]))


def _clade_toy():
    # clade A: 100 genomes, focal and target each absent in half, independently (25 of each combination);
    # clade B: 300 genomes with everything present. All genomes miss 0 universal families.
    f = np.r_[np.repeat([1, 1, 0, 0], 25), np.zeros(300)].astype(bool)
    t = np.r_[np.repeat([1, 0, 1, 0], 25), np.zeros(300)].astype(bool)
    grp = np.r_[np.repeat("A", 100), np.repeat("B", 300)]
    return f, t, grp, np.zeros(400)


def test_clade_wide_trait_looks_like_coabsence_unless_stratified_by_clade():
    f, t, grp, nm = _clade_toy()
    naive = expected_coabsence(f, t, nm)
    strat = expected_coabsence(f, t, nm, group=grp)
    assert naive["obs"] == 25 and naive["exp"] < 5 and naive["p_two"] < 1e-6      # spurious signal
    assert np.isclose(strat["exp"], 25.0) and strat["p_two"] > 0.5                 # explained by the clade


def test_true_coabsence_survives_stratification():
    # inside clade A the target is absent whenever the focal is, and in 5 of the 50 focal-present genomes
    f = np.r_[np.repeat([1, 0], 50), np.zeros(300)].astype(bool)
    t = f.copy()
    t[50:55] = True
    grp = np.r_[np.repeat("A", 100), np.repeat("B", 300)]
    r = expected_coabsence(f, t, np.zeros(400), group=grp)
    assert r["obs"] == 50 and np.isclose(r["exp"], 5.0) and r["p_two"] < 1e-6


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
