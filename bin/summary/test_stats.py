"""Known-answer tests for stats.py."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stats import group_tests, pair_regressions  # noqa: E402


def _frame(seed=0):
    rng = np.random.default_rng(seed)
    g = np.repeat(["a", "b", "c"], 40)
    shift = np.repeat([0.0, 3.0, 6.0], 40)
    return pd.DataFrame({"grp": g, "y": rng.normal(0, 1, 120) + shift, "const": 5.0})


def test_group_shift_is_significant():
    r = group_tests(_frame(), ["y"]).iloc[0]
    assert r["p_value"] < 1e-6 and r["eta2"] > 0.3 and r["effect_size"] > 0.3
    assert r["test"] == "Kruskal-Wallis" and r["n"] == 120 and r["n_groups"] == 3
    assert {"stat", "statistic", "anova_F"} <= set(group_tests(_frame(), ["y"]).columns)


def test_constant_stat_gives_nan_effect_size():
    r = group_tests(_frame(), ["const"]).iloc[0]
    assert np.isnan(r["effect_size"]) and np.isnan(r["eta2"])


def test_group_of_size_one_dropped():
    d = pd.concat([_frame(), pd.DataFrame({"grp": ["solo"], "y": [100.0], "const": 5.0})], ignore_index=True)
    r = group_tests(d, ["y"]).iloc[0]
    assert r["n_groups"] == 3 and r["n"] == 120


def test_missing_stat_column_skipped():
    assert len(group_tests(_frame(), ["y", "nope"])) == 1


def test_regression_known_slope():
    x = np.arange(1.0, 101.0)
    d = pd.DataFrame({"x": x, "y": 2 * x})
    r = pair_regressions(d, [("x", "y")], log=False).iloc[0]
    assert abs(r["slope"] - 2) < 1e-9 and abs(r["r2"] - 1) < 1e-9 and r["n"] == 100
    rl = pair_regressions(d, [("x", "y")], log=True).iloc[0]   # log-log: slope 1
    assert abs(rl["slope"] - 1) < 1e-9


def test_regression_log_drops_nonpositive_and_handles_constant():
    d = pd.DataFrame({"x": [0.0, 1, 2, 3, 4], "y": [1.0, 1, 1, 1, 1]})
    r = pair_regressions(d, [("x", "y"), ("x", "missing")]).iloc[0]
    assert r["n"] == 4 and np.isnan(r["r2"])


def test_known_answer_two_groups():
    # [1,2,3] vs [4,5,6], no ties. Rank sums 6 and 15, n=6:
    #   H = 12/(6*7)*(6^2/3 + 15^2/3) - 3*7 = 3.857142857; eps2 = H*(n+1)/(n^2-1) = H*7/35 = 0.771428571
    #   ANOVA: SSB = 13.5, SSW = 4, F = 13.5/(4/4) = 13.5; eta2 = 13.5/17.5 = 0.771428571
    d = pd.DataFrame({"grp": list("aaabbb"), "y": [1, 2, 3, 4, 5, 6]})
    r = group_tests(d, ["y"]).iloc[0]
    assert abs(r["statistic"] - 3.857142857) < 1e-6
    assert abs(r["effect_size"] - 0.771428571) < 1e-6
    assert abs(r["anova_F"] - 13.5) < 1e-6
    assert abs(r["eta2"] - 0.771428571) < 1e-6


def test_known_answer_one_zero_variance_group():
    # [5,5,5], [1,2,3], [4,6,8]; n=9. Pooled ranks: 1,2,3,4 -> 1..4; the three 5s -> 6 each; 6 -> 8; 8 -> 9.
    # Rank sums 18, 6, 21. H uncorrected = 12/(9*10)*(18^2/3+6^2/3+21^2/3) - 3*10 = 5.6.
    # Tie correction 1 - (3^3-3)/(9^3-9) = 0.966667, so H = 5.793103448; eps2 = H*10/80 = 0.724137931.
    # ANOVA: means 5,2,6, grand mean 13/3; SSB = 26, SSW = 0+2+8 = 10, F = (26/2)/(10/6) = 7.8, eta2 = 26/36.
    d = pd.DataFrame({"grp": list("aaabbbccc"), "y": [5, 5, 5, 1, 2, 3, 4, 6, 8]})
    r = group_tests(d, ["y"]).iloc[0]
    assert abs(r["statistic"] - 5.793103448) < 1e-6
    assert abs(r["effect_size"] - 0.724137931) < 1e-6
    assert abs(r["anova_F"] - 7.8) < 1e-6
    assert abs(r["eta2"] - 26 / 36) < 1e-9
    assert r["n"] == 9 and r["n_groups"] == 3


def test_regression_small_n_returns_nan_row_with_n():
    d = pd.DataFrame({"x": [1.0, 2.0], "y": [2.0, 4.0]})
    r = pair_regressions(d, [("x", "y")])
    assert len(r) == 1 and r.iloc[0]["n"] == 2
    assert r[["slope", "intercept", "r2", "p_value"]].isna().all(axis=None)


def test_bh_adjustment_known_answer():
    from stats import bh_adjust
    # p = [0.01, 0.04, 0.03], m = 3. Sorted 0.01, 0.03, 0.04 -> raw 0.03, 0.045, 0.04; running min from the top: 0.03, 0.04, 0.04.
    assert np.allclose(bh_adjust([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])


def test_bh_adjust_nan_safe_capped_and_order_preserving():
    from stats import bh_adjust
    r = bh_adjust(np.array([0.01, np.nan, 0.04]))     # m = 2: 0.01*2/1 = 0.02 and 0.04*2/2 = 0.04
    assert np.allclose(r[[0, 2]], [0.02, 0.04]) and np.isnan(r[1])
    assert np.allclose(bh_adjust([0.9, 0.8]), [0.9, 0.9])    # never above 1, monotone
    assert bh_adjust([]).shape == (0,)
    assert np.isnan(bh_adjust([np.nan])).all()
