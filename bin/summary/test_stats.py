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
