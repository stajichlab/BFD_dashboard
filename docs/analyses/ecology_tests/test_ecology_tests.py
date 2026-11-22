"""Fixture tests for 02_tests.py (synthetic data with known structure)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("et", Path(__file__).with_name("02_tests.py"))
et = importlib.util.module_from_spec(spec)
spec.loader.exec_module(et)


def synth(n_orders=30, per=12, seed=0, eco_effect=0.0, tax_sd=1.0):
    rng = np.random.default_rng(seed)
    rows = []
    for o in range(n_orders):
        base = rng.normal(0, tax_sd)
        for i in range(per):
            eco = rng.choice(["A", "B", "C"])
            y = base + eco_effect * (eco == "A") + rng.normal(0, 0.3)
            rows.append({"order": f"O{o}", "family": f"F{o}", "eco": eco, "y": y})
    return pd.DataFrame(rows)


def test_taxonomy_only_gives_no_unique_ecology():
    d = synth(eco_effect=0.0)
    r = et.partition(d.y, d.order, d.eco, n_perm=199)
    assert r["R2adj_taxonomy"] > 0.8
    assert abs(r["unique_ecology"]) < 0.02
    assert r["p"] > 0.05


def test_ecology_effect_inside_orders_is_found():
    d = synth(eco_effect=1.0)
    r = et.partition(d.y, d.order, d.eco, n_perm=199)
    assert r["unique_ecology"] > 0.05
    assert r["p"] <= 0.01


def test_shared_variance_is_ecology_that_follows_taxonomy():
    rng = np.random.default_rng(1)
    rows = []
    for o in range(40):
        eco = "A" if o % 2 == 0 else "B"              # ecology fixed within order: fully confounded with taxonomy
        for _ in range(10):
            rows.append({"order": f"O{o}", "eco": eco, "y": (1.0 if eco == "A" else 0.0) + rng.normal(0, 0.3)})
    d = pd.DataFrame(rows)
    r = et.partition(d.y, d.order, d.eco, n_perm=99)
    assert r["R2adj_ecology"] > 0.5
    assert abs(r["unique_ecology"]) < 0.02             # nothing left once order is known
    assert r["shared"] > 0.5


def test_matched_stat_equals_weighted_mean_difference_and_detects_effect():
    rng = np.random.default_rng(2)
    clade, is_a, y = [], [], []
    for c in range(8):
        for i in range(10):
            a = i < 5
            clade.append(f"F{c}")
            is_a.append(a)
            y.append(rng.normal(0, 0.2) + c * 0.5 + (0.8 if a else 0.0))
    clade, is_a, Y = np.array(clade), np.array(is_a), np.array(y)[:, None]
    obs, p, eff, agree, clades = et.matched_stat_test(Y, clade, is_a, n_perm=999)
    assert abs(obs[0] - 0.8) < 0.2
    assert p[0] <= 0.002
    assert agree[0] == 1.0
    # clade offsets must not create an effect
    Y0 = np.array([rng.normal(0, 0.2) + int(c[1:]) * 0.5 for c in clade])[:, None]
    obs0, p0, *_ = et.matched_stat_test(Y0, clade, is_a, n_perm=999)
    assert p0[0] > 0.01


def test_matched_null_p_values_are_roughly_uniform():
    rng = np.random.default_rng(3)
    clade = np.repeat([f"F{c}" for c in range(6)], 10)
    is_a = np.tile([True] * 5 + [False] * 5, 6)
    Y = rng.normal(size=(60, 400))
    _, p, *_ = et.matched_stat_test(Y, clade, is_a, n_perm=999)
    assert 0.02 < (p < 0.05).mean() < 0.09


def test_clade_design_needs_min_per_group():
    df = pd.DataFrame({"FAMILY": ["F1"] * 6 + ["F2"] * 4, "g": ["a", "a", "a", "b", "b", "b", "a", "a", "b", "b"]})
    sub, clade, is_a = et.clade_design(df, "g", "a", "b")
    assert set(clade) == {"F1"}


def test_verdict_rule():
    assert et.verdict(0.01, 0.8, 5) == "supported"
    assert et.verdict(0.01, 0.6, 5) == "mixed"
    assert et.verdict(0.2, 0.9, 5) == "no evidence"
    assert et.verdict(0.001, 1.0, 3) == "not evaluable"
