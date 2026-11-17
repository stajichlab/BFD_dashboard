"""Fixture tests for the carbon-storage tree scripts (run on the HPCC: python3.12 -m pytest -q test_carbon_trees.py)."""
import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, HERE / fname)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


# 01 imports bfd_connect lazily (inside main), so it can be loaded without the database
ext = load("ct_extract", "01_extract_domains.py")
pi = load("ct_pi", "04_paralog_index.py")


def test_trim_window_clips_to_protein():
    assert ext.trim_window(30, 400, 450) == (20, 410)
    assert ext.trim_window(5, 100, 105) == (1, 105)


def test_short_sequences_dropped():
    assert ext.keep_lengths([400, 410, 390, 100, 395]) == [0, 1, 2, 4]


def rows(*specs):
    return [(p, loc, genus, "Cls") for p, loc, genus in specs]


def test_recent_duplications_give_low_index():
    # genus A copies and genus B copies are sister pairs; different genera are far apart
    nwk = "((a1|LA:0.01,a2|LA:0.01):1.0,(b1|LB:0.01,b2|LB:0.01):1.0);"
    root = pi.fam_lib.parse_newick(nwk)
    res = pi.paralog_index(root, rows(("a1", "LA", "A"), ("a2", "LA", "A"), ("b1", "LB", "B"), ("b2", "LB", "B")),
                           np.random.default_rng(1))["Cls"]
    assert abs(res["median_within"] - 0.02) < 1e-9
    assert res["pi"] < 0.02


def test_ancient_paralogs_give_index_near_one():
    # each genome has one copy in each of two deep paralog clades
    nwk = "((a1|LA:1,b1|LB:1):0.1,(a2|LA:1,b2|LB:1):0.1);"
    root = pi.fam_lib.parse_newick(nwk)
    res = pi.paralog_index(root, rows(("a1", "LA", "A"), ("a2", "LA", "A"), ("b1", "LB", "B"), ("b2", "LB", "B")),
                           np.random.default_rng(1))["Cls"]
    assert abs(res["median_within"] - 2.2) < 1e-9
    assert 0.9 <= res["pi"] <= 1.2   # the reference median is 2.0 to 2.2 depending on the draw


def test_tp1_verdict_rule():
    def cell(v, n=20):
        return {"n_multi": n, "pi": v}
    ok = {(f, c): cell(0.9) for f in ("TPS", "TPP") for c in ("Mucoromycetes", "Mortierellomycetes")}
    assert pi.tp1_verdict(ok).startswith("supported")
    low = dict(ok)
    low[("TPS", "Mucoromycetes")] = cell(0.3)
    assert pi.tp1_verdict(low).startswith("killed")
    mid = dict(ok)
    mid[("TPP", "Mortierellomycetes")] = cell(0.65)
    assert pi.tp1_verdict(mid).startswith("mixed")
    few = dict(ok)
    few[("TPS", "Mucoromycetes")] = cell(0.9, n=2)
    assert pi.tp1_verdict(few).startswith("not evaluable")
