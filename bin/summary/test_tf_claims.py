"""Tests for the claims table (tf_claims.py). Fixtures are synthetic: no real parquet is read.

Run on the HPCC:  python3.12 -m pytest -q -W error::UserWarning bin/summary/test_tf_claims.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tf_claims as tc  # noqa: E402

CLAIM_COLS = ["claim_id", "text", "family", "metric", "level", "group_a", "group_b", "paper_direction"]
OUT_COLS = ["claim_id", "text", "family", "metric", "group_a", "group_b", "n_genera_a", "n_genera_b", "value_a", "value_b",
            "effect", "ci_lo", "ci_hi", "p", "p_adj", "observed_direction", "paper_direction", "agrees"]
FMAP = pd.DataFrame({"pfam_acc": ["PF1", "PF2", "PF3"], "family": ["Fam One", "Fam One b", "Fam Two"],
                     "paper_family": ["P1", "P1", "P2"]})
CONFIG = Path(__file__).resolve().parents[2] / "config"


def claims(rows):
    return pd.DataFrame([dict(claim_id=r[0], text="t " + r[0], family=r[1], metric=r[2], level="genus", group_a=r[3],
                              group_b=r[4], paper_direction=r[5]) for r in rows], columns=CLAIM_COLS)


def genus_df(spec, seed=0):
    """spec: {phylum: (n_genera, mean fraction of P1)}. P1 = fam_Fam One + fam_Fam One b; n_proteins_total 1000."""
    rng = np.random.default_rng(seed)
    rows = []
    for ph, (n, mu) in spec.items():
        for i in range(n):
            tot = 1000
            c = int(max(0, rng.normal(mu * tot, 0.2 * mu * tot)))
            rows.append((f"{ph[:3]}{i}", ph, tot, c // 2, c - c // 2, int(rng.integers(0, 4)), float(rng.uniform(0.5, 1.5))))
    return pd.DataFrame(rows, columns=["genus", "PHYLUM", "n_proteins_total", "fam_Fam One", "fam_Fam One b", "fam_Fam Two",
                                       "shannon"])


# ---------------------------------------------------------------- statistics
def test_rank_biserial_known_answer():
    assert tc.rank_biserial([1, 2, 3], [4, 5, 6]) == -1.0
    assert tc.rank_biserial([4, 5, 6], [1, 2, 3]) == 1.0
    # [1,2,3] vs [3,4,5]: pairs with a > b: 0; one tie (3,3) counts 0.5. U = 0.5; r = 2U/(n_a n_b) - 1 = 1/9 - 1 = -8/9
    assert abs(tc.rank_biserial([1, 2, 3], [3, 4, 5]) - (-8 / 9)) < 1e-12


def test_bootstrap_is_deterministic():
    g = genus_df({"Ascomycota": (30, 0.02), "Mucoromycota": (12, 0.03)})
    c = claims([("c1", "P1", "fraction", "Ascomycota", "Mucoromycota", "a<b")])
    r1 = tc.evaluate_claims(g, c, n_boot=300, seed=5, fmap=FMAP).iloc[0]
    r2 = tc.evaluate_claims(g, c, n_boot=300, seed=5, fmap=FMAP).iloc[0]
    r3 = tc.evaluate_claims(g, c, n_boot=300, seed=6, fmap=FMAP).iloc[0]
    assert (r1["ci_lo"], r1["ci_hi"]) == (r2["ci_lo"], r2["ci_hi"])
    assert (r1["ci_lo"], r1["ci_hi"]) != (r3["ci_lo"], r3["ci_hi"])
    assert r1["effect"] == r3["effect"]       # the point estimate does not depend on the seed


def test_agrees_logic():
    assert tc.agreement(0.5, 0.2, 0.8, "a>b") == "yes"
    assert tc.agreement(0.5, -0.1, 0.9, "a>b") == "mixed"
    assert tc.agreement(-0.5, -0.8, -0.2, "a>b") == "no"
    assert tc.agreement(-0.5, -0.9, 0.1, "a>b") == "no"          # opposite sign is "no" even when the CI includes 0
    assert tc.agreement(-0.5, -0.8, -0.2, "a<b") == "yes"
    assert tc.agreement(-0.5, -0.8, -0.2, "absent_in_a") == "yes"   # absent in a: expects a below b
    assert tc.agreement(0.5, 0.2, 0.8, "absent_in_a") == "no"
    assert tc.agreement(np.nan, np.nan, np.nan, "a>b") == "not tested"
    assert tc.observed_direction(0.3) == "a>b" and tc.observed_direction(-0.3) == "a<b" and tc.observed_direction(0.0) == "none"


def test_evaluate_columns_and_yes_no():
    g = genus_df({"Ascomycota": (30, 0.04), "Mucoromycota": (12, 0.01)})
    c = claims([("up", "P1", "fraction", "Ascomycota", "Mucoromycota", "a>b"),
                ("wrong", "P1", "fraction", "Ascomycota", "Mucoromycota", "a<b")])
    r = tc.evaluate_claims(g, c, n_boot=200, fmap=FMAP).set_index("claim_id")
    assert set(OUT_COLS) <= set(tc.evaluate_claims(g, c, n_boot=50, fmap=FMAP).columns)
    assert r.loc["up", "agrees"] == "yes" and r.loc["wrong", "agrees"] == "no"
    assert r.loc["up", "n_genera_a"] == 30 and r.loc["up", "n_genera_b"] == 12
    assert r.loc["up", "observed_direction"] == "a>b" and r.loc["up", "effect"] > 0.9
    assert r.loc["up", "value_a"] > r.loc["up", "value_b"]
    assert r.loc["up", "source"] == "BFD"


def test_family_bridge_sums_bfd_families_and_divides_by_proteome():
    # P1 = fam_Fam One + fam_Fam One b; fraction = that sum / n_proteins_total. Medians: A = 0.020, B = 0.010.
    g = pd.DataFrame({"genus": [f"a{i}" for i in range(5)] + [f"b{i}" for i in range(5)],
                      "PHYLUM": ["Ascomycota"] * 5 + ["Mucoromycota"] * 5, "n_proteins_total": 1000,
                      "fam_Fam One": [10, 10, 10, 10, 10] + [5] * 5, "fam_Fam One b": [10, 10, 10, 10, 10] + [5] * 5,
                      "fam_Fam Two": 0, "shannon": 1.0})
    c = claims([("c", "P1", "fraction", "Ascomycota", "Mucoromycota", "a>b"), ("n", "P1", "count", "Ascomycota", "Mucoromycota", "a>b")])
    r = tc.evaluate_claims(g, c, n_boot=100, fmap=FMAP).set_index("claim_id")
    assert abs(r.loc["c", "value_a"] - 0.020) < 1e-12 and abs(r.loc["c", "value_b"] - 0.010) < 1e-12
    assert r.loc["n", "value_a"] == 20 and r.loc["n", "value_b"] == 10 and r.loc["c", "effect"] == 1.0


def test_group_with_fewer_than_5_genera_not_tested():
    g = genus_df({"Ascomycota": (30, 0.04), "Mucoromycota": (4, 0.01), "Chytridiomycota": (5, 0.01)})
    c = claims([("few", "P1", "fraction", "Ascomycota", "Mucoromycota", "a>b"),
                ("five", "P1", "fraction", "Ascomycota", "Chytridiomycota", "a>b"),
                ("none", "P1", "fraction", "Ascomycota", "Glomeromycota", "a>b")])
    r = tc.evaluate_claims(g, c, n_boot=100, fmap=FMAP).set_index("claim_id")
    assert r.loc["few", "agrees"] == "not tested" and np.isnan(r.loc["few", "effect"]) and np.isnan(r.loc["few", "value_b"])
    assert "fewer than 5" in r.loc["few", "note"] and r.loc["few", "n_genera_b"] == 4
    assert r.loc["five", "agrees"] != "not tested"
    assert r.loc["none", "agrees"] == "not tested" and r.loc["none", "n_genera_b"] == 0
    assert r["p_adj"].notna().sum() == 1      # BH over the tested claims only


def test_prevalence_claim_uses_fisher():
    # A: 1 of 5 genera has the family; B: 4 of 5. Fisher table [[1,4],[4,1]]: two-sided p = 2 * (25 + 1) / 252 = 0.2063492.
    g = pd.DataFrame({"genus": list("abcdefghij"), "PHYLUM": ["Ascomycota"] * 5 + ["Mucoromycota"] * 5,
                      "n_proteins_total": 1000, "fam_Fam Two": [1, 0, 0, 0, 0, 3, 2, 1, 5, 0], "shannon": 1.0})
    c = claims([("prev", "P2", "prevalence", "Ascomycota", "Mucoromycota", "a<b")])
    r = tc.evaluate_claims(g, c, n_boot=200, fmap=FMAP).iloc[0]
    assert abs(r["value_a"] - 0.2) < 1e-12 and abs(r["value_b"] - 0.8) < 1e-12
    assert abs(r["effect"] - (-0.6)) < 1e-12
    assert abs(r["p"] - 52 / 252) < 1e-9
    assert r["agrees"] in ("yes", "mixed")         # sign matches; CI decides


def test_edf_and_dikarya_sets_resolved():
    ph = pd.Series(["Ascomycota", "Basidiomycota", "Microsporidia", "Mucoromycota", "Chytridiomycota", None])
    m = lambda s: list(tc.resolve_group(s, ph))
    assert m("Dikarya") == [True, True, False, False, False, False]
    assert m("EDF") == [False, False, False, True, True, False]            # not Microsporidia, not unknown
    assert m("EDF-Chytridiomycota") == [False, False, False, True, False, False]
    assert m("ALL-Mucoromycota") == [True, True, True, False, True, False]
    assert m("Mucoromycota") == [False, False, False, True, False, False]
    assert m("Nonesuch") == [False] * 6


def test_missing_family_column_reported_not_raised():
    g = genus_df({"Ascomycota": (10, 0.04), "Mucoromycota": (10, 0.01)})
    c = claims([("ghost", "Not In Map", "fraction", "Ascomycota", "Mucoromycota", "a>b"),
                ("unbacked", "P2x", "fraction", "Ascomycota", "Mucoromycota", "a>b"),
                ("ok", "P1", "fraction", "Ascomycota", "Mucoromycota", "a>b")])
    fm = pd.concat([FMAP, pd.DataFrame({"pfam_acc": ["PF9"], "family": ["Fam Absent"], "paper_family": ["P2x"]})], ignore_index=True)
    r = tc.evaluate_claims(g, c, n_boot=100, fmap=fm).set_index("claim_id")
    assert r.loc["ghost", "agrees"] == "not tested" and "family" in r.loc["ghost", "note"]
    assert r.loc["unbacked", "agrees"] == "not tested" and "fam_Fam Absent" in r.loc["unbacked", "note"]
    assert r.loc["ok", "agrees"] == "yes"


def test_shannon_metric_and_bh_over_claims():
    g = genus_df({"Ascomycota": (20, 0.02), "Mucoromycota": (10, 0.02)})
    g.loc[g["PHYLUM"] == "Mucoromycota", "shannon"] += 1.0
    c = claims([("sh", "", "shannon", "Mucoromycota", "Ascomycota", "a>b"),
                ("fr", "P1", "fraction", "Ascomycota", "Mucoromycota", "a>b")])
    r = tc.evaluate_claims(g, c, n_boot=100, fmap=FMAP).set_index("claim_id")
    assert r.loc["sh", "agrees"] == "yes" and r.loc["sh", "value_a"] > r.loc["sh", "value_b"]
    p = r["p"].to_numpy()
    from stats import bh_adjust
    assert np.allclose(r["p_adj"].to_numpy(), bh_adjust(p))


def test_claims_file_is_consistent_with_family_map():
    cl = tc.load_claims(CONFIG / "tf_claims.tsv")
    assert list(cl.columns) == CLAIM_COLS and cl["claim_id"].is_unique and len(cl) >= 15
    assert set(cl["metric"]) <= {"fraction", "count", "prevalence", "shannon"}
    assert set(cl["paper_direction"]) <= {"a>b", "a<b", "absent_in_a"}
    fmap = pd.read_csv(CONFIG / "tf_family_map.tsv", sep="\t", dtype=str, comment="#", keep_default_na=False)
    pf = set(fmap["paper_family"]) - {""}
    fam = set(cl.loc[cl["metric"] != "shannon", "family"])
    assert fam <= pf, fam - pf
    head = (CONFIG / "tf_claims.tsv").read_text().splitlines()[0]
    assert head.startswith("#") and "EDF" in head and "controller" in head


# ---------------------------------------------------------------- page
def test_page_writes_tables_and_forest(tmp_path):
    import make_figures as mf
    g = genus_df({"Ascomycota": (30, 0.04), "Mucoromycota": (12, 0.01), "Chytridiomycota": (6, 0.02)})
    c = claims([("up", "P1", "fraction", "Ascomycota", "Mucoromycota", "a>b"),
                ("wrong", "P1", "fraction", "Ascomycota", "Chytridiomycota", "a<b"),
                ("few", "P1", "fraction", "Ascomycota", "Glomeromycota", "a>b")])
    paper = tc.evaluate_claims(g.iloc[:30], c, n_boot=50, fmap=FMAP).assign(source="paper")
    o = mf.Out(tmp_path, "tf")
    tc.page_tf_claims(g, c, o, paper_rows=paper, genus_df_nofilter=g.iloc[:35], fmap=FMAP, n_boot=100)
    t = pd.read_csv(tmp_path / "tables" / "tf" / "T10_tf_claims.tsv", sep="\t", keep_default_na=False)
    assert set(t["source"]) == {"BFD", "paper"} and len(t) == 6
    n = pd.read_csv(tmp_path / "tables" / "tf" / "T10_tf_claims_nofilter.tsv", sep="\t", keep_default_na=False)
    assert set(n["source"]) == {"BFD"} and len(n) == 3
    for ext in ("png", "svg"):
        assert (tmp_path / "figures" / "tf" / f"T10_tf_claims_forest.{ext}").stat().st_size > 1000
    assert any(m == "T10_tf_claims_forest" and m.kind == "figure" for m in o.made)
