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
    assert set(cl["metric"]) <= {"fraction", "count", "prevalence", "shannon", "rank", "max_other_prevalence"}
    assert {"rank", "max_other_prevalence"} <= set(cl["metric"])
    assert cl.set_index("claim_id").loc[["shannon_mucoromycota", "wrky_mucoromycota", "e2f_mucoromycota"], "text"].str.contains("(pooled comparison)", regex=False).all()
    assert set(cl["paper_direction"]) <= {"a>b", "a<b", "absent_in_a"}
    fmap = pd.read_csv(CONFIG / "tf_family_map.tsv", sep="\t", dtype=str, comment="#", keep_default_na=False)
    pf = set(fmap["paper_family"]) - {""}
    fam = set(cl.loc[~cl["metric"].isin(["shannon", "rank"]), "family"])
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
    assert set(t.loc[t["source"] == "BFD", "filter"]) == {"BUSCO>=80"}
    n = pd.read_csv(tmp_path / "tables" / "tf" / "T10_tf_claims_nofilter.tsv", sep="\t", keep_default_na=False)
    assert set(n["source"]) == {"BFD"} and len(n) == 3 and set(n["filter"]) == {"none"}
    assert n["note"].str.contains("sensitivity only; low-BUSCO genomes have incomplete annotation; an a<b result here is not "
                                  "evidence against the paper", regex=False).all()
    for ext in ("png", "svg"):
        assert (tmp_path / "figures" / "tf" / f"T10_tf_claims_forest.{ext}").stat().st_size > 1000
    assert any(m == "T10_tf_claims_forest" and m.kind == "figure" for m in o.made)


# ---------------------------------------------------------------- fix round
def test_not_comparable_family_when_bfd_detects_it_in_fewer_than_3_genera():
    g = genus_df({"Ascomycota": (10, 0.04), "Mucoromycota": (10, 0.01)})
    g["fam_Fam Two"] = 0
    g.loc[[0, 15], "fam_Fam Two"] = 2                      # non-zero in 2 genera
    c = claims([("rare", "P2", "prevalence", "Ascomycota", "Mucoromycota", "a>b"),
                ("zero", "P1", "count", "Ascomycota", "Mucoromycota", "a>b")])
    r = tc.evaluate_claims(g, c, n_boot=50, fmap=FMAP).set_index("claim_id")
    assert r.loc["rare", "agrees"] == "not comparable"
    assert r.loc["rare", "note"].startswith("BFD Pfam definition does not detect this family (non-zero in 2 genera); "
                                            "not comparable with the paper")
    assert r.loc["zero", "agrees"] == "yes" and "not comparable" not in r.loc["zero", "note"]
    g["fam_Fam Two"] = 0
    assert tc.evaluate_claims(g, c, n_boot=50, fmap=FMAP).set_index("claim_id").loc["rare", "note"].find("non-zero in 0 genera") > 0
    assert tc.MIN_NONZERO_GENERA == 3


def test_paper_family_columns_are_used_without_double_counting():
    g = genus_df({"Ascomycota": (10, 0.04), "Mucoromycota": (10, 0.01)})
    g["pf_P1"] = (g["fam_Fam One"] + g["fam_Fam One b"]) // 2      # distinct proteins: half the sum
    c = claims([("f", "P1", "fraction", "Ascomycota", "Mucoromycota", "a>b")])
    r = tc.evaluate_claims(g, c, n_boot=50, fmap=FMAP).iloc[0]
    exp = float(np.median(g.loc[g["PHYLUM"] == "Ascomycota", "pf_P1"] / 1000))
    assert abs(r["value_a"] - exp) < 1e-12 and "double counting" not in r["note"]
    r0 = tc.evaluate_claims(g.drop(columns="pf_P1"), c, n_boot=50, fmap=FMAP).iloc[0]
    assert "paper-family counts summed over Pfam families (double counting possible)" in r0["note"]
    assert r0["value_a"] > r["value_a"]


def test_no_with_ci_including_zero_gets_note():
    g = genus_df({"Ascomycota": (8, 0.02), "Mucoromycota": (8, 0.02)}, seed=3)
    c = claims([("w", "P1", "fraction", "Ascomycota", "Mucoromycota", "a>b"), ("w2", "P1", "fraction", "Mucoromycota", "Ascomycota", "a>b")])
    r = tc.evaluate_claims(g, c, n_boot=200, fmap=FMAP).set_index("claim_id")
    bad = r[(r["agrees"] == "no") & (r["ci_lo"] <= 0) & (r["ci_hi"] >= 0)]
    assert len(bad) >= 1 and bad["note"].str.contains("CI includes 0").all()


def _rank_df(shannon):
    rows = []
    for ph, (n, v) in shannon.items():
        for i in range(n):
            rows.append((f"{ph[:3]}{i}", ph, 1000, v + 0.001 * i, 0, 0, v + 0.001 * i))
    return pd.DataFrame(rows, columns=["genus", "PHYLUM", "n_proteins_total", "shannon", "fam_Fam One", "fam_Fam One b", "x"]
                        ).assign(**{"fam_Fam Two": 0})


def test_rank_claim_yes_when_first_and_all_pairs_yes_small_phyla_ignored():
    g = _rank_df({"Mucoromycota": (6, 3.0), "Ascomycota": (8, 2.0), "Basidiomycota": (8, 2.5), "Zoopagomycota": (4, 9.0)})
    c = claims([("hi", "", "rank", "Mucoromycota", "ALL-Mucoromycota", "a>b")])
    r = tc.evaluate_claims(g, c, n_boot=100, fmap=FMAP).iloc[0]
    assert r["agrees"] == "yes" and "rank 1 of 3" in r["note"] and r["n_genera_a"] == 6
    assert r["value_a"] > r["value_b"] and r["effect"] == 1.0


def test_rank_claim_no_when_another_phylum_is_higher_and_names_it():
    g = _rank_df({"Mucoromycota": (6, 3.0), "Ascomycota": (8, 2.0), "Basidiomycota": (8, 3.5)})
    c = claims([("hi", "", "rank", "Mucoromycota", "ALL-Mucoromycota", "a>b")])
    r = tc.evaluate_claims(g, c, n_boot=100, fmap=FMAP).iloc[0]
    assert r["agrees"] == "no" and "rank 2 of 3" in r["note"] and "Basidiomycota" in r["note"]


def _prev_df(prev):
    rows = []
    for ph, p in prev.items():
        for i in range(10):
            rows.append((f"{ph[:3]}{i}", ph, 1000, 1.0, 1 if i < round(p * 10) else 0, 0, 0))
    return pd.DataFrame(rows, columns=["genus", "PHYLUM", "n_proteins_total", "shannon", "fam_Fam Two", "fam_Fam One", "fam_Fam One b"])


def test_max_other_prevalence_margin():
    assert tc.EXCLUSIVE_MARGIN == 0.3
    c = claims([("ex", "P2", "max_other_prevalence", "Mucoromycota", "ALL-Mucoromycota", "a>b")])
    r = tc.evaluate_claims(_prev_df({"Mucoromycota": 1.0, "Ascomycota": 0.2, "Basidiomycota": 0.6}), c, n_boot=50, fmap=FMAP).iloc[0]
    assert r["agrees"] == "yes" and abs(r["value_a"] - 1.0) < 1e-12 and abs(r["value_b"] - 0.6) < 1e-12
    assert abs(r["effect"] - 0.4) < 1e-12 and "Basidiomycota" in r["note"]
    r = tc.evaluate_claims(_prev_df({"Mucoromycota": 1.0, "Ascomycota": 0.2, "Basidiomycota": 0.8}), c, n_boot=50, fmap=FMAP).iloc[0]
    assert r["agrees"] == "mixed" and abs(r["effect"] - 0.2) < 1e-12
    r = tc.evaluate_claims(_prev_df({"Mucoromycota": 0.5, "Ascomycota": 0.2, "Basidiomycota": 0.8}), c, n_boot=50, fmap=FMAP).iloc[0]
    assert r["agrees"] == "no"


def test_forest_draws_not_comparable_and_open_markers(tmp_path):
    import make_figures as mf
    t = pd.DataFrame({"claim_id": ["a", "b", "c", "d"], "group_a": "X", "group_b": "Y", "source": "BFD",
                      "effect": [0.5, -0.1, 0.2, np.nan], "ci_lo": [0.3, -0.4, 0.1, np.nan], "ci_hi": [0.7, 0.2, 0.3, np.nan],
                      "agrees": ["yes", "no", "not comparable", "not tested"]})
    fig = tc._forest(t)
    txt = " ".join(x.get_text() for x in fig.axes[0].get_legend().get_texts())
    assert "not comparable (BFD Pfam definition does not detect this family)" in txt
    assert "BUSCO >= 80" in fig.axes[0].get_title()
    open_markers = [l for l in fig.axes[0].lines if l.get_marker() == "o" and l.get_markerfacecolor() in ("white", "w")]
    assert len(open_markers) == 1


def test_genus_table_adds_paper_family_columns():
    import tf_data
    ga = pd.DataFrame([("A1", "G1", "P", "C", "O", 90.0, 1), ("B1", "G2", "P", "C", "O", 90.0, 1)],
                      columns=["species_prefix", "GENUS", "PHYLUM", "CLASS", "ORDER", "busco_complete", "N50_bp"])
    fam = pd.DataFrame({"species_prefix": ["A1"], "family": ["F"], "n_proteins": [4]})
    summ = pd.DataFrame({"species_prefix": ["A1", "B1"], "n_proteins_total": 100, "n_tf_proteins": [3, 0],
                         "tf_fraction": [0.03, 0.0], "n_families": [1, 0], "shannon": [0.0, np.nan]})
    pf = pd.DataFrame({"species_prefix": ["A1"], "paper_family": ["PF one"], "n_proteins": [3]})
    t = tf_data.genus_table(ga, fam, summ, paper_fam=pf)
    assert dict(zip(t["genus"], t["pf_PF one"])) == {"G1": 3, "G2": 0}
    assert "pf_PF one" not in tf_data.genus_table(ga, fam, summ).columns


# ---------------------------------------------------------------- residual fix
def test_max_other_tie_is_mixed_and_lists_all_tied_phyla():
    c = claims([("ex", "P2", "max_other_prevalence", "Mucoromycota", "ALL-Mucoromycota", "a>b")])
    r = tc.evaluate_claims(_prev_df({"Mucoromycota": 1.0, "Ascomycota": 1.0, "Basidiomycota": 1.0, "Chytridiomycota": 0.3}), c,
                           n_boot=50, fmap=FMAP).iloc[0]
    assert r["effect"] == 0.0 and r["agrees"] == "mixed"
    assert "tie at the maximum prevalence with: Ascomycota, Basidiomycota" in r["note"]
    r = tc.evaluate_claims(_prev_df({"Mucoromycota": 0.5, "Ascomycota": 0.9, "Basidiomycota": 0.2}), c, n_boot=50, fmap=FMAP).iloc[0]
    assert r["agrees"] == "no" and r["effect"] < 0 and "tie" not in r["note"]


def test_max_other_and_rank_notes_list_eligible_and_excluded_phyla():
    df = _prev_df({"Mucoromycota": 1.0, "Ascomycota": 0.2, "Basidiomycota": 0.6})
    small = pd.DataFrame({"genus": ["m1", "m2", "m3"], "PHYLUM": "Mortierellomycota", "n_proteins_total": 1000, "shannon": 1.0,
                          "fam_Fam Two": 1, "fam_Fam One": 0, "fam_Fam One b": 0})
    df = pd.concat([df, small], ignore_index=True)
    c = claims([("ex", "P2", "max_other_prevalence", "Mucoromycota", "ALL-Mucoromycota", "a>b"),
                ("rk", "", "rank", "Mucoromycota", "ALL-Mucoromycota", "a>b")])
    r = tc.evaluate_claims(df, c, n_boot=50, fmap=FMAP).set_index("claim_id")
    for k in ("ex", "rk"):
        assert "eligible phyla: Ascomycota, Basidiomycota; excluded (< 5 genera): Mortierellomycota (3)" in r.loc[k, "note"]
    assert r.loc["ex", "agrees"] == "yes"            # the small phylum does not count: the note says so


def test_mucoromyceta_group_is_a_named_set():
    ph = pd.Series(["Mucoromycota", "Mortierellomycota", "Glomeromycota", "Ascomycota", None])
    assert list(tc.resolve_group("Mucoromyceta", ph)) == [True, True, True, False, False]
    assert list(tc.resolve_group("ALL-Mucoromyceta", ph)) == [False, False, False, True, False]
    cl = tc.load_claims(CONFIG / "tf_claims.tsv").set_index("claim_id")
    for k in ("wrky_mucoromyceta", "e2f_mucoromyceta", "shannon_mucoromyceta"):
        assert cl.loc[k, "group_a"] == "Mucoromyceta" and cl.loc[k, "group_b"] == "ALL-Mucoromyceta"
    assert "alternative reading" in (CONFIG / "tf_claims.tsv").read_text()


def test_not_comparable_note_singular_plural_and_basis():
    assert "non-zero in 1 genus)" in tc.not_comparable_note(1) and "non-zero in 2 genera)" in tc.not_comparable_note(2)
    assert "non-zero in 0 genera)" in tc.not_comparable_note(0)
    assert tc.not_comparable_note(1).endswith("(genus table: BUSCO>=80)")


def test_one_nonzero_definition_can_be_passed_in():
    g = genus_df({"Ascomycota": (10, 0.04), "Mucoromycota": (10, 0.01)})
    c = claims([("x", "P1", "count", "Ascomycota", "Mucoromycota", "a>b")])
    assert tc.evaluate_claims(g, c, n_boot=50, fmap=FMAP).iloc[0]["agrees"] == "yes"
    r = tc.evaluate_claims(g, c, n_boot=50, fmap=FMAP, nonzero={"P1": 2}).iloc[0]
    assert r["agrees"] == "not comparable" and "non-zero in 2 genera" in r["note"]
