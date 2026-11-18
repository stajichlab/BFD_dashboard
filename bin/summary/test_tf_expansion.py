"""Tests for the TF expansion and contraction page (tf_expansion.py). Fixtures are synthetic: no real parquet is read.

Run on the HPCC:  python3.12 -m pytest -q -W error::UserWarning bin/summary/test_tf_expansion.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tf_expansion as te  # noqa: E402
from stats import bh_adjust  # noqa: E402

FAM = "pf_Fam A"
FAM2 = "pf_Fam B"
TEST_COLS = ["family", "lineage", "parent", "n_genera", "n_rest", "median_resid", "effect", "p", "p_adj", "direction",
             "n_agree", "flagged", "mean_count", "mean_count_rest", "raw_agree", "flagged_brief"]


def make_df(spec, seed=0, slope=0.8, sd=0.3, shift=None, fams=(FAM,)):
    """spec: {phylum: {order: n_genera}}. count = expm1(0.5 + slope * log(N) + noise); shift: {order: SD units added}."""
    rng = np.random.default_rng(seed)
    rows = []
    for ph, orders in spec.items():
        for od, n in orders.items():
            for i in range(n):
                tot = float(np.exp(rng.normal(9.0, 0.5)))
                row = {"genus": f"{od}_{i}", "PHYLUM": ph, "CLASS": f"C_{od}", "ORDER": od, "busco_complete": 95.0,
                       "n_proteins_total": tot}
                for f in fams:
                    y = 0.5 + slope * np.log(tot) + rng.normal(0, sd) + (shift or {}).get(od, 0.0) * sd
                    row[f] = float(np.expm1(y))
                rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- scaling
def test_fit_scaling_recovers_known_slope():
    rng = np.random.default_rng(1)
    n = rng.uniform(5000, 15000, 60)
    df = pd.DataFrame({"genus": [f"g{i}" for i in range(60)], "PHYLUM": "P", "ORDER": "O", "n_proteins_total": n,
                       FAM: np.expm1(np.log(2.0) + 0.8 * np.log(n))})
    out = te.fit_scaling(df, FAM)
    assert out.attrs["slope"] == pytest.approx(0.8, abs=1e-9)
    assert out.attrs["intercept"] == pytest.approx(np.log(2.0), abs=1e-9)
    assert np.abs(out["resid"]).max() < 1e-9
    assert "resid" not in df.columns          # the input is not modified


def test_fit_scaling_drops_nan_proteome_and_counts_it():
    df = make_df({"P": {"O": 12}})
    df.loc[0, "n_proteins_total"] = np.nan
    df.loc[1, "n_proteins_total"] = 0
    out = te.fit_scaling(df, FAM)
    assert out["resid"].isna().sum() == 2
    assert out.attrs["n_dropped_no_size"] == 2
    assert np.isfinite(out.attrs["slope"])


def test_constant_family_gives_nan_not_error():
    df = make_df({"P": {"O": 10, "Q": 10}})
    df[FAM] = 0
    out = te.fit_scaling(df, FAM)
    assert out["resid"].isna().all()
    assert np.isnan(out.attrs["slope"])
    t = te.lineage_tests(df, [FAM])
    assert list(t.columns) == TEST_COLS
    assert len(t) == 0
    assert te.label(FAM) in " ".join(t.attrs["skipped_families"])


# ---------------------------------------------------------------- lineage tests
def planted():
    return make_df({"A": {"O1": 8, "O2": 8, "O3": 8, "O4": 8, "O5": 8}, "B": {"P1": 10, "P2": 10}}, shift={"O1": 1.5})


def test_lineage_test_flags_planted_expansion():
    t = te.lineage_tests(planted(), [FAM])
    assert list(t.columns) == TEST_COLS
    hit = t[t["lineage"] == "O1"].iloc[0]
    assert hit["parent"] == "A"
    assert hit["n_genera"] == 8 and hit["n_rest"] == 32
    assert hit["direction"] == "expanded"
    assert hit["effect"] > 0.5 and hit["p_adj"] < 0.05 and hit["n_agree"] >= 3
    assert bool(hit["flagged"])
    assert not t[t["lineage"] != "O1"]["flagged"].any()


def test_contraction_direction():
    df = make_df({"A": {"O1": 8, "O2": 8, "O3": 8, "O4": 8, "O5": 8}}, shift={"O1": -1.5})
    hit = te.lineage_tests(df, [FAM]).query("lineage == 'O1'").iloc[0]
    assert hit["direction"] == "contracted" and hit["effect"] < -0.5 and bool(hit["flagged"])


def test_rest_is_the_same_phylum_only():
    # phylum B has a very different level; it must not enter the comparison of a phylum A order
    df = make_df({"A": {"O1": 8, "O2": 8, "O3": 8}, "B": {"P1": 30}})
    df.loc[df["PHYLUM"] == "B", FAM] = df.loc[df["PHYLUM"] == "B", FAM] * 50
    t = te.lineage_tests(df, [FAM])
    assert (t[t["parent"] == "A"]["n_rest"] == 16).all()


def test_no_flag_below_min_genera():
    df = make_df({"A": {"O1": 4, "O2": 12, "O3": 12}}, shift={"O1": 3.0})
    t = te.lineage_tests(df, [FAM])
    assert "O1" not in set(t["lineage"])
    sk = t.attrs["skipped"]
    row = sk[sk["lineage"] == "O1"].iloc[0]
    assert row["n_genera"] == 4 and "fewer than 5" in row["note"]


def test_bh_across_all_tests():
    df = make_df({"A": {"O1": 8, "O2": 8, "O3": 8}}, shift={"O1": 1.5}, fams=(FAM, FAM2), seed=3)
    t = te.lineage_tests(df, [FAM, FAM2])
    assert len(t) == 6
    assert np.allclose(t["p_adj"], bh_adjust(t["p"].to_numpy()))
    s = t.sort_values("p")
    assert (np.diff(s["p_adj"].to_numpy()) >= -1e-12).all()


def test_deterministic_output_order():
    df = make_df({"A": {"O1": 8, "O2": 8, "O3": 8}}, shift={"O1": 1.5}, fams=(FAM, FAM2), seed=3)
    a = te.lineage_tests(df, [FAM, FAM2])
    b = te.lineage_tests(df.sample(frac=1, random_state=5).reset_index(drop=True), [FAM2, FAM])
    pd.testing.assert_frame_equal(a, b)
    key = list(zip(a["p_adj"], a["family"], a["lineage"]))
    assert key == sorted(key)


def test_class_level():
    df = make_df({"A": {"O1": 8, "O2": 8, "O3": 8}}, shift={"O1": 1.5})
    t = te.lineage_tests(df, [FAM], level="CLASS")
    assert set(t["lineage"]) == {"C_O1", "C_O2", "C_O3"}


def test_null_lineage_genera_stay_in_the_rest():
    df = make_df({"A": {"O1": 8, "O2": 8}})
    df.loc[df.index[-3:], "ORDER"] = None
    t = te.lineage_tests(df, [FAM])
    assert t[t["lineage"] == "O1"].iloc[0]["n_rest"] == 8


# ---------------------------------------------------------------- presence
def presence_df():
    rows = []
    for ph, od, n, a, b in [("A", "L1", 8, 0, 1), ("A", "L2", 40, 1, 1), ("A", "L3", 8, 1, 0), ("A", "S", 3, 0, 1)]:
        for i in range(n):
            rows.append({"genus": f"{od}{i}", "PHYLUM": ph, "CLASS": "C", "ORDER": od, "n_proteins_total": 9000.0,
                         FAM: float(a * 3), FAM2: float(b * 2)})
    return pd.DataFrame(rows)


def test_presence_flags_lost_and_gained():
    # FAM: absent in L1 (8), present in L2 (40) and L3 (8), S too small. Global prevalence = 48/59 = 0.81
    t = te.presence_flags(presence_df(), [FAM])
    assert list(t.columns) == ["family", "lineage", "parent", "global_prevalence", "lineage_prevalence", "flag"]
    r = t.set_index("lineage")
    assert r.loc["L1", "flag"] == "lost" and r.loc["L1", "lineage_prevalence"] == 0.0
    assert r.loc["L2", "flag"] == ""
    assert "S" not in r.index
    # FAM2: present in L1, L2 only
    g = te.presence_flags(presence_df().assign(**{FAM2: lambda d: np.where(d["ORDER"] == "L3", 2.0, 0.0)}), [FAM2])
    g = g.set_index("lineage")
    assert g.loc["L3", "flag"] == "gained" and g.loc["L3", "global_prevalence"] < 0.2


def test_presence_flag_thresholds():
    df = presence_df()
    df.loc[df["ORDER"] == "L1", FAM] = [3.0] * 2 + [0.0] * 6            # lineage prevalence 0.25: not lost
    t = te.presence_flags(df, [FAM]).set_index("lineage")
    assert t.loc["L1", "flag"] == ""


# ---------------------------------------------------------------- page
class FakeOut:
    """Collects what the page writes (no matplotlib files)."""
    def __init__(self, root):
        from make_figures import Out
        self.o = Out(root, "tf")


def page_inputs():
    full = make_df({"A": {"O1": 8, "O2": 8, "O3": 8, "O4": 8, "O5": 8}, "B": {"P1": 10, "P2": 10}}, shift={"O1": 1.5})
    low = make_df({"A": {"L1": 8}, "B": {"L2": 8}}, seed=9, shift={"L1": -3.0, "L2": 3.0})
    low["busco_complete"] = 40.0
    low["genus"] = "low_" + low["genus"]
    return full, pd.concat([full, low], ignore_index=True)


def test_page_writes_tables_and_figures(tmp_path):
    from make_figures import Out
    g80, g0 = page_inputs()
    o = Out(tmp_path, "tf")
    te.page_tf_expansion(g80, g0, o, families=[FAM])
    tab = tmp_path / "tables" / "tf"
    hits = pd.read_csv(tab / "T11_tf_expansion_hits.tsv", sep="\t")
    nof = pd.read_csv(tab / "T11_tf_expansion_nofilter.tsv", sep="\t")
    assert (tab / "T11_tf_presence_flags.tsv").exists()
    assert set(hits["filter"]) == {"BUSCO>=80"} and set(nof["filter"]) == {"none"}
    assert nof["note"].str.contains("sensitivity only; low-BUSCO genomes have incomplete annotation").all()
    assert hits["flagged"].any()
    for n in ("T11_tf_effect_heatmap", "T12_tf_hit_boxes"):
        assert (tmp_path / "figures" / "tf" / f"{n}.png").exists()
    side = json.loads((tmp_path / "figures" / "tf" / "T11_tf_effect_heatmap.json").read_text())
    assert side["signal_kind"] == "rank_biserial_max_flagged" and 0.5 <= side["signal"] <= 1.0


def test_sensitivity_run_differs_when_low_busco_genera_planted(tmp_path):
    from make_figures import Out
    g80, g0 = page_inputs()
    o = Out(tmp_path, "tf")
    te.page_tf_expansion(g80, g0, o, families=[FAM])
    tab = tmp_path / "tables" / "tf"
    hits = pd.read_csv(tab / "T11_tf_expansion_hits.tsv", sep="\t")
    nof = pd.read_csv(tab / "T11_tf_expansion_nofilter.tsv", sep="\t")
    assert len(nof) > len(hits)
    assert set(nof["lineage"]) - set(hits["lineage"]) >= {"L1", "L2"}


def test_page_without_flagged_hits_skips_boxes_without_error(tmp_path):
    from make_figures import Out
    df = make_df({"A": {"O1": 8, "O2": 8, "O3": 8}})
    o = Out(tmp_path, "tf")
    te.page_tf_expansion(df, df, o, families=[FAM])
    assert not (tmp_path / "figures" / "tf" / "T12_tf_hit_boxes.png").exists()
    assert (tmp_path / "tables" / "tf" / "T11_tf_expansion_hits.tsv").exists()


def test_default_families_leave_out_not_comparable(tmp_path):
    g80, g0 = page_inputs()
    g80["pf_Absent"] = 0
    g0["pf_Absent"] = 0
    g0["pf_Rare"] = 0.0
    g0.loc[:1, "pf_Rare"] = 4.0           # non-zero in 2 genera: below the limit of 3
    g80["pf_Rare"] = g0["pf_Rare"].iloc[:len(g80)].to_numpy()
    assert te.comparable_families(g0) == [FAM]


def test_direction_and_agreement_are_relative_to_the_rest_of_the_phylum():
    # the whole phylum sits above the global fit; O1 is above 0 but clearly BELOW the other orders of its phylum
    df = make_df({"A": {"O1": 8, "O2": 8, "O3": 8, "O4": 8}, "B": {"P1": 120}}, shift={"O1": 1.5, "O2": 5.0, "O3": 5.0, "O4": 5.0})
    hit = te.lineage_tests(df, [FAM]).query("lineage == 'O1'").iloc[0]
    assert hit["median_resid"] > 0
    assert hit["effect"] < -0.5 and hit["direction"] == "contracted"
    assert hit["n_agree"] >= 6


def test_zero_floor_artefact_is_not_flagged():
    # a rare family (count 0 nearly everywhere) and a lineage with a small proteome: the residuals of its zeros are higher than
    # the residuals of the zeros of large proteomes. That is the log1p floor, not an expansion: the raw mean is not above.
    rng = np.random.default_rng(4)
    rows = []
    for i in range(40):
        rows.append(("R", f"R{i}", float(np.exp(rng.normal(9.5, 0.2))), 1.0 if i < 4 else 0.0))
    for i in range(8):
        rows.append(("O1", f"S{i}", float(np.exp(rng.normal(7.5, 0.2))), 0.0))
    for i in range(8):
        rows.append(("O2", f"T{i}", float(np.exp(rng.normal(9.5, 0.2))), 0.0))
    df = pd.DataFrame(rows, columns=["ORDER", "genus", "n_proteins_total", FAM])
    df["PHYLUM"], df["CLASS"] = "A", "C"
    hit = te.lineage_tests(df, [FAM]).query("lineage == 'O1'").iloc[0]
    assert hit["effect"] > 0.5 and hit["direction"] == "expanded"
    assert hit["mean_count"] == 0.0 and hit["mean_count_rest"] > 0
    assert not hit["raw_agree"] and not hit["flagged"]
    assert hit["p_adj"] < 0.05 and bool(hit["flagged_brief"])
