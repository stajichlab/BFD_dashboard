"""Ecology step: stratified PERMANOVA, BH, pooling. Expectations are built by hand."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import summary_ecology as se  # noqa: E402
import summary_ordination as so  # noqa: E402


def two_cluster_dist(n=40, seed=0):
    rng = np.random.default_rng(seed)
    pts = np.vstack([rng.normal(0, 1, (n // 2, 3)), rng.normal(6, 1, (n // 2, 3))])
    d = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1))
    idx = [f"g{i}" for i in range(n)]
    return pd.DataFrame(d, index=idx, columns=idx)


def test_unstratified_matches_existing_permanova():
    dist = two_cluster_dist()
    groups = pd.Series(["x"] * 20 + ["y"] * 20, index=dist.index)
    r = se.stratified_permanova(dist, groups)
    ref = so._permanova(dist, groups)
    assert r["F"] == pytest.approx(ref["F"], rel=1e-9)
    assert r["R2"] == pytest.approx(ref["R2"], rel=1e-9)
    assert r["p"] < 0.01


def test_permute_within_keeps_each_stratum_multiset():
    labels = np.array(list("aabbabab") + list("cdcdcdcd"))
    strata = np.array(["s1"] * 8 + ["s2"] * 8)
    rng = np.random.default_rng(1)
    for _ in range(50):
        out = se.permute_within(labels, strata, rng)
        for s in ("s1", "s2"):
            m = strata == s
            assert sorted(out[m]) == sorted(labels[m])


def test_confounded_labels_give_p_one_when_stratified():
    dist = two_cluster_dist()
    # label equals stratum: constant within every stratum
    strata = pd.Series(["s1"] * 20 + ["s2"] * 20, index=dist.index)
    groups = strata.map({"s1": "x", "s2": "y"})
    assert se.stratified_permanova(dist, groups)["p"] < 0.01
    r = se.stratified_permanova(dist, groups, strata=strata)
    assert r["p"] == 1.0


def test_fewer_than_two_informative_strata_is_na():
    dist = two_cluster_dist()
    groups = pd.Series(["x", "y"] * 20, index=dist.index)
    strata = pd.Series(["only"] + [f"u{i}" for i in range(39)], index=dist.index)
    r = se.stratified_permanova(dist, groups, strata=strata)
    assert np.isnan(r["p"])
    assert "strata" in r["caveat"]


def test_null_stratum_rows_counted_and_kept_out():
    dist = two_cluster_dist()
    groups = pd.Series(["x", "y"] * 20, index=dist.index)
    strata = pd.Series(["s1"] * 20 + ["s2"] * 20, index=dist.index, dtype=object)
    strata.iloc[:3] = None
    r = se.stratified_permanova(dist, groups, strata=strata)
    assert r["n_dropped_no_stratum"] == 3
    assert r["n"] == 37
    assert se.stratified_permanova(dist, groups)["n"] == 40


def test_single_group_is_na_not_error():
    dist = two_cluster_dist()
    groups = pd.Series(["x"] * 40, index=dist.index)
    r = se.stratified_permanova(dist, groups)
    assert np.isnan(r["p"]) and r["caveat"]


def test_bh_adjust_known_values():
    q = se.bh_adjust(pd.Series([0.01, 0.02, 0.03]))
    assert q.tolist() == pytest.approx([0.03, 0.03, 0.03])
    q = se.bh_adjust(pd.Series([0.01, np.nan, 0.03]))
    assert np.isnan(q.iloc[1]) and q.iloc[0] == pytest.approx(0.02)


def test_pool_groups():
    labels = pd.Series(["a"] * 6 + ["b"] * 4 + ["c"] * 2 + [None])
    out = se.pool_groups(labels, min_group=5)
    assert (out.iloc[:6] == "a").all()
    assert (out.iloc[6:12] == "other (2 groups)").all()
    assert pd.isna(out.iloc[12])


def _reps_for_clades():
    rows, lab = [], {}

    def add(clade, n, styles):
        for i in range(n):
            sid = f"{clade}{i}"
            rows.append({"species_prefix": sid, "PHYLUM": "P", "CLASS": clade, "ORDER": f"O{clade}", "FAMILY": f"F{clade}"})
            if styles[i % len(styles)] is not None and not (clade == "D" and i >= 2):
                lab[sid] = styles[i % len(styles)]
    add("A", 40, ["x", "y"])
    add("B", 29, ["x", "y"])
    add("C", 40, ["x"])
    add("D", 31, ["x", "y"])
    rows.append({"species_prefix": "nul", "PHYLUM": "P", "CLASS": None, "ORDER": None, "FAMILY": None})
    return pd.DataFrame(rows), pd.Series(lab)


def test_select_clades_threshold_and_reasons():
    reps, labels = _reps_for_clades()
    clades, skipped = se.select_clades(reps, labels, levels=("class",))
    assert [c.name for c in clades] == ["A"]
    assert len(clades[0].ids) == 40 and clades[0].level == "class"
    s = skipped.set_index("clade")
    assert "genera" in s.loc["B", "reason"]
    assert "lifestyles" in s.loc["C", "reason"]
    assert s.loc["D", "n_labeled"] == 2 and "genera" in s.loc["D", "reason"]
    assert list(skipped.columns) == ["level", "clade", "n_labeled", "n_groups_ge_min", "reason"]
    assert None not in set(skipped["clade"]) and skipped["clade"].notna().all()


def test_select_clades_drops_null_labels_before_counting():
    reps, labels = _reps_for_clades()
    labels = labels.copy()
    labels.loc[[f"A{i}" for i in range(12)]] = np.nan   # 28 labeled left in A
    clades, skipped = se.select_clades(reps, labels, levels=("class",))
    assert clades == []
    assert skipped.set_index("clade").loc["A", "n_labeled"] == 28


def _make_db(tmp_path, cols):
    import duckdb
    p = str(tmp_path / "t.duckdb")
    con = duckdb.connect(p)
    con.sql(f"CREATE TABLE v_genome_labels AS SELECT 'g1' AS LOCUSTAG, {cols}")
    con.close()
    return p


def test_load_traits_from_duckdb(tmp_path):
    p = _make_db(tmp_path, "'plant_pathogen' AS primary_lifestyle, 'filamentous_mycelium' AS growth_form")
    df = se.load_traits(p)
    assert list(df.columns) == ["species_prefix", "primary_lifestyle", "growth_form"]
    assert df.iloc[0].tolist() == ["g1", "plant_pathogen", "filamentous_mycelium"]


def test_load_traits_names_missing_column(tmp_path):
    p = _make_db(tmp_path, "'filamentous_mycelium' AS growth_form")
    with pytest.raises(RuntimeError, match=r"primary_lifestyle.*build_bfd_db\.sh"):
        se.load_traits(p)
