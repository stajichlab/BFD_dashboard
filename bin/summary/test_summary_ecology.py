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


# ---- global figures and tables ----
import make_figures as mf  # noqa: E402,F401


def global_fixture(n=20):
    rng = np.random.default_rng(0)
    rows, meta, tr = [], [], []
    styles = ["wood_saprotroph", "soil_saprotroph", "plant_pathogen"]
    for k, ph in enumerate(["P1", "P2"]):
        for i in range(n):
            v = np.full(10, 1.0)
            v[k * 5:(k + 1) * 5] = 10.0
            v = v + rng.uniform(0, 0.5, 10)
            sid = f"{ph}_{i}"
            rows.append(v / v.sum())
            meta.append({"species_prefix": sid, "GENUS": f"{ph}g{i}", "PHYLUM": ph, "CLASS": f"{ph}c{i % 2}",
                         "ORDER": f"{ph}o{i % 4}", "FAMILY": f"{ph}f{i % 8}", "busco_complete": 90.0, "N50_bp": 1})
            tr.append({"species_prefix": sid, "primary_lifestyle": None if i < 3 else styles[(i + k) % 3],
                       "growth_form": "filamentous_mycelium" if i % 2 else "yeast"})
    cols = [f"codon_{j}" for j in range(10)]
    codon = pd.DataFrame(rows, columns=cols).assign(species_prefix=[m["species_prefix"] for m in meta])
    return pd.DataFrame(meta), codon, pd.DataFrame(tr)


def test_global_coords_match_existing_ordination_up_to_sign(tmp_path):
    reps, codon, traits = global_fixture()
    so.run_ordination(reps, {"codon": codon}, tmp_path / "a")
    coords = se.run_global(reps, {"codon": codon}, traits, tmp_path / "b")["codon"]
    ref = pd.read_csv(tmp_path / "a" / "tables" / "ordination" / "codon_coords.tsv", sep="\t").set_index("species_prefix")
    for ax in ("PCo1", "PCo2"):
        assert np.allclose(np.abs(coords[ax].reindex(ref.index).to_numpy()), np.abs(ref[ax].to_numpy()), atol=1e-9)


def test_coverage_reconciles(tmp_path):
    reps, codon, traits = global_fixture()
    se.run_global(reps, {"codon": codon}, traits, tmp_path)
    cov = pd.read_csv(tmp_path / "tables" / "ecology" / "ecology_traits_coverage.tsv", sep="\t").set_index("trait")
    assert (cov["n_labeled"] + cov["n_unlabeled"] == cov["n_genera"]).all()
    assert cov.loc["primary_lifestyle", "n_unlabeled"] == 6 and cov.loc["primary_lifestyle", "n_genera"] == 40
    assert (tmp_path / "figures" / "ecology" / "E01_codon_pcoa_lifestyle.png").exists()
    assert (tmp_path / "figures" / "ecology" / "E02_codon_pcoa_growthform.svg").exists()


def test_global_permanova_has_both_tests(tmp_path):
    reps, codon, traits = global_fixture()
    se.run_global(reps, {"codon": codon}, traits, tmp_path)
    t = pd.read_csv(tmp_path / "tables" / "ecology" / "ecology_permanova_global.tsv", sep="\t")
    assert list(t.columns) == ["matrix", "trait", "test", "n_genera", "n_groups", "n_dropped_rare_group", "n_dropped_no_stratum",
                               "n_strata_used", "n_unique_permutations_log10", "pseudo_F", "R2", "R2_adj", "p_value", "caveat"]
    life = t[t["trait"] == "primary_lifestyle"]
    assert set(life["test"]) == {"unrestricted", "within_order"}
    assert (life["n_genera"] == 34).all()


def _many_groups():
    rng = np.random.default_rng(3)
    vals = [f"s{i}" for i in range(4) for _ in range(10)] + [f"rare{i}" for i in range(24)]
    labels = pd.Series(vals)
    coords = pd.DataFrame({"PCo1": rng.normal(size=len(vals)), "PCo2": rng.normal(size=len(vals))})
    pooled = se.pool_groups(labels, min_group=5)
    pooled[3] = None   # one unlabeled genus
    return coords, pooled


def test_scatter_legend_outside_and_png_not_clipped(tmp_path):
    import matplotlib.pyplot as plt
    coords, pooled = _many_groups()
    real = [g for g in pooled.dropna().unique() if not g.startswith("other")]
    fig = se.scatter_labels(coords, pooled, se.lifestyle_palette(real), [10.0, 5.0], "t")
    ax = fig.axes[0]
    leg = ax.get_legend()
    fig.canvas.draw()
    assert leg.get_window_extent().x0 >= ax.get_window_extent().x1
    assert leg._ncols == 2 and max(t.get_fontsize() for t in leg.get_texts()) <= 6.5
    texts = [t.get_text() for t in leg.get_texts()]
    assert any(t.startswith("no label (1)") for t in texts)
    assert sum(int(t.rsplit("(", 1)[1].rstrip(")")) for t in texts) == len(coords)   # coloured + unlabeled = plotted
    o = mf.Out(tmp_path, "ecology")
    o.save(fig, "e_test")
    img = plt.imread(next(tmp_path.rglob("e_test.png")))[..., :3]
    nonwhite = np.where((img < 0.95).any(axis=2).any(axis=0))[0]
    assert nonwhite.max() < img.shape[1] - 3


def test_palette_distinct_and_reserved_greys_unused():
    names = [f"g{i}" for i in range(40)]
    pal = se.lifestyle_palette(names)
    assert set(pal) == set(names) and len(set(pal.values())) == 40
    assert se.NO_LABEL not in pal.values() and se.OTHER_GREY not in pal.values()


def test_missing_matrix_is_skipped(tmp_path):
    reps, codon, traits = global_fixture()
    coords = se.run_global(reps, {"codon": codon}, traits, tmp_path)
    assert set(coords) == {"codon"}
    assert not list((tmp_path / "figures" / "ecology").glob("E01_aa_*"))


# ---- clade zoom-ins ----
def clade_fixture(confounded=False, n=40):
    rng = np.random.default_rng(5)
    rows, meta, tr = [], [], []
    for i in range(n):
        order = i % 4
        style = ("x" if order < 2 else "y") if confounded else ("x" if i < n // 2 else "y")
        k = 0 if style == "x" else 1
        v = np.full(10, 1.0)
        v[k * 5:(k + 1) * 5] = 10.0
        v = v + rng.uniform(0, 0.5, 10)
        sid = f"K_{i}"
        rows.append(v / v.sum())
        meta.append({"species_prefix": sid, "GENUS": f"Kg{i}", "PHYLUM": "P", "CLASS": "K", "ORDER": f"Ko{order}",
                     "FAMILY": f"Kf{i % 8}"})
        tr.append({"species_prefix": sid, "primary_lifestyle": style, "growth_form": "yeast"})
    codon = pd.DataFrame(rows, columns=[f"codon_{j}" for j in range(10)]).assign(species_prefix=[m["species_prefix"] for m in meta])
    return pd.DataFrame(meta), codon, pd.DataFrame(tr)


def _clade(reps, level="class", name="K"):
    return se.Clade(level, name, reps["species_prefix"].tolist())


_NOSKIP = pd.DataFrame(columns=["level", "clade", "n_labeled", "n_groups_ge_min", "reason"])


def test_clade_slug_safe_and_unique():
    a, b, c = se.clade_slug("class", "Nectar/tap sp."), se.clade_slug("class", "Ünï"), se.clade_slug("order", "Ünï")
    assert a == "class_nectar_tap_sp"
    assert all(ch.isalnum() or ch == "_" for ch in a + b + c) and b != c and b.startswith("class_")


def test_clade_with_separated_lifestyles_is_significant(tmp_path):
    reps, codon, traits = clade_fixture()
    t = se.run_clades(reps, {"codon": codon}, traits, [_clade(reps)], _NOSKIP, tmp_path)
    u = t[t["test"] == "unrestricted"].iloc[0]
    assert u["p_value"] < 0.01 and u["q_BH"] >= u["p_value"]
    assert set(t["test"]) == {"unrestricted", "stratified"}
    assert list(t.columns) == ["level", "clade", "matrix", "test", "n_genera", "n_labeled", "n_groups", "n_dropped_rare_group",
                               "n_dropped_no_stratum", "n_strata_used", "n_unique_permutations_log10", "pseudo_F", "R2", "R2_adj", "p_value", "q_BH", "caveat"]


def test_confounded_clade_stratified_p_is_one(tmp_path):
    reps, codon, traits = clade_fixture(confounded=True)
    t = se.run_clades(reps, {"codon": codon}, traits, [_clade(reps)], _NOSKIP, tmp_path).set_index("test")
    assert t.loc["unrestricted", "p_value"] < 0.05
    assert t.loc["stratified", "p_value"] == 1.0


def test_q_bh_is_per_level_matrix_test(tmp_path):
    reps, codon, traits = clade_fixture()
    # same data twice at two levels: q within a level equals p (one clade per level/test group)
    clades = [_clade(reps, "class", "K"), _clade(reps, "phylum", "P")]
    t = se.run_clades(reps, {"codon": codon}, traits, clades, _NOSKIP, tmp_path)
    u = t[t["test"] == "unrestricted"]
    assert len(u) == 2 and np.allclose(u["q_BH"], u["p_value"])


def test_below_threshold_clade_in_skipped_table(tmp_path):
    reps, codon, traits = clade_fixture()
    skipped = pd.DataFrame([{"level": "class", "clade": "Tiny", "n_labeled": 3, "n_groups_ge_min": 0, "reason": "fewer than 30 labeled genera"}])
    se.run_clades(reps, {"codon": codon}, traits, [_clade(reps)], skipped, tmp_path)
    s = pd.read_csv(tmp_path / "tables" / "ecology" / "ecology_clades_skipped.tsv", sep="\t")
    assert s["clade"].tolist() == ["Tiny"] and "genera" in s["reason"][0]


def test_e03_files_exist_with_expected_names(tmp_path):
    reps, codon, traits = clade_fixture()
    se.run_clades(reps, {"codon": codon}, traits, [_clade(reps)], _NOSKIP, tmp_path)
    d = tmp_path / "figures" / "ecology"
    assert (d / "E03_class_k_codon_pcoa_lifestyle.png").exists() and (d / "E03_class_k_codon_pcoa_lifestyle.svg").exists()
    assert (tmp_path / "tables" / "ecology" / "ecology_permanova_clades.tsv").exists()


def test_structure_scaling_uses_clade_only():
    rng = np.random.default_rng(2)
    n = 30
    tab = pd.DataFrame({"species_prefix": [f"s{i}" for i in range(n)], "a": rng.normal(size=n), "b": rng.normal(size=n)})
    tab.loc[:14, "a"] += 5   # first half differs strongly
    reps = pd.DataFrame({"species_prefix": tab["species_prefix"][:15]})
    m = se.clade_matrix("structure", tab, reps)
    assert np.allclose(m.mean(), 0, atol=1e-9) and len(m) == 15
    whole = so._feature_matrix("structure", tab, pd.DataFrame({"species_prefix": tab["species_prefix"]}))[0].loc[m.index]
    assert not np.allclose(whole.mean(), 0, atol=1e-3)   # the global scaling would not centre the clade


# ---- command line, script, page ----
REPO = Path(__file__).resolve().parents[2]


def _write_inputs(tmp_path, with_aa=False, db_cols="primary_lifestyle, growth_form"):
    import duckdb
    reps, codon, traits = global_fixture()
    core = reps.assign(busco_complete=90.0, N50_bp=1, n_genes=100)
    d = tmp_path / "tables"
    d.mkdir()
    core.to_parquet(d / "genome_all.parquet")
    codon.to_parquet(d / "codon_genome.parquet")
    if with_aa:
        codon.rename(columns=lambda c: c.replace("codon_", "aa_")).to_parquet(d / "aa_genome.parquet")
    db = tmp_path / "t.duckdb"
    con = duckdb.connect(str(db))
    t = traits.rename(columns={"species_prefix": "LOCUSTAG"})
    if "primary_lifestyle" not in db_cols:
        t = t.drop(columns="primary_lifestyle")
    con.register("t_df", t)
    con.sql("CREATE TABLE v_genome_labels AS SELECT * FROM t_df")
    con.close()
    return d, db


def _run(tmp_path, d, db, extra=()):
    return se.main(["--tables", str(d), "--db", str(db), "--out", str(tmp_path / "out"), "--min-genera", "10", *extra])


def test_main_end_to_end_on_synthetic_inputs(tmp_path):
    d, db = _write_inputs(tmp_path)
    assert _run(tmp_path, d, db) == 0
    for n in se.GLOBAL_TABLES + se.CLADE_TABLES:
        assert (tmp_path / "out" / "tables" / "ecology" / f"{n}.tsv").exists(), n
    figs = list((tmp_path / "out" / "figures" / "ecology").glob("*.png"))
    assert any(f.name.startswith("E01_") for f in figs) and any(f.name.startswith("E03_") for f in figs)


def test_main_missing_matrix_skips_it(tmp_path, capsys):
    d, db = _write_inputs(tmp_path)
    assert _run(tmp_path, d, db) == 0
    assert "aa_genome" in capsys.readouterr().out
    assert not list((tmp_path / "out" / "figures" / "ecology").glob("*_aa_*"))


def test_main_old_db_exits_1_and_names_column(tmp_path, capsys):
    d, db = _write_inputs(tmp_path, db_cols="growth_form")
    assert _run(tmp_path, d, db) == 1
    assert "primary_lifestyle" in capsys.readouterr().err


def test_workflow_script_has_no_bash_source_and_parses():
    import subprocess
    p = REPO / "workflows" / "summary_ecology.sh"
    assert subprocess.run(["bash", "-n", str(p)]).returncode == 0
    txt = p.read_text()
    assert "BASH_SOURCE" not in txt and "SCRATCH" in txt
    for path in ("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_tables",
                 "/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/db/BFD.duckdb",
                 "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_ecology",
                 "atlas_env/venv/bin/python"):
        assert path in txt


def test_page_mentions_every_output_table():
    txt = (REPO / "docs" / "quarto_summary" / "17_ecology.qmd").read_text()
    for n in se.GLOBAL_TABLES + se.CLADE_TABLES:
        assert n in txt, n


# ---- final-review fixes ----
def test_few_unique_permutations_carry_a_caveat():
    d = two_cluster_dist(6)
    groups = pd.Series(["x", "x", "x", "y", "y", "y"], index=d.index)
    r = se.stratified_permanova(d, groups)
    assert r["n_unique_permutations_log10"] == pytest.approx(np.log10(20))
    assert "distinct" in r["caveat"]
    strata = pd.Series(["a", "a", "b", "b", "c", "c"], index=d.index)
    g2 = pd.Series(["x", "y", "x", "y", "x", "y"], index=d.index)
    r2 = se.stratified_permanova(d, g2, strata=strata)
    assert r2["n_unique_permutations_log10"] == pytest.approx(np.log10(8)) and "distinct" in r2["caveat"]


def test_many_unique_permutations_no_caveat():
    d = two_cluster_dist()
    groups = pd.Series(["x"] * 20 + ["y"] * 20, index=d.index)
    r = se.stratified_permanova(d, groups)
    assert r["caveat"] == "" and r["n_unique_permutations_log10"] > 3


def test_infinite_pseudo_f_gives_na_p():
    d = pd.DataFrame(np.array([[0, 0, 1, 1], [0, 0, 1, 1], [1, 1, 0, 0], [1, 1, 0, 0]], float), index=list("abcd"), columns=list("abcd"))
    r = se.stratified_permanova(d, pd.Series(["x", "x", "y", "y"], index=d.index))
    assert np.isnan(r["p"]) and "not finite" in r["caveat"]
