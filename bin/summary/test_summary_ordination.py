"""Genus representatives and PCoA. Expectations are constructed by hand (ties, separated clusters)."""
import sys
from pathlib import Path

import numpy as np
import pytest
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_figures as mf  # noqa: E402
import summary_ordination as so  # noqa: E402


def core_fixture():
    return pd.DataFrame({
        "species_prefix": ["a1", "a2", "a3", "b1", "b2", "c1", "n1", "n2"],
        "GENUS": ["A", "A", "A", "B", "B", "C", None, None],
        "busco_complete": [90.0, 95.0, 95.0, 80.0, 80.0, 70.0, 99.0, 98.0],
        "N50_bp": [10, 500, 900, 300, 300, 5, 1, 1],
    })


def test_pick_genus_reps_best_busco_then_n50_then_prefix():
    reps = so.pick_genus_reps(core_fixture())
    assert reps["GENUS"].is_unique and sorted(reps["GENUS"]) == ["A", "B", "C"]
    by = reps.set_index("GENUS")["species_prefix"]
    assert by["A"] == "a3"      # busco tie 95 (a2, a3): higher N50 wins
    assert by["B"] == "b1"      # busco and N50 tie: lowest species_prefix wins
    assert by["C"] == "c1"


def test_pick_genus_reps_counts_null_genus():
    reps = so.pick_genus_reps(core_fixture())
    assert reps.attrs["n_no_genus"] == 2
    assert not reps["species_prefix"].isin(["n1", "n2"]).any()


def clusters(n=20, seed=0):
    rng = np.random.default_rng(seed)
    rows, meta = [], []
    for k, ph in enumerate(["P1", "P2"]):
        for i in range(n):
            v = np.full(10, 1.0)
            v[k * 5:(k + 1) * 5] = 10.0
            v = v + rng.uniform(0, 0.5, 10)
            rows.append(v / v.sum())
            meta.append({"species_prefix": f"{ph}_{i}", "GENUS": f"{ph}g{i}", "PHYLUM": ph,
                         "CLASS": f"{ph}c{i % 2}", "ORDER": f"{ph}o{i % 4}", "busco_complete": 90.0, "N50_bp": 1})
    cols = [f"codon_{j}" for j in range(10)]
    return pd.DataFrame(meta), pd.DataFrame(rows, columns=cols).assign(species_prefix=[m["species_prefix"] for m in meta])


def test_separated_clusters_permanova_significant(tmp_path):
    reps, codon = clusters()
    res = so.run_ordination(reps, {"codon": codon}, tmp_path)
    p = res["codon"]["permanova"]["PHYLUM"]
    assert p["p"] < 0.01 and p["F"] > 10
    assert res["codon"]["n_genera"] == 40 and res["codon"]["n_dropped_nan"] == 0
    # PCo1 carries almost all the variance between two tight clusters
    assert res["codon"]["pct_variance"][0] > 50
    coords = pd.read_csv(tmp_path / "tables" / "ordination" / "codon_coords.tsv", sep="\t")
    assert {"species_prefix", "GENUS", "PCo1", "PCo2"} <= set(coords.columns)
    # the two phyla separate along PCo1
    m = coords.merge(reps[["species_prefix", "PHYLUM"]])
    a, b = (m.loc[m.PHYLUM == p_, "PCo1"] for p_ in ("P1", "P2"))
    assert (a.max() < b.min()) or (b.max() < a.min())


def test_nan_codon_row_dropped_and_reported(tmp_path):
    reps, codon = clusters()
    codon.loc[3, "codon_2"] = np.nan
    res = so.run_ordination(reps, {"codon": codon}, tmp_path)
    assert res["codon"]["n_dropped_nan"] == 1 and res["codon"]["n_genera"] == 39
    coords = pd.read_csv(tmp_path / "tables" / "ordination" / "codon_coords.tsv", sep="\t")
    assert reps.loc[3, "species_prefix"] not in set(coords["species_prefix"])


def test_one_phylum_and_singleton_groups_give_nan_not_crash(tmp_path):
    reps, codon = clusters()
    reps["PHYLUM"] = "P1"
    reps["ORDER"] = [f"o{i}" for i in range(len(reps))]       # every order has one member
    res = so.run_ordination(reps, {"codon": codon}, tmp_path)
    perm = res["codon"]["permanova"]
    assert np.isnan(perm["PHYLUM"]["p"])
    assert np.isnan(perm["ORDER"]["p"])   # singleton groups: no within-group term


def test_negative_eigenvalues_reported(tmp_path):
    # Bray-Curtis is not Euclidean: a constructed three-point set breaks the metric embedding
    reps = pd.DataFrame({"species_prefix": list("wxyz"), "GENUS": list("WXYZ"), "PHYLUM": ["A", "A", "B", "B"],
                         "CLASS": "c", "ORDER": ["o1", "o1", "o2", "o2"], "busco_complete": 1.0, "N50_bp": 1})
    m = pd.DataFrame({"species_prefix": list("wxyz"), "aa_A": [1, 0, .5, .2], "aa_B": [0, 1, .5, .3], "aa_C": [0, 0, 0, .5]})
    m[["aa_A", "aa_B", "aa_C"]] = m[["aa_A", "aa_B", "aa_C"]].div(m[["aa_A", "aa_B", "aa_C"]].sum(axis=1), axis=0)
    res = so.run_ordination(reps, {"aa": m}, tmp_path)
    assert res["aa"]["n_negative_eigen"] == 1 and 0 < res["aa"]["negative_eigen_share"] < 0.01   # eigenvalue -0.0024 of 4 axes


def test_structure_matrix_logs_positive_standardizes_drops_nan():
    reps, _ = clusters(5)
    st = pd.DataFrame({"species_prefix": reps["species_prefix"], "x": np.arange(1, 11) * 10.0,
                       "neg": np.linspace(-1, 1, 10), "allnan": np.nan,
                       "y": np.arange(10, 0, -1) * 1.0})
    st.loc[2, "y"] = np.nan
    mat = so.structure_matrix(st)
    assert list(mat.columns) == ["x", "neg", "y"] and len(mat) == 9
    assert abs(mat.mean().abs().max()) < 1e-9 and np.allclose(mat.std(ddof=0), 1.0)
    z = mat["x"].values       # x was log-transformed, then standardised
    raw = np.log(st["x"].drop(2).values)
    assert np.allclose(z, (raw - raw.mean()) / raw.std())
    # neg column has non-positive values: not logged
    r = st["neg"].drop(2).values
    assert np.allclose(mat["neg"].values, (r - r.mean()) / r.std())


def test_main_writes_ordination_outputs(tmp_path):
    from test_make_figures import synth
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    # make genera distinct per species so there are many reps
    ga = pd.read_parquet(data / "genome_all.parquet")
    ga["GENUS"] = [f"G{i % 60}" for i in range(len(ga))]
    ga.to_parquet(data / "genome_all.parquet")
    made = mf.main(data, tmp_path / "out", ("phylum",))
    assert "ordination" in made
    names = {str(m) for m in made["ordination"]}
    for n in ("codon", "aa", "structure"):
        for stem in (f"O01_{n}_pcoa_phylum", f"O02_{n}_pcoa_class"):
            assert stem in names
            for ext in ("png", "svg"):
                assert (tmp_path / "out" / "figures" / "ordination" / f"{stem}.{ext}").stat().st_size > 1000
        for t in (f"{n}_coords", f"{n}_permanova"):
            assert (tmp_path / "out" / "tables" / "ordination" / f"{t}.tsv").stat().st_size > 10
    assert not any(m.kind == "skipped" for m in made["ordination"])
    assert {m.topic for m in made["ordination"] if m.kind == "figure"} == {"ordination"}


def test_main_without_composition_tables_skips_ordination(tmp_path, capsys):
    from test_make_figures import synth
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    for n in ("gene_structure_genome", "codon_genome", "aa_genome"):
        (data / f"{n}.parquet").unlink()
    made = mf.main(data, tmp_path / "out", ("phylum",))
    assert "ordination" not in made or not made["ordination"]


def _pool_fixture():
    core = pd.DataFrame({
        "species_prefix": ["a1", "a2", "b1", "c1", "c2", "n1"],
        "GENUS": ["A", "A", "B", "C", "C", None],
        "busco_complete": [99.0, 90.0, 80.0, 70.0, 60.0, 99.0],
        "N50_bp": 1})
    cols = [f"codon_{j}" for j in range(3)]
    codon = pd.DataFrame([[.2, .3, .5]] * 6, columns=cols).assign(species_prefix=core["species_prefix"].values)
    codon = codon[codon.species_prefix != "a1"]                    # best genome of A has no codon data
    aa = pd.DataFrame([[.5, .5]] * 6, columns=["aa_A", "aa_C"]).assign(species_prefix=core["species_prefix"].values)
    st = pd.DataFrame({"species_prefix": core["species_prefix"], "x": [1., 2, 3, 4, 5, 6], "y": [2., 1, 4, 3, 6, 5]})
    st.loc[st.species_prefix == "b1", "y"] = np.nan                # the only genome of B has NaN structure stats
    return core, {"codon": codon, "aa": aa, "structure": st}


def test_reps_come_from_data_complete_genomes():
    core, mats = _pool_fixture()
    reps = so.pick_complete_reps(core, mats)
    by = reps.set_index("GENUS")["species_prefix"]
    assert by["A"] == "a2"                     # best BUSCO genome a1 lacks codon data: second best is used
    assert by["C"] == "c1"
    assert "B" not in by.index                 # no genome of B is complete
    assert reps.attrs["n_genera_no_complete_data"] == 1 and reps.attrs["n_no_genus"] == 1


def test_run_ordination_reports_pool_counts(tmp_path):
    reps, codon = clusters()
    reps.attrs.update(n_no_genus=7, n_genera_no_complete_data=3)
    res = so.run_ordination(reps, {"codon": codon}, tmp_path)
    assert res["_summary"] == {"n_no_genus": 7, "n_genera_no_complete_data": 3}
    s = pd.read_csv(tmp_path / "tables" / "ordination" / "ordination_summary.tsv", sep="\t")
    assert s.loc[0, "n_no_genus"] == 7 and s.loc[0, "n_genera_no_complete_data"] == 3


def test_permanova_table_has_real_permutation_count_and_neg_columns(tmp_path):
    reps, codon = clusters(3)       # 6 genera in groups 3+3: 20 unique permutations, so 19 run, not 999
    so.run_ordination(reps, {"codon": codon}, tmp_path)
    t = pd.read_csv(tmp_path / "tables" / "ordination" / "codon_permanova.tsv", sep="\t")
    r = t[t.grouping == "PHYLUM"].iloc[0]
    assert r.n_permutations_run == 19 and r.n_unique_permutations_log10 == pytest.approx(np.log10(20))
    assert "n_unique_permutations" not in t.columns
    assert "screening" in r.caveat
    assert {"neg_eigen_over_abs_total", "neg_eigen_over_pos_total"} <= set(t.columns)


def test_composition_matrix_uses_only_prefixed_columns():
    reps, codon = clusters(3)
    codon["extra_count"] = 5.0
    mat, _, _ = so._feature_matrix("codon", codon, reps)
    assert list(mat.columns) == [f"codon_{j}" for j in range(10)]


def test_main_survives_ordination_failure(tmp_path, capsys):
    from test_make_figures import synth
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    ga = pd.read_parquet(data / "genome_all.parquet").drop(columns=["busco_complete"])
    ga.to_parquet(data / "genome_all.parquet")
    made = mf.main(data, tmp_path / "out", ("phylum",))
    assert "F01_genomes_by_group" in made["phylum"] and "all" in made
    assert "ordination failed" in capsys.readouterr().out


def _ga_without_busco(tmp_path):
    from test_make_figures import synth
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    ga = pd.read_parquet(data / "genome_all.parquet").drop(columns=["busco_complete"])
    ga.to_parquet(data / "genome_all.parquet")
    return data


def test_run_cli_exits_1_when_ordination_fails(tmp_path, capsys):
    data = _ga_without_busco(tmp_path)
    rc = mf.run_cli(["make_figures.py", str(data), str(tmp_path / "out"), "phylum"])
    out = capsys.readouterr().out
    assert rc == 1
    assert "FAILED: ordination" in out


def test_run_cli_exits_1_when_index_fails(tmp_path, capsys, monkeypatch):
    from test_make_figures import synth
    import figure_index
    data = tmp_path / "data"
    data.mkdir()
    synth(data)

    def boom(_root):
        raise RuntimeError("index boom")
    monkeypatch.setattr(figure_index, "build_index", boom)
    rc = mf.run_cli(["make_figures.py", str(data), str(tmp_path / "out"), "phylum"])
    out = capsys.readouterr().out
    assert rc == 1 and "FAILED: figure index" in out


def test_run_cli_skipped_figure_listed_but_exit_0(tmp_path, capsys, monkeypatch):
    from test_make_figures import synth
    import figs_composition as fc
    data = tmp_path / "data"
    data.mkdir()
    synth(data)

    def boom(*a, **k):
        raise ValueError("one figure broke")
    monkeypatch.setattr(fc, "_gc3_box", boom)
    rc = mf.run_cli(["make_figures.py", str(data), str(tmp_path / "out"), "phylum"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "skipped figures:" in out and "C04_gc3" in out.split("skipped figures:")[-1]
    assert "FAILED" not in out


def test_huge_unique_permutation_count_is_written_as_log10(tmp_path, monkeypatch):
    from types import SimpleNamespace
    reps, codon = clusters(3)
    real = so.lib_ordination.permanova

    def fake(dist, groups, n_permutations, seed):
        r = real(dist, groups, n_permutations=n_permutations, seed=seed)
        return SimpleNamespace(pseudo_f=r.pseudo_f, p_value=r.p_value, n_permutations_run=r.n_permutations_run,
                               n_unique_permutations=10 ** 400, caveat=r.caveat)
    monkeypatch.setattr(so.lib_ordination, "permanova", fake)
    so.run_ordination(reps, {"codon": codon}, tmp_path)
    t = pd.read_csv(tmp_path / "tables" / "ordination" / "codon_permanova.tsv", sep="\t")
    assert np.isfinite(t["n_unique_permutations_log10"]).all()
    assert t.loc[t.grouping == "PHYLUM", "n_unique_permutations_log10"].iloc[0] == pytest.approx(400)


def test_log10_count_helper():
    assert so._log10_count(0) == 0.0 and so._log10_count(None) == 0.0
    assert so._log10_count(10 ** 400) == pytest.approx(400)


def test_run_ordination_many_genera_many_groups(tmp_path):
    rng = np.random.default_rng(5)
    n, k = 1700, 30
    ids = [f"g{i}" for i in range(n)]
    reps = pd.DataFrame({"species_prefix": ids, "GENUS": ids, "PHYLUM": [f"P{i % 5}" for i in range(n)],
                         "CLASS": [f"C{i % k}" for i in range(n)], "ORDER": [f"O{i % k}" for i in range(n)],
                         "busco_complete": 90.0, "N50_bp": 1})
    v = rng.uniform(0.5, 1.5, (n, 10))
    codon = pd.DataFrame(v / v.sum(1, keepdims=True), columns=[f"codon_{j}" for j in range(10)])
    codon["species_prefix"] = ids
    so.run_ordination(reps, {"codon": codon}, tmp_path)
    t = pd.read_csv(tmp_path / "tables" / "ordination" / "codon_permanova.tsv", sep="\t")
    assert t["n_unique_permutations_log10"].gt(100).all()
