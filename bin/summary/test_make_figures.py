"""Fixture test: synthetic aggregated tables -> every figure/table is produced and sane.

Run on the HPCC:  python3.12 -m pytest -q bin/summary/test_make_figures.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_figures as mf  # noqa: E402
from panel import PANEL, PANEL_ACCESSIONS  # noqa: E402


def synth(root, n=240, seed=0):
    rng = np.random.default_rng(seed)
    phyla = ["Ascomycota", "Basidiomycota", "Mucoromycota", "Chytridiomycota", "Zoopagomycota"]
    ids = [f"G{i:04d}" for i in range(n)]
    ph = rng.choice(phyla, n, p=[0.5, 0.3, 0.1, 0.06, 0.04])
    # Ascomycota get real class names (so the subphylum mapping and the yeast/filamentous split are exercised):
    # yeast classes have small genomes and few genes, as in the real data; SUBPHYLUM is NULL as in the database
    asco = ["Sordariomycetes", "Eurotiomycetes", "Saccharomycetes", "Pichiomycetes", "Taphrinomycetes"]
    cls = [rng.choice(asco, p=[0.35, 0.25, 0.2, 0.15, 0.05]) if p == "Ascomycota" else f"{p[:4]}_c{rng.integers(0, 3)}"
           for p in ph]
    yeast = np.isin(cls, asco[2:])
    core = pd.DataFrame({"species_prefix": ids, "PHYLUM": ph, "SUBPHYLUM": None,
                         "CLASS": cls,
                         "ORDER": "o", "FAMILY": "f", "GENUS": "g", "SPECIES": [f"s{i}" for i in range(n)],
                         "total_length_bp": rng.integers(10_000_000, 60_000_000, n), "gc_pct": rng.uniform(35, 60, n),
                         "busco_complete": rng.uniform(80, 99, n), "N50_bp": 100000, "contig_count": 50,
                         "n_genes": rng.integers(5000, 15000, n), "mean_len": rng.uniform(400, 550, n)})
    pez = (ph == "Ascomycota") & ~yeast
    core.loc[pez, "total_length_bp"] = rng.integers(25_000_000, 60_000_000, pez.sum())
    core.loc[pez, "n_genes"] = rng.integers(9000, 15000, pez.sum())
    core.loc[yeast, "total_length_bp"] = rng.integers(9_000_000, 15_000_000, yeast.sum())
    core.loc[yeast, "n_genes"] = rng.integers(4000, 6500, yeast.sum())
    core.loc[yeast, "busco_complete"] = rng.uniform(70, 95, yeast.sum())
    core.loc[ph == "Ascomycota", "N50_bp"] = rng.integers(20_000, 2_000_000, (ph == "Ascomycota").sum())
    for c in ("n_signalp", "n_signalp_09", "n_secreted_noTM", "n_secreted_strict", "n_cazy", "n_cazy_secreted",
              "n_merops", "n_merops_secreted", "n_tf_proteins"):
        core[c] = (core["n_genes"] * rng.uniform(0.01, 0.1, n)).astype(int)
    core.to_parquet(root / "genome_core.parquet")

    accs = PANEL_ACCESSIONS + [f"PF9{i:04d}" for i in range(60)]
    rows = [(i, a, int(rng.integers(1, 20)), int(rng.integers(1, 30))) for i in ids for a in accs if rng.random() < 0.7]
    pd.DataFrame(rows, columns=["species_prefix", "pfam_acc", "n_proteins", "n_domains"]).to_parquet(root / "pfam_genome.parquet")
    pd.DataFrame({"pfam_acc": accs, "pfam_id": [f"name_{a}" for a in accs], "n_genomes": 100}).to_parquet(root / "pfam_names.parquet")
    arch = [(i, d, int(rng.integers(100, 4000))) for i in ids for d in range(1, 6)]
    pd.DataFrame(arch, columns=["species_prefix", "n_domains", "n_proteins"]).to_parquet(root / "arch_genome.parquet")
    labs = list(mf.LOC_COLORS)
    pd.DataFrame([(i, l, int(rng.integers(50, 4000))) for i in ids for l in labs],
                 columns=["species_prefix", "label", "n_proteins"]).to_parquet(root / "deeploc_genome.parquet")
    wl = ["cyto", "nucl", "extr", "mito", "plas", "E.R."]
    pd.DataFrame([(i, l, int(rng.integers(50, 4000))) for i in ids for l in wl],
                 columns=["species_prefix", "localization", "n_proteins"]).to_parquet(root / "wolfpsort_genome.parquet")
    pd.DataFrame({"species_prefix": ids, "n_proteins_idp_summary": 5000, "mean_idp_fraction": rng.uniform(0.2, 0.5, n),
                  "frac_idp_gt_0_3": rng.uniform(0.2, 0.6, n)}).to_parquet(root / "idp_genome.parquet")
    cz = [(i, c, f"{c}{k}", int(rng.integers(1, 30))) for i in ids for c in ["GH", "GT", "PL", "CE", "AA", "CBM"] for k in (1, 2, 3)]
    pd.DataFrame(cz, columns=["species_prefix", "cazy_class", "cazy_family", "n_proteins"]).to_parquet(root / "cazy_genome.parquet")
    pd.DataFrame([(i, f, int(rng.integers(1, 80))) for i in ids for f in ["Zn2Cys6", "bHLH", "bZIP", "Homeobox"]],
                 columns=["species_prefix", "tf_family", "n_proteins"]).to_parquet(root / "tf_genome.parquet")
    pd.DataFrame([(i, c, int(rng.integers(10, 300))) for i in ids for c in "ACMST"],
                 columns=["species_prefix", "merops_class", "n_proteins"]).to_parquet(root / "merops_genome.parquet")


def test_all_outputs(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    made = mf.main(data, tmp_path / "out")
    assert set(made) == set(mf.LEVELS) | {"global"}
    for level in mf.LEVELS:
        names = made[level]
        for name in ["F01_genomes_by_group", "F03_genome_size_vs_genes", "F06_domains_per_protein", "F07_pfam_prevalence",
                     "F07b_widespread_fold", "F08_pfam_variance", "F10_deeploc", "F12_wolfpsort", "F15_feature_heatmap",
                     "F17_feature_pca", "F18_cazy_classes", "F20_merops_classes", "F21_tf_families", "F23_idp_fraction"]:
            assert name in names, (level, name)
            assert (tmp_path / "out" / "figures" / level / f"{name}.png").stat().st_size > 2000
        assert sum(n.startswith("F22_") for n in names) == len(PANEL)
    t = pd.read_csv(tmp_path / "out" / "tables" / "phylum" / "T01_overview.tsv", sep="\t")
    assert t["genomes"].sum() == 240
    key = pd.read_csv(tmp_path / "out" / "tables" / "class" / "T00_group_key.tsv", sep="\t")
    assert key["genomes"].sum() == 240 and (key["phylum"] != "").all()


def test_groups_order_and_colour_by_parent():
    core = pd.DataFrame({"PHYLUM": ["A"] * 30 + ["B"] * 20, "SUBPHYLUM": ["A1"] * 15 + ["A2"] * 15 + [None] * 20,
                         "CLASS": ["a1x"] * 15 + ["a2x"] * 15 + ["bx"] * 20})
    d, levels, colors = mf.assign_groups(core, "class")
    assert levels == ["a1x", "a2x", "bx"]                      # ordered by phylum, then subphylum
    assert colors["a1x"] == colors["a2x"] != colors["bx"]      # coloured by phylum
    d, levels, _ = mf.assign_groups(core, "subphylum")
    assert levels == ["A1", "A2", "B"]                          # subphylum where defined, else phylum


def test_small_groups_pool_to_other():
    core = pd.DataFrame({"PHYLUM": [f"P{i}" for i in range(12)] * 2, "CLASS": "c"})
    d, levels, colors = mf.assign_groups(core, "phylum")
    assert levels[-1] == mf.OTHER and len(levels) == mf.TOP_PHYLA + 1
    assert colors[mf.OTHER] == mf.GREY


def test_subphylum_derived_from_class_when_null():
    core = pd.DataFrame({"PHYLUM": ["Ascomycota"] * 30 + ["Basidiomycota"] * 12,
                         "SUBPHYLUM": [None] * 42,
                         "CLASS": ["Sordariomycetes"] * 15 + ["Saccharomycetes"] * 15 + ["Agaricomycetes"] * 12})
    d, levels, _ = mf.assign_groups(core, "subphylum")
    assert set(levels) == {"Pezizomycotina", "Saccharomycotina", "Agaricomycotina"}
    core.loc[core.CLASS == "Agaricomycetes", "CLASS"] = "Notaclassomycetes"   # unmapped class keeps the phylum label
    d, levels, _ = mf.assign_groups(core, "subphylum")
    assert "Basidiomycota" in levels


def test_feature_space_outputs(tmp_path):
    """The level-free feature-space set: biplot, QC overlays, one ordination per phylum panel, Ascomycota split."""
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    made = mf.main(data, tmp_path / "out", levels=("phylum",))
    g = made["global"]
    figs = ["F17b_pca_biplot", "F17c_pca_qc", "F17d_pca_ascomycota", "F17d_pca_basidiomycota",
            "F17d_pca_mucoromycota", "F17d_pca_other_phyla", "F17e_asco_pc1", "F17f_asco_features"]
    for name in figs:
        assert name in g, name
        assert (tmp_path / "out" / "figures" / "global" / f"{name}.png").stat().st_size > 2000
    tab = tmp_path / "out" / "tables" / "global"
    sep = pd.read_csv(tab / "T06f_asco_feature_separation.tsv", sep="\t")
    # the synthetic yeasts differ by genome size and gene number: these must separate them almost perfectly
    assert sep.set_index("feature").loc[["genome_Mb", "n_genes"], "balanced_accuracy"].min() > 0.95
    sub = pd.read_csv(tab / "T06d_asco_cluster_subphylum.tsv", sep="\t").set_index("subphylum")
    assert sub["genomes"].sum() == (pd.read_parquet(data / "genome_core.parquet").PHYLUM == "Ascomycota").sum()
    b90 = pd.read_csv(tab / "T06h_asco_busco90.tsv", sep="\t")
    assert len(b90) == 4 and b90["best_cut_balanced_accuracy"].iloc[0] > 0.9
    wp = pd.read_csv(tab / "T06c_within_phylum_pca.tsv", sep="\t")
    assert set(wp.panel) == {"Ascomycota", "Basidiomycota", "Mucoromycota", mf.OTHER_PHYLA}


def test_best_threshold_exact():
    v = np.array([1, 2, 3, 10, 11, 12, np.nan])
    y = np.array([1, 1, 1, 0, 0, 0, 1], bool)
    bal, acc, thr, d = mf.best_threshold(v, y)
    assert bal == acc == 1.0 and thr == 6.5 and d == "<"
    bal, acc, thr, d = mf.best_threshold(v, ~y)
    assert bal == 1.0 and d == ">="
    assert mf.best_threshold(np.array([1, 1, 2, 2]), np.array([1, 0, 1, 0], bool))[0] == 0.5   # ties: no cut helps
    assert np.isnan(mf.best_threshold(np.array([1, 2]), np.array([1, 1], bool))[0])            # one class only
    # rare class: balanced accuracy picks the cut that finds it, not the one that ignores it
    v = np.array([0.0, 1, 2, 3] + list(range(10, 106)) + [5.0])
    y = np.array([True] * 4 + [False] * 96 + [True])
    bal, acc, thr, d = mf.best_threshold(v, y)
    assert bal == 1.0 and 5 < thr < 10 and d == "<"


def test_pca_orientation_and_reconstruction():
    rng = np.random.default_rng(3)
    x = pd.DataFrame(rng.normal(size=(200, 4)) @ rng.normal(size=(4, 4)), columns=list("abcd"))
    x = x - x.mean()
    scores, ev, load = mf.pca(x, 4)
    assert np.allclose(scores @ load.values.T, x.values)            # sign flips keep scores and loadings consistent
    for j in range(4):
        assert load.iloc[:, j].abs().idxmax() == load.iloc[:, j].idxmax()   # largest |loading| is positive
    assert abs(ev.sum() - 1) < 1e-9 and (np.diff(ev) <= 0).all()
    s2, _, l2 = mf.pca(-x, 2)                                          # orientation does not depend on input sign
    assert np.allclose(l2.abs(), load.iloc[:, :2].abs())


def test_density_valley():
    rng = np.random.default_rng(4)
    v = np.concatenate([rng.normal(-3, 0.5, 400), rng.normal(2, 0.7, 900)])
    valley, m1, m2 = mf.density_valley(v)
    assert -3.5 < m1 < -2.5 and 1.5 < m2 < 2.5 and -2 < valley < 1
    assert np.isnan(mf.density_valley(rng.normal(0, 1, 2000))[0])


def test_panel_groups_pool_and_label():
    sub = pd.DataFrame({"PHYLUM": ["X"] * 40, "CLASS": ["big"] * 20 + ["mid"] * 10 + ["tiny"] * 3 + [None] * 7,
                        "subphylum": ["S1"] * 20 + ["S2"] * 10 + ["S2"] * 3 + [None] * 7})
    lab, sections, color, marker = mf.panel_groups(sub, "subphylum")
    assert set(lab) == {"big", "mid", "X (class not set)", f"Other classes (<{mf.PANEL_MIN_N} genomes each)"}
    assert sections[0] == ("S1", ["big"]) and sections[-1][0] == "pooled"
    assert marker["big"] != marker["mid"] and len({color[g] for _, gs in sections[:-1] for g in gs}) == 3
