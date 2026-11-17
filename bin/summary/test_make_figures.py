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
    core = pd.DataFrame({"species_prefix": ids, "PHYLUM": ph, "SUBPHYLUM": [f"{p[:4]}_sub{rng.integers(0, 2)}" if p == "Ascomycota" else None for p in ph],
                         "CLASS": [f"{p[:4]}_c{rng.integers(0, 3)}" for p in ph],
                         "ORDER": [None if rng.random() < 0.1 else f"{c[:4]}_o{rng.integers(0, 3)}" for c in ph],
                         "FAMILY": "f", "GENUS": "g", "SPECIES": [f"s{i}" for i in range(n)],
                         "total_length_bp": rng.integers(10_000_000, 60_000_000, n), "gc_pct": rng.uniform(35, 60, n),
                         "busco_complete": rng.uniform(80, 99, n), "N50_bp": 100000, "contig_count": 50,
                         "n_genes": rng.integers(5000, 15000, n), "mean_len": rng.uniform(400, 550, n)})
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
    # genome_all: every genome (the representatives plus non-representatives), the table of the A01-A12 scatter plots
    m = n // 2
    ext = [f"X{i:04d}" for i in range(m)]
    eph = rng.choice(phyla, m, p=[0.5, 0.3, 0.1, 0.06, 0.04])
    allids = ids + ext
    allph = np.concatenate([ph, eph])
    na = len(allids)
    size = rng.integers(10_000_000, 60_000_000, na)
    ngn = (size / 2500 + rng.normal(0, 800, na)).astype(int)
    bc = rng.uniform(60, 99, na)
    pd.DataFrame({"species_prefix": allids, "PHYLUM": allph, "CLASS": "c", "ORDER": "o", "FAMILY": "f", "GENUS": "g",
                  "SPECIES": [f"sp{i}" for i in range(na)], "is_representative": [True] * n + [False] * m,
                  "total_length_bp": size, "gc_pct": rng.uniform(35, 60, na), "busco_complete": bc,
                  "busco_single": bc - 2, "busco_duplicated": 2.0, "busco_fragmented": rng.uniform(0, 5, na),
                  "busco_missing": 100 - bc, "N50_bp": rng.integers(20_000, 5_000_000, na),
                  "contig_count": rng.integers(5, 3000, na), "n_genes": ngn, "mean_len": rng.uniform(400, 550, na),
                  "genes_per_mb": ngn / (size / 1e6), "coding_fraction": rng.uniform(0.3, 0.75, na),
                  "median_intron_len_bp": rng.uniform(50, 300, na),
                  "median_intergenic_bp": rng.uniform(300, 3000, na)}).to_parquet(root / "genome_all.parquet")

    sp = ids
    pd.DataFrame({"species_prefix": sp, "n_genes": core["n_genes"].values,
                  "median_gene_len_bp": rng.uniform(1200, 2200, n), "median_cds_len_bp": rng.uniform(1000, 1800, n),
                  "median_protein_len_aa": rng.uniform(330, 600, n), "genes_per_mb": rng.uniform(200, 500, n),
                  "coding_fraction": rng.uniform(0.3, 0.75, n), "median_exons_per_gene": rng.uniform(1.5, 4, n),
                  "intron_frac_genes_with_intron": rng.uniform(0.2, 0.9, n),
                  "median_intron_len_bp": [np.nan if (p == "Zoopagomycota" or rng.random() < 0.05) else rng.uniform(50, 300) for p in ph],
                  "median_intergenic_bp": rng.uniform(300, 3000, n)}).to_parquet(root / "gene_structure_genome.parquet")
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    stops = {"TAA", "TAG", "TGA"}
    cod = [c for c in codons if c not in stops]
    f = rng.dirichlet(np.ones(61) * 8, n)
    pd.DataFrame(f, columns=[f"codon_{c}" for c in cod]).assign(species_prefix=sp)[["species_prefix"] + [f"codon_{c}" for c in cod]] \
        .to_parquet(root / "codon_genome.parquet")
    aas = "ACDEFGHIKLMNPQRSTVWY"
    pd.DataFrame(rng.dirichlet(np.ones(20) * 20, n), columns=[f"aa_{a}" for a in aas]).assign(species_prefix=sp) \
        [["species_prefix"] + [f"aa_{a}" for a in aas]].to_parquet(root / "aa_genome.parquet")


NEW_NAMES = ["S01_gene_len", "S02_gene_density", "S03_coding_fraction", "S04_exons", "S05_intron_len", "S06_intron_frac",
             "S07_intergenic", "C01_codon_heatmap", "C02_rscu_boxes", "C03_aa_composition", "C04_gc3", "C05_aa_vs_gc"]


def test_new_pages_make_figures(tmp_path):
    import re
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    made = mf.main(data, tmp_path / "out", ("phylum", "order"))
    for level in ("phylum", "order"):
        for name in NEW_NAMES:
            for ext in ("png", "svg"):
                assert (tmp_path / "out" / "figures" / level / f"{name}.{ext}").stat().st_size > 1000, (level, name, ext)
        assert all(re.fullmatch(r"[A-Za-z0-9_]+", n) for n in made[level])
        assert len(set(made[level])) == len(made[level])
        assert not any(getattr(m, "kind", "") == "skipped" for m in made[level])
        assert {m.topic for m in made[level] if m.startswith(("S0", "C0"))} == {"structure", "composition"}
    for p in (tmp_path / "out" / "figures").rglob("*.*"):
        assert re.fullmatch(r"[A-Za-z0-9_]+\.(png|svg)", p.name), p


def test_main_tolerates_missing_new_tables(tmp_path, capsys):
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    for n in ("gene_structure_genome", "codon_genome", "aa_genome"):
        (data / f"{n}.parquet").unlink()
    made = mf.main(data, tmp_path / "out", ("phylum",))
    assert "F01_genomes_by_group" in made["phylum"] and not any(m.startswith(("S0", "C0")) for m in made["phylum"])
    assert "skipping" in capsys.readouterr().out


def test_all_outputs(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    made = mf.main(data, tmp_path / "out")
    assert set(made) == set(mf.LEVELS) | {"all", "ordination"}
    for level, names in made.items():
        if level in ("all", "ordination"):
            continue
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


def _hue(hexc):
    import colorsys
    h = hexc.lstrip("#")
    return colorsys.rgb_to_hsv(*(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)))[0]


def _order_core():
    n = 60
    return pd.DataFrame({"PHYLUM": ["A"] * 40 + ["B"] * 20, "CLASS": "c",
                         "ORDER": ["o1"] * 15 + ["o2"] * 15 + [None] * 10 + ["o3"] * 20})


def test_order_level_groups_exist():
    d, levels, colors = mf.assign_groups(_order_core(), "order")
    assert "o1" in levels and "o3" in levels and None not in levels and not d["grp"].isna().any()
    assert set(levels) <= set(colors)


def test_shades_keep_parent_hue():
    d, levels, colors = mf.assign_groups(_order_core(), "order")
    assert colors["o1"] != colors["o2"]
    assert abs(_hue(colors["o1"]) - _hue(colors["o2"])) < 0.02
    assert abs(_hue(colors["o1"]) - _hue(colors["o3"])) > 0.02
    assert mf.shade("#0072B2", 0, 3) != mf.shade("#0072B2", 2, 3)


def test_unknown_order_kept():
    core = _order_core()
    d, levels, colors = mf.assign_groups(core, "order")
    assert len(d) == len(core)
    assert (d.loc[core["ORDER"].isna(), "grp"] == "Unknown").all()
    assert levels[-1] == "Unknown" or levels[-2:] == ["Unknown", mf.OTHER]
    assert colors["Unknown"] == mf.UNKNOWN_GREY != mf.GREY
    assert colors["Unknown"] not in (colors["o1"], colors["o2"], colors["o3"])


def test_small_phylum_visible_at_order_level():
    rows = [("Big", f"b{i}") for i in range(40) for _ in range(30)]          # 40 orders x 30 genomes
    rows += [("Small", f"s{i}") for i in range(2) for _ in range(11)]       # 2 orders x 11, below the 36th largest
    core = pd.DataFrame(rows, columns=["PHYLUM", "ORDER"])
    core["CLASS"] = "c"
    d, levels, colors = mf.assign_groups(core, "order")
    assert "s0" in levels   # the phylum quota reserves its largest group (one per phylum)
    assert d.attrs["parents"]["s0"] == "Small"


def test_adjacent_shades_separable():
    import colorsys
    core = pd.DataFrame({"PHYLUM": "A", "CLASS": "c", "ORDER": [f"o{i:02d}" for i in range(20) for _ in range(12)]})
    d, levels, colors = mf.assign_groups(core, "order")
    lum = [colorsys.rgb_to_hls(*(int(colors[g][k:k + 2], 16) / 255 for k in (1, 3, 5)))[1] for g in levels]
    assert min(abs(a - b) for a, b in zip(lum, lum[1:])) >= 0.05
    assert max(abs(_hue(colors[g]) - _hue(colors[levels[0]])) for g in levels) < 0.02


def test_groups_order_and_colour_by_parent():
    core = pd.DataFrame({"PHYLUM": ["A"] * 30 + ["B"] * 20, "SUBPHYLUM": ["A1"] * 15 + ["A2"] * 15 + [None] * 20,
                         "CLASS": ["a1x"] * 15 + ["a2x"] * 15 + ["bx"] * 20})
    d, levels, colors = mf.assign_groups(core, "class")
    assert levels == ["a1x", "a2x", "bx"]                      # ordered by phylum, then subphylum
    assert colors["a1x"] != colors["a2x"] and colors["bx"] not in (colors["a1x"], colors["a2x"])  # shades of the phylum hue
    assert abs(_hue(colors["a1x"]) - _hue(colors["a2x"])) < 0.02
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


def _all_frame(n=300, n_rep=40, seed=1):
    rng = np.random.default_rng(seed)
    size = rng.integers(10_000_000, 60_000_000, n)
    return pd.DataFrame({"PHYLUM": rng.choice(["Ascomycota", "Basidiomycota", "Mucoromycota"], n, p=[.5, .4, .1]),
                         "is_representative": [True] * n_rep + [False] * (n - n_rep), "total_length_bp": size,
                         "n_genes": (size / 2500 + rng.normal(0, 500, n)).astype(int)})


def test_scatter_all_uses_every_genome(tmp_path):
    from figs_scatter import scatter_all
    df = _all_frame()
    o = mf.Out(tmp_path, "all")
    assert scatter_all(df, "total_length_bp", "n_genes", {"Ascomycota": "#0072B2"}, o, "A07_genes_vs_genome_size", logx=True, logy=True)
    m = [x for x in o.made if x == "A07_genes_vs_genome_size"][0]
    assert m.kind == "figure" and m.scope == "all_genomes"
    assert m.stats["n"] == 300 and -1 <= m.stats["pearson_r"] <= 1 and m.stats["pearson_r"] > 0.5
    assert (tmp_path / "figures" / "all" / "A07_genes_vs_genome_size.png").stat().st_size > 1000


def test_scatter_all_missing_column_is_skipped(tmp_path):
    from figs_scatter import scatter_all
    o = mf.Out(tmp_path, "all")
    assert not scatter_all(_all_frame(), "total_length_bp", "busco_missing", {}, o, "A11_x")
    assert o.made[-1].kind == "skipped"


def test_all_and_regression_outputs_via_main(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    made = mf.main(data, tmp_path / "out", ("phylum",))
    names = {m: m for m in made["all"]}
    for n in ("A01_busco_vs_genome_size", "A07_genes_vs_genome_size", "A12_n50_vs_n_genes"):
        assert names[n].kind == "figure" and names[n].stats["n"] == 360
    assert not any(m.kind == "skipped" for m in made["all"])
    assert {"R01_size_vs_genes", "R02_density_vs_intron_len"} <= set(made["phylum"])
    for t in ("T_anova", "T_regression"):
        d = pd.read_csv(tmp_path / "out" / "tables" / "phylum" / f"{t}.tsv", sep="\t")
        assert len(d) > 3
    assert (tmp_path / "out" / "figures" / "all" / "A01_busco_vs_genome_size.png").exists()


def test_phylum_colors_cover_phyla_that_get_points():
    from figs_scatter import phylum_colors, PALETTE, GREY, POINTS_BELOW
    ph = ["Big"] * 500 + [f"Small{i}" for i in range(10) for _ in range(10 + i)]
    cols = phylum_colors(pd.DataFrame({"PHYLUM": ph}))
    counts = pd.Series(ph).value_counts()
    drawn = [p for p in counts.index if counts[p] < POINTS_BELOW]
    assert all(p in cols for p in drawn)
    top = drawn[:8]
    assert all(cols[p] in PALETTE for p in top) and len({cols[p] for p in top}) == 8
    assert all(cols[p] == GREY for p in drawn[8:]) and cols["Big"] == GREY


def test_scatter_all_title_marks_log10(tmp_path):
    from figs_scatter import scatter_all
    o = mf.Out(tmp_path, "all")
    scatter_all(_all_frame(), "total_length_bp", "n_genes", {}, o, "A07_x", logx=True, logy=True)
    assert [m for m in o.made if m == "A07_x"][0].stats["title"].count("(log10)") == 1
