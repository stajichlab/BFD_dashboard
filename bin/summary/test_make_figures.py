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
    core = pd.DataFrame({"species_prefix": ids, "PHYLUM": ph, "CLASS": [f"{p[:4]}_c{rng.integers(0, 3)}" for p in ph],
                         "ORDER": "o", "FAMILY": "f", "GENUS": "g", "SPECIES": [f"s{i}" for i in range(n)],
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


def test_all_outputs(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    made = mf.main(data, tmp_path / "out")
    for name in ["F01_genomes_by_phylum", "F03_genome_size_vs_genes", "F06_domains_per_protein", "F07_pfam_prevalence",
                 "F08_pfam_phylum_variance", "F10_deeploc_by_phylum", "F12_wolfpsort_by_phylum", "F15_phylum_feature_heatmap",
                 "F17_feature_pca", "F18_cazy_classes", "F20_merops_classes", "F21_tf_families", "F23_idp_fraction"]:
        assert name in made, name
        assert (tmp_path / "out" / "figures" / f"{name}.png").stat().st_size > 2000
    assert sum(n.startswith("F22_") for n in made) == len(PANEL)
    t = pd.read_csv(tmp_path / "out" / "tables" / "T01_phylum_overview.tsv", sep="\t")
    assert t["genomes"].sum() == 240


def test_phylum_grouping_collapses_rare_to_other():
    core = pd.DataFrame({"PHYLUM": [f"P{i}" for i in range(12)] * 2})
    out, levels, colors = mf.group_phyla(core)
    assert levels[-1] == "Other" and len(levels) == mf.TOP_PHYLA + 1
    assert colors["Other"] == mf.GREY
