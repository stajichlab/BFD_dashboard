"""Gene structure figures (topic 'structure') from gene_structure_genome. One value per genome; groups are
described by the distribution over genomes."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from figs_common import join_groups, positive, safe
from make_figures import box_by_group, heat, parent_legend, width

TOPIC = "structure"
FAMILIES = [
    ("S01_gene_len", "median_gene_len_bp", "Median gene length (bp)", True, "Gene length (transcript span)"),
    ("S01b_cds_len", "median_cds_len_bp", "Median CDS length (bp)", True, "CDS length"),
    ("S01c_protein_len", "median_protein_len_aa", "Median protein length (aa)", True, "Protein length"),
    ("S02_gene_density", "genes_per_mb", "Genes per Mb", False, "Gene density"),
    ("S03_coding_fraction", "coding_fraction", "Coding fraction of genome", False, "Coding fraction"),
    ("S04_exons", "median_exons_per_gene", "Median exons per gene", False, "Exons per gene"),
    ("S05_intron_len", "median_intron_len_bp", "Median intron length (bp)", True, "Intron length (genomes with introns)"),
    ("S06_intron_frac", "intron_frac_genes_with_intron", "Fraction of genes with an intron", False, "Intron-containing genes"),
    ("S07_intergenic", "median_intergenic_bp", "Median intergenic distance (bp)", True, "Intergenic distance"),
]


def _box_fig(d, col, ylabel, title, log, levels, colors, core):
    if d[col].notna().sum() == 0:
        return None
    fig, ax = plt.subplots(figsize=(width(levels), 3.8))
    box_by_group(ax, d, col, levels, colors, log=log and positive(d, col), ylabel=ylabel)
    ax.set_title(title)
    parent_legend(ax, core, loc="best")
    return fig


def _heat_fig(d, levels, colors):
    cols = [f[1] for f in FAMILIES if f[1] in d and f[1] != "median_protein_len_aa"]
    med = d.groupby("grp")[cols].median().reindex(levels)
    sd = med.std(ddof=0).replace(0, np.nan)
    z = ((med - med.mean()) / sd).fillna(0).T
    z.index = [next(f[4] for f in FAMILIES if f[1] == c) for c in z.index]
    fig, ax = plt.subplots(figsize=(width(levels, 3.2), 3.6))
    heat(ax, z, cmap="RdBu_r", vmin=-2.5, vmax=2.5, label="z-score of group median", colors=colors)
    ax.set_title("Gene structure, group medians (z-score across groups)")
    return fig


def _size_fig(d, core, levels, colors):
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    x = d["median_gene_len_bp"]
    y = d["median_intergenic_bp"]
    ok = (x > 0) & (y > 0)
    for g in levels:
        m = ok & (d.grp == g)
        ax.scatter(x[m], y[m], s=6, alpha=0.5, color=colors[g], linewidths=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Median gene length (bp)")
    ax.set_ylabel("Median intergenic distance (bp)")
    ax.set_title("Gene length vs intergenic distance")
    parent_legend(ax, core, loc="best")
    return fig


def page_gene_structure(core, levels, colors, structure, o):
    d = join_groups(core, structure)
    for name, col, ylabel, log, title in FAMILIES:
        if col not in d:
            o.skip(name, TOPIC)
            continue
        safe(o, name, TOPIC, lambda col=col, ylabel=ylabel, title=title, log=log: _box_fig(d, col, ylabel, title, log, levels, colors, core))
    safe(o, "S08_structure_heatmap", TOPIC, lambda: _heat_fig(d, levels, colors))
    safe(o, "S09_gene_vs_intergenic", TOPIC, lambda: _size_fig(d, core, levels, colors))
