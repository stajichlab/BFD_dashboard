"""Scatter figures. A01-A12: all genomes (genome_all, not only representatives), coloured by phylum, written
once under figures/all/. R01-R04: regressions of per-genome structure values at the level's groups."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as st

from figs_common import safe
from make_figures import GREY, PALETTE, Out, parent_legend
from stats import group_tests, pair_regressions

TOPIC_ALL = "scatter"
TOPIC_REG = "regression"
POINTS_BELOW = 200     # phyla with fewer genomes than this get individual points; larger phyla only the density layer
TOP_PHYLA = 8
# name, x column, y column, logx, logy, x label, y label
A_FAMILY = [
    ("A01_busco_vs_genome_size", "total_length_bp", "busco_complete", True, False, "Genome size (bp)", "BUSCO complete (%)"),
    ("A02_busco_vs_n_genes", "n_genes", "busco_complete", False, False, "Protein-coding genes", "BUSCO complete (%)"),
    ("A03_busco_vs_n50", "N50_bp", "busco_complete", True, False, "N50 (bp)", "BUSCO complete (%)"),
    ("A04_busco_vs_contigs", "contig_count", "busco_complete", True, False, "Contigs", "BUSCO complete (%)"),
    ("A05_busco_vs_gc", "gc_pct", "busco_complete", False, False, "GC (%)", "BUSCO complete (%)"),
    ("A06_busco_vs_gene_density", "genes_per_mb", "busco_complete", False, False, "Genes per Mb", "BUSCO complete (%)"),
    ("A07_genes_vs_genome_size", "total_length_bp", "n_genes", True, True, "Genome size (bp)", "Protein-coding genes"),
    ("A08_intron_len_vs_gene_density", "genes_per_mb", "median_intron_len_bp", False, True, "Genes per Mb", "Median intron length (bp)"),
    ("A09_gc_vs_genome_size", "total_length_bp", "gc_pct", True, False, "Genome size (bp)", "GC (%)"),
    ("A10_busco_dup_vs_genome_size", "total_length_bp", "busco_duplicated", True, False, "Genome size (bp)", "BUSCO duplicated (%)"),
    ("A11_busco_missing_vs_contigs", "contig_count", "busco_missing", True, False, "Contigs", "BUSCO missing (%)"),
    ("A12_n50_vs_n_genes", "N50_bp", "n_genes", True, True, "N50 (bp)", "Protein-coding genes"),
]


def phylum_colors(df):
    """Colour per phylum. The palette goes by rank (genome count) to the phyla that get individual points
    (fewer than POINTS_BELOW genomes), up to TOP_PHYLA of them. Other phyla are grey: the large ones are shown
    only in the density layer."""
    counts = df["PHYLUM"].fillna("Unknown").value_counts()
    cols = {p: GREY for p in counts.index}
    small = [p for p in counts.index if counts[p] < POINTS_BELOW and p != "Unknown"]
    for i, p in enumerate(small[:TOP_PHYLA]):
        cols[p] = PALETTE[i]
    return cols


def _all_positive(v):
    return len(v) > 0 and bool((v > 0).all())


def scatter_all(df, x, y, colors_by_phylum, o, name, logx=False, logy=False, xlabel=None, ylabel=None):
    """Density (hexbin) of every genome in df plus individual points for phyla with fewer than POINTS_BELOW
    genomes, in the phylum hue. Records n, Pearson r and Spearman r (on the plotted scale) in the Made entry's
    .stats, with scope all_genomes. Returns True when saved, False when skipped (missing column, no data)."""
    def build():
        if x not in df.columns or y not in df.columns:
            return None
        d = df[[x, y, "PHYLUM"]].copy()
        d[[x, y]] = d[[x, y]].astype(float).replace([np.inf, -np.inf], np.nan)
        d = d.dropna(subset=[x, y])
        d["PHYLUM"] = d["PHYLUM"].fillna("Unknown")
        if len(d) < 3:
            return None
        lx = logx and _all_positive(d[x])
        ly = logy and _all_positive(d[y])
        vx, vy = (np.log10(d[x]) if lx else d[x]), (np.log10(d[y]) if ly else d[y])
        const = vx.nunique() < 2 or vy.nunique() < 2
        o.meta[name] = {"n": int(len(d)), "pearson_r": float("nan") if const else float(st.pearsonr(vx, vy)[0]),
                        "spearman_r": float("nan") if const else float(st.spearmanr(vx, vy)[0]),
                        "logx": bool(lx), "logy": bool(ly)}
        fig, ax = plt.subplots(figsize=(5.6, 4.4))
        if not const:
            hb = ax.hexbin(d[x], d[y], gridsize=45, bins="log", mincnt=1, cmap="Greys",
                           xscale="log" if lx else "linear", yscale="log" if ly else "linear", linewidths=0, zorder=1)
            fig.colorbar(hb, ax=ax, label="Genomes per bin", pad=0.02, fraction=0.05)
        counts = d["PHYLUM"].value_counts()
        for ph in counts.index[::-1]:
            if counts[ph] >= POINTS_BELOW:
                continue
            m = d["PHYLUM"] == ph
            ax.scatter(d.loc[m, x], d.loc[m, y], s=7, alpha=0.8, color=colors_by_phylum.get(ph, GREY), linewidths=0,
                       rasterized=True, zorder=2, label=f"{ph} (n={counts[ph]})")
        if lx:
            ax.set_xscale("log")
        if ly:
            ax.set_yscale("log")
        ax.set_xlabel(xlabel or x)
        ax.set_ylabel(ylabel or y)
        s = o.meta[name]
        tag = " (log10)" if (lx or ly) else ""
        s["title"] = f"n={s['n']} genomes (all)   Pearson r={s['pearson_r']:.2f}   Spearman r={s['spearman_r']:.2f}{tag}"
        ax.set_title(s["title"], fontsize=8.5)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=6.5, frameon=False, loc="best", title=f"Phyla with <{POINTS_BELOW} genomes", title_fontsize=6.5)
        return fig
    return safe(o, name, TOPIC_ALL, build)


def page_all_genomes(genome_all, out):
    """A01-A12 from genome_all (every genome), written under figures/all/. Returns the Made list."""
    o = Out(out, "all")
    colors = phylum_colors(genome_all)
    for name, x, y, lx, ly, xl, yl in A_FAMILY:
        scatter_all(genome_all, x, y, colors, o, name, logx=lx, logy=ly, xlabel=xl, ylabel=yl)
    print(f"all: {len(genome_all)} genomes, {len(o.made)} outputs", flush=True)
    return o.made


# ---------------------------------------------------------------- group tests and regressions (per level)
REG_LOG = [("total_length_bp", "n_genes"), ("genes_per_mb", "median_intron_len_bp"),
           ("total_length_bp", "median_intergenic_bp"), ("median_gene_len_bp", "median_intergenic_bp")]
REG_LINEAR = [("gc_pct", "genes_per_mb")]
R_FIGS = [("R01_size_vs_genes", "total_length_bp", "n_genes", True, "Genome size (bp)", "Protein-coding genes"),
          ("R02_density_vs_intron_len", "genes_per_mb", "median_intron_len_bp", True, "Genes per Mb", "Median intron length (bp)"),
          ("R03_size_vs_intergenic", "total_length_bp", "median_intergenic_bp", True, "Genome size (bp)", "Median intergenic distance (bp)"),
          ("R04_gc_vs_gene_density", "gc_pct", "genes_per_mb", False, "GC (%)", "Genes per Mb")]
TEST_COLS = ["median_gene_len_bp", "median_cds_len_bp", "median_protein_len_aa", "genes_per_mb", "coding_fraction",
             "median_exons_per_gene", "intron_frac_genes_with_intron", "median_intron_len_bp", "median_intergenic_bp",
             "gc_pct", "busco_complete", "total_length_bp", "n_genes"]


def _reg_fig(d, x, y, log, xl, yl, core, levels, colors):
    r = pair_regressions(d, [(x, y)], log=log)
    if r.empty:
        return None
    r = r.iloc[0]
    ok = (d[x] > 0) & (d[y] > 0) if log else d[x].notna() & d[y].notna()
    ok = ok & d[x].notna() & d[y].notna()
    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    for g in levels:
        m = ok & (d.grp == g)
        ax.scatter(d.loc[m, x], d.loc[m, y], s=6, alpha=0.5, color=colors[g], linewidths=0, rasterized=True)
    xs = np.linspace(d.loc[ok, x].min(), d.loc[ok, x].max(), 100)
    if log:
        ax.plot(xs, 10 ** (r["intercept"] + r["slope"] * np.log10(xs)), color="black", lw=1.2)
        ax.set_xscale("log")
        ax.set_yscale("log")
    elif not np.isnan(r["slope"]):
        ax.plot(xs, r["intercept"] + r["slope"] * xs, color="black", lw=1.2)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title(f"slope={r['slope']:.2f}{' (log-log)' if log else ''}  R$^2$={r['r2']:.2f}  n={int(r['n'])}", fontsize=9)
    parent_legend(ax, core, loc="best")
    return fig


def page_regressions(core, levels, colors, structure, o):
    """T_anova.tsv (Kruskal-Wallis + ANOVA per stat across the level's groups), T_regression.tsv and R01-R04."""
    d = structure.merge(core[["species_prefix", "grp", "gc_pct", "busco_complete", "total_length_bp"]], on="species_prefix", how="inner")
    if d.empty:
        return
    o.table(group_tests(d, TEST_COLS, "grp"), "T_anova")
    reg = pd.concat([pair_regressions(d, REG_LOG, log=True).assign(scale="log10"),
                     pair_regressions(d, REG_LINEAR, log=False).assign(scale="linear")], ignore_index=True)
    o.table(reg, "T_regression")
    for name, x, y, log, xl, yl in R_FIGS:
        if x not in d or y not in d:
            o.skip(name, TOPIC_REG)
            continue
        safe(o, name, TOPIC_REG, lambda x=x, y=y, log=log, xl=xl, yl=yl: _reg_fig(d, x, y, log, xl, yl, core, levels, colors))
