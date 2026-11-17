#!/usr/bin/env python3
"""Figures and tables for the Quarto summary, from the parquet files of build_summary_tables.py.

    python make_figures.py DATA_DIR OUT_DIR
      -> OUT_DIR/figures/<level>/*.png|svg, OUT_DIR/tables/<level>/*.tsv   for level in phylum, subphylum, class

Levels: 'phylum'; 'subphylum' (the subphylum where defined, else the phylum); 'class'. Groups are ordered
and coloured by their parent: subphyla by phylum, classes by phylum and ordered within a phylum by
subphylum. Groups with fewer than MIN_N genomes are pooled as 'Other (small groups)' at subphylum/class level.

Unit: one representative genome per species. Group means/medians are taken over genomes, never over
proteins, so large genomes do not dominate. Nothing here is a statistical test: the figures describe the
database. The .qmd pages state the unit for every figure.
"""
import colorsys
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from panel import PANEL  # noqa: E402
from taxonomy import CLASS_TO_SUBPHYLUM  # noqa: E402

LEVELS = ("phylum", "subphylum", "class", "order")
TOP_PHYLA = 8      # phylum level: this many phyla shown separately, the rest pooled
MIN_N = 10         # subphylum/class level: smaller groups are pooled
MAX_GROUPS = 36
OTHER = "Other (small groups)"
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#F0E442", "#8C564B"]
GREY = "#9a9a9a"
LOC_COLORS = {"cytoplasm": "#4C78A8", "nucleus": "#F58518", "extracellular": "#54A24B", "cell_membrane": "#E45756",
              "mitochondrion": "#B279A2", "plastid": "#72B7B2", "endoplasmic_reticulum": "#EECA3B",
              "lysosome_vacuole": "#9D755D", "golgi": "#BAB0AC", "peroxisome": "#FF9DA6"}

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 10, "axes.titleweight": "bold", "figure.dpi": 110,
                     "savefig.bbox": "tight", "svg.fonttype": "none"})


class Out:
    def __init__(self, root, level):
        self.fig = Path(root) / "figures" / level
        self.tab = Path(root) / "tables" / level
        self.fig.mkdir(parents=True, exist_ok=True)
        self.tab.mkdir(parents=True, exist_ok=True)
        self.made = []

    def save(self, fig, name):
        for ext in ("png", "svg"):
            fig.savefig(self.fig / f"{name}.{ext}", dpi=200 if ext == "png" else None)
        plt.close(fig)
        self.made.append(name)

    def table(self, df, name):
        df.to_csv(self.tab / f"{name}.tsv", sep="\t", index=False)
        self.made.append(name)


def width(levels, base=6.8, per=0.3):
    """Figure width that grows with the number of groups."""
    return max(base, per * len(levels) + 2.6)


# ---------------------------------------------------------------- grouping
def shade(hex_color, i, n):
    """Shade i of n of a colour: same hue and saturation, lightness spread from dark (i=0) to light (i=n-1)."""
    h, l, s = colorsys.rgb_to_hls(*(int(hex_color.lstrip("#")[k:k + 2], 16) / 255 for k in (0, 2, 4)))
    lo, hi = max(0.25, l - 0.18), min(0.78, l + 0.22)
    li = lo if n <= 1 else lo + (hi - lo) * i / (n - 1)
    if n <= 1:
        li = l
    return "#%02x%02x%02x" % tuple(round(v * 255) for v in colorsys.hls_to_rgb(h, li, s))


def assign_groups(core, level):
    """Add 'grp' (the plotted group) and 'parent' (phylum used for colour/order); return core, levels, colors."""
    d = core.copy()
    ph = d["PHYLUM"].fillna("Unknown")
    sub = d["SUBPHYLUM"] if "SUBPHYLUM" in d.columns else pd.Series(np.nan, index=d.index)
    # SUBPHYLUM is NULL in the database for all representatives: derive it from the class (taxonomy.py);
    # classes that are not in the mapping keep the phylum label
    sub = sub.fillna(d["CLASS"].map(CLASS_TO_SUBPHYLUM))
    d["PHYLUM"] = ph
    d["parent"] = ph
    d["subphylum_or_phylum"] = sub.fillna(ph)
    if level == "phylum":
        d["grp"] = ph
    elif level == "subphylum":
        d["grp"] = d["subphylum_or_phylum"]
    elif level == "order":
        d["grp"] = d["ORDER"].fillna("Unknown")
    else:
        d["grp"] = d["CLASS"].fillna("Unknown")
    # rank phyla by genome count: colours and ordering
    ph_rank = ph.value_counts()
    ph_order = list(ph_rank.index)
    n = d["grp"].value_counts()
    if level == "phylum":
        keep = list(n.index[:TOP_PHYLA])
    else:
        keep = [g for g in n.index if n[g] >= MIN_N][:MAX_GROUPS]
    pooled = (d.loc[~d["grp"].isin(keep)].groupby("grp")
              .agg(genomes=("parent", "size"), phylum=("parent", lambda x: x.value_counts().index[0]),
                   subphylum=("subphylum_or_phylum", lambda x: x.value_counts().index[0]))
              .reset_index().rename(columns={"grp": "pooled_group"}).sort_values("genomes", ascending=False))
    d.loc[~d["grp"].isin(keep), "grp"] = OTHER
    # parent of each kept group: the phylum with most genomes in it
    par = d[d.grp != OTHER].groupby("grp")["parent"].agg(lambda x: x.value_counts().index[0])
    subs = d[d.grp != OTHER].groupby("grp")["subphylum_or_phylum"].agg(lambda x: x.value_counts().index[0])
    top_ph = ph_order[:len(PALETTE)]
    pcolor = {p: (PALETTE[i] if p in top_ph else GREY) for i, p in enumerate(ph_order[:len(PALETTE)])}

    def sort_key(g):
        return (ph_order.index(par[g]), str(subs[g]), -int(n[g]))
    levels = sorted(par.index, key=sort_key)
    if (d.grp == OTHER).any():
        levels.append(OTHER)
    colors = {}
    for p_, members in par.groupby(par):
        ranked = sorted(members.index, key=lambda g: -int(n[g]))   # largest group in the phylum gets shade 0
        for i, g in enumerate(ranked):
            base = pcolor.get(p_, GREY)
            colors[g] = base if (level == "phylum" or base == GREY) else shade(base, i, len(ranked))
    colors[OTHER] = GREY
    d.attrs["parents"] = {g: (par[g] if g in par.index else "Other") for g in levels}
    d.attrs["parent_colors"] = {p: pcolor.get(p, GREY) for p in dict.fromkeys(d.attrs["parents"].values())}
    d.attrs["subphylum_of"] = {g: (subs[g] if g in subs.index else "") for g in levels}
    d.attrs["other_members"] = pooled
    return d, levels, colors


def parent_legend(ax, core, loc="upper right", fontsize=7):
    pc = core.attrs.get("parent_colors", {})
    if len(pc) > 1:
        handles = [Patch(color=c, label=p) for p, c in pc.items() if p != "Other"]
        if "Other" in pc:
            handles.append(Patch(color=GREY, label="Other = pooled small groups"))
        ax.legend(handles=handles, frameon=False, fontsize=fontsize, title="Phylum", title_fontsize=fontsize, loc=loc)


def color_ticks(ax, levels, colors, axis="x"):
    labs = ax.get_xticklabels() if axis == "x" else ax.get_yticklabels()
    for t, g in zip(labs, levels):
        t.set_color(colors.get(g, "black"))


def tick_font(levels):
    return 8 if len(levels) <= 12 else (7 if len(levels) <= 24 else 6)


def box_by_group(ax, core, col, levels, colors, log=False, ylabel=None):
    data = [core.loc[core.grp == g, col].dropna().values for g in levels]
    bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False, medianprops={"color": "black", "lw": 1.2})
    rng = np.random.default_rng(1)
    for i, (d, g) in enumerate(zip(data, levels), start=1):
        bp["boxes"][i - 1].set(facecolor=colors[g], alpha=0.35, edgecolor=colors[g])
        if len(d):
            ax.scatter(i + rng.uniform(-0.18, 0.18, len(d)), d, s=3, color=colors[g], alpha=0.35, linewidths=0)
    ax.set_xticks(range(1, len(levels) + 1))
    ax.set_xticklabels([f"{g}\n(n={len(d)})" for g, d in zip(levels, data)], rotation=60 if len(levels) > 12 else 35,
                       ha="right", fontsize=tick_font(levels))
    color_ticks(ax, levels, colors)
    if log:
        ax.set_yscale("log")
    if ylabel:
        ax.set_ylabel(ylabel)


def heat(ax, mat, cmap="viridis", vmin=None, vmax=None, label=None, colors=None):
    im = ax.imshow(mat.values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=8)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=60, ha="right", fontsize=tick_font(mat.columns))
    if colors:
        color_ticks(ax, list(mat.columns), colors)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    if label:
        cb.set_label(label, fontsize=8)
    return im


def stacked_by_group(ax, frac, colors_map, levels, ylabel, group_colors=None):
    bottom = np.zeros(len(levels))
    for lab in frac.columns:
        v = frac.loc[levels, lab].values
        ax.bar(range(len(levels)), v, bottom=bottom, color=colors_map.get(lab, GREY), label=lab, width=0.8)
        bottom += v
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels(levels, rotation=60 if len(levels) > 12 else 35, ha="right", fontsize=tick_font(levels))
    if group_colors:
        color_ticks(ax, levels, group_colors)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 1)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False, fontsize=8)


def per_genome_fraction(df, label_col, count_col="n_proteins"):
    piv = df.pivot_table(index="species_prefix", columns=label_col, values=count_col, aggfunc="sum", fill_value=0)
    return piv.div(piv.sum(axis=1), axis=0)


def group_mean(frac, core, levels):
    j = frac.join(core.set_index("species_prefix")["grp"], how="inner")
    return j.groupby("grp").mean(numeric_only=True).reindex(levels).fillna(0)


def group_key_table(core, levels):
    d = core.groupby("grp").agg(genomes=("species_prefix", "size")).reindex(levels).reset_index()
    d["phylum"] = [core.attrs["parents"].get(g, "") for g in levels]
    d["subphylum"] = [core.attrs["subphylum_of"].get(g, "") for g in levels]
    return d.rename(columns={"grp": "group"})[["group", "subphylum", "phylum", "genomes"]]


# ---------------------------------------------------------------- page 1: gene content
def page_gene_content(core, levels, colors, o):
    n = core.groupby("grp").size().reindex(levels)
    fig, ax = plt.subplots(figsize=(width(levels, 6.5), 3.4))
    ax.bar(range(len(levels)), n.values, color=[colors[g] for g in levels])
    for i, v in enumerate(n.values):
        ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=tick_font(levels))
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels(levels, rotation=60 if len(levels) > 12 else 35, ha="right", fontsize=tick_font(levels))
    color_ticks(ax, levels, colors)
    ax.set_ylabel("Representative genomes")
    ax.set_title("One genome per species")
    parent_legend(ax, core)
    o.save(fig, "F01_genomes_by_group")
    o.table(group_key_table(core, levels), "T00_group_key")
    o.table(core.attrs["other_members"], "T00b_other_members")

    for col, name, ylab, log in [("n_genes", "F02_gene_count", "Protein-coding genes", True),
                                 ("total_length_bp", "F02b_genome_size", "Genome size (bp)", True),
                                 ("mean_len", "F04_protein_length", "Mean protein length (aa)", False),
                                 ("busco_complete", "F05_busco", "BUSCO complete (%)", False),
                                 ("gc_pct", "F05b_gc", "Genome GC (%)", False)]:
        fig, ax = plt.subplots(figsize=(width(levels), 3.8))
        box_by_group(ax, core, col, levels, colors, log=log, ylabel=ylab)
        ax.set_title(ylab)
        parent_legend(ax, core, loc="best")
        o.save(fig, name)

    fig, ax = plt.subplots(figsize=(5.8, 4.6))
    for g in levels:
        d = core[core.grp == g]
        ax.scatter(d.total_length_bp / 1e6, d.n_genes, s=6, alpha=0.5, color=colors[g], linewidths=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Genome size (Mb)")
    ax.set_ylabel("Protein-coding genes")
    ax.set_title("Genome size vs gene number")
    if core.attrs["parent_colors"] and len(levels) > 10:
        parent_legend(ax, core, loc="upper left")
    else:
        ax.legend(handles=[Patch(color=colors[g], label=g) for g in levels], frameon=False, fontsize=7)
    o.save(fig, "F03_genome_size_vs_genes")

    t = core.groupby("grp").agg(genomes=("species_prefix", "size"), genome_Mb=("total_length_bp", lambda x: x.median() / 1e6),
                                genes=("n_genes", "median"), protein_len=("mean_len", "median"),
                                gc_pct=("gc_pct", "median"), busco=("busco_complete", "median")).reindex(levels).round(2)
    o.table(t.reset_index().rename(columns={"grp": "group"}), "T01_overview")


# ---------------------------------------------------------------- page 2: domains
def page_domains(core, levels, colors, pfam, names, arch, o):
    n_genes = core.set_index("species_prefix")["n_genes"]
    a = arch.pivot_table(index="species_prefix", columns="n_domains", values="n_proteins", aggfunc="sum", fill_value=0)
    a = a.reindex(core.species_prefix).fillna(0)
    a.insert(0, 0, (n_genes.reindex(a.index) - a.sum(axis=1)).clip(lower=0))
    frac = a.div(a.sum(axis=1), axis=0)
    frac.columns = ["0" if c == 0 else ("5+" if c >= 5 else str(int(c))) for c in frac.columns]
    gm = group_mean(frac, core, levels)
    fig, ax = plt.subplots(figsize=(width(levels), 3.8))
    cmap = plt.get_cmap("Blues")
    stacked_by_group(ax, gm, {c: cmap(0.2 + 0.16 * i) for i, c in enumerate(gm.columns)}, levels,
                     "Share of proteins (genome mean)", colors)
    ax.set_title("Pfam domains per protein (0 = no Pfam hit)")
    o.save(fig, "F06_domains_per_protein")

    wide = pfam.pivot_table(index="species_prefix", columns="pfam_acc", values="n_proteins", aggfunc="sum", fill_value=0)
    wide = wide.reindex(core.species_prefix).fillna(0)
    grp = core.set_index("species_prefix")["grp"].reindex(wide.index)
    present = (wide > 0).groupby(grp).mean().reindex(levels)
    label = names.set_index("pfam_acc")["pfam_id"]
    top = present.mean().sort_values(ascending=False).index[:30]
    # copies per 10,000 genes (mean over genomes): removes the genome-size effect (reduced genomes look
    # low on every domain otherwise) and avoids integer-median ties (median 1 vs 2 copies)
    per10k = wide.div(n_genes.reindex(wide.index), axis=0) * 1e4
    cnt = per10k.groupby(grp).mean().reindex(levels)

    # F07: the widespread domains are present almost everywhere (a heatmap of presence is flat), so show
    # presence only for domains that are neither near-universal nor rare overall (5-95% of genomes),
    # ranked by how much their prevalence differs between groups
    overall = (wide > 0).mean()
    mid = overall[(overall >= 0.05) & (overall <= 0.95)].index
    vp = present[mid].var().sort_values(ascending=False).index[:40]
    m = present[vp].T
    m.index = [f"{a} {label.get(a, '')}" for a in m.index]
    fig, ax = plt.subplots(figsize=(width(levels), 8.6))
    heat(ax, m, cmap="magma", vmin=0, vmax=1, label="Fraction of genomes with the domain", colors=colors)
    ax.set_title("Domains whose presence varies between groups")
    o.save(fig, "F07_pfam_prevalence")
    o.table(pd.DataFrame({"pfam_acc": vp, "pfam_id": [label.get(a, "") for a in vp],
                          "overall_prevalence": [float(overall[a]) for a in vp],
                          "between_group_variance": [float(present[a].var()) for a in vp]}), "T02b_variable_presence")

    # F07b: the 30 most widespread domains differ in copy number, not presence: log2 fold of each group's
    # mean copies per 10,000 genes over the all-genome mean
    base = per10k.mean()
    fold = np.log2((cnt[top] + 0.05) / (base[top] + 0.05))
    m3 = fold.T
    m3.index = [f"{a} {label.get(a, '')}" for a in m3.index]
    fig, ax = plt.subplots(figsize=(width(levels), 7.2))
    heat(ax, m3, cmap="RdBu_r", vmin=-1.5, vmax=1.5, label="log2(group mean / all-genome mean), copies per 10,000 genes", colors=colors)
    ax.set_title("Most widespread Pfam domains: copy-number deviation")
    o.save(fig, "F07b_widespread_fold")
    keep = [c for c in wide.columns if (present[c].max() >= 0.3) and per10k[c].mean() >= 0.5]
    z = np.log1p(cnt[keep])
    z = (z - z.mean()) / z.std(ddof=0).replace(0, np.nan)
    var = z.var().sort_values(ascending=False).index[:40]
    m2 = z[var].T
    m2.index = [f"{a} {label.get(a, '')}" for a in m2.index]
    fig, ax = plt.subplots(figsize=(width(levels), 8.6))
    heat(ax, m2, cmap="RdBu_r", vmin=-2.5, vmax=2.5, label="z-score of log mean copies per 10,000 genes", colors=colors)
    ax.set_title("Domains that differ most between groups")
    o.save(fig, "F08_pfam_variance")

    nfam = (wide > 0).sum(axis=1).rename("n_pfam_families_in_top_set")
    d = core.set_index("species_prefix").join(nfam)
    fig, ax = plt.subplots(figsize=(width(levels), 3.8))
    box_by_group(ax, d.reset_index(), "n_pfam_families_in_top_set", levels, colors,
                 ylabel="Families present (of the retained set)")
    ax.set_title("Pfam repertoire breadth")
    o.save(fig, "F09_pfam_breadth")

    t = pd.DataFrame({"pfam_acc": top, "pfam_id": [label.get(a, "") for a in top],
                      "genomes_with_domain": [int((wide[a] > 0).sum()) for a in top],
                      "mean_copies": [float(wide[a].mean()) for a in top]})
    o.table(t, "T02_top_pfam")
    return wide, label


# ---------------------------------------------------------------- page 3: localization
def page_localization(core, levels, colors, deeploc, wolf, o):
    fr = per_genome_fraction(deeploc, "label")
    gm = group_mean(fr, core, levels)
    gm = gm[[c for c in LOC_COLORS if c in gm.columns]]
    fig, ax = plt.subplots(figsize=(width(levels, 6.8), 4.0))
    stacked_by_group(ax, gm, LOC_COLORS, levels, "Share of proteins (genome mean)", colors)
    ax.set_title("DeepLoc primary localization (argmax of class probabilities)")
    o.save(fig, "F10_deeploc")
    o.table(gm.round(4).reset_index().rename(columns={"grp": "group"}), "T03_deeploc")

    wf = per_genome_fraction(wolf, "localization")
    wg = group_mean(wf, core, levels)
    wg = wg[wg.mean().sort_values(ascending=False).index[:10]]
    cmap = plt.get_cmap("tab10")
    fig, ax = plt.subplots(figsize=(width(levels, 6.8), 4.0))
    stacked_by_group(ax, wg, {c: cmap(i) for i, c in enumerate(wg.columns)}, levels, "Share of proteins (genome mean)", colors)
    ax.set_title("WoLF PSORT localization (top 10 classes)")
    o.save(fig, "F12_wolfpsort")

    core = core.copy()
    core["secreted_frac"] = core["n_secreted_strict"] / core["n_genes"]
    core["signalp_frac"] = core["n_signalp_09"] / core["n_genes"]
    for col, name, ylab in [("secreted_frac", "F13_secretome_fraction", "Strict secretome / genes"),
                            ("signalp_frac", "F13b_signal_peptide_fraction", "SignalP (p>=0.9) / genes")]:
        fig, ax = plt.subplots(figsize=(width(levels), 3.8))
        box_by_group(ax, core, col, levels, colors, ylabel=ylab)
        ax.set_title(ylab)
        o.save(fig, name)


# ---------------------------------------------------------------- page 4: patterns
def genome_features(core, deeploc, idp):
    d = core.set_index("species_prefix").copy()
    d["genome_Mb"] = d["total_length_bp"] / 1e6
    d["secreted_frac"] = d["n_secreted_strict"] / d["n_genes"]
    d["cazy_per_gene"] = d["n_cazy"] / d["n_genes"]
    d["merops_per_gene"] = d["n_merops"] / d["n_genes"]
    d["tf_per_gene"] = d["n_tf_proteins"] / d["n_genes"]
    fr = per_genome_fraction(deeploc, "label")
    for c in ("nucleus", "mitochondrion", "extracellular"):
        d[f"loc_{c}"] = fr[c] if c in fr.columns else np.nan
    d = d.join(idp.set_index("species_prefix")[["mean_idp_fraction"]])
    return d


FEATURES = ["genome_Mb", "n_genes", "mean_len", "gc_pct", "secreted_frac", "cazy_per_gene", "merops_per_gene",
            "tf_per_gene", "loc_nucleus", "loc_mitochondrion", "loc_extracellular", "mean_idp_fraction"]


def page_patterns(core, levels, colors, deeploc, idp, o):
    d = genome_features(core, deeploc, idp)
    d.attrs = core.attrs
    feats = [f for f in FEATURES if f in d.columns]
    med = d.groupby("grp")[feats].median().reindex(levels)
    z = (med - med.mean()) / med.std(ddof=0).replace(0, np.nan)
    fig, ax = plt.subplots(figsize=(width(levels, 7), 4.4))
    heat(ax, z.T, cmap="RdBu_r", vmin=-2.5, vmax=2.5, label="z-score of group median", colors=colors)
    ax.set_title("Group profiles across genome-level features")
    o.save(fig, "F15_feature_heatmap")
    o.table(med.round(4).reset_index().rename(columns={"grp": "group"}), "T04_feature_medians")

    x = d[feats].copy()
    x = x.fillna(x.median())
    x = np.log1p(x.clip(lower=0))
    x = (x - x.mean()) / x.std(ddof=0).replace(0, np.nan)
    x = x.fillna(0)
    u, s, vt = np.linalg.svd(x.values, full_matrices=False)
    pc = u[:, :2] * s[:2]
    ev = (s ** 2) / (s ** 2).sum()
    fig, ax = plt.subplots(figsize=(5.8, 4.8))
    for g in levels:
        mk = (d["grp"] == g).values
        ax.scatter(pc[mk, 0], pc[mk, 1], s=7, alpha=0.55, color=colors[g], linewidths=0)
    ax.set_xlabel(f"PC1 ({ev[0]:.0%})")
    ax.set_ylabel(f"PC2 ({ev[1]:.0%})")
    ax.set_title("Genomes in feature space (log, z-scored)")
    if len(levels) > 10:
        parent_legend(ax, d)
    else:
        ax.legend(handles=[Patch(color=colors[g], label=g) for g in levels], frameon=False, fontsize=7)
    o.save(fig, "F17_feature_pca")
    load = pd.DataFrame(vt[:2].T, index=feats, columns=["PC1", "PC2"]).round(3)
    o.table(load.reset_index().rename(columns={"index": "feature"}), "T06_pca_loadings")


# ---------------------------------------------------------------- page 5: functional
def page_functional(core, levels, colors, wide, label, cazy, tf, merops, idp, o):
    ng = core.set_index("species_prefix")["n_genes"]
    grp = core.set_index("species_prefix")["grp"]

    cz = cazy.pivot_table(index="species_prefix", columns="cazy_class", values="n_proteins", aggfunc="sum", fill_value=0)
    cz = cz.reindex(core.species_prefix).fillna(0)
    classes = [c for c in ["GH", "GT", "PL", "CE", "AA", "CBM"] if c in cz.columns]
    if classes:
        per_gene = cz[classes].div(ng.reindex(cz.index), axis=0)
        mean_pg = per_gene.groupby(grp.reindex(per_gene.index)).mean().reindex(levels)
        # horizontal stacked bars: groups on the y axis (readable labels at any number of groups), one colour per
        # CAZy class with the full class names in the legend
        names = {"GH": "GH  glycoside hydrolases", "GT": "GT  glycosyltransferases", "PL": "PL  polysaccharide lyases",
                 "CE": "CE  carbohydrate esterases", "AA": "AA  auxiliary activities (redox)",
                 "CBM": "CBM  carbohydrate-binding modules"}
        ccol = {"GH": "#0072B2", "GT": "#E69F00", "PL": "#009E73", "CE": "#CC79A7", "AA": "#D55E00", "CBM": "#56B4E9"}
        fig, ax = plt.subplots(figsize=(7.4, max(3.0, 0.3 * len(levels) + 1.6)))
        left = np.zeros(len(levels))
        ypos = np.arange(len(levels))
        for c in classes:
            v = mean_pg.loc[levels, c].values
            ax.barh(ypos, v, left=left, color=ccol[c], label=names[c], height=0.78)
            left += v
        ax.set_yticks(ypos)
        ax.set_yticklabels(levels, fontsize=tick_font(levels) + 1)
        ax.invert_yaxis()
        color_ticks(ax, levels, colors, axis="y")
        ax.set_xlabel("CAZy proteins per gene (genome mean)")
        ax.legend(title="CAZy class", bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False, fontsize=8)
        ax.set_title("CAZy classes per gene")
        o.save(fig, "F18_cazy_classes")

        fam = cazy.pivot_table(index="species_prefix", columns="cazy_family", values="n_proteins", aggfunc="sum", fill_value=0)
        fam = fam.reindex(core.species_prefix).fillna(0)
        mean = fam.groupby(grp.reindex(fam.index)).mean().reindex(levels)
        top = fam.sum().sort_values(ascending=False).index[:30]
        fig, ax = plt.subplots(figsize=(width(levels), 7))
        heat(ax, np.log2(mean[top].T + 1), cmap="viridis", label="log2(mean proteins per genome + 1)", colors=colors)
        ax.set_title("Top 30 CAZy families")
        o.save(fig, "F19_cazy_families")
        o.table(mean[top].T.reset_index().rename(columns={"cazy_family": "family"}), "T07_cazy_family_means")

    mer = merops.pivot_table(index="species_prefix", columns="merops_class", values="n_proteins", aggfunc="sum", fill_value=0)
    mer = mer.reindex(core.species_prefix).fillna(0)
    if mer.shape[1]:
        mf = mer.div(ng.reindex(mer.index), axis=0)
        gm = mf.groupby(grp.reindex(mf.index)).mean().reindex(levels)
        fig, ax = plt.subplots(figsize=(width(levels), 3.8))
        bottom = np.zeros(len(levels))
        cmap = plt.get_cmap("Set2")
        for i, c in enumerate(gm.columns):
            ax.bar(range(len(levels)), gm[c].values, bottom=bottom, label=c, color=cmap(i % 8), width=0.8)
            bottom += gm[c].values
        ax.set_xticks(range(len(levels)))
        ax.set_xticklabels(levels, rotation=60 if len(levels) > 12 else 35, ha="right", fontsize=tick_font(levels))
        color_ticks(ax, levels, colors)
        ax.set_ylabel("MEROPS proteins per gene (genome mean)")
        ax.legend(title="class", bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
        ax.set_title("Peptidase classes (MEROPS first letter)")
        o.save(fig, "F20_merops_classes")

    tfw = tf.pivot_table(index="species_prefix", columns="tf_family", values="n_proteins", aggfunc="sum", fill_value=0)
    tfw = tfw.reindex(core.species_prefix).fillna(0)
    if tfw.shape[1]:
        mean = tfw.groupby(grp.reindex(tfw.index)).mean().reindex(levels)
        top = tfw.sum().sort_values(ascending=False).index[:25]
        fig, ax = plt.subplots(figsize=(width(levels), 6.2))
        heat(ax, np.log2(mean[top].T + 1), cmap="viridis", label="log2(mean proteins per genome + 1)", colors=colors)
        ax.set_title("Transcription-factor families")
        o.save(fig, "F21_tf_families")
        o.table(mean[top].T.reset_index().rename(columns={"tf_family": "family"}), "T08_tf_family_means")

    pan = []
    for k, (cat, accs) in enumerate(PANEL.items(), start=1):
        have = [a for a in accs if a in wide.columns]
        if not have:
            continue
        mean = wide[have].groupby(grp.reindex(wide.index)).mean().reindex(levels)  # mean proteins per genome
        prev = (wide[have] > 0).groupby(grp.reindex(wide.index)).mean().reindex(levels)
        # left: absolute abundance; right: deviation from the all-genome mean per 10,000 genes, which shows
        # variation even when a domain is present in nearly every genome (prevalence is in the table)
        p10 = wide[have].div(ng.reindex(wide.index), axis=0) * 1e4
        gp10 = p10.groupby(grp.reindex(wide.index)).mean().reindex(levels)
        fold = np.log2((gp10 + 0.05) / (p10.mean() + 0.05))
        rows = [f"{a} {label.get(a, accs[a])}" for a in have]
        m = np.log2(mean.T + 1)
        m.index = rows
        f2 = fold.T
        f2.index = rows
        w = width(levels, 5.2, 0.28)
        rowh = max(1.6, 0.34 * len(have) + 0.9)
        stack = len(levels) > 16  # many groups: stack the two heatmaps instead of putting them side by side
        if stack:
            fig, (ax, ax2) = plt.subplots(2, 1, figsize=(w + 2.6, 2 * rowh + 1.2), gridspec_kw={"hspace": 0.95})
        else:
            fig, (ax, ax2) = plt.subplots(1, 2, figsize=(2 * w + 2.6, rowh + 0.8), gridspec_kw={"wspace": 0.08})
        heat(ax, m, cmap="viridis", label="log2(mean proteins per genome + 1)", colors=colors)
        ax.set_title("Abundance")
        heat(ax2, f2, cmap="RdBu_r", vmin=-2, vmax=2, label="log2(group / all-genome mean per 10,000 genes)", colors=colors)
        if not stack:
            ax2.set_yticks([])
        ax2.set_title("Deviation from the all-genome mean")
        fig.suptitle(cat, y=1.02, fontweight="bold")
        slug = "".join(ch if ch.isalnum() else "_" for ch in cat.lower()).strip("_")
        o.save(fig, f"F22_{k:02d}_panel_{slug}")
        for a in have:
            for g in levels:
                pan.append({"category": cat, "pfam_acc": a, "pfam_name": label.get(a, ""), "description": accs[a], "group": g,
                            "mean_proteins": float(mean.loc[g, a]), "fraction_genomes": float(prev.loc[g, a])})
    if pan:
        o.table(pd.DataFrame(pan), "T09_functional_panel")

    d = core.set_index("species_prefix").join(idp.set_index("species_prefix")[["mean_idp_fraction", "frac_idp_gt_0_3"]]).reset_index()
    for col, name, ylab in [("mean_idp_fraction", "F23_idp_fraction", "Mean IDR fraction per protein"),
                            ("frac_idp_gt_0_3", "F23b_idp_gt03", "Share of proteins with IDR fraction > 0.3")]:
        fig, ax = plt.subplots(figsize=(width(levels), 3.8))
        box_by_group(ax, d, col, levels, colors, ylabel=ylab)
        ax.set_title(ylab)
        o.save(fig, name)


def run_level(level, core0, tables, out):
    core, levels, colors = assign_groups(core0, level)
    o = Out(out, level)
    rd = tables
    page_gene_content(core, levels, colors, o)
    wide, label = page_domains(core, levels, colors, rd("pfam_genome"), rd("pfam_names"), rd("arch_genome"), o)
    deeploc = rd("deeploc_genome")
    page_localization(core, levels, colors, deeploc, rd("wolfpsort_genome"), o)
    idp = rd("idp_genome")
    page_patterns(core, levels, colors, deeploc, idp, o)
    page_functional(core, levels, colors, wide, label, rd("cazy_genome"), rd("tf_genome"), rd("merops_genome"), idp, o)
    print(f"{level}: {len(levels)} groups, {len(o.made)} outputs", flush=True)
    return o.made


def main(data, out, levels=LEVELS):
    data = Path(data)
    cache = {}

    def rd(n):
        if n not in cache:
            cache[n] = pd.read_parquet(data / f"{n}.parquet")
        return cache[n]
    core0 = rd("genome_core")
    core0 = core0[core0["n_genes"].notna()]  # representatives without a proteome are not plotted
    made = {}
    for level in levels:
        made[level] = run_level(level, core0, rd, out)
    return made


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], tuple(sys.argv[3].split(",")) if len(sys.argv) > 3 else LEVELS)
