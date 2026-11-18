#!/usr/bin/env python3
"""Figures and tables for the Quarto summary, from the parquet files of build_summary_tables.py.

    python make_figures.py DATA_DIR OUT_DIR
      -> OUT_DIR/figures/<level>/*.png|svg, OUT_DIR/tables/<level>/*.tsv   for level in phylum, subphylum, class
      -> OUT_DIR/figures/global/, OUT_DIR/tables/global/   level-free feature-space set (biplot, QC overlays,
         one ordination per major phylum coloured by class, the Ascomycota split); see page_feature_space

Levels: 'phylum'; 'subphylum' (the subphylum where defined, else the phylum); 'class'. Groups are ordered
and coloured by their parent: subphyla by phylum, classes by phylum and ordered within a phylum by
subphylum. Groups with fewer than MIN_N genomes are pooled as 'Other (small groups)' at subphylum/class level.

Unit: one representative genome per species. Group means/medians are taken over genomes, never over
proteins, so large genomes do not dominate. Nothing here is a statistical test: the figures describe the
database. The .qmd pages state the unit for every figure.
"""
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

LEVELS = ("phylum", "subphylum", "class")
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
    colors = {g: pcolor.get(par[g], GREY) for g in par.index}
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


FEATURE_LABELS = {"genome_Mb": "genome size", "n_genes": "genes", "mean_len": "protein length", "gc_pct": "GC",
                  "secreted_frac": "secreted / gene", "cazy_per_gene": "CAZy / gene", "merops_per_gene": "MEROPS / gene",
                  "tf_per_gene": "TF / gene", "loc_nucleus": "nuclear share", "loc_mitochondrion": "mitochondrial share",
                  "loc_extracellular": "extracellular share", "mean_idp_fraction": "IDR fraction",
                  "busco_complete": "BUSCO complete (%)", "N50_bp": "N50 (bp)", "contig_count": "contigs"}


def standardize(df):
    """PCA input: missing values -> column median, log1p (features are counts, sizes and fractions, all >= 0),
    then z-score; a constant column becomes 0."""
    x = df.copy()
    x = x.fillna(x.median())
    x = np.log1p(x.clip(lower=0))
    x = (x - x.mean()) / x.std(ddof=0).replace(0, np.nan)
    return x.fillna(0)


def pca(x, k=2):
    """PCA by SVD of a centred matrix (a DataFrame). Returns scores (n x k), explained-variance ratios (k) and
    loadings (features x k). The SVD sign is arbitrary, so each axis is flipped to make its largest-|loading|
    feature positive: the orientation is then fixed for a given input and stable across reruns."""
    u, s, vt = np.linalg.svd(np.asarray(x, dtype=float), full_matrices=False)
    k = min(k, len(s))
    load, scores = vt[:k].T.copy(), u[:, :k] * s[:k]
    for j in range(k):
        if load[np.argmax(np.abs(load[:, j])), j] < 0:
            load[:, j] *= -1
            scores[:, j] *= -1
    ev = (s ** 2) / (s ** 2).sum()
    return scores, ev[:k], pd.DataFrame(load, index=list(x.columns), columns=[f"PC{j + 1}" for j in range(k)])


def halo(lw=2.5):
    import matplotlib.patheffects as pe
    return [pe.withStroke(linewidth=lw, foreground="white")]


def label_centroids(ax, pc, grp, levels, colors, fontsize=7):
    """Mark the median of each group with a numbered, black-edged dot and list number -> group in a legend on the
    right (group names coloured by phylum). Numbers overlap far less than names when there are 30+ groups."""
    from matplotlib.lines import Line2D
    handles, labels = [], []
    k = 0
    for g in levels:
        mk = (grp == g)
        if g == OTHER or mk.sum() == 0:
            continue
        k += 1
        cx, cy = np.median(pc[mk, 0]), np.median(pc[mk, 1])
        ax.scatter([cx], [cy], s=70, color=colors[g], edgecolor="black", linewidths=0.8, zorder=5)
        ax.text(cx, cy, str(k), fontsize=6, ha="center", va="center", zorder=6, color="white", fontweight="bold")
        handles.append(Line2D([], [], ls="none", marker="o", markersize=7, color=colors[g], markeredgecolor="black"))
        labels.append(f"{k:>2}  {g} ({int(mk.sum()):,})")
    leg = ax.legend(handles, labels, title="Group median (n)", fontsize=fontsize, title_fontsize=8, frameon=False,
                    bbox_to_anchor=(1.01, 1), loc="upper left", labelspacing=0.3, ncol=1 if k <= 24 else 2)
    for t, g in zip(leg.get_texts(), [g for g in levels if g != OTHER and (grp == g).sum()]):
        t.set_color(colors[g])
    return leg


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

    pc, ev, load = pca(standardize(d[feats]), 2)
    grp = d["grp"].values
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    for g in levels:
        mk = grp == g
        ax.scatter(pc[mk, 0], pc[mk, 1], s=7, alpha=0.45, color=colors[g], linewidths=0)
    ax.set_xlabel(f"PC1 ({ev[0]:.0%})")
    ax.set_ylabel(f"PC2 ({ev[1]:.0%})")
    ax.set_title("Genomes in feature space (log, z-scored)")
    if len(levels) > 10:
        # many groups share a phylum colour: number each group at its median and list every group on the right
        label_centroids(ax, pc, grp, levels, colors)
    else:
        ax.legend(handles=[Patch(color=colors[g], label=g) for g in levels], frameon=False, fontsize=8,
                  bbox_to_anchor=(1.01, 1), loc="upper left")
    o.save(fig, "F17_feature_pca")
    load = load.round(3)
    o.table(load.reset_index().rename(columns={"index": "feature"}), "T06_pca_loadings")


# ---------------------------------------------------------------- page 4b: feature space (one set, no level tabs)
# Written once to figures/global/ and tables/global/: the biplot, QC overlays, one ordination per major phylum
# coloured by class, and the Ascomycota split. Unit: one representative genome per species, as everywhere.
PANEL_PHYLA = ["Ascomycota", "Basidiomycota", "Mucoromycota"]
OTHER_PHYLA = "Other phyla"
PANEL_MIN_N = 5           # classes with fewer genomes are pooled within a phylum panel
PANEL_MAX_NAMED = 19      # at most this many named classes per panel (plus the pooled entry)
# 20 colours of maximum contrast (Kelly 1965), most distinct first; grey is kept for pooled entries
KELLY = ["#0067A5", "#BE0032", "#008856", "#F38400", "#875692", "#A1CAF1", "#E68FAC", "#882D17", "#F3C300",
         "#2B3D26", "#604E97", "#F99379", "#8DB600", "#B3446C", "#C2B280", "#654522", "#DCD300", "#E25822",
         "#F6A600", "#848482"]
MARKERS = ["o", "^", "s", "D", "v", "P", "X", "<", ">", "h"]
ASCO_SPLIT = -1.0         # global PC1 cut between the two Ascomycota modes (yeasts below), see page_feature_space
SUBPH_COLORS = {"Pezizomycotina": "#0072B2", "Saccharomycotina": "#D55E00", "Taphrinomycotina": "#009E73",
                "unassigned": GREY}
LOG_FEATS = {"genome_Mb", "n_genes", "N50_bp", "contig_count"}


def best_threshold(v, y):
    """Best single cut on v for boolean labels y (NaN in v dropped), chosen by balanced accuracy (mean of the
    true-positive and true-negative rates), so a rare class cannot be ignored: with 132 yeasts against 3,479
    filamentous genomes, calling everything filamentous already has 96% plain accuracy but 50% balanced.
    Returns (balanced_accuracy, accuracy, threshold, direction): direction '<' means 'v < threshold' predicts
    True, '>=' means 'v >= threshold' predicts True. Every cut between distinct sorted values (midpoints) is
    tried, so the result is exact, not a grid search."""
    v = np.asarray(v, dtype=float)
    y = np.asarray(y, dtype=bool)
    ok = ~np.isnan(v)
    v, y = v[ok], y[ok]
    if len(v) < 2 or y.all() or not y.any():
        return np.nan, np.nan, np.nan, ""
    o = np.argsort(v, kind="mergesort")
    v, y = v[o], y[o]
    n, npos = len(v), int(y.sum())
    nneg = n - npos
    pos_below = np.cumsum(y)[:-1]                    # cut after position i (i = 0..n-2): below = v[:i+1]
    neg_below = np.arange(1, n) - pos_below
    # rule '<' (below predicts True): TPR = pos_below/npos, TNR = (nneg - neg_below)/nneg; rule '>=' is its mirror
    bal_lt = 0.5 * (pos_below / npos + (nneg - neg_below) / nneg)
    acc_lt = (pos_below + nneg - neg_below) / n
    valid = v[1:] > v[:-1]
    if not valid.any():
        return np.nan, np.nan, np.nan, ""
    score = np.where(valid, np.maximum(bal_lt, 1 - bal_lt), -1)
    i = int(np.argmax(score))
    lt = bal_lt[i] >= 1 - bal_lt[i]
    return float(score[i]), float(acc_lt[i] if lt else 1 - acc_lt[i]), float((v[i] + v[i + 1]) / 2), ("<" if lt else ">=")


def density_valley(v):
    """Lowest point of a Gaussian kernel density (Silverman bandwidth) between its two highest modes.
    Used as a data-driven check on the fixed ASCO_SPLIT cut.
    Returns (valley, low_mode, high_mode); NaNs when the density has fewer than two modes."""
    v = np.asarray(v, dtype=float)
    v = v[~np.isnan(v)]
    if len(v) < 10:
        return np.nan, np.nan, np.nan
    grid = np.linspace(v.min(), v.max(), 400)
    bw = 0.9 * min(v.std(), np.subtract(*np.percentile(v, [75, 25])) / 1.34) * len(v) ** -0.2
    dens = np.zeros_like(grid)
    for chunk in np.array_split(v, max(1, len(v) // 2000)):
        dens += np.exp(-0.5 * ((grid[:, None] - chunk[None, :]) / bw) ** 2).sum(axis=1)
    # modes lower than 10% of the highest are tail wiggles, not modes
    peaks = [i for i in range(1, len(grid) - 1) if dens[i] > dens[i - 1] and dens[i] >= dens[i + 1]
             and dens[i] >= 0.1 * dens.max()]
    if len(peaks) < 2:
        return np.nan, np.nan, np.nan
    a, b = sorted(sorted(peaks, key=lambda i: dens[i])[-2:])
    valley = a + int(np.argmin(dens[a:b + 1]))
    return float(grid[valley]), float(grid[a]), float(grid[b])


def panel_groups(sub, by):
    """Legend groups for one phylum panel: the class (or '<phylum> (class not set)'), small classes pooled.
    `by` is the column that sets the marker (subphylum within a phylum, phylum in the pooled panel).
    Returns labels (Series), ordered [(section, [labels])], colour and marker per label."""
    lab = sub["CLASS"].where(sub["CLASS"].notna(), sub["PHYLUM"] + " (class not set)")
    n = lab.value_counts()
    named = [g for g in n.index if n[g] >= PANEL_MIN_N][:PANEL_MAX_NAMED]
    pooled = f"Other classes (<{PANEL_MIN_N} genomes each)"
    lab = lab.where(lab.isin(named), pooled)
    sec_of = pd.DataFrame({"lab": lab, "sec": sub[by].fillna("unassigned")}).groupby("lab")["sec"].agg(
        lambda x: x.value_counts().index[0])
    sec_n = sub[by].fillna("unassigned").value_counts()
    sections = []
    for sec in sorted(set(sec_of[named]), key=lambda s_: -sec_n.get(s_, 0)):
        sections.append((sec, sorted([g for g in named if sec_of[g] == sec], key=lambda g: -n[g])))
    order = [g for _, gs in sections for g in gs]
    color = {g: KELLY[i % len(KELLY)] for i, g in enumerate(order)}
    marker = {g: MARKERS[i % len(MARKERS)] for i, (_, gs) in enumerate(sections) for g in gs}
    if (lab == pooled).any():
        sections.append(("pooled", [pooled]))
        color[pooled], marker[pooled] = "#bbbbbb", "."
    return lab, sections, color, marker


def sectioned_legend(fig, ax, lab, sections, color, marker, title):
    """Legend outside the axes on the right: one bold header per section (subphylum or phylum), then its
    classes with their genome count. Every class is listed, so the plot is readable without the phylum tabs."""
    from matplotlib.lines import Line2D
    n = lab.value_counts()
    handles, labels, bold = [], [], []
    for sec, gs in sections:
        if len(sections) > 1:
            handles.append(Line2D([], [], ls="none"))
            labels.append(sec)
            bold.append(len(labels) - 1)
        for g in gs:
            handles.append(Line2D([], [], ls="none", marker=marker[g], color=color[g], markersize=6))
            labels.append(f"{g} ({n[g]:,})")
    leg = ax.legend(handles, labels, title=title, title_fontsize=8.5, fontsize=8, frameon=False,
                    bbox_to_anchor=(1.02, 1), loc="upper left", handletextpad=0.4, labelspacing=0.35,
                    ncol=1 if len(labels) <= 26 else 2)
    for i in bold:
        leg.get_texts()[i].set_fontweight("bold")
    return leg


def scatter_groups(ax, xy, lab, sections, color, marker, s=10, alpha=0.6):
    n = lab.value_counts()
    order = sorted([g for _, gs in sections for g in gs], key=lambda g: -n[g])   # small groups drawn on top
    for g in order:
        mk = (lab == g).values
        ax.scatter(xy[mk, 0], xy[mk, 1], s=s, alpha=alpha, color=color[g], marker=marker[g], linewidths=0)


def page_feature_space(core0, deeploc, idp, o):
    core, _, pcolors = assign_groups(core0, "phylum")
    d = genome_features(core, deeploc, idp)
    d["subphylum"] = d["CLASS"].map(CLASS_TO_SUBPHYLUM)
    feats = [f for f in FEATURES if f in d.columns]
    pc, ev, load = pca(standardize(d[feats]), 3)
    d["PC1"], d["PC2"], d["PC3"] = pc[:, 0], pc[:, 1], pc[:, 2]
    xlab, ylab = f"global PC1 ({ev[0]:.0%})", f"global PC2 ({ev[1]:.0%})"
    phyla = [g for g in pcolors if g != OTHER]

    # ---- biplot: points by phylum, one arrow per feature (loading x scale)
    fig, ax = plt.subplots(figsize=(7.6, 6.0))
    for g in sorted(pcolors, key=lambda g_: -(d.grp == g_).sum()):
        mk = (d.grp == g).values
        ax.scatter(pc[mk, 0], pc[mk, 1], s=5, alpha=0.3, color=pcolors[g], linewidths=0)
    # arrows scaled so the longest one reaches 90% of the 99th percentile of the scores: they stay inside the axes
    scale = 0.9 * np.percentile(np.abs(pc[:, :2]), 99, axis=0).min() / np.hypot(load["PC1"], load["PC2"]).max()
    for f in feats:
        lx, ly = load.loc[f, "PC1"] * scale, load.loc[f, "PC2"] * scale
        ax.annotate("", xy=(lx, ly), xytext=(0, 0), arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 1.1},
                    zorder=5)
        ax.text(lx * 1.08, ly * 1.08, FEATURE_LABELS.get(f, f), fontsize=8, ha="left" if lx >= 0 else "right",
                va="bottom" if ly >= 0 else "top", path_effects=halo(3), zorder=6)
    ax.axhline(0, color="#dddddd", lw=0.6, zorder=0)
    ax.axvline(0, color="#dddddd", lw=0.6, zorder=0)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title("Biplot: genomes and feature loadings")
    hs = [Patch(color=pcolors[g], label=f"{g} ({(d.grp == g).sum():,})") for g in phyla]
    if OTHER in pcolors:
        hs.append(Patch(color=GREY, label=f"other phyla ({(d.grp == OTHER).sum():,})"))
    ax.legend(handles=hs, title="Phylum", fontsize=8, title_fontsize=8.5, frameon=False, bbox_to_anchor=(1.01, 1),
              loc="upper left")
    o.save(fig, "F17b_pca_biplot")
    t = load.copy()
    t.loc["explained_variance"] = ev
    o.table(t.round(3).reset_index().rename(columns={"index": "feature"}), "T06a_global_pca")

    # ---- QC overlays: assembly and annotation quality are not PCA inputs; colour the same points by them
    # (column, colour-bar label, transform, True if low values are the poor-quality end)
    qc = [("busco_complete", "BUSCO complete (%), fungi_odb12", None, True), ("N50_bp", "log10 N50 (bp)", np.log10, True),
          ("contig_count", "log10 contigs", np.log10, False)]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4), sharex=True, sharey=True)
    for ax, (col, title, fn, low_bad) in zip(axs, qc):
        v = d[col].astype(float).values
        v = fn(np.where(v > 0, v, np.nan)) if fn else v
        # poor-quality genomes drawn last (on top); colour map runs from poor (dark) to good (yellow)
        key = np.nan_to_num(v, nan=np.nanmedian(v))
        o_ = np.argsort(-key if low_bad else key, kind="mergesort")
        lo, hi = np.nanpercentile(v, [2, 98])
        sc = ax.scatter(pc[o_, 0], pc[o_, 1], c=v[o_], s=5, cmap="viridis" if low_bad else "viridis_r",
                        vmin=lo, vmax=hi, linewidths=0, alpha=0.7)
        plt.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
        ax.set_title(title)
        ax.set_xlabel(xlab)
    axs[0].set_ylabel(ylab)
    fig.suptitle("Quality overlays on the global PCA (not PCA inputs)", fontweight="bold", y=1.02)
    o.save(fig, "F17c_pca_qc")
    r = d[["PC1", "PC2", "busco_complete", "N50_bp", "contig_count", "PHYLUM", "subphylum"]].copy()
    r["N50_bp"] = np.log10(r["N50_bp"].where(r["N50_bp"] > 0))
    r["contig_count"] = np.log10(r["contig_count"].where(r["contig_count"] > 0))
    qrows = []
    for name, mk in [("all fungi", r.index == r.index), ("Pezizomycotina", r.subphylum == "Pezizomycotina"),
                     ("Saccharomycotina", r.subphylum == "Saccharomycotina"), ("Basidiomycota", r.PHYLUM == "Basidiomycota")]:
        if mk.sum() < 10:
            continue
        cm = r.loc[mk, ["PC1", "PC2", "busco_complete", "N50_bp", "contig_count"]].corr(method="spearman")
        for ax_ in ("PC1", "PC2"):
            qrows.append({"subset": name, "genomes": int(mk.sum()), "axis": ax_, "busco_complete": cm.loc[ax_, "busco_complete"],
                          "log10_N50": cm.loc[ax_, "N50_bp"], "log10_contigs": cm.loc[ax_, "contig_count"]})
    o.table(pd.DataFrame(qrows).round(3),
            "T06b_qc_correlation")

    # ---- one ordination per major phylum, coloured by class: within-phylum PCA | the same genomes on global axes
    rows = []
    sets = [(p, d[d.PHYLUM == p], "subphylum") for p in PANEL_PHYLA]
    sets.append((OTHER_PHYLA, d[~d.PHYLUM.isin(PANEL_PHYLA)], "PHYLUM"))
    for name, sub, by in sets:
        if len(sub) < 5:
            continue
        lab, sections, color, marker = panel_groups(sub, by)
        wpc, wev, wload = pca(standardize(sub[feats]), 2)
        for f in feats:
            rows.append({"panel": name, "feature": f, "PC1": wload.loc[f, "PC1"], "PC2": wload.loc[f, "PC2"]})
        rows.append({"panel": name, "feature": "explained_variance", "PC1": wev[0], "PC2": wev[1]})
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.9), gridspec_kw={"wspace": 0.25})
        scatter_groups(a1, wpc, lab, sections, color, marker)
        a1.set_xlabel(f"{name} PC1 ({wev[0]:.0%})")
        a1.set_ylabel(f"{name} PC2 ({wev[1]:.0%})")
        top1 = wload["PC1"].abs().sort_values(ascending=False).index[:3]
        a1.set_title(f"PCA within {name} (n={len(sub):,})")
        a1.text(0.01, 0.01, "PC1 top loadings: " + ", ".join(f"{FEATURE_LABELS[f]} {wload.loc[f, 'PC1']:+.2f}" for f in top1),
                transform=a1.transAxes, fontsize=7, color="#444444", path_effects=halo())
        a2.scatter(pc[:, 0], pc[:, 1], s=3, color="#e3e3e3", linewidths=0)
        scatter_groups(a2, sub[["PC1", "PC2"]].values, lab, sections, color, marker)
        if name == "Ascomycota":
            a2.axvline(ASCO_SPLIT, color="black", ls="--", lw=0.8)
        a2.set_xlabel(xlab)
        a2.set_ylabel(ylab)
        a2.set_title("Same genomes on the global axes (grey: all fungi)")
        sectioned_legend(fig, a2, lab, sections, color, marker,
                         "Class (n)" + (", by subphylum" if by == "subphylum" else ", by phylum"))
        slug = name.lower().replace(" ", "_")
        o.save(fig, f"F17d_pca_{slug}")
    o.table(pd.DataFrame(rows).round(3), "T06c_within_phylum_pca")

    asco_split(d, feats, ev, o)
    return d


def asco_split(d, feats, ev, o):
    """Two Ascomycota modes on global PC1: yeasts (Saccharomycotina, Taphrinomycotina) below ASCO_SPLIT,
    Pezizomycotina above. Tables: cluster x subphylum and x class, single-feature separation, misfits, and the
    same split restricted to BUSCO >= 90 (on the global axes, and with the PCA refitted on that subset)."""
    a = d[d.PHYLUM == "Ascomycota"].copy()
    a["subphylum"] = a["subphylum"].fillna("unassigned")
    if not {"Saccharomycotina", "Pezizomycotina"} <= set(a.subphylum):
        return
    lowname, highname = f"low PC1 (< {ASCO_SPLIT:g})", f"high PC1 (>= {ASCO_SPLIT:g})"
    a["cluster"] = np.where(a.PC1 < ASCO_SPLIT, lowname, highname)
    hi90 = a.busco_complete >= 90

    ct = pd.crosstab(a.subphylum, a.cluster).reindex(columns=[lowname, highname], fill_value=0)
    ct90 = pd.crosstab(a.subphylum[hi90], a.cluster[hi90]).reindex(index=ct.index, columns=[lowname, highname], fill_value=0)
    t = ct.copy()
    t.columns = [f"{c} all" for c in ct.columns]
    for c in ct90.columns:
        t[f"{c} BUSCO>=90"] = ct90[c]
    t.insert(0, "genomes", ct.sum(axis=1))
    t["share low PC1"] = (ct[lowname] / ct.sum(axis=1)).round(3)
    o.table(t.reset_index(), "T06d_asco_cluster_subphylum")

    cls = a["CLASS"].fillna("(class not set)")
    cc = pd.crosstab([cls.rename("class"), a.subphylum], a.cluster).reindex(columns=[lowname, highname], fill_value=0)
    cc["genomes"] = cc.sum(axis=1)
    cc["share low PC1"] = (cc[lowname] / cc["genomes"]).round(3)
    cc["median BUSCO"] = a.groupby([cls.rename("class"), a.subphylum])["busco_complete"].median().round(1)
    o.table(cc.sort_values(["subphylum", "genomes"], ascending=[True, False]).reset_index(), "T06e_asco_cluster_class")

    # single-feature separation, Saccharomycotina vs Pezizomycotina (Taphrinomycotina and unassigned left out)
    cols = ["PC1"] + feats + ["busco_complete", "N50_bp", "contig_count"]
    sp = a[a.subphylum.isin(["Saccharomycotina", "Pezizomycotina"])]
    y = (sp.subphylum == "Saccharomycotina").values
    rows = []
    for c in cols:
        bal, acc, thr, dirn = best_threshold(sp[c].values, y)
        m90 = sp.busco_complete.values >= 90
        bal90, _, _, _ = best_threshold(sp[c].values[m90], y[m90])
        med = a.groupby("subphylum")[c].median()
        rows.append({"feature": c, "label": FEATURE_LABELS.get(c, c), "pca_input": c in feats,
                     "balanced_accuracy": bal, "accuracy": acc, "yeast_side": f"{dirn} {thr:.4g}" if dirn else "",
                     "balanced_accuracy_BUSCO>=90": bal90,
                     "median_Saccharomycotina": med.get("Saccharomycotina", np.nan),
                     "median_Pezizomycotina": med.get("Pezizomycotina", np.nan),
                     "median_Taphrinomycotina": med.get("Taphrinomycotina", np.nan)})
    sep = pd.DataFrame(rows).sort_values("balanced_accuracy", ascending=False)
    o.table(sep.round(4), "T06f_asco_feature_separation")

    # misfits: Pezizomycotina on the yeast side, Saccharomycotina on the filamentous side
    mis = a[((a.subphylum == "Pezizomycotina") & (a.PC1 < ASCO_SPLIT)) |
            ((a.subphylum == "Saccharomycotina") & (a.PC1 >= ASCO_SPLIT))].copy()
    mis = mis.reset_index().rename(columns={"index": "species_prefix"})
    keep = ["species_prefix", "subphylum", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES", "PC1", "PC2",
            "busco_complete", "genome_Mb", "n_genes", "N50_bp", "contig_count", "secreted_frac", "loc_mitochondrion"]
    mis = mis[[c for c in keep if c in mis.columns]].sort_values(["subphylum", "CLASS", "PC1"])
    o.table(mis.round(3), "T06g_asco_misfits")

    # BUSCO >= 90: same axes, and the PCA refitted on all fungi with BUSCO >= 90
    keep90 = d.busco_complete >= 90
    rpc, rev, rload = pca(standardize(d.loc[keep90, feats]), 2)
    a90 = a[hi90].copy()
    for j in (0, 1):
        a90[f"PC{j + 1}_refit"] = pd.Series(rpc[:, j], index=d.index[keep90]).reindex(a90.index)
    summ = []
    for label, frame, col in [("all Ascomycota, global PC1", a, "PC1"), ("BUSCO >= 90, global PC1", a90, "PC1"),
                              ("BUSCO >= 90, PCA refitted on BUSCO >= 90 fungi: PC1", a90, "PC1_refit"),
                              ("BUSCO >= 90, PCA refitted on BUSCO >= 90 fungi: PC2", a90, "PC2_refit")]:
        s_ = frame[frame.subphylum.isin(["Saccharomycotina", "Pezizomycotina"])]
        yy = (s_.subphylum == "Saccharomycotina").values
        bal, acc, thr, dirn = best_threshold(s_[col].values, yy)
        valley, m1, m2 = density_valley(frame[col].values)
        low = (s_[col] < ASCO_SPLIT).values
        summ.append({"subset": label, "Ascomycota": len(frame),
                     "Saccharomycotina": int(yy.sum()), "Pezizomycotina": int((~yy).sum()),
                     "Taphrinomycotina": int((frame.subphylum == "Taphrinomycotina").sum()),
                     "explained_variance": {"PC1": ev[0], "PC1_refit": rev[0], "PC2_refit": rev[1]}[col],
                     "cut_-1_balanced_accuracy": 0.5 * (low[yy].mean() + (~low[~yy]).mean()) if col == "PC1" else np.nan,
                     "best_cut": f"{dirn} {thr:.3g}" if dirn else "", "best_cut_balanced_accuracy": bal,
                     "density_valley": valley, "low_mode": m1, "high_mode": m2})
    rl = rload.copy()
    rl.columns = ["PC1_refit", "PC2_refit"]
    o.table(pd.DataFrame(summ).round(3), "T06h_asco_busco90")
    o.table(rl.round(3).reset_index().rename(columns={"index": "feature"}), "T06i_busco90_refit_loadings")

    # ---- figure: PC1 distribution by subphylum (all, BUSCO >= 90, refit)
    subs = [s_ for s_ in ["Pezizomycotina", "Saccharomycotina", "Taphrinomycotina", "unassigned"] if (a.subphylum == s_).any()]
    # third panel: the refitted axis that separates the subphyla best (balanced accuracy, table T06h)
    rcol = max(["PC1_refit", "PC2_refit"], key=lambda c: next(r_["best_cut_balanced_accuracy"] for r_ in summ
                                                                if "refitted" in r_["subset"] and r_["subset"].endswith(c[:3])))
    fig, axs = plt.subplots(1, 3, figsize=(15, 3.9))
    for ax, (title, frame, col) in zip(axs, [("All Ascomycota", a, "PC1"), ("BUSCO >= 90 only", a90, "PC1"),
                                              (f"BUSCO >= 90, PCA refitted ({rcol[:3]})", a90, rcol)]):
        v = frame[col].dropna()
        bins = np.linspace(v.min(), v.max(), 60) if len(v) > 1 else 10
        for s_ in subs:
            vv = frame.loc[frame.subphylum == s_, col].dropna()
            if len(vv):
                ax.hist(vv, bins=bins, histtype="stepfilled", alpha=0.45, color=SUBPH_COLORS[s_],
                        label=f"{s_} ({len(vv):,})")
                ax.hist(vv, bins=bins, histtype="step", color=SUBPH_COLORS[s_], lw=1)
        if col == "PC1":
            ax.axvline(ASCO_SPLIT, color="black", ls="--", lw=0.8)
        # the density valley only on the full set: with BUSCO >= 90 the yeast mode is too small to show
        valley = density_valley(v.values)[0] if frame is a else np.nan
        if not np.isnan(valley):
            ax.axvline(valley, color="#888888", ls=":", lw=1)
        ax.set_title(title)
        ax.set_xlabel(f"{col[:3]} refitted on BUSCO >= 90 fungi" if col.endswith("_refit") else xlab_ev(ev))
        ax.set_ylim(0, ax.get_ylim()[1] * 1.45)               # headroom: the legend sits above the bars
        ax.legend(frameon=False, fontsize=8, loc="upper left")
    axs[0].set_ylabel("Ascomycota genomes")
    fig.suptitle("Ascomycota on PC1 (dashed: cut at -1; dotted, left: density valley between the two modes)", fontweight="bold", y=1.03)
    o.save(fig, "F17e_asco_pc1")

    # ---- figure: the eight best single features, by subphylum, with the best cut
    top = [r_ for r_ in sep.feature if r_ != "PC1"][:8]
    fig, axs = plt.subplots(2, 4, figsize=(15, 6.4))
    for ax, c in zip(axs.flat, top):
        v = a[c].astype(float)
        logx = c in LOG_FEATS
        vv = np.log10(v.where(v > 0)) if logx else v
        lo, hi = np.nanpercentile(vv, [0.5, 99.5])           # a few extreme genomes would squash the x axis
        bins = np.linspace(lo, hi, 50)
        for s_ in subs:
            x_ = vv[a.subphylum == s_].dropna()
            if len(x_):
                ax.hist(x_, bins=bins, histtype="stepfilled", alpha=0.45, color=SUBPH_COLORS[s_], label=s_)
        rr = sep.set_index("feature").loc[c]
        thr = float(rr["yeast_side"].split()[-1]) if rr["yeast_side"] else np.nan
        if not np.isnan(thr):
            ax.axvline(np.log10(thr) if logx else thr, color="black", ls="--", lw=0.8)
        ax.set_title(f"{FEATURE_LABELS.get(c, c)}: {rr['balanced_accuracy']:.1%}", fontsize=9)
        ax.set_xlabel(("log10 " if logx else "") + FEATURE_LABELS.get(c, c), fontsize=8)
    axs.flat[0].legend(frameon=False, fontsize=7.5)
    fig.suptitle("Best single features separating Saccharomycotina from Pezizomycotina (title: balanced accuracy of the best cut)",
                 fontweight="bold", y=1.01)
    fig.tight_layout()
    o.save(fig, "F17f_asco_features")


def xlab_ev(ev):
    return f"global PC1 ({ev[0]:.0%})"


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
    o = Out(out, "global")
    page_feature_space(core0, rd("deeploc_genome"), rd("idp_genome"), o)
    print(f"global: {len(o.made)} outputs", flush=True)
    made["global"] = o.made
    return made


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], tuple(sys.argv[3].split(",")) if len(sys.argv) > 3 else LEVELS)
