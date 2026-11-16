#!/usr/bin/env python3
"""Figures and tables for the Quarto summary, from the parquet files of build_summary_tables.py.

    python make_figures.py DATA_DIR OUT_DIR   ->  OUT_DIR/figures/*.png|svg, OUT_DIR/tables/*.tsv

Unit: one representative genome per species. Group means/medians are taken over genomes, never over
proteins, so large genomes do not dominate. Nothing here is a statistical test: the figures describe
the database. The .qmd pages state n for every figure.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from panel import PANEL  # noqa: E402

TOP_PHYLA = 8
# Okabe-Ito first, then muted extras; 'Other' is grey
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#F0E442", "#8C564B"]
GREY = "#9a9a9a"
LOC_COLORS = {"cytoplasm": "#4C78A8", "nucleus": "#F58518", "extracellular": "#54A24B", "cell_membrane": "#E45756",
              "mitochondrion": "#B279A2", "plastid": "#72B7B2", "endoplasmic_reticulum": "#EECA3B",
              "lysosome_vacuole": "#9D755D", "golgi": "#BAB0AC", "peroxisome": "#FF9DA6"}

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 10, "axes.titleweight": "bold", "figure.dpi": 110,
                     "savefig.bbox": "tight", "svg.fonttype": "none"})


class Out:
    def __init__(self, root):
        self.fig = Path(root) / "figures"
        self.tab = Path(root) / "tables"
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


def group_phyla(core):
    order = core["PHYLUM"].fillna("Unknown").value_counts()
    top = list(order.index[:TOP_PHYLA])
    core = core.copy()
    core["grp"] = core["PHYLUM"].fillna("Unknown").where(core["PHYLUM"].fillna("Unknown").isin(top), "Other")
    levels = top + (["Other"] if (core["grp"] == "Other").any() else [])
    colors = {p: (PALETTE[i % len(PALETTE)] if p != "Other" else GREY) for i, p in enumerate(levels)}
    return core, levels, colors


def box_by_group(ax, core, col, levels, colors, log=False, ylabel=None):
    data = [core.loc[core.grp == g, col].dropna().values for g in levels]
    bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False,
                    medianprops={"color": "black", "lw": 1.4})
    rng = np.random.default_rng(1)
    for i, (d, g) in enumerate(zip(data, levels), start=1):
        bp["boxes"][i - 1].set(facecolor=colors[g], alpha=0.35, edgecolor=colors[g])
        if len(d):
            ax.scatter(i + rng.uniform(-0.18, 0.18, len(d)), d, s=3, color=colors[g], alpha=0.35, linewidths=0)
    ax.set_xticks(range(1, len(levels) + 1))
    ax.set_xticklabels([f"{g}\n(n={len(d)})" for g, d in zip(levels, data)], rotation=35, ha="right")
    if log:
        ax.set_yscale("log")
    if ylabel:
        ax.set_ylabel(ylabel)


def heat(ax, mat, cmap="viridis", vmin=None, vmax=None, label=None, fmt=None):
    im = ax.imshow(mat.values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=8)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=60, ha="right", fontsize=8)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    if label:
        cb.set_label(label, fontsize=8)
    if fmt:
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat.values[i, j]
                if np.isfinite(v):
                    ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=6, color="white")
    return im


def stacked_by_group(ax, frac, colors_map, levels, ylabel):
    bottom = np.zeros(len(levels))
    for lab in frac.columns:
        v = frac.loc[levels, lab].values
        ax.bar(range(len(levels)), v, bottom=bottom, color=colors_map.get(lab, GREY), label=lab, width=0.75)
        bottom += v
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels(levels, rotation=35, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 1)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False, fontsize=8)


def per_genome_fraction(df, label_col, count_col="n_proteins"):
    piv = df.pivot_table(index="species_prefix", columns=label_col, values=count_col, aggfunc="sum", fill_value=0)
    return piv.div(piv.sum(axis=1), axis=0)


def group_mean(frac, core, levels):
    j = frac.join(core.set_index("species_prefix")["grp"], how="inner")
    return j.groupby("grp").mean(numeric_only=True).reindex(levels).fillna(0)


# ---------------------------------------------------------------- page 1: gene content
def page_gene_content(core, levels, colors, o):
    n = core.groupby("grp").size().reindex(levels)
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    ax.bar(range(len(levels)), n.values, color=[colors[g] for g in levels])
    for i, v in enumerate(n.values):
        ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels(levels, rotation=35, ha="right")
    ax.set_ylabel("Representative genomes")
    ax.set_title("One genome per species, by phylum")
    o.save(fig, "F01_genomes_by_phylum")

    for col, name, ylab, log in [("n_genes", "F02_gene_count", "Protein-coding genes", True),
                                 ("mean_len", "F04_protein_length", "Mean protein length (aa)", False),
                                 ("busco_complete", "F05_busco", "BUSCO complete (%)", False),
                                 ("gc_pct", "F05b_gc", "Genome GC (%)", False)]:
        fig, ax = plt.subplots(figsize=(6.8, 3.6))
        box_by_group(ax, core, col, levels, colors, log=log, ylabel=ylab)
        ax.set_title(ylab + " by phylum")
        o.save(fig, name)

    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    for g in levels:
        d = core[core.grp == g]
        ax.scatter(d.total_length_bp / 1e6, d.n_genes, s=6, alpha=0.5, color=colors[g], label=g, linewidths=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Genome size (Mb)")
    ax.set_ylabel("Protein-coding genes")
    ax.set_title("Genome size vs gene number")
    ax.legend(frameon=False, fontsize=7, markerscale=2)
    o.save(fig, "F03_genome_size_vs_genes")

    t = core.groupby("grp").agg(genomes=("species_prefix", "size"), genome_Mb=("total_length_bp", lambda x: x.median() / 1e6),
                                genes=("n_genes", "median"), protein_len=("mean_len", "median"),
                                gc_pct=("gc_pct", "median"), busco=("busco_complete", "median")).reindex(levels).round(2)
    o.table(t.reset_index().rename(columns={"grp": "phylum"}), "T01_phylum_overview")


# ---------------------------------------------------------------- page 2: domains
def page_domains(core, levels, colors, pfam, names, arch, o):
    n_genes = core.set_index("species_prefix")["n_genes"]
    a = arch.pivot_table(index="species_prefix", columns="n_domains", values="n_proteins", aggfunc="sum", fill_value=0)
    a = a.reindex(core.species_prefix).fillna(0)
    a.insert(0, 0, (n_genes.reindex(a.index) - a.sum(axis=1)).clip(lower=0))
    frac = a.div(a.sum(axis=1), axis=0)
    frac.columns = ["0" if c == 0 else ("5+" if c >= 5 else str(int(c))) for c in frac.columns]
    gm = group_mean(frac, core, levels)
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    cmap = plt.get_cmap("Blues")
    stacked_by_group(ax, gm, {c: cmap(0.2 + 0.16 * i) for i, c in enumerate(gm.columns)}, levels,
                     "Share of proteins (genome mean)")
    ax.set_title("Pfam domains per protein (0 = no Pfam hit)")
    o.save(fig, "F06_domains_per_protein")

    wide = pfam.pivot_table(index="species_prefix", columns="pfam_acc", values="n_proteins", aggfunc="sum", fill_value=0)
    wide = wide.reindex(core.species_prefix).fillna(0)
    grp = core.set_index("species_prefix")["grp"].reindex(wide.index)
    present = (wide > 0).groupby(grp).mean().reindex(levels)
    label = names.set_index("pfam_acc")["pfam_id"]
    top = present.mean().sort_values(ascending=False).index[:30]
    m = present[top].T
    m.index = [f"{a} {label.get(a, '')}" for a in m.index]
    fig, ax = plt.subplots(figsize=(6.8, 7.2))
    heat(ax, m, cmap="magma", vmin=0, vmax=1, label="Fraction of genomes with the domain")
    ax.set_title("Most widespread Pfam domains")
    o.save(fig, "F07_pfam_prevalence")

    # copies per 10,000 genes (mean over genomes): removes the genome-size effect (reduced genomes look
    # low on every domain otherwise) and avoids integer-median ties (median 1 vs 2 copies)
    per10k = wide.div(n_genes.reindex(wide.index), axis=0) * 1e4
    cnt = per10k.groupby(grp).mean().reindex(levels)
    keep = [c for c in wide.columns if (present[c].max() >= 0.3) and per10k[c].mean() >= 0.5]
    z = np.log1p(cnt[keep])
    z = (z - z.mean()) / z.std(ddof=0).replace(0, np.nan)
    var = z.var().sort_values(ascending=False).index[:40]
    m2 = z[var].T
    m2.index = [f"{a} {label.get(a, '')}" for a in m2.index]
    fig, ax = plt.subplots(figsize=(6.8, 8.6))
    heat(ax, m2, cmap="RdBu_r", vmin=-2.5, vmax=2.5, label="z-score of log mean copies per 10,000 genes")
    ax.set_title("Domains that differ most between phyla")
    o.save(fig, "F08_pfam_phylum_variance")

    nfam = (wide > 0).sum(axis=1).rename("n_pfam_families_in_top_set")
    d = core.set_index("species_prefix").join(nfam)
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    box_by_group(ax, d.reset_index(), "n_pfam_families_in_top_set", levels, colors,
                 ylabel="Families present (of the top set)")
    ax.set_title("Pfam repertoire breadth")
    o.save(fig, "F09_pfam_breadth")

    t = pd.DataFrame({"pfam_acc": top, "pfam_id": [label.get(a, "") for a in top],
                      "genomes_with_domain": [int((wide[a] > 0).sum()) for a in top],
                      "median_copies": [float(wide[a].median()) for a in top]})
    o.table(t, "T02_top_pfam")
    return wide, label


# ---------------------------------------------------------------- page 3: localization
def page_localization(core, levels, colors, deeploc, wolf, o):
    fr = per_genome_fraction(deeploc, "label")
    gm = group_mean(fr, core, levels)
    gm = gm[[c for c in LOC_COLORS if c in gm.columns]]
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    stacked_by_group(ax, gm, LOC_COLORS, levels, "Share of proteins (genome mean)")
    ax.set_title("DeepLoc primary localization (argmax of class probabilities)")
    o.save(fig, "F10_deeploc_by_phylum")
    o.table(gm.round(4).reset_index().rename(columns={"grp": "phylum"}), "T03_deeploc_by_phylum")

    cls = core.groupby("CLASS").size()
    big = cls[cls >= 20].index
    j = fr.join(core.set_index("species_prefix")["CLASS"], how="inner")
    cm = j[j.CLASS.isin(big)].groupby("CLASS").mean(numeric_only=True)
    cm = cm[[c for c in LOC_COLORS if c in cm.columns]]
    if len(cm):
        fig, ax = plt.subplots(figsize=(6.5, max(3, 0.25 * len(cm) + 1.5)))
        heat(ax, cm, cmap="YlGnBu", vmin=0, label="Share of proteins")
        ax.set_title("DeepLoc localization by class (classes with 20+ genomes)")
        o.save(fig, "F11_deeploc_by_class")

    wf = per_genome_fraction(wolf, "localization")
    wg = group_mean(wf, core, levels)
    keepc = wg.mean().sort_values(ascending=False).index[:10]
    wg = wg[keepc]
    cmap = plt.get_cmap("tab10")
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    stacked_by_group(ax, wg, {c: cmap(i) for i, c in enumerate(wg.columns)}, levels, "Share of proteins (genome mean)")
    ax.set_title("WoLF PSORT localization (top 10 classes)")
    o.save(fig, "F12_wolfpsort_by_phylum")

    core = core.copy()
    core["secreted_frac"] = core["n_secreted_strict"] / core["n_genes"]
    core["signalp_frac"] = core["n_signalp_09"] / core["n_genes"]
    for col, name, ylab in [("secreted_frac", "F13_secretome_fraction", "Strict secretome / genes"),
                            ("signalp_frac", "F13b_signal_peptide_fraction", "SignalP (p>=0.9) / genes")]:
        fig, ax = plt.subplots(figsize=(6.8, 3.6))
        box_by_group(ax, core, col, levels, colors, ylabel=ylab)
        ax.set_title(ylab + " by phylum")
        o.save(fig, name)


# ---------------------------------------------------------------- page 4: phylum patterns
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
    feats = [f for f in FEATURES if f in d.columns]
    med = d.groupby("grp")[feats].median().reindex(levels)
    z = (med - med.mean()) / med.std(ddof=0).replace(0, np.nan)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    heat(ax, z.T, cmap="RdBu_r", vmin=-2.5, vmax=2.5, label="z-score of phylum median", fmt=None)
    ax.set_title("Phylum profiles across genome-level features")
    o.save(fig, "F15_phylum_feature_heatmap")
    o.table(med.round(4).reset_index().rename(columns={"grp": "phylum"}), "T04_phylum_feature_medians")

    cls = d.groupby("CLASS").size()
    big = cls[cls >= 20].index
    cm = d[d.CLASS.isin(big)].groupby("CLASS")[feats].median()
    if len(cm) > 1:
        zc = (cm - cm.mean()) / cm.std(ddof=0).replace(0, np.nan)
        fig, ax = plt.subplots(figsize=(7.2, max(3.5, 0.28 * len(zc) + 1.6)))
        heat(ax, zc, cmap="RdBu_r", vmin=-2.5, vmax=2.5, label="z-score of class median")
        ax.set_title("Class profiles (classes with 20+ genomes)")
        o.save(fig, "F16_class_feature_heatmap")
        o.table(cm.round(4).reset_index(), "T05_class_feature_medians")

    x = d[feats].copy()
    x = x.fillna(x.median())
    x = np.log1p(x.clip(lower=0))
    x = (x - x.mean()) / x.std(ddof=0).replace(0, np.nan)
    x = x.fillna(0)
    u, s, vt = np.linalg.svd(x.values, full_matrices=False)
    pc = u[:, :2] * s[:2]
    ev = (s ** 2) / (s ** 2).sum()
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    for g in levels:
        mk = (d["grp"] == g).values
        ax.scatter(pc[mk, 0], pc[mk, 1], s=7, alpha=0.55, color=colors[g], label=g, linewidths=0)
    ax.set_xlabel(f"PC1 ({ev[0]:.0%})")
    ax.set_ylabel(f"PC2 ({ev[1]:.0%})")
    ax.set_title("Genomes in feature space (log, z-scored)")
    ax.legend(frameon=False, fontsize=7, markerscale=2)
    o.save(fig, "F17_feature_pca")
    load = pd.DataFrame(vt[:2].T, index=feats, columns=["PC1", "PC2"]).round(3)
    o.table(load.reset_index().rename(columns={"index": "feature"}), "T06_pca_loadings")


# ---------------------------------------------------------------- page 5: functional domains
def page_functional(core, levels, colors, wide, label, cazy, tf, merops, idp, o):
    ng = core.set_index("species_prefix")["n_genes"]
    grp = core.set_index("species_prefix")["grp"]

    cz = cazy.pivot_table(index="species_prefix", columns="cazy_class", values="n_proteins", aggfunc="sum", fill_value=0)
    cz = cz.reindex(core.species_prefix).fillna(0)
    classes = [c for c in ["GH", "GT", "PL", "CE", "AA", "CBM"] if c in cz.columns]
    if classes:
        fig, axs = plt.subplots(1, len(classes), figsize=(2.2 * len(classes), 3.6), sharey=False)
        axs = np.atleast_1d(axs)
        for ax, c in zip(axs, classes):
            tmp = core.set_index("species_prefix").assign(v=cz[c] / ng).reset_index()
            box_by_group(ax, tmp, "v", levels, colors)
            ax.set_title(c)
            ax.set_xticklabels([])
        axs[0].set_ylabel("Proteins per gene")
        axs[0].set_xticklabels([f"{g}" for g in levels], rotation=60, ha="right", fontsize=7)
        fig.suptitle("CAZy classes per gene, by phylum", y=1.02, fontweight="bold")
        o.save(fig, "F18_cazy_classes")

        fam = cazy.pivot_table(index="species_prefix", columns="cazy_family", values="n_proteins", aggfunc="sum", fill_value=0)
        fam = fam.reindex(core.species_prefix).fillna(0)
        med = fam.groupby(grp.reindex(fam.index)).median().reindex(levels)
        top = fam.sum().sort_values(ascending=False).index[:30]
        fig, ax = plt.subplots(figsize=(6.8, 7))
        heat(ax, np.log2(med[top].T + 1), cmap="viridis", label="log2(median proteins + 1)")
        ax.set_title("Top 30 CAZy families")
        o.save(fig, "F19_cazy_families")
        o.table(med[top].T.reset_index().rename(columns={"cazy_family": "family"}), "T07_cazy_family_medians")

    mer = merops.pivot_table(index="species_prefix", columns="merops_class", values="n_proteins", aggfunc="sum", fill_value=0)
    mer = mer.reindex(core.species_prefix).fillna(0)
    if mer.shape[1]:
        mf = mer.div(ng.reindex(mer.index), axis=0)
        gm = mf.groupby(grp.reindex(mf.index)).mean().reindex(levels)
        fig, ax = plt.subplots(figsize=(6.8, 3.6))
        bottom = np.zeros(len(levels))
        cmap = plt.get_cmap("Set2")
        for i, c in enumerate(gm.columns):
            ax.bar(range(len(levels)), gm[c].values, bottom=bottom, label=c, color=cmap(i % 8), width=0.75)
            bottom += gm[c].values
        ax.set_xticks(range(len(levels)))
        ax.set_xticklabels(levels, rotation=35, ha="right")
        ax.set_ylabel("MEROPS proteins per gene (genome mean)")
        ax.legend(title="class", bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
        ax.set_title("Peptidase classes (MEROPS first letter)")
        o.save(fig, "F20_merops_classes")

    tfw = tf.pivot_table(index="species_prefix", columns="tf_family", values="n_proteins", aggfunc="sum", fill_value=0)
    tfw = tfw.reindex(core.species_prefix).fillna(0)
    if tfw.shape[1]:
        med = tfw.groupby(grp.reindex(tfw.index)).median().reindex(levels)
        top = tfw.sum().sort_values(ascending=False).index[:25]
        fig, ax = plt.subplots(figsize=(6.8, 6.2))
        heat(ax, np.log2(med[top].T + 1), cmap="viridis", label="log2(median proteins + 1)")
        ax.set_title("Transcription-factor families")
        o.save(fig, "F21_tf_families")
        o.table(med[top].T.reset_index().rename(columns={"tf_family": "family"}), "T08_tf_family_medians")

    pan = []
    for k, (cat, accs) in enumerate(PANEL.items(), start=1):
        have = [a for a in accs if a in wide.columns]
        if not have:
            continue
        med = wide[have].groupby(grp.reindex(wide.index)).mean().reindex(levels)  # mean proteins per genome
        prev = (wide[have] > 0).groupby(grp.reindex(wide.index)).mean().reindex(levels)
        m = np.log2(med.T + 1)
        m.index = [f"{a} {accs[a]}" for a in m.index]
        fig, ax = plt.subplots(figsize=(6.8, max(2.2, 0.34 * len(have) + 1.4)))
        heat(ax, m, cmap="viridis", label="log2(mean proteins per genome + 1)")
        pv = prev.T.values
        for i in range(pv.shape[0]):
            for j in range(pv.shape[1]):
                ax.text(j, i, f"{pv[i, j]:.0%}", ha="center", va="center", fontsize=6,
                        color="white" if m.values[i, j] < m.values.max() * 0.6 else "black")
        ax.set_title(cat + " (cell text: % of genomes with the domain)")
        slug = "".join(ch if ch.isalnum() else "_" for ch in cat.lower()).strip("_")
        o.save(fig, f"F22_{k:02d}_panel_{slug}")
        for a in have:
            for g in levels:
                pan.append({"category": cat, "pfam_acc": a, "name": accs[a], "group": g,
                            "mean_proteins": float(med.loc[g, a]), "fraction_genomes": float(prev.loc[g, a])})
    if pan:
        o.table(pd.DataFrame(pan), "T09_functional_panel")

    d = core.set_index("species_prefix").join(idp.set_index("species_prefix")[["mean_idp_fraction", "frac_idp_gt_0_3"]]).reset_index()
    for col, name, ylab in [("mean_idp_fraction", "F23_idp_fraction", "Mean IDR fraction per protein"),
                            ("frac_idp_gt_0_3", "F23b_idp_gt03", "Share of proteins with IDR fraction > 0.3")]:
        fig, ax = plt.subplots(figsize=(6.8, 3.6))
        box_by_group(ax, d, col, levels, colors, ylabel=ylab)
        ax.set_title(ylab + " by phylum")
        o.save(fig, name)


def main(data, out):
    data = Path(data)
    o = Out(out)
    rd = lambda n: pd.read_parquet(data / f"{n}.parquet")  # noqa: E731
    core = rd("genome_core")
    core = core[core["n_genes"].notna()]  # representatives without a proteome are not plotted
    core, levels, colors = group_phyla(core)
    page_gene_content(core, levels, colors, o)
    wide, label = page_domains(core, levels, colors, rd("pfam_genome"), rd("pfam_names"), rd("arch_genome"), o)
    deeploc = rd("deeploc_genome")
    page_localization(core, levels, colors, deeploc, rd("wolfpsort_genome"), o)
    idp = rd("idp_genome")
    page_patterns(core, levels, colors, deeploc, idp, o)
    page_functional(core, levels, colors, wide, label, rd("cazy_genome"), rd("tf_genome"), rd("merops_genome"), idp, o)
    print(f"{len(o.made)} outputs")
    return o.made


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
