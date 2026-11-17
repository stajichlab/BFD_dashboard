"""Codon and amino-acid composition figures (topic 'composition') from codon_genome and aa_genome.
Frequencies are per genome (rows sum to 1, stop codons dropped). Groups are summarised over genomes."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from figs_common import join_groups, safe
from make_figures import box_by_group, heat, parent_legend, stacked_by_group, width

TOPIC = "composition"
BASES = "TCAG"
AA_STRING = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
# Standard genetic code only (table 1). CTG-clade yeasts (CTG = Ser) and mitochondrial codes would be mislabelled.
CODON_AA = {a + b + c: AA_STRING[i] for i, (a, b, c) in enumerate((a, b, c) for a in BASES for b in BASES for c in BASES)}
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_COLORS = {a: plt.get_cmap("tab20")(i / 20 + 0.01) for i, a in enumerate(AA_ORDER)}


def _codon_cols(codon):
    cols = [c for c in codon.columns if c.startswith("codon_") and c[6:] in CODON_AA and CODON_AA[c[6:]] != "*"]
    return sorted(cols, key=lambda c: (CODON_AA[c[6:]], c))


def _rscu(f):
    """Relative synonymous codon usage per genome: freq * family size / family sum."""
    out = f.copy()
    fam = {}
    for c in f.columns:
        fam.setdefault(CODON_AA[c[6:]], []).append(c)
    for aa, cs in fam.items():
        s = f[cs].sum(axis=1).replace(0, np.nan)
        for c in cs:
            out[c] = f[c] / s * len(cs)
    return out


def _gc3(codon):
    cols = _codon_cols(codon)
    third = [c[-1] in "GC" for c in cols]
    f = codon[cols]
    return (f.loc[:, third].sum(axis=1) / f.sum(axis=1).replace(0, np.nan)).rename("gc3")


def _codon_heat(codon, core, levels, colors):
    cols = _codon_cols(codon)
    d = join_groups(core, codon[["species_prefix"] + cols])
    m = d.groupby("grp")[cols].mean().reindex(levels).T
    m.index = [f"{c[6:]} {CODON_AA[c[6:]]}" for c in m.index]
    fig, ax = plt.subplots(figsize=(width(levels, 3.5), 11))
    heat(ax, m, cmap="viridis", label="Mean codon frequency", colors=colors)
    ax.set_yticklabels(m.index, fontsize=6, family="monospace")
    ax.set_title("Codon usage by group")
    return fig


def _rscu_heat(codon, core, levels, colors):
    cols = _codon_cols(codon)
    r = _rscu(codon[cols])
    d = join_groups(core, r.assign(species_prefix=codon["species_prefix"].values))
    m = d.groupby("grp")[cols].mean().reindex(levels).T
    m.index = [f"{c[6:]} {CODON_AA[c[6:]]}" for c in m.index]
    fig, ax = plt.subplots(figsize=(width(levels, 3.5), 11))
    heat(ax, m, cmap="PuOr_r", vmin=0, vmax=max(2.0, float(np.nanmax(m.values))), label="Mean RSCU", colors=colors)
    ax.set_yticklabels(m.index, fontsize=6, family="monospace")
    ax.set_title("Relative synonymous codon usage by group")
    return fig


def _rscu_boxes(codon):
    cols = _codon_cols(codon)
    r = _rscu(codon[cols])
    fig, ax = plt.subplots(figsize=(11, 3.8))
    data = [r[c].dropna().values for c in cols]
    bp = ax.boxplot(data, patch_artist=True, widths=0.7, showfliers=False, medianprops={"color": "black", "lw": 1})
    aas = [CODON_AA[c[6:]] for c in cols]
    for patch, a in zip(bp["boxes"], aas):
        patch.set(facecolor=AA_COLORS[a], alpha=0.6, edgecolor="grey")
    ax.axhline(1, color="grey", lw=0.6, ls=":")
    ax.set_xticks(range(1, len(cols) + 1))
    ax.set_xticklabels([f"{c[6:]} {a}" for c, a in zip(cols, aas)], rotation=90, fontsize=6, family="monospace")
    ax.set_ylabel("RSCU per genome")
    ax.set_title("Relative synonymous codon usage, all genomes (colour = amino acid)")
    return fig


def _aa_cols(aa):
    return [f"aa_{a}" for a in AA_ORDER if f"aa_{a}" in aa.columns]


def _aa_stack(aa, core, levels, colors):
    cols = _aa_cols(aa)
    d = join_groups(core, aa[["species_prefix"] + cols])
    m = d.groupby("grp")[cols].mean().reindex(levels).fillna(0)
    m = m.div(m.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    m.columns = [c[3:] for c in cols]
    fig, ax = plt.subplots(figsize=(width(levels, 6.5), 4.2))
    stacked_by_group(ax, m, AA_COLORS, levels, "Mean amino acid frequency", group_colors=colors)
    ax.set_ylim(0, 1.0)
    ax.set_title("Amino acid composition by group")
    return fig


def _aa_heat(aa, core, levels, colors):
    cols = _aa_cols(aa)
    d = join_groups(core, aa[["species_prefix"] + cols])
    m = d.groupby("grp")[cols].mean().reindex(levels)
    sd = m.std(ddof=0).replace(0, np.nan)
    z = ((m - m.mean()) / sd).fillna(0).T
    z.index = [c[3:] for c in z.index]
    fig, ax = plt.subplots(figsize=(width(levels, 3.2), 5))
    heat(ax, z, cmap="RdBu_r", vmin=-2.5, vmax=2.5, label="z-score of group mean", colors=colors)
    ax.set_title("Amino acid composition, z-score across groups")
    return fig


def _gc3_box(g3, core, levels, colors):
    d = join_groups(core, g3.reset_index())
    fig, ax = plt.subplots(figsize=(width(levels), 3.8))
    box_by_group(ax, d, "gc3", levels, colors, ylabel="GC3 (G or C at codon position 3)")
    ax.set_title("GC3 per genome")
    parent_legend(ax, core, loc="best")
    return fig


def _gc3_vs_gc(g3, core, levels, colors):
    d = join_groups(core, g3.reset_index())
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    for g in levels:
        m = d[d.grp == g]
        ax.scatter(m["gc_pct"], m["gc3"] * 100, s=6, alpha=0.5, color=colors[g], linewidths=0)
    ax.set_xlabel("Genome GC (%)")
    ax.set_ylabel("GC3 (%)")
    ax.set_title("GC3 vs genome GC")
    parent_legend(ax, core, loc="best")
    return fig


def _aa_vs_gc(aa, g3, core, levels, colors):
    cols = _aa_cols(aa)
    d = join_groups(core, aa[["species_prefix"] + cols]).merge(g3.reset_index(), on="species_prefix")
    fig, axes = plt.subplots(4, 5, figsize=(11, 8), sharex=True)
    for ax, a in zip(axes.flat, AA_ORDER):
        c = f"aa_{a}"
        if c not in d:
            ax.axis("off")
            continue
        for g in levels:
            m = d[d.grp == g]
            ax.scatter(m["gc3"] * 100, m[c] * 100, s=3, alpha=0.4, color=colors[g], linewidths=0)
        ax.set_title(a, fontsize=8)
        ax.tick_params(labelsize=6)
    for ax in axes[-1]:
        ax.set_xlabel("GC3 (%)", fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel("aa (%)", fontsize=7)
    fig.suptitle("Amino acid frequency vs GC3, each point a genome", fontsize=10, fontweight="bold")
    return fig


def _aa_class_box(aa, core, levels, colors, aas, name_label):
    cols = [f"aa_{a}" for a in aas if f"aa_{a}" in aa.columns]
    d = join_groups(core, aa[["species_prefix"] + cols])
    d["v"] = d[cols].sum(axis=1)
    fig, ax = plt.subplots(figsize=(width(levels), 3.8))
    box_by_group(ax, d, "v", levels, colors, ylabel=f"{name_label} (fraction of residues)")
    ax.set_title(name_label)
    return fig


def page_composition(core, levels, colors, codon, aa, o):
    g3 = _gc3(codon).set_axis(codon["species_prefix"].values)
    g3.index.name = "species_prefix"
    safe(o, "C01_codon_heatmap", TOPIC, lambda: _codon_heat(codon, core, levels, colors))
    safe(o, "C01b_rscu_heatmap", TOPIC, lambda: _rscu_heat(codon, core, levels, colors))
    safe(o, "C02_rscu_boxes", TOPIC, lambda: _rscu_boxes(codon))
    safe(o, "C03_aa_composition", TOPIC, lambda: _aa_stack(aa, core, levels, colors))
    safe(o, "C03b_aa_heatmap", TOPIC, lambda: _aa_heat(aa, core, levels, colors))
    safe(o, "C04_gc3", TOPIC, lambda: _gc3_box(g3, core, levels, colors))
    safe(o, "C04b_gc3_vs_gc", TOPIC, lambda: _gc3_vs_gc(g3, core, levels, colors))
    safe(o, "C05_aa_vs_gc", TOPIC, lambda: _aa_vs_gc(aa, g3, core, levels, colors))
    safe(o, "C06_charged_aa", TOPIC, lambda: _aa_class_box(aa, core, levels, colors, "DEKR", "Charged residues (D,E,K,R)"))
    safe(o, "C07_aromatic_aa", TOPIC, lambda: _aa_class_box(aa, core, levels, colors, "FWY", "Aromatic residues (F,W,Y)"))
