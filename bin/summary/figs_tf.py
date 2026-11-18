"""TF landscape figures (topic 'tf'), T01-T06, from tf_family_genome and tf_summary_genome.

Unit: one representative genome per species (the level pages). Group values are means or medians over genomes.
T06 is the exception: it plots every genome (scope all_genomes). Family counts are counts of proteins per
family; a protein with domains of two families counts once in each, so the family counts of a genome can add up
to more than its n_tf_proteins. No figure here says anything about regulation or binding.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figs_common import safe
from figs_scatter import scatter_all
from make_figures import GREY, box_by_group, group_mean, heat, parent_legend, per_genome_fraction, phylum_color_map, \
    phylum_order, stacked_by_group, width
from tf_data import rare_family_labels

TOPIC = "tf"
FIGURES = ["T01_tf_fraction", "T02_tf_composition", "T03_tf_family_heatmap", "T04_tf_rare_presence", "T05_tf_shannon",
           "T06_tf_vs_proteome"]
# E2F_TDP, Tub, RFX_DNA_binding, WRKY, CP2, STAT (three Pfams, one label). Labels come from config/tf_family_map.tsv.
RARE_ACCESSIONS = ["PF02319", "PF01167", "PF02257", "PF03106", "PF04516", "PF01017", "PF02864", "PF21354"]
TOP_FAMILIES = 12
OTHER_FAMILIES = "Other TF families"
LABEL_MAX = 44
CAT = [plt.get_cmap("tab20")(i) for i in (0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 1, 3)]


def _short(s):
    s = str(s).replace("Fungal-specific transcription factor domain", "Fungal TF domain")
    return s if len(s) <= LABEL_MAX else s[:LABEL_MAX - 1] + "…"


def _short_unique(labels):
    """label -> short label; labels whose short forms would collide keep the full label (no merged bars or rows)."""
    sh = {l: _short(l) for l in labels}
    dup = pd.Series(list(sh.values())).duplicated(keep=False).values
    return {l: (l if d else sh[l]) for l, d in zip(sh, dup)}


def _group_table(d, levels):
    rows = []
    for g in levels:
        x = d[d["grp"] == g]
        has = x["n_proteins_total"] > 0
        rows.append({"group": g, "n_genomes": int(len(x)), "n_with_proteome": int(has.sum()),
                     "n_with_tf": int((x["n_tf_proteins"] > 0).sum()),
                     "median_tf_fraction": float(x["tf_fraction"].median()) if x["tf_fraction"].notna().any() else np.nan,
                     "median_shannon": float(x["shannon"].median()) if x["shannon"].notna().any() else np.nan})
    return pd.DataFrame(rows)


def _box(d, col, ylabel, title, levels, colors, core):
    if d[col].notna().sum() == 0:
        return None
    fig, ax = plt.subplots(figsize=(width(levels), 3.8), layout="constrained")
    box_by_group(ax, d, col, levels, colors, ylabel=ylabel)
    ax.set_title(title, fontsize=9)
    parent_legend(ax, core, loc="upper left")
    if ax.get_legend() is not None:     # outside the axes: inside, it covers the boxes of the upper-right groups
        ax.get_legend().set_bbox_to_anchor((1.01, 1.0))
    return fig


def _composition(d, fam, levels, colors):
    f = fam[fam["species_prefix"].isin(set(d["species_prefix"]))]
    if f.empty:
        return None
    frac = per_genome_fraction(f, "family")          # rows sum to 1; genomes with no TF family count are absent
    top = list(frac.mean().sort_values(ascending=False).index[:TOP_FAMILIES])
    rest = [c for c in frac.columns if c not in top]
    out = frac[top].copy()
    if rest:
        out[OTHER_FAMILIES] = frac[rest].sum(axis=1)
    gm = group_mean(out, d, levels)
    cmap = {c: CAT[i] for i, c in enumerate(top)}
    if rest:
        cmap[OTHER_FAMILIES] = GREY
    fig, ax = plt.subplots(figsize=(width(levels) + 3.2, 4.2), layout="constrained")
    sm = _short_unique(list(gm.columns))
    stacked_by_group(ax, gm.rename(columns=sm), {sm[k]: v for k, v in cmap.items()}, levels,
                     "Share of TF family counts (genome mean)", colors)
    ax.set_title("TF family composition: share of TF protein counts per family, genomes with at least one TF protein",
                 fontsize=9)
    return fig


def family_z(d, fam, levels):
    """Families (rows) x groups (columns): group mean of n_proteins / n_proteins_total per genome, z-scored per family
    across groups. A genome with an empty proteome is left out of the means. A family with the same mean in every
    group is 0; a group with no genome with a proteome stays NaN. None when there is no family data."""
    f = fam[fam["species_prefix"].isin(set(d["species_prefix"]))]
    if f.empty:
        return None
    wide = f.pivot_table(index="species_prefix", columns="family", values="n_proteins", aggfunc="sum", fill_value=0)
    tot = d.set_index("species_prefix")["n_proteins_total"]
    wide = wide.reindex(tot.index).fillna(0)
    frac = wide.div(tot.where(tot > 0), axis=0)
    m = frac.groupby(d.set_index("species_prefix")["grp"]).mean().reindex(levels)
    m = m[m.mean().sort_values(ascending=False).index]
    sd = m.std(ddof=0)
    z = (m - m.mean()) / sd.replace(0, np.nan)
    return z.mask(m.notna() & sd.eq(0), 0.0).T      # & aligns sd (indexed by family) with the columns of m


def _family_heat(d, fam, levels, colors):
    z = family_z(d, fam, levels)
    if z is None:
        return None
    sm = _short_unique(list(z.index))
    z.index = [sm[i] for i in z.index]
    fig, ax = plt.subplots(figsize=(width(levels, 4.5), 1.6 + 0.22 * len(z)), layout="constrained")
    heat(ax, z, cmap="RdBu_r", vmin=-2.5, vmax=2.5, label="z-score across groups", colors=colors)
    ax.set_title("TF family proteins / all proteins, group mean, z-score per family", fontsize=9)
    return fig


def _rare_heat(d, fam, levels, colors):
    labels = rare_family_labels(RARE_ACCESSIONS)
    ok = d[d["n_proteins_total"] > 0]
    if ok.empty:
        return None
    f = fam[fam["species_prefix"].isin(set(ok["species_prefix"]))]
    has = f[f["n_proteins"] >= 1].groupby("family")["species_prefix"].agg(set)
    grp = ok.set_index("species_prefix")["grp"]
    rows = {}
    for lab in labels:
        ids = has.get(lab, set())
        rows[lab] = (grp.index.isin(ids)).astype(float)
    pres = pd.DataFrame(rows, index=grp.index).groupby(grp).mean().reindex(levels).T
    fig, ax = plt.subplots(figsize=(width(levels, 4.0), 1.8 + 0.3 * len(pres)), layout="constrained")
    heat(ax, pres, cmap="magma", vmin=0, vmax=1, label="Fraction of genomes with at least one protein", colors=colors)
    if len(levels) <= 16:
        for i in range(pres.shape[0]):
            for j in range(pres.shape[1]):
                v = pres.iat[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}".rstrip("0").rstrip(".") if v not in (0, 1) else f"{int(v)}", ha="center",
                            va="center", fontsize=6, color="white" if v < 0.6 else "black")
    ax.set_title("Rare TF families: presence by group (genomes with a proteome)", fontsize=9)
    return fig


def _scatter_all_genomes(o, genome_all, summ, name):
    """T06: every genome with a proteome and at least one TF protein, coloured by the site-wide phylum map."""
    d = genome_all[["species_prefix", "PHYLUM"]].merge(summ[["species_prefix", "n_proteins_total", "n_tf_proteins"]],
                                                        on="species_prefix", how="inner")
    d = d[(d["n_proteins_total"] > 0) & (d["n_tf_proteins"] > 0)]
    colors = phylum_color_map(phylum_order(genome_all))
    scope, info = o.scope, o.info
    o.scope, o.info = "all_genomes", {"n_genomes": None, "n_groups": None, "pooled_share": None}
    try:    # this figure is not a level figure: the sidecar must not carry the level's genome and group counts
        return scatter_all(d, "n_proteins_total", "n_tf_proteins", colors, o, name, logx=True, logy=True,
                           xlabel="Proteome size (proteins, log scale); genomes with no TF protein not shown",
                           ylabel="Proteins with a TF domain (log scale)", topic=TOPIC)
    finally:
        o.scope, o.info = scope, info


def page_tf_landscape(core, levels, colors, fam, summ, o, genome_all=None):
    """T01-T06 and the table T_tf_groups. core: representatives with 'grp' (assign_groups); fam: tf_family_genome;
    summ: tf_summary_genome; genome_all: every genome with PHYLUM (for T06; None skips T06)."""
    d = core[["species_prefix", "grp"]].merge(summ, on="species_prefix", how="inner")
    o.table(_group_table(d, levels), "T_tf_groups")
    o.stat_of["T01_tf_fraction"] = "tf_fraction"
    o.stat_of["T05_tf_shannon"] = "shannon"
    safe(o, FIGURES[0], TOPIC, lambda: _box(d, "tf_fraction", "Proteins with a TF domain / all proteins",
                                            "TF fraction: proteins with a TF domain / all proteins (genome values)",
                                            levels, colors, core))
    safe(o, FIGURES[1], TOPIC, lambda: _composition(d, fam, levels, colors))
    safe(o, FIGURES[2], TOPIC, lambda: _family_heat(d, fam, levels, colors))
    safe(o, FIGURES[3], TOPIC, lambda: _rare_heat(d, fam, levels, colors))
    safe(o, FIGURES[4], TOPIC, lambda: _box(d, "shannon", "Shannon index of TF family counts (nats)",
                                            "TF family diversity per genome (0 = one family; genomes with no TF omitted)",
                                            levels, colors, core))
    if genome_all is None:
        o.skip(FIGURES[5], TOPIC)
    else:
        _scatter_all_genomes(o, genome_all, summ, FIGURES[5])
