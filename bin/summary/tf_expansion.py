"""TF expansions and contractions (broader search than the paper). Every result here is a hypothesis.

Unit: one representative genome per genus (tf_data.genus_table). Families are the paper families (columns pf_<name>:
distinct proteins per genome, so no protein is counted twice). Families that BFD cannot detect (non-zero in fewer than
tf_claims.MIN_NONZERO_GENERA genera) are left out.

Method:
1. Scaling model per family: log1p(count) ~ log(n_proteins_total), ordinary least squares (numpy.polyfit) over all genera
   with a proteome size. The residual of a genus is its count relative to what its proteome size predicts.
2. Lineage test per ORDER (and CLASS) with at least 5 genera: residuals of the lineage against the other genera of the
   SAME phylum (Mann-Whitney U; effect = rank-biserial, positive = lineage above the rest; direction = sign of the effect).
   Benjamini-Hochberg over ALL family-by-lineage tests of the run. A hit is flagged when p_adj < 0.05, |effect| >= 0.5 and at
   least 3 genera of the lineage lie on the side of the effect of the median of the rest of the phylum (n_agree). The brief
   said 'sign of the lineage median'; the fit is global, so a whole phylum can sit above 0 and that rule gave a direction
   opposite to the effect.
   Extra guard (raw_agree): the mean raw count of the lineage must lie on the side of the direction (above the rest for
   'expanded', below for 'contracted'). Without it, a rare family (zero in most genera) looks 'expanded' in a lineage with a
   small proteome, because the residual of a zero rises as the fitted value falls (log1p floor). `flagged_brief` is the flag
   without this guard.
3. Presence flags: a family is 'lost' in a lineage when it is present in >= 80% of all genera and in <= 20% of the lineage;
   'gained' for the reverse (global <= 20%, lineage >= 80%).

Limits: genera of one lineage are related, and no phylogenetic correction is made, so p-values overstate the evidence.
Annotation gaps (fragmented or low-BUSCO assemblies, missed Pfam hits) look like contractions. The sensitivity table uses
genera of every BUSCO level and is NOT evidence: low-BUSCO genomes have incomplete annotation.
"""
import numpy as np
import pandas as pd
from scipy import stats as st

from stats import bh_adjust
from tf_claims import MIN_NONZERO_GENERA

PF_PREFIX = "pf_"
FAM_PREFIX = "fam_"
MIN_GENERA = 5
MIN_AGREE = 3
P_FLAG = 0.05
EFFECT_FLAG = 0.5
PREV_HIGH = 0.8
PREV_LOW = 0.2
MAX_BOX_HITS = 12
TEST_COLS = ["family", "lineage", "parent", "n_genera", "n_rest", "median_resid", "effect", "p", "p_adj", "direction",
             "n_agree", "flagged", "mean_count", "mean_count_rest", "raw_agree", "flagged_brief"]
PRES_COLS = ["family", "lineage", "parent", "global_prevalence", "lineage_prevalence", "flag"]
LEVEL_COL = {"ORDER": "ORDER", "CLASS": "CLASS"}
HYP = "hypothesis; genus values; no phylogenetic control"
SENS_NOTE = "sensitivity only; low-BUSCO genomes have incomplete annotation"
SIGNAL_KIND = "rank_biserial_max_flagged"


def label(col):
    """Family name of a table column: the pf_ / fam_ prefix is removed."""
    for p in (PF_PREFIX, FAM_PREFIX):
        if str(col).startswith(p):
            return str(col)[len(p):]
    return str(col)


def comparable_families(genus_df, min_nonzero=MIN_NONZERO_GENERA):
    """pf_ columns of the table with a non-zero count in at least min_nonzero genera (the rule of the claims page: a
    family that BFD does not detect is not comparable). Sorted by name."""
    out = []
    for c in sorted(genus_df.columns):
        if c.startswith(PF_PREFIX) and int((genus_df[c].astype(float) > 0).sum()) >= min_nonzero:
            out.append(c)
    return out


def fit_scaling(genus_df, family):
    """Copy of genus_df with `resid`: residual of log1p(count) regressed on log(n_proteins_total), fitted over the genera
    with a count and a proteome size > 0. NaN where either is missing, and NaN everywhere (slope NaN) when the family is
    constant or fewer than 3 genera can be fitted. attrs: slope, intercept, n_fit, n_dropped_no_size (genera with no usable
    proteome size)."""
    d = genus_df.copy()
    n = pd.to_numeric(d["n_proteins_total"], errors="coerce").to_numpy(dtype=float)
    cnt = pd.to_numeric(d[family], errors="coerce").to_numpy(dtype=float)
    size_ok = np.isfinite(n) & (n > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        x = np.where(size_ok, np.log(np.where(size_ok, n, 1.0)), np.nan)
        y = np.where(np.isfinite(cnt) & (cnt >= 0), np.log1p(np.where(np.isfinite(cnt), cnt, 0.0)), np.nan)
    ok = np.isfinite(x) & np.isfinite(y)
    resid = np.full(len(d), np.nan)
    slope = intercept = np.nan
    if ok.sum() >= 3 and np.ptp(y[ok]) > 0 and np.ptp(x[ok]) > 0:
        slope, intercept = np.polyfit(x[ok], y[ok], 1)
        resid[ok] = y[ok] - (intercept + slope * x[ok])
    d["resid"] = resid
    d.attrs.update({"slope": float(slope), "intercept": float(intercept), "n_fit": int(ok.sum()),
                    "n_dropped_no_size": int((~size_ok).sum())})
    return d


def _skipped(genus_df, col, min_genera):
    """Lineages that are not tested, once per lineage (not per family): too few genera, or a phylum rest that is too small."""
    d = genus_df[[col, "PHYLUM", "n_proteins_total"]].copy()
    d["n_proteins_total"] = pd.to_numeric(d["n_proteins_total"], errors="coerce")
    d = d[d[col].notna() & d["PHYLUM"].notna() & (d["n_proteins_total"] > 0)]
    rows = []
    for (ph, lin), g in d.groupby(["PHYLUM", col]):
        n = len(g)
        n_rest = int((d["PHYLUM"] == ph).sum() - n)
        if n < min_genera:
            rows.append((lin, ph, n, f"not tested: fewer than {min_genera} genera"))
        elif n_rest < min_genera:
            rows.append((lin, ph, n, f"not tested: fewer than {min_genera} other genera in the phylum"))
    return pd.DataFrame(rows, columns=["lineage", "parent", "n_genera", "note"])


def lineage_tests(genus_df, families, level="ORDER", min_genera=MIN_GENERA, min_agree=MIN_AGREE):
    """One row per family-by-lineage test (see the module docstring). families: table columns (pf_<name>); the family
    column of the result is the name without prefix. Lineages with fewer than min_genera genera (or a phylum rest with fewer)
    are not tested: they are in attrs['skipped'] with a note. Constant families give no rows and are named in
    attrs['skipped_families']. Sorted by p_adj, family, lineage. attrs: n_tests, n_dropped_no_size."""
    col = LEVEL_COL[level]
    rows, skipped_fams = [], []
    n_drop = 0
    for fam in families:
        d = fit_scaling(genus_df, fam)
        n_drop = d.attrs["n_dropped_no_size"]
        if not np.isfinite(d.attrs["slope"]):
            skipped_fams.append(f"{label(fam)}: constant or too few genera, no scaling fit")
            continue
        d = d[d["resid"].notna() & d["PHYLUM"].notna()]
        for ph, gp in d.groupby("PHYLUM"):
            for lin, g in gp[gp[col].notna()].groupby(col):
                rest = gp[gp[col] != lin]
                rest = pd.concat([rest, gp[gp[col].isna()]]).drop_duplicates()   # genera with no lineage stay in the rest
                a, b = g["resid"].to_numpy(), rest["resid"].to_numpy()
                if len(a) < min_genera or len(b) < min_genera:
                    continue
                res = st.mannwhitneyu(a, b, alternative="two-sided")
                eff = 2 * float(res.statistic) / (len(a) * len(b)) - 1
                med = float(np.median(a))
                # direction = sign of the effect; agreement = genera on that side of the median of the rest of the phylum
                # (the fit is global, so the residuals of a whole phylum can sit above or below 0)
                side = np.sign(eff)
                agree = int((np.sign(a - np.median(b)) == side).sum()) if side != 0 else 0
                ca, cb = float(g[fam].astype(float).mean()), float(rest[fam].astype(float).mean())
                raw = bool(ca > cb) if eff > 0 else bool(ca < cb)     # the raw mean must not contradict the direction
                rows.append((label(fam), lin, ph, len(a), len(b), med, eff, float(res.pvalue), np.nan,
                             "expanded" if eff > 0 else "contracted", agree, False, ca, cb, raw, False))
    t = pd.DataFrame(rows, columns=TEST_COLS)
    t["p_adj"] = bh_adjust(t["p"].to_numpy(dtype=float))
    t["flagged_brief"] = ((t["p_adj"] < P_FLAG) & (t["effect"].abs() >= EFFECT_FLAG) & (t["n_agree"] >= min_agree)).astype(bool)
    t["flagged"] = (t["flagged_brief"] & t["raw_agree"]).astype(bool)
    t = t.sort_values(["p_adj", "family", "lineage", "parent"], kind="stable").reset_index(drop=True)
    t.attrs.update({"skipped": _skipped(genus_df, col, min_genera), "skipped_families": skipped_fams,
                    "n_tests": int(len(t)), "n_dropped_no_size": int(n_drop)})
    return t


def presence_flags(genus_df, families, level="ORDER", min_genera=MIN_GENERA):
    """Presence (count > 0) prevalence per family and lineage with at least min_genera genera, against the prevalence over all
    genera. flag: 'lost' (global >= 0.8, lineage <= 0.2), 'gained' (global <= 0.2, lineage >= 0.8), else ''."""
    col = LEVEL_COL[level]
    rows = []
    for fam in families:
        cnt = pd.to_numeric(genus_df[fam], errors="coerce")
        pres = (cnt > 0).where(cnt.notna())
        glob = float(pres.mean()) if pres.notna().any() else np.nan
        d = pd.DataFrame({"lin": genus_df[col], "ph": genus_df["PHYLUM"], "p": pres}).dropna()
        for (ph, lin), g in d.groupby(["ph", "lin"]):
            if len(g) < min_genera:
                continue
            lp = float(g["p"].mean())
            flag = "lost" if (glob >= PREV_HIGH and lp <= PREV_LOW) else \
                   "gained" if (glob <= PREV_LOW and lp >= PREV_HIGH) else ""
            rows.append((label(fam), lin, ph, glob, lp, flag))
    t = pd.DataFrame(rows, columns=PRES_COLS)
    return t.sort_values(["family", "parent", "lineage"], kind="stable").reset_index(drop=True)


# ---------------------------------------------------------------- figures
def _phylum_colors(porder, phyla):
    from make_figures import GREY, phylum_color_map
    order = list(porder) if porder is not None else sorted(phyla)
    cm = phylum_color_map(order)
    return {p: cm.get(p, GREY) for p in phyla}, order


def _heatmap(tests, porder, n_boot_note=""):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch, Rectangle
    if tests.empty:
        return None
    colors, order = _phylum_colors(porder, set(tests["parent"]))
    rank = {p: i for i, p in enumerate(order)}
    lin = tests[["parent", "lineage"]].drop_duplicates()
    lin = lin.assign(r=lin["parent"].map(lambda p: rank.get(p, len(rank)))).sort_values(["r", "parent", "lineage"])
    cols = list(zip(lin["parent"], lin["lineage"]))
    fams = sorted(tests["family"].unique())
    eff = tests.pivot_table(index="family", columns=["parent", "lineage"], values="effect", aggfunc="first")
    eff = eff.reindex(index=fams, columns=pd.MultiIndex.from_tuples(cols))
    flg = tests.pivot_table(index="family", columns=["parent", "lineage"], values="flagged", aggfunc="first")
    flg = flg.reindex(index=fams, columns=pd.MultiIndex.from_tuples(cols)).fillna(False).astype(bool)
    w = max(7.5, 0.26 * len(cols) + 3.2)
    h = 2.6 + 0.26 * len(fams)
    fig, (ax0, ax) = plt.subplots(2, 1, figsize=(w, h), gridspec_kw={"height_ratios": [0.35, h - 1.4]}, layout="constrained")
    cmap = mpl.colormaps["RdBu_r"].with_extremes(bad="#bdbdbd")
    im = ax.imshow(eff.to_numpy(dtype=float), aspect="auto", cmap=cmap, vmin=-1, vmax=1)
    ax.set_yticks(range(len(fams)))
    ax.set_yticklabels(fams, fontsize=7)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c[1] for c in cols], rotation=70, ha="right", fontsize=6)
    for tl, (ph, _) in zip(ax.get_xticklabels(), cols):
        tl.set_color(tuple(v * 0.6 for v in mpl.colors.to_rgb(colors[ph])))
    for i in range(len(fams)):
        for j in range(len(cols)):
            if flg.iat[i, j]:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="black", lw=1.6))
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("Effect (rank-biserial of scaling residuals; positive = more than expected)", fontsize=7)
    strip = np.array([[mpl.colors.to_rgb(colors[ph]) for ph, _ in cols]])
    ax0.imshow(strip, aspect="auto")
    ax0.set_axis_off()
    hs = [Patch(color=colors[p], label=p) for p in order if p in colors]
    hs.append(Line2D([], [], marker="s", mfc="none", mec="black", ls="", mew=1.6,
                     label="flagged: p_adj < 0.05, |effect| >= 0.5, >= 3 genera agree"))
    hs.append(Patch(color="#bdbdbd", label="not tested (< 5 genera)"))
    ax.legend(handles=hs, fontsize=6, ncol=min(4, len(hs)), loc="upper center", bbox_to_anchor=(0.5, -0.42 if len(cols) < 25 else -0.3),
              frameon=False)
    ax0.set_title("Hypothesis: TF family expansions (red) and contractions (blue) by lineage. Residual of log1p(count) on log(proteome "
                  "size), genus values,\nlineage vs other genera of the same phylum, BH over all tests. No phylogenetic control; "
                  "annotation gaps can look like contractions.", fontsize=8)
    return fig


def _hit_boxes(g80, tests, porder, max_hits=MAX_BOX_HITS):
    import matplotlib.pyplot as plt
    fl = tests[tests["flagged"]].copy()
    if fl.empty:
        return None
    fl["_abs"] = -fl["effect"].abs()
    fl = fl.sort_values(["p_adj", "_abs", "family", "lineage"], kind="stable").head(max_hits)
    colors, _ = _phylum_colors(porder, set(fl["parent"]))
    cols_of = {label(c): c for c in g80.columns if c.startswith(PF_PREFIX) or c.startswith(FAM_PREFIX)}
    n = len(fl)
    ncol = 4 if n > 6 else 3 if n > 3 else n
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.9 * ncol, 2.6 * nrow + 0.8), squeeze=False, layout="constrained")
    rng = np.random.default_rng(1)
    cache = {}
    for ax, (_, r) in zip(axes.ravel(), fl.iterrows()):
        if r["family"] not in cache:
            cache[r["family"]] = fit_scaling(g80, cols_of[r["family"]])
        d = cache[r["family"]]
        d = d[d["resid"].notna() & (d["PHYLUM"] == r["parent"])]
        col = LEVEL_COL.get(r.get("level", "ORDER"), "ORDER")
        a = d[d[col] == r["lineage"]]["resid"].to_numpy()
        b = d[d[col] != r["lineage"]]["resid"].to_numpy()
        c = colors[r["parent"]]
        bp = ax.boxplot([a, b], patch_artist=True, widths=0.6, showfliers=False, medianprops={"color": "black"})
        for box, fc in zip(bp["boxes"], (c, "#bdbdbd")):
            box.set(facecolor=fc, alpha=0.45, edgecolor=fc)
        for i, v in ((1, a), (2, b)):
            ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=6, color=c if i == 1 else "#7f7f7f", alpha=0.6, linewidths=0)
        ax.axhline(0, color="black", lw=0.5, ls=":")
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"{r['lineage']}\n(n={len(a)})", f"rest of {r['parent']}\n(n={len(b)})"], fontsize=6)
        ax.tick_params(axis="y", labelsize=6)
        ax.set_title(f"{r['family']}\n{r['direction']}, effect {r['effect']:.2f}, p_adj {r['p_adj']:.1e}", fontsize=6.5)
    for ax in axes.ravel()[n:]:
        ax.set_axis_off()
    fig.supylabel("Scaling residual (log1p count vs log proteome size)", fontsize=7)
    fig.suptitle("Hypothesis: top flagged TF family hits (residuals, lineage vs the rest of its phylum; genus values, no "
                 "phylogenetic control)", fontsize=8)
    return fig


# ---------------------------------------------------------------- page
def _labelled(t, filter_label, note):
    t = t.copy()
    t.attrs = {}     # attrs hold DataFrames; pandas compares them (and fails) when it concatenates
    t.insert(0, "filter", filter_label)
    t["note"] = note
    return t


def page_tf_expansion(genus_df_busco, genus_df_all, o, families=None, porder=None, levels=("ORDER", "CLASS")):
    """Tables T11_tf_expansion_hits (genus_df_busco: BUSCO >= 80; every tested family-by-lineage row, flagged column),
    T11_tf_expansion_nofilter (genus_df_all, no BUSCO filter: sensitivity only), T11_tf_expansion_skipped (lineages not
    tested), T11_tf_presence_flags (flagged rows only), and figures T11_tf_effect_heatmap, T12_tf_hit_boxes. o: an Out for 'tf'.
    families: table columns; default = comparable pf_ columns of genus_df_all. levels: lineage levels tested (rows carry a
    `level` column; the heat map and boxes use ORDER rows). Returns the BUSCO>=80 test table."""
    from figs_common import safe
    fams = list(families) if families is not None else comparable_families(genus_df_all)
    o.scope = "genus_reps"
    main, nof, skip, pres = [], [], [], []
    for lv in levels:
        t = lineage_tests(genus_df_busco, fams, level=lv)
        t.insert(0, "level", lv)
        s = t.attrs["skipped"].copy()
        s.insert(0, "level", lv)
        skip.append(s)
        main.append(_labelled(t, "BUSCO>=80", HYP))
        n = lineage_tests(genus_df_all, fams, level=lv)
        n.insert(0, "level", lv)
        nof.append(_labelled(n, "none", f"{HYP}; {SENS_NOTE}"))
        p = presence_flags(genus_df_busco, fams, level=lv)
        p = p[p["flag"] != ""]
        p.insert(0, "level", lv)
        pres.append(_labelled(p, "BUSCO>=80", HYP))
    srt = ["p_adj", "family", "lineage"]
    hits = pd.concat(main, ignore_index=True).sort_values(["level"] + srt, kind="stable").reset_index(drop=True)
    sens = pd.concat(nof, ignore_index=True).sort_values(["level"] + srt, kind="stable").reset_index(drop=True)
    o.table(hits, "T11_tf_expansion_hits")
    o.table(sens, "T11_tf_expansion_nofilter")
    o.table(pd.concat(skip, ignore_index=True), "T11_tf_expansion_skipped")
    o.table(pd.concat(pres, ignore_index=True), "T11_tf_presence_flags")
    order_tests = hits[hits["level"] == "ORDER"] if "ORDER" in levels else hits.iloc[0:0]
    o.info = {"n_genomes": int(len(genus_df_busco)), "n_groups": int(order_tests["lineage"].nunique()), "pooled_share": 0.0}
    fl = order_tests[order_tests["flagged"]]
    # The heat map has no single effect size; its signal is the largest |effect| among flagged cells (rank-biserial), used
    # only as the second sort key of figures_index.tsv. A run with no flagged cell is left unscored.
    if len(fl):
        sig = {"signal": float(fl["effect"].abs().max()), "signal_kind": SIGNAL_KIND}
        o.meta["T11_tf_effect_heatmap"] = dict(sig)
        o.meta["T12_tf_hit_boxes"] = dict(sig)
    safe(o, "T11_tf_effect_heatmap", "tf", lambda: _heatmap(order_tests, porder))
    safe(o, "T12_tf_hit_boxes", "tf", lambda: _hit_boxes(genus_df_busco, order_tests, porder))
    return hits
