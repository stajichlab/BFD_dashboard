"""TF expansions and contractions (broader search than the paper). Every result here is a hypothesis.

Unit: one representative genome per genus (tf_data.genus_table). Families are the paper families (columns pf_<name>:
distinct proteins per genome, so no protein is counted twice). Families that BFD cannot detect (non-zero in fewer than
tf_claims.MIN_NONZERO_GENERA genera) are left out.

Method:
1. Scaling model per family: log1p(count) ~ log(n_proteins_total), ordinary least squares (numpy.polyfit). Default
   scaling="phylum": one fit per phylum, so an order is compared with the rest of its phylum on the same footing. A phylum
   with fewer than MIN_FIT_GENERA (30) eligible genera uses the global fit (all phyla) instead; the table column `scaling`
   and attrs['fallback_phyla'] record this. scaling="global" fits all genera together. The residual of a genus is its count
   relative to what its proteome size predicts. A family that is constant within a phylum has no residual there (NaN).
2. Lineage test per ORDER (and CLASS) with at least 5 genera and at least 5 other genera in the phylum: residuals of the
   lineage against the other genera of the SAME phylum (Mann-Whitney U; effect = rank-biserial, positive = lineage above the
   rest; direction = sign of the effect, empty for an effect of exactly 0). Benjamini-Hochberg is applied PER LEVEL (ORDER,
   CLASS) and PER FILTER (BUSCO>=80, none), over all family-by-lineage tests of that level and filter (the levels are nested
   and dependent, so pooling them adds little). Column bh_scope says so.
3. `flagged` = p_adj < 0.05 AND |effect| >= 0.5 AND at least 3 genera of the lineage lie on the side of the effect of the
   median of the rest of the phylum (n_agree) AND raw_agree (the mean raw count of the lineage lies on the side of the
   direction: above the rest for 'expanded', below for 'contracted'). raw_agree removes a log1p-floor artefact: for a rare
   family (zero in most genera) a lineage with a small proteome looks 'expanded' although its count is 0.
   `flagged_brief` is the same without raw_agree (the thresholds of the brief, direction against the rest of the phylum). The
   literal rule of the brief (n_agree counted against the sign of the LINEAGE median) is not reproducible: with a global fit
   a whole phylum can sit above or below 0, so that sign could contradict the effect.
4. `lineage_wide`: the lineage is flagged for MORE THAN HALF (LINEAGE_WIDE_FRACTION = 0.5, a judgment, not a tested value) of the
   families tested for it (and at least LINEAGE_WIDE_MIN_FAMILIES = 4 were tested). That pattern suggests a proteome-size or annotation effect, not a family-specific change.
5. Presence flags: a family is 'lost' in a lineage when it is present in >= 80% of all genera and in <= 20% of the lineage;
   'gained' for the reverse (global <= 20%, lineage >= 80%). A phylum with <= 2 lineages in the table gets the note
   'phylum-level statement'.

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
MIN_FIT_GENERA = 30          # a phylum with fewer eligible genera uses the global scaling fit
LINEAGE_WIDE_FRACTION = 0.5  # lineage_wide when flagged for MORE THAN this share of its families (a judgment)
LINEAGE_WIDE_MIN_FAMILIES = 4   # a share of fewer tested families means nothing (1 of 1 is 100%): never lineage-wide (a judgment)
SCALINGS = ("phylum", "global")
P_FLAG = 0.05
EFFECT_FLAG = 0.5
PREV_HIGH = 0.8
PREV_LOW = 0.2
MAX_BOX_HITS = 12
TEST_COLS = ["family", "lineage", "parent", "n_genera", "n_rest", "median_resid", "effect", "p", "p_adj", "direction",
             "n_agree", "flagged", "mean_count", "mean_count_rest", "raw_agree", "flagged_brief", "lineage_wide", "scaling"]
PRES_COLS = ["family", "lineage", "parent", "global_prevalence", "lineage_prevalence", "flag", "n_genera", "n_orders_in_phylum",
             "note"]
LEVEL_COL = {"ORDER": "ORDER", "CLASS": "CLASS"}
SENS_NOTE = "sensitivity only; low-BUSCO genomes have incomplete annotation"
SIGNAL_KIND = "n_flagged_cells"
BH_SCOPE = "per level and filter"
FLAG_DEF = ("flagged = p_adj < 0.05 and |effect| >= 0.5 and >= 3 genera on the effect side of the rest-of-phylum median and "
            "mean raw count on the side of the direction")
HYP = "hypothesis; genus values; no phylogenetic control; " + FLAG_DEF
WIDE_NOTE = "lineage-wide shift: likely a proteome-size or annotation effect, not family-specific"


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


def direction_of(effect):
    """'expanded' for an effect > 0, 'contracted' for < 0, '' for exactly 0 or NaN."""
    if not np.isfinite(effect) or effect == 0:
        return ""
    return "expanded" if effect > 0 else "contracted"


def _ols(x, y):
    """(slope, intercept) of y on x, or None when there are fewer than 3 points or x or y is constant."""
    if len(x) < 3 or np.ptp(y) == 0 or np.ptp(x) == 0:
        return None
    sl, ic = np.polyfit(x, y, 1)
    return float(sl), float(ic)


def fit_scaling(genus_df, family, scaling="phylum", min_fit=MIN_FIT_GENERA):
    """Copy of genus_df with `resid`: residual of log1p(count) regressed on log(n_proteins_total) over the genera with a count
    and a proteome size > 0. NaN where either is missing. scaling='phylum' (default): one fit per phylum; a phylum with fewer
    than min_fit eligible genera uses the global fit; a phylum whose family is constant has NaN residuals. scaling='global':
    one fit over all genera. The whole family is NaN (slope NaN) when the global fit is impossible (constant family or fewer
    than 3 genera). attrs: slope, intercept (global fit), n_fit, n_dropped_no_size (genera with no usable proteome size),
    scaling, fallback_phyla (sorted list), n_fallback_phyla, phylum_scaling ({phylum: 'phylum' | 'global (fallback: ...)'})."""
    if scaling not in SCALINGS:
        raise ValueError(f"scaling must be one of {SCALINGS}, not {scaling!r}")
    d = genus_df.copy()
    n = pd.to_numeric(d["n_proteins_total"], errors="coerce").to_numpy(dtype=float)
    cnt = pd.to_numeric(d[family], errors="coerce").to_numpy(dtype=float)
    size_ok = np.isfinite(n) & (n > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        x = np.where(size_ok, np.log(np.where(size_ok, n, 1.0)), np.nan)
        y = np.where(np.isfinite(cnt) & (cnt >= 0), np.log1p(np.where(np.isfinite(cnt), cnt, 0.0)), np.nan)
    ok = np.isfinite(x) & np.isfinite(y)
    resid = np.full(len(d), np.nan)
    gfit = _ols(x[ok], y[ok])
    fallback, how = [], {}
    if gfit is not None:
        if scaling == "global":
            resid[ok] = y[ok] - (gfit[1] + gfit[0] * x[ok])
        else:
            ph = d["PHYLUM"].to_numpy()
            for p in pd.unique(pd.Series(ph).dropna()):
                m = ok & (ph == p)
                if m.sum() < min_fit:
                    fallback.append(str(p))
                    how[str(p)] = f"global (fallback: fewer than {min_fit} genera)"
                    resid[m] = y[m] - (gfit[1] + gfit[0] * x[m])
                    continue
                how[str(p)] = "phylum"
                f = _ols(x[m], y[m])
                if f is not None:
                    resid[m] = y[m] - (f[1] + f[0] * x[m])
    if scaling == "global":
        how = {str(p): "global" for p in pd.unique(d["PHYLUM"].dropna())}
    d["resid"] = resid
    d.attrs.update({"slope": gfit[0] if gfit else np.nan, "intercept": gfit[1] if gfit else np.nan, "n_fit": int(ok.sum()),
                    "n_dropped_no_size": int((~size_ok).sum()), "scaling": scaling, "fallback_phyla": sorted(fallback),
                    "n_fallback_phyla": len(fallback), "phylum_scaling": how})
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


def lineage_tests(genus_df, families, level="ORDER", min_genera=MIN_GENERA, min_agree=MIN_AGREE, scaling="phylum"):
    """One row per family-by-lineage test (see the module docstring). families: table columns (pf_<name>); the family
    column of the result is the name without prefix. Lineages with fewer than min_genera genera (or a phylum rest with fewer)
    are not tested: they are in attrs['skipped'] with a note. Constant families give no rows and are named in
    attrs['skipped_families']. BH is over all tests of this call (one level, one table). Sorted by p_adj, family, lineage.
    attrs: skipped, skipped_families, n_tests, n_dropped_no_size, scaling, n_fallback_phyla, fallback_phyla."""
    col = LEVEL_COL[level]
    rows, skipped_fams = [], []
    n_drop, fb = 0, set()
    for fam in families:
        d = fit_scaling(genus_df, fam, scaling=scaling)
        n_drop = d.attrs["n_dropped_no_size"]
        fb |= set(d.attrs["fallback_phyla"])
        if not np.isfinite(d.attrs["slope"]):
            skipped_fams.append(f"{label(fam)}: constant or too few genera, no scaling fit")
            continue
        how = d.attrs["phylum_scaling"]
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
                side = np.sign(eff)
                agree = int((np.sign(a - np.median(b)) == side).sum()) if side != 0 else 0
                ca, cb = float(g[fam].astype(float).mean()), float(rest[fam].astype(float).mean())
                direction = direction_of(eff)
                raw = (ca > cb) if direction == "expanded" else (ca < cb) if direction == "contracted" else False
                rows.append((label(fam), lin, ph, len(a), len(b), med, eff, float(res.pvalue), np.nan, direction, agree, False,
                             ca, cb, bool(raw), False, False, how.get(str(ph), scaling)))
    t = pd.DataFrame(rows, columns=TEST_COLS)
    t["p_adj"] = bh_adjust(t["p"].to_numpy(dtype=float))
    t["flagged_brief"] = ((t["p_adj"] < P_FLAG) & (t["effect"].abs() >= EFFECT_FLAG) & (t["n_agree"] >= min_agree)).astype(bool)
    t["flagged"] = (t["flagged_brief"] & t["raw_agree"]).astype(bool)
    if len(t):
        g = t.groupby(["parent", "lineage"])["flagged"]
        t["lineage_wide"] = ((g.transform("mean") > LINEAGE_WIDE_FRACTION) & (g.transform("size") >= LINEAGE_WIDE_MIN_FAMILIES)
                             ).astype(bool)
    t["lineage_wide"] = t["lineage_wide"].astype(bool)
    t = t.sort_values(["p_adj", "family", "lineage", "parent"], kind="stable").reset_index(drop=True)
    t.attrs.update({"skipped": _skipped(genus_df, col, min_genera), "skipped_families": skipped_fams,
                    "n_tests": int(len(t)), "n_dropped_no_size": int(n_drop), "scaling": scaling,
                    "fallback_phyla": sorted(fb), "n_fallback_phyla": len(fb)})
    return t


def presence_flags(genus_df, families, level="ORDER", min_genera=MIN_GENERA):
    """Presence (count > 0) prevalence per family and lineage with at least min_genera genera, against the prevalence over all
    genera. flag: 'lost' (global >= 0.8, lineage <= 0.2), 'gained' (global <= 0.2, lineage >= 0.8), else ''. n_genera =
    genera of the lineage; n_orders_in_phylum = lineages (of this level) with at least one genus in that phylum in the table;
    note = 'phylum-level statement (the phylum has k order(s))' when k <= 2."""
    col = LEVEL_COL[level]
    word = "order(s)" if level == "ORDER" else "class(es)"
    n_lin = genus_df.dropna(subset=["PHYLUM", col]).groupby("PHYLUM")[col].nunique().to_dict()
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
            k = int(n_lin.get(ph, 0))
            note = f"phylum-level statement (the phylum has {k} {word})" if k <= 2 else ""
            rows.append((label(fam), lin, ph, glob, lp, flag, len(g), k, note))
    t = pd.DataFrame(rows, columns=PRES_COLS)
    return t.sort_values(["family", "parent", "lineage"], kind="stable").reset_index(drop=True)


# ---------------------------------------------------------------- figures
def _phylum_colors(porder, phyla):
    from make_figures import GREY, phylum_color_map
    order = list(porder) if porder is not None else sorted(phyla)
    cm = phylum_color_map(order)
    return {p: cm.get(p, GREY) for p in phyla}, order


def _heatmap(tests, porder):
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
    wide = set(zip(tests.loc[tests["lineage_wide"], "parent"], tests.loc[tests["lineage_wide"], "lineage"]))
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
    for j, c in enumerate(cols):
        if c in wide:      # lineage-wide shift: grey hatched column; the colours stay faintly visible below it
            ax.add_patch(Rectangle((j - 0.5, -0.5), 1, len(fams), facecolor="#d0d0d0", alpha=0.72, hatch="////", ec="#8a8a8a", lw=0))
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
                     label="flagged: p_adj < 0.05, |effect| >= 0.5, >= 3 genera on the effect side of the rest-of-phylum median, "
                           "mean raw count on the side of the direction"))
    hs.append(Patch(facecolor="#d0d0d0", hatch="////", edgecolor="#8a8a8a",
                    label="lineage-wide shift (flagged in > half of its families; likely proteome-size or annotation effect)"))
    hs.append(Patch(color="#bdbdbd", label="not tested (< 5 genera, or family constant in the phylum)"))
    fig.legend(handles=hs, fontsize=6, ncol=2, loc="outside lower center", frameon=False)
    ax0.set_title("Hypothesis: TF family expansions (red) and contractions (blue) by lineage. Residual of log1p(count) on log(proteome "
                  "size), genus values,\nlineage vs other genera of the same phylum, BH over all tests of the level. No phylogenetic "
                  "control; annotation gaps can look like contractions.", fontsize=8)
    return fig


def _hit_boxes(g80, tests, porder, max_hits=MAX_BOX_HITS, scaling="phylum"):
    """Up to max_hits flagged hits that are not lineage-wide (best p_adj first; lineage-wide hits are skipped and the next
    hits fill the panels). The jitter of the points uses a fixed seed (1), so the figure is the same on every run."""
    import matplotlib.pyplot as plt
    fl = tests[tests["flagged"] & ~tests["lineage_wide"]].copy()
    if fl.empty:
        return None
    fl["_abs"] = -fl["effect"].abs()
    fl = fl.sort_values(["p_adj", "_abs", "family", "lineage"], kind="stable").head(max_hits)
    colors, _ = _phylum_colors(porder, set(fl["parent"]))
    cols_of = {label(c): c for c in g80.columns if c.startswith(PF_PREFIX) or c.startswith(FAM_PREFIX)}
    n = len(fl)
    ncol = 4 if n > 6 else 3 if n > 3 else n
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.9 * ncol, 2.6 * nrow + 1.0), squeeze=False, layout="constrained")
    rng = np.random.default_rng(1)
    cache = {}
    for ax, (_, r) in zip(axes.ravel(), fl.iterrows()):
        if r["family"] not in cache:
            cache[r["family"]] = fit_scaling(g80, cols_of[r["family"]], scaling=scaling)
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
    short = f" Only {n} flagged hits are not lineage-wide." if n < max_hits else ""
    fig.suptitle("Hypothesis: top flagged TF family hits that are not lineage-wide (residuals, lineage vs the rest of its phylum; "
                 "genus values, no phylogenetic control)." + short + "\nflagged = p_adj < 0.05, |effect| >= 0.5, >= 3 genera on the "
                 "effect side of the rest-of-phylum median, mean raw count on the side of the direction.", fontsize=7)
    return fig


# ---------------------------------------------------------------- page
def _labelled(t, filter_label, note):
    t = t.copy()
    t.attrs = {}     # attrs hold DataFrames; pandas compares them (and fails) when it concatenates
    t.insert(0, "filter", filter_label)
    t["bh_scope"] = BH_SCOPE
    base = pd.Series(note, index=t.index)
    if "note" in t.columns:     # a note already on the rows (presence flags) goes first
        base = t["note"].where(t["note"].astype(str) == "", t["note"].astype(str) + "; ") + note
    if "lineage_wide" in t.columns:
        base = base.where(~t["lineage_wide"], base + "; " + WIDE_NOTE)
    t["note"] = base
    return t


def page_tf_expansion(genus_df_busco, genus_df_all, o, families=None, porder=None, levels=("ORDER", "CLASS"), scaling="phylum"):
    """Tables T11_tf_expansion_hits (genus_df_busco: BUSCO >= 80; every tested family-by-lineage row, flagged column),
    T11_tf_expansion_nofilter (genus_df_all, no BUSCO filter: sensitivity only), T11_tf_expansion_skipped (lineages not
    tested), T11_tf_presence_flags (flagged rows only), and figures T11_tf_effect_heatmap, T12_tf_hit_boxes. o: an Out for 'tf'.
    families: table columns; default = comparable pf_ columns of genus_df_all. levels: lineage levels tested (rows carry a
    `level` column; the heat map and boxes use ORDER rows). BH is per level and per filter (column bh_scope). The sidecar signal
    of the heat map is the number of flagged ORDER cells, of the boxes the number of panels (signal_kind n_flagged_cells; only
    comparable with itself). Returns the BUSCO>=80 test table."""
    from figs_common import safe
    fams = list(families) if families is not None else comparable_families(genus_df_all)
    o.scope = "genus_reps"
    main, nof, skip, pres = [], [], [], []
    for lv in levels:
        t = lineage_tests(genus_df_busco, fams, level=lv, scaling=scaling)
        t.insert(0, "level", lv)
        s = t.attrs["skipped"].copy()
        s.insert(0, "level", lv)
        skip.append(s)
        main.append(_labelled(t, "BUSCO>=80", HYP))
        n = lineage_tests(genus_df_all, fams, level=lv, scaling=scaling)
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
    n_cells = int(order_tests["flagged"].sum())
    n_panels = int((order_tests["flagged"] & ~order_tests["lineage_wide"]).sum())
    # The heat map has no single effect size; its signal is the number of flagged cells (the box figure: its panels), used only
    # as the second sort key of figures_index.tsv. A figure with none is left unscored.
    if n_cells:
        o.meta["T11_tf_effect_heatmap"] = {"signal": n_cells, "signal_kind": SIGNAL_KIND}
    if n_panels:
        o.meta["T12_tf_hit_boxes"] = {"signal": min(n_panels, MAX_BOX_HITS), "signal_kind": SIGNAL_KIND}
    safe(o, "T11_tf_effect_heatmap", "tf", lambda: _heatmap(order_tests, porder))
    safe(o, "T12_tf_hit_boxes", "tf", lambda: _hit_boxes(genus_df_busco, order_tests, porder, scaling=scaling))
    return hits
