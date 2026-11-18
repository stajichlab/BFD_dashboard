"""The paper's lineage claims (Lax et al. 2025) tested on BFD genus values, and the forest figure (T10).

Unit: one value per genus (genus representatives from tf_data.genus_table). Claims come from config/tf_claims.tsv. Each claim
compares group a with group b on one metric of one paper family:
  fraction    family proteins / all proteins of the genome   (rank-biserial effect, Mann-Whitney p)
  count       family proteins                                (same)
  shannon     Shannon index of the genome's TF family counts (same; family is empty)
  prevalence  share of genera with at least one protein      (effect = difference in prevalence, Fisher exact p)
The BFD family labels are bridged to the paper family through config/tf_family_map.tsv (several BFD families can map to one
paper family; their protein counts are added, so a protein with two domains of one paper family counts twice).
Groups are phyla or the sets Dikarya (Ascomycota + Basidiomycota), EDF (all phyla except Ascomycota, Basidiomycota and
Microsporidia), ALL (every phylum); 'X-Y' removes the phylum Y from the set X (for example ALL-Mucoromycota).
Confidence intervals are percentile bootstraps over genera (n_boot resamples within each group, fixed seed). p-values are
Benjamini-Hochberg adjusted over all tested claims. Phylogeny is not controlled: related genera are not independent.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

from stats import bh_adjust

CONFIG = Path(__file__).resolve().parents[2] / "config"
CLAIMS_FILE = CONFIG / "tf_claims.tsv"
FAMILY_MAP = CONFIG / "tf_family_map.tsv"
MIN_GENERA = 5
DIKARYA = ("Ascomycota", "Basidiomycota")
NOT_EDF = DIKARYA + ("Microsporidia",)
FAM_PREFIX = "fam_"
UNIT = {"genus": "genera", "species": "species"}
CLAIM_COLS = ["claim_id", "text", "family", "metric", "level", "group_a", "group_b", "paper_direction"]
OUT_COLS = ["claim_id", "text", "family", "metric", "level", "group_a", "group_b", "n_genera_a", "n_genera_b", "value_a",
            "value_b", "effect", "ci_lo", "ci_hi", "p", "p_adj", "observed_direction", "paper_direction", "agrees", "note",
            "source"]
AGREE_COLOR = {"yes": "#009E73", "mixed": "#E69F00", "no": "#D55E00", "not tested": "#8c8c8c"}   # Okabe-Ito


def load_claims(path=None):
    """config/tf_claims.tsv (lines starting with # are comments) with the columns CLAIM_COLS."""
    c = pd.read_csv(path or CLAIMS_FILE, sep="\t", dtype=str, comment="#", keep_default_na=False)
    missing = [k for k in CLAIM_COLS if k not in c.columns]
    if missing:
        raise ValueError(f"claims file lacks columns {missing}")
    return c[CLAIM_COLS]


def load_fmap(path=None):
    return pd.read_csv(path or FAMILY_MAP, sep="\t", dtype=str, comment="#", keep_default_na=False)


def resolve_group(spec, phylum):
    """Boolean numpy mask of the rows of `phylum` (Series) in the group `spec`. Unknown phylum (NaN) is in no group."""
    ph = pd.Series(phylum).reset_index(drop=True)
    known = ph.notna()
    base, *minus = [s.strip() for s in str(spec).split("-")]
    if base == "Dikarya":
        m = ph.isin(DIKARYA)
    elif base == "EDF":
        m = known & ~ph.isin(NOT_EDF)
    elif base == "ALL":
        m = known
    else:
        m = ph == base
    for s in minus:
        m = m & (ph != s)
    return m.to_numpy(dtype=bool)


def _ranks_u(a, b):
    """Rows: resamples. U of a = (pairs with a > b) + 0.5 * ties, from average ranks. a: (k, na), b: (k, nb)."""
    na = a.shape[1]
    r = st.rankdata(np.concatenate([a, b], axis=1), axis=1)
    return r[:, :na].sum(axis=1) - na * (na + 1) / 2


def rank_biserial(a, b):
    """2 U_a / (n_a n_b) - 1: +1 when every a is above every b, -1 when every a is below, 0 for no difference."""
    a = np.asarray(a, dtype=float)[None, :]
    b = np.asarray(b, dtype=float)[None, :]
    return float(2 * _ranks_u(a, b)[0] / (a.shape[1] * b.shape[1]) - 1)


def observed_direction(effect):
    if not np.isfinite(effect):
        return ""
    return "a>b" if effect > 0 else "a<b" if effect < 0 else "none"


def agreement(effect, lo, hi, direction):
    """yes: the sign matches the paper's direction and the CI excludes 0; mixed: sign matches, CI includes 0 (or effect 0);
    no: the sign is opposite. absent_in_a expects a below b. NaN effect: 'not tested'."""
    if not (np.isfinite(effect) and np.isfinite(lo) and np.isfinite(hi)):
        return "not tested"
    want = 1 if direction == "a>b" else -1
    if effect == 0:
        return "mixed"
    if np.sign(effect) != want:
        return "no"
    return "yes" if (lo > 0 or hi < 0) else "mixed"


def _boot_ci(x, y, binary, n_boot, rng):
    """Percentile 95% CI of the effect over resamples of the genera of each group (with replacement)."""
    ia = rng.integers(0, len(x), (n_boot, len(x)))
    ib = rng.integers(0, len(y), (n_boot, len(y)))
    a, b = x[ia], y[ib]
    if binary:
        e = a.mean(axis=1) - b.mean(axis=1)
    else:
        e = 2 * _ranks_u(a, b) / (len(x) * len(y)) - 1
    lo, hi = np.percentile(e, [2.5, 97.5])
    return float(lo), float(hi)


def _empty(claim, source, note):
    r = {k: claim[k] for k in CLAIM_COLS}
    r.update(n_genera_a=0, n_genera_b=0, value_a=np.nan, value_b=np.nan, effect=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan,
             p_adj=np.nan, observed_direction="", agrees="not tested", note=note, source=source)
    return r


def evaluate_table(df, claims, value_fn, n_boot=2000, seed=1, source="BFD", level="genus"):
    """Core of the claim tests. df: one row per unit (genus or species) with a PHYLUM column. value_fn(family, metric) ->
    (values aligned to df.index, note): metric is 'fraction', 'count' or 'shannon'; a prevalence claim uses 'count' and a value
    of at least 1. values None = not available (the note says why). Returns OUT_COLS, p_adj by Benjamini-Hochberg."""
    rows = []
    for i, c in claims.reset_index(drop=True).iterrows():
        c = c.to_dict()
        c["level"] = level
        metric = c["metric"]
        vals, note = value_fn(c["family"], "count" if metric == "prevalence" else metric)
        ma, mb = resolve_group(c["group_a"], df["PHYLUM"]), resolve_group(c["group_b"], df["PHYLUM"])
        if vals is None:
            r = _empty(c, source, note)
            r["n_genera_a"], r["n_genera_b"] = int(ma.sum()), int(mb.sum())
            rows.append(r)
            continue
        v = np.asarray(vals, dtype=float)
        x, y = v[ma][np.isfinite(v[ma])], v[mb][np.isfinite(v[mb])]
        r = _empty(c, source, note or "")
        r["n_genera_a"], r["n_genera_b"] = len(x), len(y)
        small = [g for g, a in (("a", x), ("b", y)) if len(a) < MIN_GENERA]
        if small:
            r["note"] = (r["note"] + "; " if r["note"] else "") + \
                f"fewer than {MIN_GENERA} {UNIT.get(level, level)} in group {' and '.join(small)}"
            rows.append(r)
            continue
        rng = np.random.default_rng([seed, i])
        if metric == "prevalence":
            xb, yb = (x >= 1).astype(float), (y >= 1).astype(float)
            ka, kb = int(xb.sum()), int(yb.sum())
            r["value_a"], r["value_b"] = ka / len(x), kb / len(y)
            r["effect"] = r["value_a"] - r["value_b"]
            r["p"] = float(st.fisher_exact([[ka, len(x) - ka], [kb, len(y) - kb]])[1])
            r["ci_lo"], r["ci_hi"] = _boot_ci(xb, yb, True, n_boot, rng)
        else:
            r["value_a"], r["value_b"] = float(np.median(x)), float(np.median(y))
            r["effect"] = rank_biserial(x, y)
            r["p"] = float(st.mannwhitneyu(x, y, alternative="two-sided").pvalue) if np.ptp(np.r_[x, y]) > 0 else 1.0
            r["ci_lo"], r["ci_hi"] = _boot_ci(x, y, False, n_boot, rng)
        r["observed_direction"] = observed_direction(r["effect"])
        r["agrees"] = agreement(r["effect"], r["ci_lo"], r["ci_hi"], c["paper_direction"])
        rows.append(r)
    out = pd.DataFrame(rows, columns=OUT_COLS)
    out["p_adj"] = bh_adjust(out["p"].to_numpy(dtype=float))
    return out


def evaluate_claims(genus_df, claims, n_boot=2000, seed=1, fmap=None):
    """Test each claim on genus values. genus_df: tf_data.genus_table output (PHYLUM, n_proteins_total, shannon, fam_<label>).
    fmap: the family map (default config/tf_family_map.tsv). Returns OUT_COLS (source = 'BFD')."""
    fm = load_fmap() if fmap is None else fmap
    bridge = fm[(fm["paper_family"] != "") & (fm["family"] != "")][["family", "paper_family"]].drop_duplicates()
    labels = bridge.groupby("paper_family")["family"].agg(list).to_dict()
    tot = genus_df["n_proteins_total"].astype(float).where(genus_df["n_proteins_total"] > 0)

    def value_fn(family, metric):
        if metric == "shannon":
            return (genus_df["shannon"].to_numpy(dtype=float), "") if "shannon" in genus_df else (None, "no shannon column")
        if family not in labels:
            return None, f"family '{family}' is not a paper family in the family map"
        cols = [FAM_PREFIX + f for f in labels[family]]
        have = [c for c in cols if c in genus_df.columns]
        if not have:
            return None, f"family column {', '.join(cols)} missing from the genus table (no genome has this family)"
        cnt = genus_df[have].astype(float).sum(axis=1)
        miss = [c for c in cols if c not in have]
        note = f"column {', '.join(miss)} missing: counted as 0" if miss else ""
        return ((cnt / tot).to_numpy() if metric == "fraction" else cnt.to_numpy()), note

    return evaluate_table(genus_df, claims, value_fn, n_boot, seed, "BFD", "genus")


def _short(g):
    return {"Dikarya": "Dikarya", "EDF": "EDF"}.get(g, g)


def _forest(t):
    import matplotlib.pyplot as plt
    t = t.reset_index(drop=True)
    ids = list(dict.fromkeys(t["claim_id"]))
    n = len(ids)
    fig, ax = plt.subplots(figsize=(8.6, 1.6 + 0.42 * n), layout="constrained")
    ypos = {c: n - 1 - k for k, c in enumerate(ids)}
    off = {"BFD": 0.14, "paper": -0.14}
    mark = {"BFD": "o", "paper": "s"}
    for _, r in t.iterrows():
        y = ypos[r["claim_id"]] + off.get(r["source"], 0)
        col = AGREE_COLOR.get(r["agrees"], "#8c8c8c")
        if np.isfinite(r["effect"]):
            ax.plot([r["ci_lo"], r["ci_hi"]], [y, y], color=col, lw=1.6, solid_capstyle="butt")
            ax.plot([r["effect"]], [y], marker=mark.get(r["source"], "o"), color=col, ms=6, mec="black", mew=0.5)
        else:
            ax.plot([0], [y], marker="x", color=AGREE_COLOR["not tested"], ms=5)
    ax.axvline(0, color="black", lw=0.8)
    first = t.drop_duplicates("claim_id").set_index("claim_id")
    ax.set_yticks([ypos[c] for c in ids])
    ax.set_yticklabels([f"{c}  ({_short(first.loc[c, 'group_a'])} vs {_short(first.loc[c, 'group_b'])})" for c in ids], fontsize=7)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Effect (rank-biserial for fraction, count, Shannon; difference in prevalence for prevalence); positive = a above b",
                  fontsize=7)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", color=AGREE_COLOR[k], ls="", mec="black", mew=0.5, label=f"agrees: {k}") for k in AGREE_COLOR]
    h += [Line2D([], [], marker=m, color="white", mec="black", ls="", label=s) for s, m in mark.items()]
    ax.legend(handles=h, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    ax.set_title("Claims of Lax et al. 2025 at BFD scale: effect with 95% bootstrap CI (bars), per claim\n"
                 "Genus values, no phylogenetic control. 'agrees' compares the sign with the paper's direction.", fontsize=8)
    return fig


def page_tf_claims(genus_df, claims_table, o, paper_rows=None, genus_df_nofilter=None, fmap=None, n_boot=2000, seed=1):
    """T10_tf_claims.tsv (BFD rows on genus_df, plus paper_rows when given), T10_tf_claims_nofilter.tsv (BFD rows on
    genus_df_nofilter, when given) and the figure T10_tf_claims_forest. o: an Out for the 'tf' level."""
    from figs_common import safe
    bfd = evaluate_claims(genus_df, claims_table, n_boot, seed, fmap)
    tab = bfd if paper_rows is None else pd.concat([bfd, paper_rows], ignore_index=True)
    o.table(tab, "T10_tf_claims")
    if genus_df_nofilter is not None:
        o.table(evaluate_claims(genus_df_nofilter, claims_table, n_boot, seed, fmap), "T10_tf_claims_nofilter")
    safe(o, "T10_tf_claims_forest", "tf", lambda: _forest(tab))
    return tab
