"""The paper's lineage claims (Lax et al. 2025) tested on BFD genus values, and the forest figure (T10).

Unit: one value per genus (genus representatives from tf_data.genus_table). Claims come from config/tf_claims.tsv. Each claim
compares group a with group b on one metric of one paper family:
  fraction    family proteins / all proteins of the genome   (rank-biserial effect, Mann-Whitney p)
  count       family proteins                                (same)
  shannon     Shannon index of the genome's TF family counts (same; family is empty)
  prevalence  share of genera with at least one protein      (effect = difference in prevalence, Fisher exact p)
  rank        (family empty) Shannon index: the median of group a is ranked among all phyla of group b with at least
              MIN_GENERA genera; group a is tested against EACH of those phyla (rank-biserial, Mann-Whitney, Benjamini-Hochberg
              over the pairwise tests). agrees = yes only if group a ranks first and every pairwise comparison is yes
              (mixed if it ranks first but a pair is not yes; no if it does not rank first). Row effect = the weakest pair.
  max_other_prevalence  prevalence in group a against the MAXIMUM prevalence among the other phyla (>= MIN_GENERA genera);
              effect = difference; agrees = yes when it is at least EXCLUSIVE_MARGIN, mixed when above 0, no otherwise.
The BFD paper-family counts come from tf_paper_family_genome (distinct proteins per paper family; genus table columns
pf_<paper family>). Without that table the counts of the BFD families of one paper family are added (a protein with two
domains of one paper family then counts twice); every such row has a note.
A paper family is "not comparable" when BFD detects it in fewer than MIN_NONZERO_GENERA genera of the table: agrees is then
'not comparable' and the note says so (the statistics are still shown).
Groups are phyla or the sets Dikarya (Ascomycota + Basidiomycota), EDF (all phyla except Ascomycota, Basidiomycota and
Microsporidia), Mucoromyceta (Mucoromycota + Mortierellomycota + Glomeromycota), ALL (every phylum); 'X-Y' removes the phylum Y from the set X (for example ALL-Mucoromycota).
Confidence intervals are percentile bootstraps over genera (n_boot resamples within each group, fixed seed). p-values are
Benjamini-Hochberg adjusted over the tested claims of ONE call, so separately for the BFD table (BUSCO >= 80), the BFD no-filter
sensitivity table and the paper-side table. Phylogeny is not controlled: related genera are not independent.
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
MIN_NONZERO_GENERA = 3      # a paper family seen in fewer genera than this is "not comparable"
EXCLUSIVE_MARGIN = 0.3      # max_other_prevalence: yes when the prevalence in a exceeds the best other phylum by this much
PF_PREFIX = "pf_"
SUM_NOTE = "paper-family counts summed over Pfam families (double counting possible)"
SENS_NOTE = ("sensitivity only; low-BUSCO genomes have incomplete annotation; an a<b result here is not evidence against "
             "the paper")
DIKARYA = ("Ascomycota", "Basidiomycota")
NOT_EDF = DIKARYA + ("Microsporidia",)
MUCOROMYCETA = ("Mucoromycota", "Mortierellomycota", "Glomeromycota")   # one reading of "Mucoromycota" in the paper
FAM_PREFIX = "fam_"
UNIT = {"genus": "genera", "species": "species"}
CLAIM_COLS = ["claim_id", "text", "family", "metric", "level", "group_a", "group_b", "paper_direction"]
OUT_COLS = ["claim_id", "text", "family", "metric", "level", "group_a", "group_b", "n_genera_a", "n_genera_b", "value_a",
            "value_b", "effect", "ci_lo", "ci_hi", "p", "p_adj", "observed_direction", "paper_direction", "agrees", "note",
            "source", "filter"]
AGREE_COLOR = {"yes": "#009E73", "mixed": "#E69F00", "no": "#D55E00", "not tested": "#8c8c8c",
               "not comparable": "#b3b3b3"}   # Okabe-Ito; the greys differ by marker (x for not tested)
NC_LEGEND = "not comparable (BFD Pfam definition does not detect this family)"


def load_claims(path=None):
    """config/tf_claims.tsv (lines starting with # are comments) with the columns CLAIM_COLS."""
    c = pd.read_csv(path or CLAIMS_FILE, sep="\t", dtype=str, comment="#", keep_default_na=False)
    missing = [k for k in CLAIM_COLS if k not in c.columns]
    if missing:
        raise ValueError(f"claims file lacks columns {missing}")
    return c[CLAIM_COLS]


def load_fmap(path=None):
    return pd.read_csv(path or FAMILY_MAP, sep="\t", dtype=str, comment="#", keep_default_na=False)


def _base_mask(name, ph, known):
    if name == "Dikarya":
        return ph.isin(DIKARYA)
    if name == "EDF":
        return known & ~ph.isin(NOT_EDF)
    if name == "Mucoromyceta":
        return ph.isin(MUCOROMYCETA)
    if name == "ALL":
        return known
    return ph == name


def resolve_group(spec, phylum):
    """Boolean numpy mask of the rows of `phylum` (Series) in the group `spec`: a phylum, or one of the sets Dikarya, EDF,
    Mucoromyceta (Mucoromycota + Mortierellomycota + Glomeromycota), ALL; 'X-Y' removes the set or phylum Y from X. Unknown
    phylum (NaN) is in no group."""
    ph = pd.Series(phylum).reset_index(drop=True)
    known = ph.notna()
    base, *minus = [s.strip() for s in str(spec).split("-")]
    m = _base_mask(base, ph, known)
    for s in minus:
        m = m & ~_base_mask(s, ph, known)
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


def _empty(claim, source, note, filter_label=""):
    r = {k: claim[k] for k in CLAIM_COLS}
    r.update(n_genera_a=0, n_genera_b=0, value_a=np.nan, value_b=np.nan, effect=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan,
             p_adj=np.nan, observed_direction="", agrees="not tested", note=note, source=source, **{"filter": filter_label})
    return r


def _join(*parts):
    return "; ".join(p for p in parts if p)


def _pair(x, y, n_boot, rng):
    """(effect, ci_lo, ci_hi, p) of a rank-biserial comparison of x with y."""
    lo, hi = _boot_ci(x, y, False, n_boot, rng)
    p = float(st.mannwhitneyu(x, y, alternative="two-sided").pvalue) if np.ptp(np.r_[x, y]) > 0 else 1.0
    return rank_biserial(x, y), lo, hi, p


def _phylum_values(v, ph, mask, min_n):
    """{phylum: finite values} for the phyla of the rows in `mask` with at least min_n finite values."""
    out = {}
    for name, idx in pd.Series(np.arange(len(v))[mask]).groupby(ph.to_numpy()[mask]):
        x = v[idx.to_numpy()]
        x = x[np.isfinite(x)]
        if len(x) >= min_n:
            out[name] = x
    return out


def _eligible_note(v, ph, mask, others, unit):
    """'eligible phyla: ...; excluded (< 5 unit): Name (n), ...' for the phyla of `mask`."""
    cnt = {}
    for name, idx in pd.Series(np.arange(len(v))[mask]).groupby(ph.to_numpy()[mask]):
        cnt[name] = int(np.isfinite(v[idx.to_numpy()]).sum())
    elig = ", ".join(sorted(others))
    exc = ", ".join(f"{k} ({n})" for k, n in sorted(cnt.items()) if k not in others and n > 0)
    return f"eligible phyla: {elig}" + (f"; excluded (< {MIN_GENERA} {unit}): {exc}" if exc else "")


def _rank_row(r, v, ph, ma, mb, n_boot, rng, unit):
    a = v[ma][np.isfinite(v[ma])]
    others = _phylum_values(v, ph, mb & ~ma, MIN_GENERA)
    r["n_genera_a"], r["n_genera_b"] = len(a), int(sum(len(x) for x in others.values()))
    if len(a) < MIN_GENERA or not others:
        r["note"] = _join(r["note"], f"fewer than {MIN_GENERA} {unit} in group a or no other phylum with {MIN_GENERA} {unit}")
        return r
    med = {k: float(np.median(x)) for k, x in others.items()}
    rank = 1 + sum(m > np.median(a) for m in med.values())
    pairs = {k: _pair(a, x, n_boot, rng) for k, x in sorted(others.items())}
    padj = dict(zip(pairs, bh_adjust([q[3] for q in pairs.values()])))
    verdict = {k: agreement(q[0], q[1], q[2], "a>b") for k, q in pairs.items()}
    weakest = min(pairs, key=lambda k: pairs[k][0])
    e, lo, hi, p = pairs[weakest]
    notyes = [k for k, a_ in verdict.items() if a_ != "yes"]
    r.update(value_a=float(np.median(a)), value_b=max(med.values()), effect=e, ci_lo=lo, ci_hi=hi, p=max(padj.values()),
             observed_direction=observed_direction(e))
    r["agrees"] = "yes" if rank == 1 and not notyes else "mixed" if rank == 1 else "no"
    r["note"] = _join(r["note"], f"rank {rank} of {len(med) + 1} phyla with >= {MIN_GENERA} {unit}; weakest pair vs {weakest} "
                      f"(effect {e:.2f}); p = largest BH-adjusted pairwise p" + (f"; not yes vs: {', '.join(notyes)}" if notyes else ""),
                      _eligible_note(v, ph, mb & ~ma, others, unit))
    return r


def _maxprev_row(r, v, ph, ma, mb, n_boot, rng, unit):
    a = v[ma][np.isfinite(v[ma])]
    others = _phylum_values(v, ph, mb & ~ma, MIN_GENERA)
    r["n_genera_a"], r["n_genera_b"] = len(a), int(sum(len(x) for x in others.values()))
    if len(a) < MIN_GENERA or not others:
        r["note"] = _join(r["note"], f"fewer than {MIN_GENERA} {unit} in group a or no other phylum with {MIN_GENERA} {unit}")
        return r
    prev = {k: float((x >= 1).mean()) for k, x in sorted(others.items())}
    top = max(prev, key=lambda k: prev[k])
    xa, xb = (a >= 1).astype(float), (others[top] >= 1).astype(float)
    ka, kb = int(xa.sum()), int(xb.sum())
    e = float(xa.mean() - xb.mean())
    lo, hi = _boot_ci(xa, xb, True, n_boot, rng)
    r.update(value_a=float(xa.mean()), value_b=prev[top], effect=e, ci_lo=lo, ci_hi=hi, observed_direction=observed_direction(e),
             p=float(st.fisher_exact([[ka, len(xa) - ka], [kb, len(xb) - kb]])[1]))
    r["agrees"] = "yes" if e >= EXCLUSIVE_MARGIN else "mixed" if e >= 0 else "no"
    tied = [k for k in prev if prev[k] == prev[top]]
    what = (f"tie at the maximum prevalence with: {', '.join(tied)} ({prev[top]:.2f})" if e == 0 else
            f"highest other prevalence: {top} ({prev[top]:.2f})")
    r["note"] = _join(r["note"], f"{what}; yes needs a margin of at least {EXCLUSIVE_MARGIN} over it",
                      _eligible_note(v, ph, mb & ~ma, others, unit))
    return r


def evaluate_table(df, claims, value_fn, n_boot=2000, seed=1, source="BFD", level="genus", flag_fn=None, filter_label="",
                   extra_note=""):
    """Core of the claim tests. df: one row per unit (genus or species) with a PHYLUM column. value_fn(family, metric) ->
    (values aligned to df.index, note): metric is 'fraction', 'count' or 'shannon'; a prevalence claim uses 'count' and a value
    of at least 1. values None = not available (the note says why). flag_fn(family) -> (not_comparable, note) or None: when
    not_comparable, agrees becomes 'not comparable' (the statistics stay); the note is added in both cases. extra_note is
    added to every row. Returns OUT_COLS, p_adj by Benjamini-Hochberg over this call."""
    unit = UNIT.get(level, level)
    ph = pd.Series(df["PHYLUM"]).reset_index(drop=True)
    rows = []
    for i, c in claims.reset_index(drop=True).iterrows():
        c = c.to_dict()
        c["level"] = level
        metric, fam = c["metric"], c["family"]
        vm = "shannon" if (metric == "rank" and not fam) else "count" if metric in ("prevalence", "max_other_prevalence") else \
            "shannon" if metric == "shannon" else metric
        vals, note = value_fn(fam, vm)
        ma, mb = resolve_group(c["group_a"], df["PHYLUM"]), resolve_group(c["group_b"], df["PHYLUM"])
        if vals is None:
            r = _empty(c, source, _join(note, extra_note), filter_label)
            r["n_genera_a"], r["n_genera_b"] = int(ma.sum()), int(mb.sum())
            rows.append(r)
            continue
        v = np.asarray(vals, dtype=float)
        flag = flag_fn(fam) if (flag_fn and fam and metric != "shannon") else None
        r = _empty(c, source, note or "", filter_label)
        rng = np.random.default_rng([seed, i])
        if metric in ("rank", "max_other_prevalence"):
            r = (_rank_row if metric == "rank" else _maxprev_row)(r, v, ph, ma, mb, n_boot, rng, unit)
        else:
            x, y = v[ma][np.isfinite(v[ma])], v[mb][np.isfinite(v[mb])]
            r["n_genera_a"], r["n_genera_b"] = len(x), len(y)
            small = [g for g, a in (("a", x), ("b", y)) if len(a) < MIN_GENERA]
            if small:
                r["note"] = _join(r["note"], f"fewer than {MIN_GENERA} {unit} in group {' and '.join(small)}")
            elif metric == "prevalence":
                xb, yb = (x >= 1).astype(float), (y >= 1).astype(float)
                ka, kb = int(xb.sum()), int(yb.sum())
                r["value_a"], r["value_b"] = ka / len(x), kb / len(y)
                r["effect"] = r["value_a"] - r["value_b"]
                r["p"] = float(st.fisher_exact([[ka, len(x) - ka], [kb, len(y) - kb]])[1])
                r["ci_lo"], r["ci_hi"] = _boot_ci(xb, yb, True, n_boot, rng)
            else:
                r["value_a"], r["value_b"] = float(np.median(x)), float(np.median(y))
                r["effect"], r["ci_lo"], r["ci_hi"], r["p"] = _pair(x, y, n_boot, rng)
            if not small:
                r["observed_direction"] = observed_direction(r["effect"])
                r["agrees"] = agreement(r["effect"], r["ci_lo"], r["ci_hi"], c["paper_direction"])
        if r["agrees"] == "no" and np.isfinite(r["ci_lo"]) and r["ci_lo"] <= 0 <= r["ci_hi"]:
            r["note"] = _join(r["note"], "CI includes 0")
        if flag and np.isfinite(r["effect"]):
            if flag[0]:
                r["agrees"] = "not comparable"
            r["note"] = _join(flag[1], r["note"]) if flag[0] else _join(r["note"], flag[1])
        r["note"] = _join(r["note"], extra_note)
        rows.append(r)
    out = pd.DataFrame(rows, columns=OUT_COLS)
    out["p_adj"] = bh_adjust(out["p"].to_numpy(dtype=float))
    return out


def nonzero_genera(genus_df, fmap=None):
    """{paper family: number of genera with at least one protein}, from the pf_ columns of the genus table (or the sum of
    the fam_ columns of its BFD families when the table has no pf_ columns)."""
    fm = load_fmap() if fmap is None else fmap
    labels = _labels(fm)
    out = {}
    for f, labs in labels.items():
        col = PF_PREFIX + f
        if col in genus_df.columns:
            cnt = genus_df[col].astype(float)
        else:
            have = [FAM_PREFIX + x for x in labs if FAM_PREFIX + x in genus_df.columns]
            cnt = genus_df[have].astype(float).sum(axis=1) if have else pd.Series(0.0, index=genus_df.index)
        out[f] = int((cnt > 0).sum())
    return out


def not_comparable_note(k, basis="BUSCO>=80"):
    """k = number of genera with a protein in the BUSCO>=80 genus table (the one definition used for T10, T09 and the paper
    side)."""
    n = f"{k} genus" if k == 1 else f"{k} genera"
    return (f"BFD Pfam definition does not detect this family (non-zero in {n}); not comparable with the paper "
            f"(genus table: {basis})")


def _labels(fm):
    bridge = fm[(fm["paper_family"] != "") & (fm["family"] != "")][["family", "paper_family"]].drop_duplicates()
    return bridge.groupby("paper_family")["family"].agg(list).to_dict()


def evaluate_claims(genus_df, claims, n_boot=2000, seed=1, fmap=None, filter_label="", extra_note="", nonzero=None):
    """Test each claim on genus values. genus_df: tf_data.genus_table output (PHYLUM, n_proteins_total, shannon, fam_<label>
    and, with the paper-family table, pf_<paper family>). fmap: the family map (default config/tf_family_map.tsv). Returns
    OUT_COLS (source = 'BFD'; filter = filter_label). nonzero: {paper family: genera with a protein} used for the
    not-comparable flag (default: counted on genus_df); pass the BUSCO>=80 table's numbers so every table uses one definition."""
    fm = load_fmap() if fmap is None else fmap
    labels = _labels(fm)
    pf_mode = any(c.startswith(PF_PREFIX) for c in genus_df.columns)
    tot = genus_df["n_proteins_total"].astype(float).where(genus_df["n_proteins_total"] > 0)
    nz = dict(nonzero) if nonzero is not None else nonzero_genera(genus_df, fm)

    def value_fn(family, metric):
        if metric == "shannon":
            return (genus_df["shannon"].to_numpy(dtype=float), "") if "shannon" in genus_df else (None, "no shannon column")
        if family not in labels:
            return None, f"family '{family}' is not a paper family in the family map"
        if pf_mode:
            col = PF_PREFIX + family
            cnt = genus_df[col].astype(float) if col in genus_df.columns else pd.Series(0.0, index=genus_df.index)
            note = "" if col in genus_df.columns else "no genome has this family"
        else:
            cols = [FAM_PREFIX + f for f in labels[family]]
            have = [c for c in cols if c in genus_df.columns]
            if not have:
                return None, f"family column {', '.join(cols)} missing from the genus table (no genome has this family)"
            cnt = genus_df[have].astype(float).sum(axis=1)
            miss = [c for c in cols if c not in have]
            note = _join(f"column {', '.join(miss)} missing: counted as 0" if miss else "",
                         SUM_NOTE if len(labels[family]) > 1 else "")
        return ((cnt / tot).to_numpy() if metric == "fraction" else cnt.to_numpy()), note

    def flag_fn(family):
        k = nz.get(family)
        if k is None:
            return None
        return (True, not_comparable_note(k)) if k < MIN_NONZERO_GENERA else (False, "")

    return evaluate_table(genus_df, claims, value_fn, n_boot, seed, "BFD", "genus", flag_fn, filter_label, extra_note)


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
            if np.isfinite(r["ci_lo"]):
                ax.plot([r["ci_lo"], r["ci_hi"]], [y, y], color=col, lw=1.6, solid_capstyle="butt")
            weak_no = r["agrees"] == "no" and np.isfinite(r["ci_lo"]) and r["ci_lo"] <= 0 <= r["ci_hi"]
            ax.plot([r["effect"]], [y], marker=mark.get(r["source"], "o"), color=col, ms=6, mec=col if weak_no else "black",
                    mew=1.2 if weak_no else 0.5, mfc="white" if weak_no else col)
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
    h = [Line2D([], [], marker="o", color=AGREE_COLOR[k], ls="", mec="black", mew=0.5,
                label=f"agrees: {NC_LEGEND if k == 'not comparable' else k}") for k in AGREE_COLOR]
    h += [Line2D([], [], marker="o", mfc="white", mec="#D55E00", ls="", mew=1.2, label="no, but CI includes 0 (open marker)")]
    h += [Line2D([], [], marker=m, color="white", mec="black", ls="", label=s) for s, m in mark.items()]
    ax.legend(handles=h, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.045), ncol=2, frameon=False)
    ax.set_title("Claims of Lax et al. 2025 at BFD scale: effect with 95% bootstrap CI\n"
                 "BFD rows: genus representatives with BUSCO >= 80. Genus values, no phylogenetic control.\n"
                 "'agrees' compares the sign with the paper's direction.",
                 fontsize=8)
    return fig


def page_tf_claims(genus_df, claims_table, o, paper_rows=None, genus_df_nofilter=None, fmap=None, n_boot=2000, seed=1,
                   nonzero=None):
    """T10_tf_claims.tsv (BFD rows on genus_df with filter BUSCO>=80, plus paper_rows when given),
    T10_tf_claims_nofilter.tsv (BFD rows on genus_df_nofilter, filter none, every note says sensitivity only) and the figure
    T10_tf_claims_forest. o: an Out for the 'tf' level."""
    from figs_common import safe
    if nonzero is None:
        nonzero = nonzero_genera(genus_df, load_fmap() if fmap is None else fmap)
    bfd = evaluate_claims(genus_df, claims_table, n_boot, seed, fmap, filter_label="BUSCO>=80", nonzero=nonzero)
    tab = bfd if paper_rows is None else pd.concat([bfd, paper_rows], ignore_index=True)
    o.table(tab, "T10_tf_claims")
    if genus_df_nofilter is not None:
        o.table(evaluate_claims(genus_df_nofilter, claims_table, n_boot, seed, fmap, filter_label="none", extra_note=SENS_NOTE,
                        nonzero=nonzero),
                "T10_tf_claims_nofilter")
    safe(o, "T10_tf_claims_forest", "tf", lambda: _forest(tab))
    return tab
