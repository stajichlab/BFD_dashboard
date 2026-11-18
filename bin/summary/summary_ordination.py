"""Genus representatives and PCoA ordination of codon usage, amino-acid composition and gene structure.

This file is named summary_ordination.py, not ordination.py, because lib/ordination.py has the PCoA and
PERMANOVA code and both folders can be on sys.path. The library module is loaded here by file path under the
name lib_ordination, so `import ordination` cannot pick the wrong file.
"""
import importlib.util
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figs_common import safe  # noqa: E402
from make_figures import Out, assign_groups, parent_legend  # noqa: E402

_spec = importlib.util.spec_from_file_location("lib_ordination", Path(__file__).resolve().parents[2] / "lib" / "ordination.py")
lib_ordination = importlib.util.module_from_spec(_spec)
sys.modules["lib_ordination"] = lib_ordination
_spec.loader.exec_module(lib_ordination)

N_AXES = 5
N_PERM = 999
SEED = 1
METRIC = {"codon": "bray-curtis", "aa": "bray-curtis", "structure": "euclidean"}
TESTS = ("PHYLUM", "ORDER")
TOPIC = "ordination"
PREFIX = {"codon": "codon_", "aa": "aa_"}


def pick_genus_reps(core):
    """One row per non-null GENUS: highest busco_complete, then highest N50_bp, then lowest species_prefix.
    Rows with a NULL genus are left out; their number is in attrs['n_no_genus']."""
    has = core["GENUS"].notna()
    d = core[has].sort_values(["GENUS", "busco_complete", "N50_bp", "species_prefix"],
                              ascending=[True, False, False, True], na_position="last")
    reps = d.drop_duplicates("GENUS", keep="first").reset_index(drop=True)
    reps.attrs["n_no_genus"] = int((~has).sum())
    return reps


def _complete_ids(name, table):
    """species_prefix values with usable data in one matrix (same rule as the ordination itself)."""
    if name == "structure":
        return set(structure_matrix(table).index)
    mat = _composition(name, table)
    return set(mat.index[mat.notna().all(axis=1)])


def _composition(name, table):
    cols = [c for c in table.columns if c.startswith(PREFIX[name])]
    return table.set_index("species_prefix")[cols].astype(float)


def pick_complete_reps(core, matrices):
    """Genus representatives chosen from genomes that have usable data in ALL the given matrices.
    The pool is built first and pick_genus_reps runs on it, so a genus survives through its best data-complete
    genome. One rep set serves all ordinations. attrs: n_no_genus (genomes with NULL genus, whole core) and
    n_genera_no_complete_data (genera with no data-complete genome)."""
    ids = None
    for name, table in matrices.items():
        have = _complete_ids(name, table)
        ids = have if ids is None else ids & have
    pool = core[core["species_prefix"].isin(ids if ids is not None else set())]
    reps = pick_genus_reps(pool)
    reps.attrs["n_no_genus"] = int(core["GENUS"].isna().sum())
    reps.attrs["n_genera_no_complete_data"] = int(core["GENUS"].dropna().nunique() - len(reps))
    return reps


CLIP_Z = 4.0


def structure_matrix(table, clip=CLIP_Z):
    """Gene-structure stats as a matrix: all-NaN columns dropped, positive columns log-transformed, every column
    standardized (zero mean, unit SD), constant columns dropped, rows with NaN dropped. Index: species_prefix.
    Each standardized value is then clipped to [-clip, clip] (clip=None: no clipping) so a few extreme genomes do
    not dominate the axes. attrs: n_clipped_values (values changed) and n_genera_clipped (rows with any)."""
    m = table.set_index("species_prefix").select_dtypes("number").astype(float)
    m = m.dropna(axis=1, how="all")
    for c in m.columns:
        v = m[c].dropna()
        if len(v) and (v > 0).all():
            m[c] = np.log(m[c])
    m = m.dropna(axis=0)
    sd = m.std(ddof=0)
    m = m.loc[:, sd > 0]
    z = (m - m.mean()) / m.std(ddof=0)
    hit = (z.abs() > clip) if clip is not None else pd.DataFrame(False, index=z.index, columns=z.columns)
    out = z.clip(-clip, clip) if clip is not None else z
    out.attrs["n_clipped_values"] = int(hit.to_numpy().sum())
    out.attrs["n_genera_clipped"] = int(hit.any(axis=1).sum())
    return out


def _feature_matrix(name, table, reps):
    """Feature matrix of the representatives: (matrix, n_missing_from_table, n_dropped_nan). With reps from
    pick_complete_reps both counts are 0; they are nonzero only if reps were chosen without that pool."""
    t = table[table["species_prefix"].isin(reps["species_prefix"])]
    n_missing = len(reps) - t["species_prefix"].nunique()
    if name == "structure":
        full = t.set_index("species_prefix").select_dtypes("number")
        mat = structure_matrix(t)
        n_nan = len(full) - len(mat)
        return mat, n_missing, n_nan
    mat = _composition(name, t)
    ok = mat.notna().all(axis=1)
    return mat[ok], n_missing, int((~ok).sum())


def _log10_count(n):
    """log10 of an exact (possibly astronomically large) Python int; 0.0 when the count is 0 or unknown.
    The raw int cannot go into a float64 DataFrame column (OverflowError for 1,600+ genera in many groups)."""
    try:
        n = int(n)
    except (TypeError, ValueError):
        return 0.0
    return math.log10(n) if n > 0 else 0.0


def _r2_adj(r2, n, k):
    """Adjusted R2 = 1 - (1 - R2) * (n - 1) / (n - k); NaN when n <= k or R2 is not finite."""
    return float(1 - (1 - r2) * (n - 1) / (n - k)) if np.isfinite(r2) and n > k else float("nan")


def _permanova(dist, groups):
    """PERMANOVA of one grouping. Rows with a NULL group are left out. Failures and degenerate designs give NaN."""
    g = groups.reindex(dist.index)
    keep = g.notna()
    nan = {"F": np.nan, "R2": np.nan, "R2_adj": np.nan, "p": np.nan, "n": int(keep.sum()), "n_groups": int(g[keep].nunique()),
           "n_permutations_run": 0, "n_unique_permutations_log10": 0.0, "caveat": "not tested: fewer than 3 genera or fewer than 2 groups"}
    if nan["n"] < 3 or nan["n_groups"] < 2:
        return nan
    sub = dist.loc[keep, keep]
    try:
        r = lib_ordination.permanova(sub, g[keep], n_permutations=N_PERM, seed=SEED)
    except Exception as e:  # a degenerate design must not stop the run
        print(f"  permanova failed: {type(e).__name__}: {e}", flush=True)
        return nan
    k, n = nan["n_groups"], nan["n"]
    f = float(r.pseudo_f)
    r2 = f * (k - 1) / (f * (k - 1) + (n - k)) if np.isfinite(f) and n > k else np.nan   # PERMANOVA R^2 from pseudo-F
    return {"F": f, "R2": float(r2), "R2_adj": _r2_adj(r2, n, k), "p": float(r.p_value), "n": nan["n"], "n_groups": nan["n_groups"],
            "n_permutations_run": int(r.n_permutations_run), "n_unique_permutations_log10": _log10_count(r.n_unique_permutations),
            "caveat": r.caveat}


def _scatter(coords, groups, colors, pct, src, title):
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    c = groups.map(colors).fillna("#9a9a9a")
    ax.scatter(coords["PCo1"], coords["PCo2"] if "PCo2" in coords else 0, s=9, c=c, alpha=0.85, linewidths=0, rasterized=True)
    ax.set_xlabel(f"PCo1 ({pct[0]:.1f}%)")
    ax.set_ylabel(f"PCo2 ({pct[1]:.1f}%)" if len(pct) > 1 else "PCo2")
    ax.set_title(title)
    leg = parent_legend(ax, src, outside=True, ncol=2, fontsize=6, title_fontsize=6.5)
    if leg is not None:   # make room on the right for the legend: savefig is not called with bbox_inches="tight"
        fig.canvas.draw()
        w_in = leg.get_window_extent().width / fig.dpi
        right = 1 - (w_in + 0.05 * fig.get_figwidth() / 8.4 + 0.1) / fig.get_figwidth()
        fig.subplots_adjust(right=max(0.5, right))
    return fig


def run_ordination(reps, matrices, out, phylum_order=None):
    """PCoA and PERMANOVA for each of codon, aa and structure that is in `matrices` (name -> per-genome table with
    species_prefix). Writes tables/ordination and figures/ordination under `out`. Returns a dict per matrix:
    n_genera, n_missing, n_dropped_nan (reps without a table row / with NaN features; 0 when reps come from
    pick_complete_reps), pct_variance (first axes), n_negative_eigen, negative_eigen_share (= sum of negative
    eigenvalues / sum of |eigenvalues|; table column neg_eigen_over_abs_total), neg_eigen_over_pos_total
    (= |sum of negative| / sum of positive), permanova. "_summary" has n_no_genus and n_genera_no_complete_data
    from reps.attrs (also written to ordination_summary.tsv)."""
    o = Out(out, TOPIC)
    o.scope = "genus_reps"
    result = {}
    summary = {"n_no_genus": int(reps.attrs.get("n_no_genus", 0)),
               "n_genera_no_complete_data": int(reps.attrs.get("n_genera_no_complete_data", 0)),
               "n_clipped_values": 0, "n_genera_clipped": 0}
    print(f"ordination: {summary['n_no_genus']} genomes with no genus excluded; {summary['n_genera_no_complete_data']} "
          f"genera have no genome with complete data; {len(reps)} genus representatives", flush=True)
    result["_summary"] = summary
    ph, _, ph_col = assign_groups(reps, "phylum", phylum_order=phylum_order)
    cl, _, cl_col = assign_groups(reps, "class", phylum_order=phylum_order)
    ph_idx = ph.set_index("species_prefix")
    cl_idx = cl.set_index("species_prefix")
    for name in ("codon", "aa", "structure"):
        if name not in matrices:
            continue
        mat, n_missing, n_nan = _feature_matrix(name, matrices[name], reps)
        info = {"n_genera": int(len(mat)), "n_missing": int(n_missing), "n_dropped_nan": int(n_nan),
                "n_clipped_values": int(mat.attrs.get("n_clipped_values", 0)), "n_genera_clipped": int(mat.attrs.get("n_genera_clipped", 0)),
                "pct_variance": [], "n_negative_eigen": 0, "negative_eigen_share": 0.0,
                "neg_eigen_over_pos_total": 0.0,
                "permanova": {t: {"F": np.nan, "R2": np.nan, "R2_adj": np.nan, "p": np.nan, "n": 0, "n_groups": 0, "n_permutations_run": 0,
                                  "n_unique_permutations_log10": 0.0, "caveat": "not tested"} for t in TESTS}}
        result[name] = info
        if name == "structure":
            summary["n_clipped_values"], summary["n_genera_clipped"] = info["n_clipped_values"], info["n_genera_clipped"]
        if len(mat) < 3:
            print(f"  {name}: {len(mat)} genera, too few for PCoA", flush=True)
            continue
        dist = lib_ordination.distance_matrix(mat, METRIC[name])
        ord_ = lib_ordination.pcoa(dist, n_axes=N_AXES)
        ev = ord_.eigenvalues
        neg = ev < -1e-10
        info["pct_variance"] = [float(x) for x in ord_.pct_variance[:ord_.coords.shape[1]]]
        info["n_negative_eigen"] = int(neg.sum())
        info["negative_eigen_share"] = float(-ev[neg].sum() / np.abs(ev).sum()) if neg.any() else 0.0
        pos = ev > 1e-10
        info["neg_eigen_over_pos_total"] = float(-ev[neg].sum() / ev[pos].sum()) if neg.any() and pos.any() else 0.0
        meta = ph_idx.reindex(mat.index)
        for t in TESTS:
            info["permanova"][t] = _permanova(dist, meta[t])
        coords = ord_.coords.copy()
        coords.index.name = "species_prefix"
        tab = coords.reset_index().merge(reps[["species_prefix", "GENUS", "PHYLUM", "CLASS", "ORDER"]], on="species_prefix", how="left")
        o.table(tab, f"{name}_coords")
        rows = [{"matrix": name, "grouping": t, "metric": METRIC[name], "n_genera": r["n"], "n_groups": r["n_groups"],
                 "pseudo_F": r["F"], "R2": r["R2"], "R2_adj": r["R2_adj"], "p_value": r["p"], "n_permutations_requested": N_PERM,
                 "n_permutations_run": r["n_permutations_run"], "n_unique_permutations_log10": r["n_unique_permutations_log10"], "seed": SEED,
                 "n_negative_eigen": info["n_negative_eigen"],
                 "neg_eigen_over_abs_total": info["negative_eigen_share"],
                 "neg_eigen_over_pos_total": info["neg_eigen_over_pos_total"], "caveat": r["caveat"],
                 "n_dropped_nan": info["n_dropped_nan"], "n_missing": info["n_missing"],
                 "n_clipped_values": info["n_clipped_values"], "n_genera_clipped": info["n_genera_clipped"]}
                for t, r in info["permanova"].items()]
        o.table(pd.DataFrame(rows), f"{name}_permanova")
        pct = info["pct_variance"]
        for stem, by, cols, src, idx in ((f"O01_{name}_pcoa_phylum", "phylum", ph_col, ph, ph_idx),
                                         (f"O02_{name}_pcoa_class", "class (phylum hue, shade by class)", cl_col, cl, cl_idx)):
            groups = idx["grp"].reindex(mat.index)
            ph_p = info["permanova"]["PHYLUM"]
            o.meta[stem] = {"n": int(len(mat)), "n_groups": int(groups.nunique()), "signal_kind": "permanova_R2",
                            "signal": ph_p["R2"] if np.isfinite(ph_p["R2"]) else None}
            title = f"{name} PCoA of {len(mat):,} genus representatives, coloured by {by}"
            safe(o, stem, TOPIC, lambda groups=groups, cols=cols, src=src, title=title: _scatter(coords, groups, cols, pct, src, title))
    o.table(pd.DataFrame([summary]), "ordination_summary")
    print(f"ordination: {len(result)} matrices, {len(o.made)} outputs", flush=True)
    result["_made"] = o.made
    return result
