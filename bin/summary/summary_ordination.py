"""Genus representatives and PCoA ordination of codon usage, amino-acid composition and gene structure.

This file is named summary_ordination.py, not ordination.py, because lib/ordination.py has the PCoA and
PERMANOVA code and both folders can be on sys.path. The library module is loaded here by file path under the
name lib_ordination, so `import ordination` cannot pick the wrong file.
"""
import importlib.util
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


def pick_genus_reps(core):
    """One row per non-null GENUS: highest busco_complete, then highest N50_bp, then lowest species_prefix.
    Rows with a NULL genus are left out; their number is in attrs['n_no_genus']."""
    has = core["GENUS"].notna()
    d = core[has].sort_values(["GENUS", "busco_complete", "N50_bp", "species_prefix"],
                              ascending=[True, False, False, True], na_position="last")
    reps = d.drop_duplicates("GENUS", keep="first").reset_index(drop=True)
    reps.attrs["n_no_genus"] = int((~has).sum())
    return reps


def structure_matrix(table):
    """Gene-structure stats as a matrix: all-NaN columns dropped, positive columns log-transformed, every column
    standardized (zero mean, unit SD), constant columns dropped, rows with NaN dropped. Index: species_prefix."""
    m = table.set_index("species_prefix").select_dtypes("number").astype(float)
    m = m.dropna(axis=1, how="all")
    for c in m.columns:
        v = m[c].dropna()
        if len(v) and (v > 0).all():
            m[c] = np.log(m[c])
    m = m.dropna(axis=0)
    sd = m.std(ddof=0)
    m = m.loc[:, sd > 0]
    return (m - m.mean()) / m.std(ddof=0)


def _feature_matrix(name, table, reps):
    """Feature matrix of the representatives: (matrix, n_missing_from_table, n_dropped_nan)."""
    t = table[table["species_prefix"].isin(reps["species_prefix"])]
    n_missing = len(reps) - t["species_prefix"].nunique()
    if name == "structure":
        full = t.set_index("species_prefix").select_dtypes("number")
        mat = structure_matrix(t)
        n_nan = len(full) - len(mat)
        return mat, n_missing, n_nan
    mat = t.set_index("species_prefix").select_dtypes("number").astype(float)
    ok = mat.notna().all(axis=1)
    return mat[ok], n_missing, int((~ok).sum())


def _permanova(dist, groups):
    """PERMANOVA of one grouping. Rows with a NULL group are left out. Failures and degenerate designs give NaN."""
    g = groups.reindex(dist.index)
    keep = g.notna()
    nan = {"F": np.nan, "p": np.nan, "n": int(keep.sum()), "n_groups": int(g[keep].nunique())}
    if nan["n"] < 3 or nan["n_groups"] < 2:
        return nan
    sub = dist.loc[keep, keep]
    try:
        r = lib_ordination.permanova(sub, g[keep], n_permutations=N_PERM, seed=SEED)
    except Exception as e:  # a degenerate design must not stop the run
        print(f"  permanova failed: {type(e).__name__}: {e}", flush=True)
        return nan
    return {"F": float(r.pseudo_f), "p": float(r.p_value), "n": nan["n"], "n_groups": nan["n_groups"]}


def _scatter(coords, groups, colors, pct, src, title):
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    c = groups.map(colors).fillna("#9a9a9a")
    ax.scatter(coords["PCo1"], coords["PCo2"] if "PCo2" in coords else 0, s=9, c=c, alpha=0.85, linewidths=0)
    ax.set_xlabel(f"PCo1 ({pct[0]:.1f}%)")
    ax.set_ylabel(f"PCo2 ({pct[1]:.1f}%)" if len(pct) > 1 else "PCo2")
    ax.set_title(title)
    parent_legend(ax, src)
    return fig


def run_ordination(reps, matrices, out):
    """PCoA and PERMANOVA for each of codon, aa and structure that is in `matrices` (name -> per-genome table with
    species_prefix). Writes tables/ordination and figures/ordination under `out`. Returns a dict per matrix:
    n_genera, n_missing, n_dropped_nan, pct_variance (first axes), n_negative_eigen, negative_eigen_share, permanova."""
    o = Out(out, TOPIC)
    o.scope = "genus_reps"
    result = {}
    ph, _, ph_col = assign_groups(reps, "phylum")
    cl, _, cl_col = assign_groups(reps, "class")
    ph_idx = ph.set_index("species_prefix")
    cl_idx = cl.set_index("species_prefix")
    for name in ("codon", "aa", "structure"):
        if name not in matrices:
            continue
        mat, n_missing, n_nan = _feature_matrix(name, matrices[name], reps)
        info = {"n_genera": int(len(mat)), "n_missing": int(n_missing), "n_dropped_nan": int(n_nan),
                "pct_variance": [], "n_negative_eigen": 0, "negative_eigen_share": 0.0,
                "permanova": {t: {"F": np.nan, "p": np.nan, "n": 0, "n_groups": 0} for t in TESTS}}
        result[name] = info
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
        meta = ph_idx.reindex(mat.index)
        for t in TESTS:
            info["permanova"][t] = _permanova(dist, meta[t])
        coords = ord_.coords.copy()
        coords.index.name = "species_prefix"
        tab = coords.reset_index().merge(reps[["species_prefix", "GENUS", "PHYLUM", "CLASS", "ORDER"]], on="species_prefix", how="left")
        o.table(tab, f"{name}_coords")
        rows = [{"matrix": name, "grouping": t, "metric": METRIC[name], "n_genera": r["n"], "n_groups": r["n_groups"],
                 "pseudo_F": r["F"], "p_value": r["p"], "n_permutations": N_PERM, "seed": SEED,
                 "n_negative_eigen": info["n_negative_eigen"], "negative_eigen_share": info["negative_eigen_share"],
                 "n_dropped_nan": info["n_dropped_nan"], "n_missing": info["n_missing"]}
                for t, r in info["permanova"].items()]
        o.table(pd.DataFrame(rows), f"{name}_permanova")
        pct = info["pct_variance"]
        for stem, by, cols, src, idx in ((f"O01_{name}_pcoa_phylum", "phylum", ph_col, ph, ph_idx),
                                         (f"O02_{name}_pcoa_class", "class (phylum hue, shade by class)", cl_col, cl, cl_idx)):
            groups = idx["grp"].reindex(mat.index)
            title = f"{name} PCoA of {len(mat):,} genus representatives, coloured by {by}"
            safe(o, stem, TOPIC, lambda groups=groups, cols=cols, src=src, title=title: _scatter(coords, groups, cols, pct, src, title))
    print(f"ordination: {len(result)} matrices, {len(o.made)} outputs", flush=True)
    result["_made"] = o.made
    return result
