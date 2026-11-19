"""Ecology (FungalTraits) in the summary ordinations: a separate step that runs after make_figures.

Reads the same per-genome tables as make_figures read-only, joins lifestyle and growth form from the DB view
v_genome_labels, and writes figures and tables under its own output directory. It changes no shared file.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import summary_ordination as so  # noqa: E402


def permute_within(labels: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Return a copy of labels with the values shuffled inside each stratum only."""
    out = np.array(labels, copy=True)
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]
        if len(idx) > 1:
            out[idx] = out[rng.permutation(idx)]
    return out


def _na(n, k, caveat, n_strata=0, n_drop=0):
    return {"F": np.nan, "R2": np.nan, "R2_adj": np.nan, "p": np.nan, "n": int(n), "n_groups": int(k),
            "n_strata_used": int(n_strata), "n_dropped_no_stratum": int(n_drop), "n_permutations_run": 0, "caveat": caveat}


def stratified_permanova(dist: pd.DataFrame, groups: pd.Series, strata: pd.Series | None = None,
                         n_permutations: int = so.N_PERM, seed: int = so.SEED) -> dict:
    """PERMANOVA of `groups`. With `strata`, labels are permuted only inside a stratum. Rows with a NULL group
    are left out; with strata, rows with a NULL stratum are also left out and counted. Strata with one row do
    not move. Degenerate designs give NaN and a caveat."""
    g = groups.reindex(dist.index)
    keep = g.notna()
    n_drop = 0
    st = None
    if strata is not None:
        st = strata.reindex(dist.index)
        n_drop = int((keep & st.isna()).sum())
        keep = keep & st.notna()
    n = int(keep.sum())
    lab = g[keep].to_numpy()
    k = len(np.unique(lab))
    if n < 3 or k < 2:
        return _na(n, k, "not tested: fewer than 3 genera or fewer than 2 groups", n_drop=n_drop)
    sv = None
    n_used = 0
    if st is not None:
        sv = st[keep].to_numpy()
        sizes = pd.Series(sv).value_counts()
        n_used = int((sizes >= 2).sum())
        if n_used < 2:
            return _na(n, k, "not tested: fewer than 2 strata with at least 2 genera", n_used, n_drop)
    d = dist.loc[keep, keep].to_numpy()
    f = float(so.lib_ordination._pseudo_f(d, lab))
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_permutations):
        perm = permute_within(lab, sv, rng) if sv is not None else lab[rng.permutation(n)]
        if so.lib_ordination._pseudo_f(d, perm) >= f:
            ge += 1
    p = (1 + ge) / (1 + n_permutations)
    r2 = f * (k - 1) / (f * (k - 1) + (n - k)) if np.isfinite(f) and n > k else np.nan
    caveat = ""
    if np.isinf(f):
        caveat = "pseudo-F is infinite: no within-group variation"
    return {"F": f, "R2": float(r2), "R2_adj": so._r2_adj(r2, n, k), "p": float(p), "n": n, "n_groups": k,
            "n_strata_used": n_used, "n_dropped_no_stratum": n_drop, "n_permutations_run": n_permutations, "caveat": caveat}


def bh_adjust(p: pd.Series) -> pd.Series:
    """Benjamini-Hochberg q values. NaN stays NaN and is not counted."""
    q = pd.Series(np.nan, index=p.index, dtype=float)
    ok = p.notna()
    v = p[ok].to_numpy(dtype=float)
    m = len(v)
    if m == 0:
        return q
    order = np.argsort(v)
    adj = v[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def pool_groups(labels: pd.Series, min_group: int, other: str = "other") -> pd.Series:
    """Groups with fewer than min_group members become 'other (k groups)'. NULL stays NULL."""
    counts = labels.value_counts()
    small = set(counts[counts < min_group].index)
    if not small:
        return labels.copy()
    name = f"{other} ({len(small)} groups)"
    return labels.where(~labels.isin(small), name)
