"""Small statistics helpers shared by th2-th4 (tested in test_stats_util.py).

bh_qvalues: Benjamini-Hochberg (1995) adjusted p-values, returned in the order of the input.
  The first version of the survey computed q-values on sorted p-values and assigned them back in the
  unsorted row order (th4_predict.py), and applied the step-up to rows ordered by effect size, not by p (th3_coloss.py).
  Both gave q-values that belong to other rows (Opus review of the branch, 2026-11).

expected_coabsence: Mantel-Haenszel-style expected co-absence of a target family in genomes missing a focal family,
  with the rate of target absence estimated among focal-present genomes within strata. Strata are bins of the number of
  universal families a genome misses (annotation propensity) and, optionally, a taxonomic group such as the class.
  Stratifying by class is needed because a trait shared within a clade (for example 72 of 78 TPS-absent genomes are
  Lecanoromycetes) looks like co-absence if the null ignores the clade.
"""
from __future__ import annotations

import numpy as np
from scipy import stats

MISSING_BINS = [1, 2, 3, 5, 8, 12, 20, 40]


def bh_qvalues(p):
    """Benjamini-Hochberg adjusted p-values in input order; NaN p-values stay NaN and are not counted."""
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return q
    order = np.argsort(pv, kind="stable")
    adj = pv[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def expected_coabsence(f, t, n_missing, group=None, min_denominator=5):
    """Observed and expected number of genomes missing both focal (f) and target (t).

    f, t: boolean arrays (True = absent). n_missing: universal families missed per genome.
    group: optional array of taxon labels; strata become (missing-count bin, group).
    Returns dict(obs, exp, p_two) or None when f has < 3 absences or the expectation is 0.
    Strata with fewer than min_denominator focal-present genomes are skipped (their focal-absent genomes add no expectation).
    """
    f, t = np.asarray(f, bool), np.asarray(t, bool)
    if f.sum() < 3:
        return None
    bins = np.digitize(np.asarray(n_missing), MISSING_BINS)
    if group is not None:
        _, gi = np.unique(np.asarray(group).astype(str), return_inverse=True)
        key = bins * (gi.max() + 1) + gi
    else:
        key = bins
    obs = int((f & t).sum())
    exp = 0.0
    for b in np.unique(key):
        inb = key == b
        denom = int((~f & inb).sum())
        if denom < min_denominator:
            continue
        rate = (t & ~f & inb).sum() / denom
        exp += rate * (f & inb).sum()
    if exp <= 0:
        return None
    p_hi = stats.poisson.sf(obs - 1, exp)
    p_lo = stats.poisson.cdf(obs, exp)
    return dict(obs=obs, exp=float(exp), p_two=float(min(1.0, 2 * min(p_hi, p_lo))))
