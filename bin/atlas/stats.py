"""Statistics helpers for the BFD summary atlas.

Units of replication are chosen by the caller (genus medians, families, pairs).
All functions return plain dicts so the renderer can print them into Markdown.
"""
import numpy as np
import pandas as pd
from scipy import stats

try:
    import statsmodels.formula.api as smf
    HAVE_SM = True
except Exception:  # statsmodels missing in /usr/bin/python3.12
    HAVE_SM = False

RNG = np.random.default_rng(20261108)


def boot_ci(x, func=np.median, n=2000, alpha=0.05):
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if len(x) < 2:
        return (np.nan, np.nan)
    idx = RNG.integers(0, len(x), size=(n, len(x)))
    vals = func(x[idx], axis=1)
    return tuple(np.quantile(vals, [alpha / 2, 1 - alpha / 2]))


def boot_diff_ci(a, b, func=np.median, n=2000):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return (np.nan, np.nan)
    ia = RNG.integers(0, len(a), size=(n, len(a))); ib = RNG.integers(0, len(b), size=(n, len(b)))
    d = func(a[ia], axis=1) - func(b[ib], axis=1)
    return tuple(np.quantile(d, [0.025, 0.975]))


def cliffs_delta(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) == 0 or len(b) == 0:
        return np.nan
    # U statistic based
    u = stats.mannwhitneyu(a, b, alternative='two-sided').statistic
    return 2 * u / (len(a) * len(b)) - 1


def two_group(a, b, label_a='yeast', label_b='filamentous'):
    """Mann-Whitney on two independent samples + median difference with bootstrap CI + Cliff's delta."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    out = {'n_' + label_a: int(len(a)), 'n_' + label_b: int(len(b)),
           'median_' + label_a: float(np.median(a)) if len(a) else np.nan,
           'median_' + label_b: float(np.median(b)) if len(b) else np.nan}
    if len(a) >= 2 and len(b) >= 2:
        mw = stats.mannwhitneyu(a, b, alternative='two-sided')
        lo, hi = boot_diff_ci(a, b)
        out.update({'median_diff': float(np.median(a) - np.median(b)), 'diff_ci_lo': lo, 'diff_ci_hi': hi,
                    'cliffs_delta': cliffs_delta(a, b), 'mw_p': float(mw.pvalue)})
    return out


def sign_test(x):
    """Fraction below zero with exact binomial test (two-sided)."""
    x = np.asarray(x, float); x = x[~np.isnan(x) & (x != 0)]
    n = len(x); k = int((x < 0).sum())
    if n == 0:
        return {'n': 0}
    p = stats.binomtest(k, n, 0.5).pvalue
    return {'n': n, 'n_below_zero': k, 'frac_below_zero': k / n, 'binom_p': float(p)}


def kruskal(groups):
    """groups: dict name -> array. Screening statistic only."""
    arrs = [np.asarray(v, float) for v in groups.values()]
    arrs = [a[~np.isnan(a)] for a in arrs if len(a) >= 2]
    if len(arrs) < 2:
        return {'kw_H': np.nan, 'kw_p': np.nan}
    h, p = stats.kruskal(*arrs)
    return {'kw_H': float(h), 'kw_p': float(p), 'k_groups': len(arrs)}


def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = ~np.isnan(x) & ~np.isnan(y)
    if m.sum() < 3:
        return {'n': int(m.sum()), 'rho': np.nan, 'p': np.nan}
    r = stats.spearmanr(x[m], y[m])
    return {'n': int(m.sum()), 'rho': float(r.statistic), 'p': float(r.pvalue)}


def loglog_slope(x, y):
    """OLS slope of log10(y) on log10(x) with 95% CI."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = (x > 0) & (y > 0) & ~np.isnan(x) & ~np.isnan(y)
    if m.sum() < 3:
        return {'n': int(m.sum()), 'slope': np.nan}
    lx, ly = np.log10(x[m]), np.log10(y[m])
    r = stats.linregress(lx, ly)
    t = stats.t.ppf(0.975, m.sum() - 2)
    return {'n': int(m.sum()), 'slope': float(r.slope), 'slope_ci_lo': float(r.slope - t * r.stderr),
            'slope_ci_hi': float(r.slope + t * r.stderr), 'intercept': float(r.intercept), 'r2': float(r.rvalue ** 2)}


def mixed_model(df, response, fixed, group='CLASS', vc='ORDER', extra_cols=None):
    """Nested random intercepts (CLASS, ORDER within CLASS) with statsmodels MixedLM.
    df rows are the unit (genus medians). Returns coefficient table for the fixed effects.
    """
    if not HAVE_SM:
        return {'available': False, 'note': 'statsmodels not available in this interpreter'}
    cols = [response, group, vc] + [f for f in fixed if f in df.columns] + [c for c in (extra_cols or []) if c in df.columns]
    d = df.dropna(subset=cols).copy()
    if d[group].nunique() < 2 or len(d) < 30:
        return {'available': False, 'note': f'too few rows ({len(d)}) or groups'}
    formula = f"{response} ~ " + ' + '.join(fixed)
    try:
        md = smf.mixedlm(formula, d, groups=d[group], vc_formula={vc: f'0 + C({vc})'}, re_formula='1')
        fit = md.fit(reml=True, method=['lbfgs'], maxiter=200)
        params = fit.params; ci = fit.conf_int(); pv = fit.pvalues
        rows = []
        for k in params.index:
            if k in ('Group Var',) or k.endswith(' Var'):
                continue
            rows.append({'term': k, 'estimate': float(params[k]), 'ci_lo': float(ci.loc[k, 0]),
                         'ci_hi': float(ci.loc[k, 1]), 'p': float(pv[k])})
        return {'available': True, 'n': int(len(d)), 'n_groups': int(d[group].nunique()),
                'converged': bool(fit.converged), 'terms': rows, 'formula': formula + f' + (1|{group}/{vc})'}
    except Exception as e:
        return {'available': False, 'note': f'fit failed: {e}'}


def permanova(D, labels, n_perm=999, strata=None):
    """Pseudo-F PERMANOVA on a distance matrix (numpy). strata: restrict permutations within strata."""
    D = np.asarray(D, float); n = D.shape[0]
    labels = np.asarray(labels);
    if strata is not None:
        strata = np.asarray(strata)
    A = -0.5 * D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    G = J @ A @ J
    sst = np.trace(G)

    def pseudo_f(lab):
        groups = pd.unique(lab)
        H = np.zeros((n, n))
        for g in groups:
            idx = np.where(lab == g)[0]
            H[np.ix_(idx, idx)] = 1.0 / len(idx)
        ssw = np.trace(G) - np.trace(H @ G)
        ssa = sst - ssw
        a = len(groups)
        return (ssa / (a - 1)) / (ssw / (n - a)), ssa / sst

    f0, r2 = pseudo_f(labels)
    cnt = 0
    for _ in range(n_perm):
        if strata is None:
            perm = RNG.permutation(labels)
        else:
            perm = labels.copy()
            for s in pd.unique(strata):
                idx = np.where(strata == s)[0]
                perm[idx] = RNG.permutation(labels[idx])
        f, _ = pseudo_f(perm)
        if f >= f0:
            cnt += 1
    return {'pseudo_F': float(f0), 'R2': float(r2), 'p': (cnt + 1) / (n_perm + 1), 'n': int(n),
            'n_groups': int(len(pd.unique(labels))), 'n_perm': n_perm}


def fmt(x, nd=3):
    if x is None:
        return 'NA'
    if isinstance(x, (int, np.integer)):
        return f'{int(x):,}'
    if isinstance(x, float):
        if np.isnan(x):
            return 'NA'
        if abs(x) < 1e-3 and x != 0:
            return f'{x:.2e}'
        return f'{x:.{nd}f}'
    return str(x)
