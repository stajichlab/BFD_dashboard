"""Group tests and pairwise regressions for the summary pages. Unit: one value per genome."""
import numpy as np
import pandas as pd
from scipy import stats as st

GT_COLS = ["stat", "test", "statistic", "p_value", "effect_size", "n", "n_groups", "anova_F", "eta2"]
REG_COLS = ["x", "y", "n", "slope", "intercept", "r2", "p_value"]


def group_tests(df, stats, group_col="grp"):
    """Kruskal-Wallis across groups for each column in `stats`, with epsilon^2 = H / ((n^2 - 1) / (n + 1)) as
    effect size, plus one-way ANOVA F and eta^2. Missing values are dropped per stat; groups with fewer than
    2 genomes are dropped. Columns that are absent, or have fewer than 2 groups, are skipped. A constant
    column gives NaN effect sizes."""
    rows = []
    for s in stats:
        if s not in df.columns:
            continue
        d = df[[group_col, s]].dropna()
        d = d[np.isfinite(d[s].astype(float))]
        groups = [v[s].astype(float).to_numpy() for _, v in d.groupby(group_col) if len(v) >= 2]
        if len(groups) < 2:
            continue
        n = sum(len(g) for g in groups)
        allv = np.concatenate(groups)
        if np.ptp(allv) == 0:
            rows.append((s, "Kruskal-Wallis", np.nan, np.nan, np.nan, n, len(groups), np.nan, np.nan))
            continue
        h, p = st.kruskal(*groups)
        eps2 = h / ((n ** 2 - 1) / (n + 1))
        f, _ = st.f_oneway(*groups)
        gm = allv.mean()
        sst = ((allv - gm) ** 2).sum()
        ssb = sum(len(g) * (g.mean() - gm) ** 2 for g in groups)
        rows.append((s, "Kruskal-Wallis", h, p, eps2, n, len(groups), f, ssb / sst if sst > 0 else np.nan))
    return pd.DataFrame(rows, columns=GT_COLS)


def pair_regressions(df, pairs, log=True):
    """Ordinary least squares (scipy.stats.linregress) of y on x for each (x, y) pair, on log10 values when
    `log` (rows with a non-positive value in either column are dropped). Pairs with a missing column or fewer
    than 3 usable genomes are skipped.

    These regressions do not control for phylogeny: genomes of related species are not independent, so the
    p-values overstate the evidence and the slopes describe the database, not an evolutionary rate."""
    rows = []
    for x, y in pairs:
        if x not in df.columns or y not in df.columns:
            continue
        d = df[[x, y]].astype(float).replace([np.inf, -np.inf], np.nan).dropna()
        if log:
            d = d[(d[x] > 0) & (d[y] > 0)]
            d = np.log10(d)
        if len(d) < 3:
            continue
        if d[x].nunique() < 2 or d[y].nunique() < 2:
            rows.append((x, y, len(d), np.nan, np.nan, np.nan, np.nan))
            continue
        r = st.linregress(d[x], d[y])
        rows.append((x, y, len(d), r.slope, r.intercept, r.rvalue ** 2, r.pvalue))
    return pd.DataFrame(rows, columns=REG_COLS)
