"""Univariate placement of Antarctolithica: percentiles / robust z within groups, Kruskal-Wallis + epsilon^2 across classes."""
import numpy as np, pandas as pd
from scipy import stats
from common import *
X = pd.concat([feat('struct'), feat('cell'), pd.read_parquet('data/feat/wgd_telo.parquet').reindex(U.index).drop(columns=['telomere_top_monomer']).astype(float)], axis=1)
tg = U.index[U.is_genus_unit & (U.role == 'train')]
cls_n = U.loc[tg].CLASS.value_counts(); cls_ok = list(cls_n[cls_n >= 5].index)
groups = {'Ascomycota': tg, 'Pezizomycotina': tg[U.loc[tg].subphylum == 'Pezizomycotina'], 'Saccharomycotina': tg[U.loc[tg].subphylum == 'Saccharomycotina'],
          'Taphrinomycotina': U.index[(U.subphylum == 'Taphrinomycotina')]}   # Taphrinomycotina: species units (only 8 genera)
for c in ['Dothideomycetes', 'Eurotiomycetes', 'Lecanoromycetes', 'Sordariomycetes', 'Leotiomycetes', 'Arthoniomycetes', 'Lichinomycetes', 'Taphrinomycetes', 'Dipodascomycetes', 'Lipomycetes']:
    groups[c] = (tg if c != 'Taphrinomycetes' else U.index)[(U.loc[tg if c != 'Taphrinomycetes' else U.index].CLASS == c)]
rows = []
for f in X.columns:
    ao_hq = X.loc[HQ, f].median(); ao_all = X.loc[AO, f]
    r = dict(feature=f, AO_HQ_median=ao_hq, AO_min6=ao_all.min(), AO_max6=ao_all.max())
    samp = [X.loc[tg[U.loc[tg].CLASS == c], f].dropna().values for c in cls_ok]; samp = [s for s in samp if len(s) >= 3]
    if np.isfinite(ao_hq) and len(samp) > 2:
        H, p = stats.kruskal(*samp); n = sum(map(len, samp)); k = len(samp)
        r.update(KW_H=H, KW_p=p, epsilon2=(H - k + 1) / (n - k))
    for g, ids in groups.items():
        v = X.loc[ids, f].dropna().values
        if len(v) < 3 or not np.isfinite(ao_hq): continue
        med = np.median(v); mad = 1.4826 * np.median(np.abs(v - med))
        r[f'pct_{g}'] = 100 * (v < ao_hq).mean() + 50 * (v == ao_hq).mean(); r[f'z_{g}'] = (ao_hq - med) / mad if mad > 0 else np.nan; r[f'med_{g}'] = med
    rows.append(r)
R = pd.DataFrame(rows); R.to_csv('results/univariate.csv', index=False)
pd.set_option('display.width', 280); pd.set_option('display.max_rows', 200)
print(R[['feature', 'AO_HQ_median', 'AO_min6', 'AO_max6', 'epsilon2', 'pct_Ascomycota', 'pct_Pezizomycotina', 'pct_Saccharomycotina', 'pct_Taphrinomycotina', 'med_Pezizomycotina', 'med_Saccharomycotina', 'med_Taphrinomycotina', 'med_Dothideomycetes','med_Lecanoromycetes']].round(3).to_string())
