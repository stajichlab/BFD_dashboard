"""Amino-acid, codon and tRNA composition: CLR deviations of A. obscura from subphylum medians (robust z)."""
import numpy as np, pandas as pd
from common import *
tg = U.index[U.is_genus_unit & (U.role == 'train')]
G = {'Pezizomycotina': tg[U.loc[tg].subphylum == 'Pezizomycotina'], 'Saccharomycotina': tg[U.loc[tg].subphylum == 'Saccharomycotina'], 'Taphrinomycotina': U.index[U.subphylum == 'Taphrinomycotina']}
out = []
for name, X in [('aa', clr(feat('aa'), 1e-6)), ('codon', clr(feat('codon'), 1e-5)), ('trna', clr(feat('trna'), 0.5))]:
    raw = feat(name); ao = X.loc[HQ].mean()
    for c in X.columns:
        r = dict(block=name, feature=c, AO_raw=raw.loc[HQ, c].mean(), AO_clr=ao[c])
        for g, ids in G.items():
            v = X.loc[ids, c]; med = v.median(); mad = 1.4826 * (v - med).abs().median(); r[f'z_{g}'] = (ao[c] - med) / mad if mad > 0 else np.nan; r[f'rawmed_{g}'] = raw.loc[ids, c].median()
        out.append(r)
O = pd.DataFrame(out); O.to_csv('results/composition_z.csv', index=False)
pd.set_option('display.width', 220)
for b in ['aa', 'codon', 'trna']:
    o = O[O.block == b]; print(b); print(o.reindex(o.z_Pezizomycotina.abs().sort_values(ascending=False).index).head(10).round(3).to_string(index=False))
t = feat('trna'); print('tRNA genes per MAG:', t.loc[AO].sum(1).to_dict(), 'distinct anticodon classes:', (t.loc[AO] > 0).sum(1).to_dict())
print('reference medians (tRNA genes; distinct classes):', {g: (t.loc[i].sum(1).median(), (t.loc[i] > 0).sum(1).median()) for g, i in G.items()})
print('classes absent in all AO MAGs but present in >=90% Pezizomycotina:', [c for c in t.columns if t.loc[AO, c].sum() == 0 and (t.loc[G['Pezizomycotina'], c] > 0).mean() >= .9])
W = pd.read_parquet('data/feat/wgd_telo.parquet'); print(W.loc[AO].to_string())
