"""Distinctive gene content: Pfam/CAZy/MEROPS/TF gains, losses, expansions and contractions in Antarctolithica."""
import numpy as np, pandas as pd
from common import *
pf = feat('pfam'); nprot = 10 ** feat('struct').log10_n_genes
tg = U.index[U.is_genus_unit & (U.role == 'train')]
sub = U.loc[tg].subphylum
G = {'Asco': tg, 'Pez': tg[sub == 'Pezizomycotina'], 'Sac': tg[sub == 'Saccharomycotina'], 'Taph': U.index[U.subphylum == 'Taphrinomycotina']}
for c in ['Dothideomycetes', 'Eurotiomycetes', 'Lecanoromycetes']: G[c[:5]] = tg[U.loc[tg].CLASS == c]
def table(M, label):
    P = M > 0; T = pd.DataFrame(index=M.columns)
    T['AO_HQ_present'] = P.loc[HQ].sum(); T['AO_all6_present'] = P.loc[AO].sum(); T['AO_HQ_median_n'] = M.loc[HQ].median()
    for g, ids in G.items():
        T[f'prev_{g}'] = P.loc[ids].mean(); T[f'med_{g}'] = M.loc[ids].median()
    # proteome-normalised expansion score vs Ascomycota genus units (log2 ratio of per-1k-protein rate, +0.5 pseudo)
    r_ao = (M.loc[HQ].median() + 0.5) / nprot.loc[HQ].median(); r_ref = ((M.loc[tg] + 0.5).div(nprot.loc[tg], axis=0))
    T['log2_rate_vs_Asco_median'] = np.log2(r_ao / r_ref.median()); T['pct_rate_in_Asco'] = 100 * (r_ref.lt(r_ao, axis=1)).mean()
    T['pct_rate_in_Pez'] = 100 * (r_ref.loc[G['Pez']].lt(r_ao, axis=1)).mean()
    T.to_csv(f'results/content_{label}.csv'); return T
T = table(pf, 'pfam')
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 400); pd.set_option('display.max_colwidth', 45)
cols = ['AO_HQ_present', 'AO_all6_present', 'AO_HQ_median_n', 'prev_Asco', 'prev_Pez', 'prev_Sac', 'prev_Taph', 'prev_Dothi', 'prev_Eurot', 'prev_Lecan', 'med_Pez', 'med_Taph']
print('n Pfam in >=1 HQ MAG:', (T.AO_HQ_present > 0).sum(), ' in all 4:', (T.AO_HQ_present == 4).sum())
print('\n== Pfam in all 4 HQ MAGs, rare in Ascomycota (prev<10%)'); print(T[(T.AO_HQ_present == 4) & (T.prev_Asco < 0.10)].sort_values('prev_Asco')[cols].round(3).to_string())
print('\n== Pfam near-universal (>=95% of Pezizomycotina AND Taphrinomycotina AND >=90% Asco) absent from all 6 MAGs')
L = T[(T.AO_all6_present == 0) & (T.prev_Pez >= 0.95) & (T.prev_Taph >= 0.90) & (T.prev_Asco >= 0.90)].sort_values('prev_Asco', ascending=False); print(len(L)); print(L[cols].round(3).to_string())
print('\n== Pfam >=95% Pezizomycotina, absent all 6 MAGs (Pezizomycotina-typical families missing):', ((T.AO_all6_present == 0) & (T.prev_Pez >= 0.95)).sum())
L2 = T[(T.AO_all6_present == 0) & (T.prev_Pez >= 0.95)].copy(); L2['also_absent_Taph'] = L2.prev_Taph < 0.2; L2['also_absent_Sac'] = L2.prev_Sac < 0.2
print(L2[cols].round(2).sort_values('prev_Taph').to_string())
L2.to_csv('results/content_pfam_pez_universal_absent.csv')
# universal counts for expected-loss calibration
for g in ['Pez', 'Taph', 'Sac']:
    uni = T.index[T[f'prev_{g}'] >= 0.95]; print(g, 'universal Pfams:', len(uni), 'absent from all6 AO:', int((T.loc[uni].AO_all6_present == 0).sum()), 'present in all4 HQ:', int((T.loc[uni].AO_HQ_present == 4).sum()))
# how many of the Pezizomycotina-universal Pfams does each reference genome lack? (calibration of "absences")
uni = T.index[T.prev_Pez >= 0.95]; miss = (pf.loc[:, uni] == 0).sum(1)
cal = miss.loc[tg].groupby(U.loc[tg].CLASS).median().sort_values(); print('median # of Pez-universal Pfams (n=%d) missing per genome, by class:' % len(uni)); print(cal.to_string())
print('AO MAGs:', miss.loc[AO].to_dict(), ' AO union of 6 MAGs missing:', int((pf.loc[AO, uni].sum() == 0).sum())); print('Taphrinomycotina spp:', miss.loc[G['Taph']].describe().round(1).to_dict())
pd.DataFrame({'missing_pez_universal': miss, 'CLASS': U.CLASS, 'role': U.role, 'subphylum': U.subphylum}).to_csv('results/pez_universal_missing_per_genome.csv')
print('\n== expansions: Pfam with >=5 copies, rate > 99th pct of Ascomycota'); E = T[(T.AO_HQ_median_n >= 5) & (T.pct_rate_in_Asco >= 99)].sort_values('log2_rate_vs_Asco_median', ascending=False)
print(E[['AO_HQ_median_n', 'med_Pez', 'med_Sac', 'med_Taph', 'log2_rate_vs_Asco_median', 'pct_rate_in_Asco']].round(2).head(50).to_string())
print('\n== contractions: Pfam with Asco median >=8 copies and AO rate < 1st pct'); K = T[(T.med_Asco >= 8) & (T.pct_rate_in_Asco <= 1)].sort_values('log2_rate_vs_Asco_median')
print(K[['AO_HQ_median_n', 'med_Pez', 'med_Sac', 'med_Taph', 'log2_rate_vs_Asco_median', 'pct_rate_in_Asco']].round(2).head(40).to_string())
# candidate adaptation sets
cand = {'PKS ketosynthase (melanin/DHN, polyketides)': ['PF00109', 'PF02801'], 'PKS dehydratase PS-DH': ['PF14765'], 'scytalone dehydratase': ['PF02982'], 'multicopper oxidase/laccase': ['PF00394', 'PF07731', 'PF07732'], 'tyrosinase': ['PF00264'],
        'NRPS condensation': ['PF00668'], 'terpene synthase': ['PF03936', 'PF01397', 'PF19086'], 'DMATS prenyltransferase': ['PF11991'], 'photolyase/cryptochrome': ['PF00875', 'PF03441'], 'UV-damage endonuclease UVDE': ['PF03851'],
        'opsin (bacteriorhodopsin-like)': ['PF01036'], 'phytoene synthase (carotenoids)': ['PF00494'], 'trehalose-P synthase': ['PF00982'], 'trehalose phosphatase': ['PF02358'], 'trehalase': ['PF01204'],
        'LEA proteins': ['PF02987', 'PF03168', 'PF03760', 'PF04927'], 'ice-binding (DUF3494)': ['PF11999'], 'cold-shock domain': ['PF00313'], 'catalase': ['PF00199'], 'Cu/Zn + Fe/Mn SOD': ['PF00080', 'PF00081', 'PF02777'],
        'aquaporin/MIP': ['PF00230'], 'HSP70': ['PF00012'], 'HSP20': ['PF00011'], 'HSP90': ['PF00183'], 'fatty acid desaturase': ['PF00487'], 'glycerol-3-P dehydrogenase': ['PF01210', 'PF07479'], 'mannitol dehydrogenase': ['PF01232', 'PF08125'],
        'chitin synthase': ['PF03142', 'PF01644', 'PF08407'], 'hydrophobin': ['PF01185', 'PF06766'], 'HET': ['PF06985'], 'NACHT/NB-ARC': ['PF05729', 'PF00931'], 'cytochrome P450': ['PF00067'], 'MFS transporter': ['PF07690', 'PF00083'], 'ABC transporter': ['PF00005'],
        'Dicer RNase III': ['PF00636'], 'Argonaute Piwi': ['PF02171'], 'RdRP': ['PF05183'], 'C5 DNA methyltransferase': ['PF00145'], 'Spo11': ['PF04406'], 'MAT alpha box': ['PF04769'], 'MAT HMG (generic HMG box)': ['PF00505'],
        'reverse transcriptase (TE)': ['PF00078', 'PF07727'], 'DDE transposase': ['PF03184', 'PF13358', 'PF01359'], 'phytanoyl-CoA dioxygenase (KanJ-like)': ['PF05721'], 'histidine kinase': ['PF00512'], 'flavin monooxygenase': ['PF00743', 'PF13738'],
        'polyketide cyclase SnoaL': ['PF07366', 'PF12680'], 'lichen-type PKS ACP PP-binding': ['PF00550'], 'glycosyl hydrolase 18 chitinase': ['PF00704'], 'cellulase GH5/7/6': ['PF00150', 'PF00840', 'PF01341'], 'LPMO AA9': ['PF03443'], 'cutinase': ['PF01083'], 'pectate lyase': ['PF00544', 'PF03211'],
        'peroxin/ Woronin HEX1': ['PF17820'], 'mycosporine DDGS (DHQ_synthase-like 2-epi-5-epi-valiolone)': ['PF01761'], 'ATP-grasp': ['PF02222', 'PF13535'], 'dynein heavy chain': ['PF03028'], 'kinesin': ['PF00225'], 'septin': ['PF00735']}
acc = pd.Series([c.split('|')[0] for c in pf.columns], index=pf.columns); rows = []
for name, accs in cand.items():
    cols_ = acc.index[acc.isin(accs)]
    if len(cols_) == 0: rows.append(dict(set=name, note='no Pfam column')); continue
    v = pf[cols_].max(1); r = dict(set=name, pfams=','.join(c.split('|')[1] for c in cols_), AO_HQ_median=v.loc[HQ].median(), AO_range6=f'{int(v.loc[AO].min())}-{int(v.loc[AO].max())}')
    for g, ids in G.items(): r['med_' + g] = v.loc[ids].median()
    rate = v / nprot * 1000; r['AO_per1k'] = rate.loc[HQ].median(); r['pct_rate_Asco'] = 100 * (rate.loc[tg] < rate.loc[HQ].median()).mean(); r['pct_rate_Taph'] = 100 * (rate.loc[G['Taph']] < rate.loc[HQ].median()).mean()
    rows.append(r)
Cn = pd.DataFrame(rows); Cn.to_csv('results/content_candidates.csv', index=False); print(Cn.round(2).to_string())
for n in ['cazy', 'merops', 'tf']:
    M = feat(n); t = table(M, n)
    print(f'\n== {n}: families present in AO (HQ median) vs medians'); print(t[(t.AO_HQ_present >= 3)].sort_values('AO_HQ_median_n', ascending=False)[['AO_HQ_median_n', 'med_Pez', 'med_Sac', 'med_Taph', 'pct_rate_in_Asco']].round(1).head(25).to_string())
    print(f'-- {n}: >=90% of Asco, absent all 6 MAGs'); print(t[(t.AO_all6_present == 0) & (t.prev_Asco >= 0.8)][['prev_Asco', 'prev_Pez', 'prev_Sac', 'prev_Taph', 'med_Pez']].round(2).to_string())
