"""Figures 5 (RF class probabilities per block and MAG) and 6 (feature importance, LOCO novelty calibration)."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from common import *
P = pd.read_csv('results/rf_predictions.csv'); perf = pd.read_csv('results/rf_performance.csv'); I = pd.read_csv('results/rf_importance.csv'); loco = pd.read_csv('results/rf_loco.csv')
blocks_ = ['struct', 'struct_all', 'cell', 'codon', 'aa', 'trna', 'tf', 'cazy', 'merops', 'pfam_pa', 'pfam_ab', 'all_nonpfam', 'all']
pc = P[(P.target == 'class') & P.LOCUSTAG.isin(AO)]; classes = [c for c in pc.columns if c not in ('LOCUSTAG', 'target', 'block') and pc[c].notna().any()]
M = pc[pc.LOCUSTAG.isin(HQ)].groupby('block')[classes].mean().reindex(blocks_); corder = M.mean().sort_values(ascending=False).index
fig, ax = plt.subplots(1, 3, figsize=(21, 7.5), gridspec_kw=dict(width_ratios=[1.55, .5, 1.1]))
a = ax[0]; im = a.imshow(M[corder].values, cmap='Blues', vmin=0, vmax=1, aspect='auto')
a.set_xticks(range(len(corder))); a.set_xticklabels(corder, rotation=55, ha='right'); ba = perf[perf.target == 'class'].set_index('block').balanced_acc
a.set_yticks(range(len(blocks_))); a.set_yticklabels([f'{b}  (CV bal. acc {ba[b]:.2f})' for b in blocks_])
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[corder].values[i, j]
        if v >= .05: a.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=8, color='white' if v > .55 else '#111111')
a.set_title('A. RF class probabilities for A. obscura (mean of 4 high-quality MAGs), by feature block', loc='left'); plt.colorbar(im, ax=a, fraction=.025, label='probability (fraction of trees)')
ps = P[(P.target == 'subphylum') & P.LOCUSTAG.isin(HQ)].groupby('block')[['Pezizomycotina', 'Saccharomycotina', 'Taphrinomycotina']].mean().reindex(blocks_)
a = ax[1]; im2 = a.imshow(ps.values, cmap='Blues', vmin=0, vmax=1, aspect='auto'); a.set_xticks(range(3)); a.set_xticklabels(ps.columns, rotation=55, ha='right'); a.set_yticks(range(len(blocks_))); a.set_yticklabels([])
for i in range(ps.shape[0]):
    for j in range(3): a.text(j, i, f'{ps.values[i, j]:.2f}', ha='center', va='center', fontsize=8, color='white' if ps.values[i, j] > .55 else '#111111')
a.set_title('B. Subphylum-level RF', loc='left')
# per-MAG top class for 'all' and 'pfam_ab'
a = ax[2]; sub = pc[pc.block == 'all'].set_index('LOCUSTAG').loc[AO, list(corder[:6])]; bottom = np.zeros(len(AO)); cols = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#6250d6']
for c, col in zip(corder[:6], cols):
    a.bar(range(len(AO)), sub[c], bottom=bottom, color=col, label=c, edgecolor='white', linewidth=1.5); bottom += sub[c].values
a.bar(range(len(AO)), 1 - bottom, bottom=bottom, color='#c9c8c0', label='all other classes', edgecolor='white', linewidth=1.5)
a.set_xticks(range(len(AO))); a.set_xticklabels([f"{MAGNAME[i]}\nBUSCO {U.busco_complete[i]:.0f}%" for i in AO], fontsize=8); a.set_ylabel('RF class probability (all blocks)'); a.set_title('C. Per-MAG prediction, all blocks', loc='left'); a.legend(frameon=False, fontsize=8, ncol=2, loc='upper center', bbox_to_anchor=(.5, -.12)); a.grid(False)
for s in ['top', 'right']: a.spines[s].set_visible(False)
plt.tight_layout(); plt.savefig(f'{FIG}/fable_C_ao_fig5_rf_probabilities.png', dpi=130); plt.close()
# ---- Fig 6
fig, ax = plt.subplots(1, 3, figsize=(21, 8))
for a, b, t in [(ax[0], 'all', 'A. Top 25 features, RF on all blocks (Gini importance)'), (ax[1], 'all_nonpfam', 'B. Top 25 features, all blocks except Pfam')]:
    d = I[(I.target == 'class') & (I.block == b)].sort_values('gini', ascending=False).head(25)[::-1]
    a.barh(range(len(d)), d.gini, color='#2a78d6', height=.7); a.set_yticks(range(len(d))); a.set_yticklabels(d.feature.str.slice(0, 48), fontsize=8); a.set_xlabel('mean decrease in Gini impurity'); a.set_title(t, loc='left', fontsize=10)
    a.spines[['top', 'right']].set_visible(False); a.grid(axis='x', color='#e6e5df', lw=.5)
a = ax[2]; l = loco[loco.block == 'all'].sort_values('median_maxprob'); oof = pd.read_csv('results/rf_oof_class_all.csv')
a.barh(range(len(l)), l.median_maxprob, color='#c9c8c0', height=.7, label='left-out class: median max probability')
a.scatter(l.q90_maxprob, range(len(l)), marker='|', c='#55554f', s=80, label='90th percentile')
aomax = pc[(pc.block == 'all') & pc.LOCUSTAG.isin(HQ)][classes].mean().max(); a.axvline(aomax, color='#111111', lw=2, label=f'A. obscura max class probability ({aomax:.2f})')
a.axvline(np.median(oof.maxprob[oof.true == oof.pred]), color='#2a78d6', lw=1.5, ls='--', label='known classes, correct in CV: median')
a.set_yticks(range(len(l))); a.set_yticklabels([f'{c} -> {t} ({f:.0%})' for c, t, f in zip(l.left_out, l.top_assigned, l.top_assigned_frac)], fontsize=8); a.set_xlim(0, 1); a.set_xlabel('max class probability'); a.set_title('C. Leave-one-class-out (label: left-out -> usual assignment)', loc='left', fontsize=9.5)
a.legend(frameon=False, fontsize=8, loc='upper center', bbox_to_anchor=(.5, -.09)); a.spines[['top', 'right']].set_visible(False); a.grid(axis='x', color='#e6e5df', lw=.5)
plt.tight_layout(); plt.savefig(f'{FIG}/fable_C_ao_fig6_rf_importance_loco.png', dpi=130); plt.close()
print(M[corder[:8]].round(2).to_string()); print(ps.round(2).to_string()); print(loco.round(2).to_string())
