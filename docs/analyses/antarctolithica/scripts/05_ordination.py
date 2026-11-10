"""PCoA + PERMANOVA per feature block; position of Antarctolithica relative to classes."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from common import *
from lib.ordination import pcoa, permanova, _pseudo_f
rng = np.random.default_rng(0)
B = blocks()
units = U.index[(U.is_genus_unit & (U.role != 'AO')) | (U.role == 'AO')]
train_g = U.index[U.is_genus_unit & (U.role == 'train')]
cls_n = U.loc[train_g].CLASS.value_counts(); cls_ok = cls_n[cls_n >= 5].index
perm_ids = []
for c in cls_ok:
    ids = list(train_g[U.loc[train_g].CLASS == c]); perm_ids += list(rng.choice(ids, min(40, len(ids)), replace=False))
perm_ids = pd.Index(perm_ids)
rows, pos_rows, nn_rows, cons_rows = [], [], [], []
fig, axes = plt.subplots(3, 4, figsize=(20, 14.5)); axes = axes.ravel()
order = ['Basidiomycota (outgroup)', 'other Pezizomycotina', 'Sordariomycetes', 'Dothideomycetes', 'Eurotiomycetes', 'Lecanoromycetes', 'Leotiomycetes', 'Saccharomycotina', 'Taphrinomycotina', 'unplaced Ascomycota']
for k, (name, (X, metric, desc)) in enumerate(B.items()):
    X = X.loc[units]
    if metric == 'z-euclidean':
        X = zscore(X, train_g).clip(-8, 8); D = squareform(pdist(X.values, 'euclidean'))
    else:
        D = squareform(pdist(X.values, metric))
    D = pd.DataFrame(D, index=units, columns=units)
    o = pcoa(D, 4); C = o.coords
    # PERMANOVA by class (genus units, classes >=5 genera, <=40 genera per class)
    Dp = D.loc[perm_ids, perm_ids]; lab = U.loc[perm_ids].CLASS
    t = permanova(Dp, lab, 999, seed=1)
    a, n = lab.nunique(), len(lab); r2 = t.pseudo_f * (a - 1) / (t.pseudo_f * (a - 1) + (n - a))
    lab2 = U.loc[perm_ids].subphylum; t2 = permanova(Dp, lab2, 999, seed=1); a2 = lab2.nunique(); r2s = t2.pseudo_f * (a2 - 1) / (t2.pseudo_f * (a2 - 1) + (n - a2))
    rows.append(dict(block=name, description=desc, n_features=X.shape[1], n_genera=n, n_classes=a, pseudoF_class=t.pseudo_f, R2_class=r2, p_class=t.p_value,
                     pseudoF_subphylum=t2.pseudo_f, R2_subphylum=r2s, p_subphylum=t2.p_value, dispersion_ratio=t.dispersion_ratio, PCo1_pct=o.pct_variance[0], PCo2_pct=o.pct_variance[1]))
    # AO position: distance from HQ-MAG mean profile to class centroids (Gower identity) and members
    d_ao = D.loc[HQ, train_g].mean(0)                       # mean distance of 4 HQ MAGs to each reference genus
    for c in cls_n[cls_n >= 2].index:
        ids = train_g[U.loc[train_g].CLASS == c]; Dc = D.loc[ids, ids].values
        within = (Dc ** 2).sum() / (2 * len(ids) ** 2)
        d2cent = max((D.loc[HQ, ids].values ** 2).mean() - within, 0) ** .5
        own = np.sqrt(np.clip((Dc ** 2).mean(1) - within, 0, None))       # members' distance to own centroid
        pos_rows.append(dict(block=name, CLASS=c, subphylum=U.loc[ids[0]].subphylum, n_genera=len(ids), dist_to_centroid=d2cent, median_member_to_centroid=np.median(own),
                             ratio=d2cent / np.median(own), pct_members_farther=100 * (own > d2cent).mean(), mean_dist_to_members=d_ao[ids].mean(), min_dist_to_members=d_ao[ids].min()))
    nn = d_ao.sort_values().head(15)
    for r, (i, v) in enumerate(nn.items()): nn_rows.append(dict(block=name, rank=r + 1, dist=v, LOCUSTAG=i, CLASS=U.CLASS[i], ORDER=U.ORDER_[i], SPECIES=U.SPECIES[i], genome_mb=U.total_length_bp[i] / 1e6))
    # pseudo-replicate consistency: AO pairwise distances vs reference nearest-neighbour distances
    Dt = D.loc[train_g, train_g].values.copy(); np.fill_diagonal(Dt, np.inf)
    hqd = D.loc[HQ, HQ].values[np.triu_indices(4, 1)]; alld = D.loc[AO, AO].values[np.triu_indices(len(AO), 1)]
    cons_rows.append(dict(block=name, AO_HQ_pairwise_median=np.median(hqd), AO_all6_pairwise_max=alld.max(), ref_nearest_neighbour_median=np.median(Dt.min(1)), AO_to_nearest_ref=d_ao.min(), ref_all_pairs_median=np.median(Dt[np.isfinite(Dt)])))
    ax = axes[k]
    for gname in order:
        ids = units[U.loc[units].pg == gname]
        ax.scatter(C.loc[ids, 'PCo1'], C.loc[ids, 'PCo2'], s=14 if gname != 'unplaced Ascomycota' else 26, c=COL[gname], alpha=.75, linewidths=.3, edgecolors='white', marker='o' if gname != 'unplaced Ascomycota' else 'D', label=gname, zorder=2)
    ax.scatter(C.loc[AO, 'PCo1'], C.loc[AO, 'PCo2'], s=190, marker='*', c=COL['Antarctolithica'], edgecolors='white', linewidths=1.2, label='Antarctolithica obscura (6 MAGs)', zorder=5)
    cx, cy = C.loc[HQ, 'PCo1'].mean(), C.loc[HQ, 'PCo2'].mean()
    ax.annotate('A. obscura', (cx, cy), xytext=(12, 12), textcoords='offset points', fontsize=10, fontweight='bold', arrowprops=dict(arrowstyle='-', lw=.6))
    ax.set_title(f'{name}: {desc}\nPERMANOVA class R²={r2:.2f}, p={t.p_value:.3f}', fontsize=9.5, loc='left')
    ax.set_xlabel(f'PCo1 ({o.pct_variance[0]:.1f}%)'); ax.set_ylabel(f'PCo2 ({o.pct_variance[1]:.1f}%)')
    ax.grid(color='#e6e5df', lw=.5, zorder=0); [s.set_visible(False) for s in (ax.spines['top'], ax.spines['right'])]
    C.join(U[['CLASS', 'pg', 'SPECIES']]).to_csv(f'results/pcoa_{name}.csv')
ax = axes[-1]; ax.axis('off'); h, l = axes[0].get_legend_handles_labels(); ax.legend(h, l, loc='center', frameon=False, fontsize=11, markerscale=1.6)
plt.tight_layout(); plt.savefig(f'{FIG}/fable_C_ao_fig2_pcoa_blocks.png', dpi=130); plt.close()
pd.DataFrame(rows).to_csv('results/permanova.csv', index=False)
P = pd.DataFrame(pos_rows); P.to_csv('results/ao_class_position.csv', index=False)
pd.DataFrame(nn_rows).to_csv('results/ao_nearest_neighbours.csv', index=False)
pd.DataFrame(cons_rows).to_csv('results/ao_consistency.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(pd.DataFrame(rows).drop(columns='description').round(3).to_string())
print(pd.DataFrame(cons_rows).round(3).to_string())
for b in B:
    p = P[(P.block == b) & (P.n_genera >= 2)].sort_values('dist_to_centroid').head(6)
    print(b); print(p[['CLASS', 'n_genera', 'dist_to_centroid', 'median_member_to_centroid', 'ratio', 'pct_members_farther', 'min_dist_to_members']].round(3).to_string(index=False))
