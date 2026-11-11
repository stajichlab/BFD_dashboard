"""Figures 1 (MAG QC), 3 (compaction), 4 (univariate percentiles), 7 (SwissProt placement signal)."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from common import *
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': '#e6e5df', 'grid.linewidth': .5, 'axes.axisbelow': True})
L = pd.read_parquet('data/labels/genome_labels.parquet'); A8 = L[L.SPECIES.str.contains('Antarctolith', na=False)].sort_values('busco_complete', ascending=False)
# ---------------- Fig 1: MAG QC
fig, ax = plt.subplots(2, 2, figsize=(14, 10))
a = ax[0, 0]; y = np.arange(len(A8))[::-1]; left = np.zeros(len(A8))
for col, c, lab in [('busco_single', '#2a78d6', 'complete, single'), ('busco_duplicated', '#eda100', 'complete, duplicated'), ('busco_fragmented', '#1baf7a', 'fragmented'), ('busco_missing', '#c9c8c0', 'missing')]:
    a.barh(y, A8[col], left=left, color=c, label=lab, edgecolor='white', linewidth=1.5, height=.7); left += A8[col].values
a.set_yticks(y); a.set_yticklabels([f"{s}{'' if q else '  (not annotated)'}" for s, q in zip(A8.STRAIN, A8.LOCUSTAG.isin(AO))]); a.set_xlabel('% of 1,122 fungi_odb12 BUSCOs (genome mode)'); a.set_title('A. BUSCO completeness of the 8 MAGs', loc='left'); a.legend(ncol=4, fontsize=8, frameon=False, loc='upper center', bbox_to_anchor=(.5, -.12)); a.grid(False)
for yy, v in zip(y, A8.busco_complete): a.text(1, yy, f'{v:.1f}', va='center', color='white', fontsize=9)
a = ax[0, 1]; x = A8.busco_complete.values; s = A8.total_length_bp.values / 1e6; b, c0 = np.polyfit(x, s, 1); xx = np.linspace(60, 100, 10)
a.plot(xx, b * xx + c0, color='#8a8a84', lw=1.2, ls='--', label=f'OLS: {b:.3f} Mb per BUSCO %; at 100% = {b*100+c0:.1f} Mb'); a.scatter(x, s, s=70, c='#111111', edgecolors='white', zorder=3)
for i, r in A8.iterrows(): a.annotate(r.STRAIN.replace('Mars-', ''), (r.busco_complete, r.total_length_bp / 1e6), xytext=(5, 4), textcoords='offset points', fontsize=8)
a.set_xlabel('BUSCO complete (%)'); a.set_ylabel('assembly size (Mb)'); a.set_title('B. MAG size scales with completeness (incompleteness, not biology)', loc='left'); a.legend(frameon=False, fontsize=9)
C = pd.read_csv('results/mag_contigs.csv'); a = ax[1, 0]
for flag, col, lab, z in [(False, '#2a78d6', 'contig', 2), (True, '#e34948', 'flagged (|GC dev| > 6 pts, >=50% prokaryote-like hits, or >=10 kb gene-free)', 3)]:
    d = C[C.flag == flag]; a.scatter(d.length / 1e3, d.gc, s=10 if not flag else 22, c=col, alpha=.5 if not flag else .9, label=lab, zorder=z, linewidths=0)
a.set_xscale('log'); a.set_xlabel('contig length (kb, log scale)'); a.set_ylabel('contig GC (%)'); a.set_title('C. Contig GC vs length, 6 annotated MAGs (1,575 contigs)', loc='left'); a.legend(frameon=False, fontsize=8, loc='lower right')
a = ax[1, 1]; S = pd.read_csv('results/mag_contam_summary.csv').set_index('strain').loc[[s for s in A8.STRAIN if s in set(C.strain)]]; yy = np.arange(len(S))[::-1]
a.barh(yy, S.flagged_pct_bp, color='#e34948', height=.6); a.set_yticks(yy); a.set_yticklabels(S.index); a.set_xlabel('% of assembly bp on flagged contigs'); a.set_title('D. Putative foreign/atypical sequence per MAG', loc='left'); a.grid(axis='y', visible=False)
for y_, v, g in zip(yy, S.flagged_pct_bp, S.genes_on_flagged): a.text(v + .05, y_, f'{v:.2f}%  ({int(g)} genes)', va='center', fontsize=9)
plt.tight_layout(); plt.savefig(f'{FIG}/fable_C_ao_fig1_mag_qc.png', dpi=130); plt.close()
print('size at 100% BUSCO (OLS over 8 MAGs):', b * 100 + c0, ' HQ4 mean size', A8.head(4).total_length_bp.mean() / 1e6)
# ---------------- Fig 3: compaction
s = feat('struct'); g = U[(U.is_genus_unit) | (U.role == 'AO')].index
order = ['Basidiomycota (outgroup)', 'other Pezizomycotina', 'Sordariomycetes', 'Dothideomycetes', 'Eurotiomycetes', 'Lecanoromycetes', 'Leotiomycetes', 'Saccharomycotina', 'unplaced Ascomycota']
tx = U.index[U.subphylum == 'Taphrinomycotina']
def scat(a, xc, yc, logx=False, logy=False):
    X = s[xc] if xc != 'genome_mb' else 10 ** s.log10_genome_mb
    for gn in order:
        ids = g[U.loc[g].pg == gn]; a.scatter(X[ids], s.loc[ids, yc], s=12, c=COL[gn], alpha=.7, linewidths=0, label=gn)
    a.scatter(X[tx], s.loc[tx, yc], s=34, c=COL['Taphrinomycotina'], edgecolors='white', linewidths=.6, label='Taphrinomycotina (all species)', zorder=4)
    a.scatter(X[AO], s.loc[AO, yc], s=170, marker='*', c='#111111', edgecolors='white', linewidths=1, label='A. obscura (6 MAGs)', zorder=5)
    a.annotate('A. obscura', (X[HQ].mean(), s.loc[HQ, yc].mean()), xytext=(14, 10), textcoords='offset points', fontweight='bold', arrowprops=dict(arrowstyle='-', lw=.6))
    if logx: a.set_xscale('log')
    if logy: a.set_yscale('log')
fig, ax = plt.subplots(2, 3, figsize=(19, 11)); ax = ax.ravel()
for a, (xc, yc, lx, ly, xl, yl, t) in zip(ax, [('genome_mb', 'genes_per_mb', True, False, 'assembly size (Mb, log)', 'protein-coding genes per Mb', 'A. Gene density'),
        ('genome_mb', 'coding_fraction', True, False, 'assembly size (Mb, log)', 'coding fraction of assembly', 'B. Coding fraction'),
        ('genome_mb', 'median_intergenic', True, True, 'assembly size (Mb, log)', 'median intergenic distance (bp, log)', 'C. Intergenic distance'),
        ('introns_per_transcript', 'median_intron_len', False, True, 'introns per transcript', 'median intron length (bp, log)', 'D. Intron number and length'),
        ('median_protein_len', 'frac_intronless', False, False, 'median protein length (aa)', 'fraction of intronless genes', 'E. Protein length vs intron-less genes'),
        ('gc_pct', 'intron_gc_minus_genome_gc', False, False, 'genome GC (%)', 'intron GC minus genome GC (pts)', 'F. GC')]):
    scat(a, xc, yc, lx, ly); a.set_xlabel(xl); a.set_ylabel(yl); a.set_title(t, loc='left')
h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, ncol=6, frameon=False, loc='lower center', fontsize=10, markerscale=1.5); plt.tight_layout(rect=(0, .04, 1, 1)); plt.savefig(f'{FIG}/fable_C_ao_fig3_compaction.png', dpi=130); plt.close()
# ---------------- Fig 4: percentile plot
R = pd.read_csv('results/univariate.csv'); R = R[~R.feature.isin(['identical_groups', 'telomere_tracts', 'telomere_scaffolds', 'sig_hets_per1k', 'sig_sesa_per1k', 'sig_te_per1k', 'partial_fraction', 'n_trna', 'q10_protein_len', 'mean_intergenic', 'frac_tm2plus', 'frac_best_hit_yeast', 'frac_best_hit_neucr_emeni', 'frac_with_swissprot'])]
R = R.sort_values('pct_Ascomycota').reset_index(drop=True); fig, a = plt.subplots(figsize=(11, 17)); y = np.arange(len(R))
a.axvspan(0, 2.5, color='#f3e3e0', zorder=0); a.axvspan(97.5, 100, color='#f3e3e0', zorder=0)
for col, c, m, lab in [('pct_Pezizomycotina', '#eb6834', 'o', 'within Pezizomycotina (genera)'), ('pct_Saccharomycotina', '#008300', 's', 'within Saccharomycotina (genera)'), ('pct_Taphrinomycotina', '#6250d6', '^', 'within Taphrinomycotina (25 species)')]:
    a.scatter(R[col], y, c=c, marker=m, s=34, label=lab, zorder=3, edgecolors='white', linewidths=.5)
a.scatter(R.pct_Ascomycota, y, c='#111111', marker='|', s=160, label='within all Ascomycota (genera)', zorder=4)
a.set_yticks(y); a.set_yticklabels([f"{f}  [{v:.3g}]" for f, v in zip(R.feature, R.AO_HQ_median)], fontsize=8); a.set_xlabel('percentile of the A. obscura value (median of 4 high-quality MAGs) within each reference group'); a.set_xlim(-1, 101)
a.set_title('Where A. obscura falls for each scalar feature (value in brackets); shaded = outer 2.5% tails', loc='left', fontsize=10); a.legend(frameon=False, loc='upper center', bbox_to_anchor=(.5, -.035), ncol=2); a.grid(axis='y', lw=.3)
plt.tight_layout(); plt.savefig(f'{FIG}/fable_C_ao_fig4_percentiles.png', dpi=130); plt.close()
# ---------------- Fig 7: swissprot signal
F = pd.read_parquet('data/feat/sprot_scalars.parquet'); fig, ax = plt.subplots(1, 2, figsize=(15, 6.5))
for a, (xc, yc, xl, yl, t) in zip(ax, [('frac_SCHPO', 'pid_NEUCR_minus_SCHPO', 'fraction of best hits to S. pombe proteins', 'mean % identity to N. crassa hits minus S. pombe hits', 'A. Best-hit organism share vs relative identity'),
                                      ('pid_SCHPO_minus_YEAST', 'pid_NEUCR_minus_YEAST', 'mean % identity: S. pombe hits minus S. cerevisiae hits', 'mean % identity: N. crassa hits minus S. cerevisiae hits', 'B. Relative identity to three anchors')]):
    for gn in order[1:]:
        ids = g[U.loc[g].pg == gn]; a.scatter(F.loc[ids, xc], F.loc[ids, yc], s=12, c=COL[gn], alpha=.7, linewidths=0, label=gn)
    t2 = [i for i in tx if U.CLASS[i] != 'Schizosaccharomycetes']
    a.scatter(F.loc[t2, xc], F.loc[t2, yc], s=34, c=COL['Taphrinomycotina'], edgecolors='white', linewidths=.6, label='Taphrinomycotina (excl. Schizosaccharomyces)', zorder=4)
    a.scatter(F.loc[AO, xc], F.loc[AO, yc], s=170, marker='*', c='#111111', edgecolors='white', linewidths=1, label='A. obscura (6 MAGs)', zorder=5)
    a.annotate('A. obscura', (F.loc[HQ, xc].mean(), F.loc[HQ, yc].mean()), xytext=(14, 12), textcoords='offset points', fontweight='bold', arrowprops=dict(arrowstyle='-', lw=.6)); a.set_xlabel(xl); a.set_ylabel(yl); a.set_title(t, loc='left')
h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, ncol=5, frameon=False, loc='lower center', fontsize=9, markerscale=1.5); plt.tight_layout(rect=(0, .07, 1, 1)); plt.savefig(f'{FIG}/fable_C_ao_fig7_swissprot.png', dpi=130); plt.close()
