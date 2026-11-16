"""Compare a rerun (AO_BASE/results) with a reference results directory (default: the committed laptop results next to this script).
Usage: python 18_compare_runs.py [reference_results_dir]   (run from AO_BASE)
Writes results/rerun_vs_reference.csv (section, item, reference, rerun, diff) and prints it."""
import os, sys, numpy as np, pandas as pd
from pathlib import Path

REF = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / 'results'
NEW = Path('results')
HQ = ['FF1FC3B3', 'F46C92EF', 'F7B4B958', 'F7E47E53']
rows = []
def add(section, item, a, b):
    d = (b - a) if isinstance(a, (int, float, np.floating, np.integer)) and isinstance(b, (int, float, np.floating, np.integer)) else ('same' if a == b else 'DIFFERENT')
    rows.append(dict(section=section, item=item, reference=a, rerun=b, diff=d))
def both(name, **kw):
    return pd.read_csv(REF / name, **kw), pd.read_csv(NEW / name, **kw)

# --- set sizes
a, b = both('taxa.tsv', sep='\t'); add('set', 'phylogeny taxa', len(a), len(b)); add('set', 'phylogeny taxa identical LOCUSTAG set', True, set(a.LOCUSTAG) == set(b.LOCUSTAG))
a, b = both('pcoa_struct.csv', index_col=0); add('set', 'genus units + AO MAGs (ordination)', len(a), len(b)); add('set', 'units in common', len(a), len(a.index.intersection(b.index)))
a, b = both('rf_oof_class_all.csv'); add('set', 'RF training species', len(a), len(b)); add('set', 'RF classes', a.true.nunique(), b.true.nunique())
# --- PERMANOVA
a, b = both('permanova.csv'); a, b = a.set_index('block'), b.set_index('block')
for blk in a.index:
    for c in ['R2_class', 'R2_subphylum', 'n_genera']: add('PERMANOVA', f'{blk} {c}', a.at[blk, c], b.at[blk, c])
# --- RF performance
a, b = both('rf_performance.csv'); a, b = a.set_index(['target', 'block']), b.set_index(['target', 'block'])
for k in a.index: add('RF balanced accuracy', f'{k[0]} {k[1]}', a.at[k, 'balanced_acc'], b.at[k, 'balanced_acc'])
# --- RF class probabilities for A. obscura (mean of 4 HQ MAGs)
a, b = both('rf_predictions.csv')
for blk in ['all', 'all_nonpfam', 'pfam_ab', 'codon', 'struct']:
    pa = a[(a.target == 'class') & (a.block == blk) & a.LOCUSTAG.isin(HQ)].drop(columns=['LOCUSTAG', 'target', 'block']).mean()
    pb = b[(b.target == 'class') & (b.block == blk) & b.LOCUSTAG.isin(HQ)].drop(columns=['LOCUSTAG', 'target', 'block']).mean()
    add('RF A. obscura', f'{blk}: top class', pa.idxmax(), pb.idxmax()); add('RF A. obscura', f'{blk}: max class probability', pa.max(), pb.max())
    for c in pa.sort_values(ascending=False).index[:3]: add('RF A. obscura', f'{blk}: P({c})', pa[c], pb.get(c, np.nan))
    sa = a[(a.target == 'subphylum') & (a.block == blk) & a.LOCUSTAG.isin(HQ)][['Pezizomycotina', 'Saccharomycotina', 'Taphrinomycotina']].mean()
    sb = b[(b.target == 'subphylum') & (b.block == blk) & b.LOCUSTAG.isin(HQ)][['Pezizomycotina', 'Saccharomycotina', 'Taphrinomycotina']].mean()
    for c in sa.index: add('RF A. obscura', f'{blk}: P({c})', sa[c], sb[c])
a, b = both('rf_loco.csv')
for blk in ['all', 'pfam_ab']:
    add('RF LOCO', f'{blk}: median of class median maxprob', a[a.block == blk].median_maxprob.median(), b[b.block == blk].median_maxprob.median())
    add('RF LOCO', f'{blk}: AO maxprob (median over folds)', a[a.block == blk].AO_maxprob.median(), b[b.block == blk].AO_maxprob.median())
# --- univariate percentiles
a, b = both('univariate.csv'); a, b = a.set_index('feature'), b.set_index('feature')
for f in ['log10_genome_mb', 'genes_per_mb', 'coding_fraction', 'median_intergenic', 'median_intron_len', 'introns_per_transcript', 'gc_pct', 'n_trna',
          'frac_signalp', 'frac_secreted_strict', 'median_idp_fraction', 'frac_zn2cys6', 'telomere_tracts', 'paralog_families_per1k']:
    if f in a.index and f in b.index:
        add('univariate', f'{f}: AO value (HQ median)', a.at[f, 'AO_HQ_median'], b.at[f, 'AO_HQ_median'])
        add('univariate', f'{f}: pct in Ascomycota', a.at[f, 'pct_Ascomycota'], b.at[f, 'pct_Ascomycota'])
# --- PCoA: sign-aligned axes, correlation over common units, A. obscura position
for blk in ['struct', 'codon', 'aa', 'trna', 'pfam_pa', 'pfam_ab', 'cazy', 'tf', 'cell']:
    a, b = both(f'pcoa_{blk}.csv', index_col=0); ix = a.index.intersection(b.index)
    for ax in ['PCo1', 'PCo2']:
        r = np.corrcoef(a.loc[ix, ax], b.loc[ix, ax])[0, 1]; s = np.sign(r) or 1
        add('PCoA', f'{blk} {ax}: |r| reference vs rerun over common units', 1.0, abs(r))
        add('PCoA', f'{blk} {ax}: A. obscura (HQ mean, sign-aligned)', a.loc[HQ, ax].mean(), s * b.loc[HQ, ax].mean())
        add('PCoA', f'{blk} {ax}: A. obscura percentile', 100 * (a[ax] < a.loc[HQ, ax].mean()).mean(), 100 * ((s * b[ax]) < s * b.loc[HQ, ax].mean()).mean())
# --- nearest classes (full-space distances)
a, b = both('ao_class_position.csv')
for blk in ['struct', 'codon', 'pfam_ab', 'pfam_pa', 'trna']:
    na = a[(a.block == blk) & (a.n_genera >= 3)].sort_values('dist_to_centroid').CLASS.head(3).tolist()
    nb = b[(b.block == blk) & (b.n_genera >= 3)].sort_values('dist_to_centroid').CLASS.head(3).tolist()
    add('position', f'{blk}: 3 nearest class centroids', ', '.join(na), ', '.join(nb))
# --- phylogeny
a, b = both('phylo_robustness.csv')
for i in range(len(a)):
    add('phylogeny', f'{a.tree[i]}: placement', a.placement[i], b.placement[i] if i < len(b) else None)
    add('phylogeny', f'{a.tree[i]}: support AO+sister', a.support_of_AO_plus_sister[i], b.support_of_AO_plus_sister[i] if i < len(b) else np.nan)
    add('phylogeny', f'{a.tree[i]}: AO monophyletic', a.AO_monophyletic[i], b.AO_monophyletic[i] if i < len(b) else None)
a, b = both('phylo_patristic.csv', index_col=0)
add('phylogeny', 'nearest 5 tips (patristic)', '; '.join(a.SPECIES.head(5)), '; '.join(b.SPECIES.head(5)))
add('phylogeny', 'nearest tip distance', a.patristic_from_AO.iloc[0], b.patristic_from_AO.iloc[0])
# --- MAG QC
a, b = both('mag_contam_summary.csv'); a, b = a.set_index('strain'), b.set_index('strain')
for s in a.index: add('MAG QC', f'{s}: flagged % bp', a.at[s, 'flagged_pct_bp'], b.at[s, 'flagged_pct_bp'] if s in b.index else np.nan)
if (NEW / 'mag_aai.csv').exists():
    a, b = both('mag_aai.csv'); a, b = a.set_index('mag'), b.set_index('mag')
    for s in a.index:
        add('MAG QC', f'{s}: median AAI vs Mars-82.9', a.at[s, 'median_AAI'], b.at[s, 'median_AAI'])
        add('MAG QC', f'{s}: frac identical vs Mars-82.9', a.at[s, 'frac_identical'], b.at[s, 'frac_identical'])
a, b = both('sprot_class_table.csv')
for c in ['frac_SCHPO', 'frac_YEAST', 'frac_NEUCR', 'frac_Fungi', 'frac_prok_like', 'pid_SCHPO', 'pidcov_SCHPO', 'pid_NEUCR']:
    add('SwissProt', f'A. obscura {c}', a[c].iloc[0], b[c].iloc[0])
R = pd.DataFrame(rows); R.to_csv(NEW / 'rerun_vs_reference.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 1000); pd.set_option('display.max_colwidth', 90)
print('reference:', REF); print(R.to_string(index=False))
