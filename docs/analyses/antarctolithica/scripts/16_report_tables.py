"""Render result CSVs into markdown table fragments (report/tables.json) used by 17_assemble_report.py."""
import json, os, numpy as np, pandas as pd
from common import *
from mdtab import md
T = {}
L = pd.read_parquet('data/labels/genome_labels.parquet'); A8 = L[L.SPECIES.str.contains('Antarctolith', na=False)].sort_values('busco_complete', ascending=False)
g = pd.read_parquet('data/A2/genome_gene_stats.parquet').set_index('genome'); cs = pd.read_csv('results/mag_contam_summary.csv').set_index('strain'); aai = pd.read_csv('results/mag_aai.csv').set_index('mag')
q = pd.DataFrame({'MAG': A8.STRAIN.values, 'LOCUSTAG': A8.LOCUSTAG.values, 'size Mb': (A8.total_length_bp / 1e6).round(2).values, 'contigs': A8.contig_count.values, 'N50 kb': (A8.N50_bp / 1e3).round(1).values, 'GC %': A8.gc_pct.values,
                  'BUSCO C %': A8.busco_complete.values, 'dup %': A8.busco_duplicated.values, 'frag %': A8.busco_fragmented.values, 'missing %': A8.busco_missing.values})
q['genes'] = [f"{int(g.n_protein_coding[i]):,}" if i in g.index else 'not annotated' for i in A8.LOCUSTAG]
q['flagged bp %'] = [f"{cs.flagged_pct_bp[s]:.2f}" if s in cs.index else '' for s in A8.STRAIN]
q['median AAI vs 82.9'] = [f"{aai.median_AAI[s]:.1f}" if s in aai.index else ('ref' if s == 'Mars-82.9' else '') for s in A8.STRAIN]
q['identical proteins vs 82.9'] = [f"{aai.frac_identical[s]:.0%}" if s in aai.index else '' for s in A8.STRAIN]
T['qc'] = md(q)
R = pd.read_csv('results/univariate.csv').set_index('feature')
def uni(feats):
    d = R.loc[feats, ['AO_HQ_median', 'AO_min6', 'AO_max6', 'med_Pezizomycotina', 'med_Saccharomycotina', 'med_Taphrinomycotina', 'med_Taphrinomycetes', 'pct_Ascomycota', 'pct_Pezizomycotina', 'pct_Saccharomycotina', 'pct_Taphrinomycotina', 'z_Pezizomycotina', 'epsilon2']].copy()
    d.columns = ['A. obscura (HQ median)', 'min of 6', 'max of 6', 'Pezizo. median', 'Saccharo. median', 'Taphrino. median', 'Taphrinomycetes median', 'pct in Asco.', 'pct in Pezizo.', 'pct in Saccharo.', 'pct in Taphrino.', 'robust z vs Pezizo.', 'KW eps2']
    return md(d.reset_index(), 3)
T['uni_struct'] = uni(['log10_genome_mb', 'log10_n_genes', 'genes_per_mb', 'coding_fraction', 'median_intergenic', 'q10_intergenic', 'q90_intergenic', 'median_gene_span', 'median_cds_len', 'median_protein_len', 'frac_long_ge1000', 'introns_per_transcript', 'frac_intronless', 'median_intron_len', 'q90_intron_len', 'mean_intron_gc', 'intron_gc_minus_genome_gc', 'frac_gc_ag', 'gc_pct', 'gc3', 'masked_pct', 'n_trna', 'trna_per_mb'])
T['uni_cell'] = uni(['frac_signalp', 'frac_secreted_strict', 'frac_gpi', 'frac_tm1plus', 'median_idp_fraction', 'frac_idp_ge50', 'wolf_nucl', 'wolf_mito', 'wolf_cyto', 'wolf_extr', 'wolf_plas', 'wolf_cysk', 'frac_cazy', 'frac_cazy_GH', 'frac_cazy_GT', 'frac_cazy_AA', 'frac_merops', 'frac_tf', 'frac_zn2cys6', 'frac_with_pfam', 'frac_multidomain', 'sig_kinase_per1k', 'sig_ras_arf_per1k', 'sig_calcium_per1k', 'sig_wd40_per1k', 'sig_het_per1k', 'sig_nacht_per1k', 'sig_ank_per1k'])
T['uni_wgd'] = uni(['paralog_nodes_per1k', 'young_nodes_per1k', 'nearzero_nodes_per1k', 'paralog_families_per1k', 'identical_groups', 'telomere_tracts'])
P = pd.read_csv('results/permanova.csv'); T['permanova'] = md(P[['block', 'n_features', 'n_genera', 'n_classes', 'pseudoF_class', 'R2_class', 'p_class', 'pseudoF_subphylum', 'R2_subphylum', 'p_subphylum', 'dispersion_ratio', 'PCo1_pct', 'PCo2_pct']], 3)
pos = pd.read_csv('results/ao_class_position.csv'); rows = []
for b in P.block:
    d = pos[(pos.block == b) & (pos.n_genera >= 3)].sort_values('dist_to_centroid'); top = d.head(3)
    rows.append(dict(block=b, **{f'nearest {i+1}': f"{r.CLASS} (ratio {r.ratio:.2f}; {r.pct_members_farther:.0f}% of members farther)" for i, r in enumerate(top.itertuples())}))
T['position'] = md(pd.DataFrame(rows))
nn = pd.read_csv('results/ao_nearest_neighbours.csv'); T['nn'] = md(pd.DataFrame([dict(block=b, nearest_genus_units='; '.join(f"{r.SPECIES} [{str(r.CLASS)[:10]}]" for r in d.head(5).itertuples())) for b, d in nn.groupby('block', sort=False)]))
c = pd.read_csv('results/ao_consistency.csv'); T['consistency'] = md(c, 3)
z = pd.read_csv('results/composition_z.csv'); rows = []
for b in ['aa', 'codon', 'trna']:
    d = z[z.block == b]; d = d.reindex(d.z_Pezizomycotina.abs().sort_values(ascending=False).index).head(6)
    for r in d.itertuples(): rows.append(dict(block=b, feature=r.feature, AO=r.AO_raw, Pezizo_median=r.rawmed_Pezizomycotina, Saccharo_median=r.rawmed_Saccharomycotina, Taphrino_median=r.rawmed_Taphrinomycotina, z_Pezizo=r.z_Pezizomycotina, z_Saccharo=r.z_Saccharomycotina, z_Taphrino=r.z_Taphrinomycotina))
T['composition'] = md(pd.DataFrame(rows), 3)
s = pd.read_csv('results/sprot_class_table.csv'); s = s.rename(columns={s.columns[0]: 'group', s.columns[1]: 'class'})
s = s[(s.n_genera >= 3) | (s.group == 'Taphrinomycotina')][['group', 'class', 'n_genera', 'frac_SCHPO', 'frac_YEAST', 'frac_NEUCR', 'frac_EMENI_ASP', 'frac_Fungi', 'frac_prok_like', 'pid_SCHPO', 'pid_YEAST', 'pid_NEUCR', 'pid_EMENI', 'pid_NEUCR_minus_SCHPO', 'pid_SCHPO_minus_YEAST']]; T['sprot'] = md(s, 3)
cm = pd.read_csv('results/mag_contam_summary.csv'); T['contam'] = md(cm, 3)
ct = pd.read_csv('results/content_pfam.csv', index_col=0)
rare = ct[(ct.AO_HQ_present == 4) & (ct.prev_Asco < 0.05)].sort_values('prev_Asco')[['AO_all6_present', 'prev_Asco', 'prev_Pez', 'prev_Sac', 'prev_Taph']]; rare.index.name = 'Pfam'; T['rare'] = md(rare.reset_index(), 2)
ex = ct[(ct.AO_HQ_median_n >= 5) & (ct.pct_rate_in_Asco >= 99)].sort_values('log2_rate_vs_Asco_median', ascending=False).head(20)[['AO_HQ_median_n', 'med_Pez', 'med_Sac', 'med_Taph', 'log2_rate_vs_Asco_median', 'pct_rate_in_Asco']]; ex.index.name = 'Pfam'; T['expansions'] = md(ex.reset_index(), 3)
cand = pd.read_csv('results/content_candidates.csv'); cand = cand[~cand['set'].str.contains('Woronin')]
T['candidates'] = md(cand[['set', 'pfams', 'AO_HQ_median', 'AO_range6', 'med_Pez', 'med_Dothi', 'med_Eurot', 'med_Lecan', 'med_Sac', 'med_Taph', 'pct_rate_Asco']], 3)
lc = pd.read_csv('results/loss_check_blastx.csv'); lc = lc.rename(columns={lc.columns[0]: 'Pfam domain | Saitoella query'}); T['losscheck'] = md(lc)
miss = pd.read_csv('results/pez_universal_missing_per_genome.csv', index_col=0); tgm = miss[(miss.role == 'train')].groupby('CLASS').missing_pez_universal.median().sort_values()
T['missing_cal'] = md(pd.concat([tgm, pd.Series({'A. obscura (4 HQ MAGs, median)': miss.loc[HQ].missing_pez_universal.median(), 'A. obscura (absent from all 6)': 312})]).rename('Pezizomycotina-universal Pfams missing (of 3,745)').rename_axis('group').reset_index().sort_values('Pezizomycotina-universal Pfams missing (of 3,745)'), 4)
for n in ['cazy', 'merops']:
    t = pd.read_csv(f'results/content_{n}.csv', index_col=0); d = t[(t.AO_HQ_present >= 3)]; d = d[(d.pct_rate_in_Asco >= 97) | (d.pct_rate_in_Asco <= 10)].sort_values('pct_rate_in_Asco', ascending=False)[['AO_HQ_median_n', 'med_Pez', 'med_Sac', 'med_Taph', 'pct_rate_in_Asco']]; d.index.name = 'family'; T[n] = md(d.reset_index(), 3)
    a = t[(t.AO_all6_present == 0) & (t.prev_Asco >= 0.8)][['prev_Asco', 'prev_Pez', 'prev_Sac', 'prev_Taph', 'med_Pez']]; a.index.name = 'family'; T[n + '_absent'] = md(a.reset_index(), 2)
tf = pd.read_csv('results/content_tf.csv', index_col=0)[['AO_HQ_median_n', 'med_Pez', 'med_Sac', 'med_Taph', 'pct_rate_in_Asco']].sort_values('AO_HQ_median_n', ascending=False); tf.index.name = 'TF family'; T['tf'] = md(tf.reset_index(), 3)
if os.path.exists('results/rf_performance.csv') and os.path.exists('results/rf_loco.csv'):
    pf = pd.read_csv('results/rf_performance.csv'); T['rf_perf'] = md(pf[['target', 'block', 'n_features', 'balanced_acc', 'macro_F1', 'acc', 'balanced_acc_thin0.93', 'balanced_acc_thin0.75']], 3)
    Pp = pd.read_csv('results/rf_predictions.csv'); pc = Pp[(Pp.target == 'class') & Pp.LOCUSTAG.isin(HQ)]; cl = [c for c in pc.columns if c not in ('LOCUSTAG', 'target', 'block') and pc[c].notna().any() and c not in ('Pezizomycotina', 'Saccharomycotina', 'Taphrinomycotina')]
    M = pc.groupby('block', sort=False)[cl].mean(); rows = []
    for b, r in M.iterrows():
        t3 = r.sort_values(ascending=False).head(3); rows.append(dict(block=b, **{f'class {i+1}': f'{k} ({v:.2f})' for i, (k, v) in enumerate(t3.items())}))
    ps = Pp[(Pp.target == 'subphylum') & Pp.LOCUSTAG.isin(HQ)].groupby('block', sort=False)[['Pezizomycotina', 'Saccharomycotina', 'Taphrinomycotina']].mean().round(2)
    d = pd.DataFrame(rows).set_index('block').join(ps); T['rf_pred'] = md(d.reset_index(), 2)
    pm = Pp[(Pp.target == 'class') & (Pp.block == 'all') & Pp.LOCUSTAG.isin(AO)].set_index('LOCUSTAG').loc[AO, cl]; top = pm.mean().sort_values(ascending=False).index[:6]
    pm = pm[top].round(2); pm.insert(0, 'BUSCO C %', U.busco_complete[AO]); pm.insert(0, 'MAG', [MAGNAME[i] for i in AO]); T['rf_permag'] = md(pm, 2)
    un = Pp[(Pp.target == 'class') & (Pp.block == 'all') & ~Pp.LOCUSTAG.isin(AO)].set_index('LOCUSTAG'); rows = []
    for i, r in un[cl].iterrows(): rows.append(dict(species=U.SPECIES[i], strain=U.STRAIN[i], top_class=r.idxmax(), prob=r.max()))
    T['rf_unplaced'] = md(pd.DataFrame(rows).sort_values('prob', ascending=False), 2)
    I = pd.read_csv('results/rf_importance.csv'); rows = []
    for b in ['struct', 'cell', 'codon', 'aa', 'trna', 'tf', 'cazy', 'merops', 'pfam_ab', 'all']:
        d = I[(I.target == 'class') & (I.block == b)].sort_values('gini', ascending=False).head(6)
        rows.append(dict(block=b, top_features_gini='; '.join(f"{f.split('|')[-1]} ({g_:.3f})" for f, g_ in zip(d.feature, d.gini)), top_by_permutation='; '.join(f"{f.split('|')[-1]} ({p:.3f})" for f, p in d.sort_values('perm_bal_acc', ascending=False).head(3)[['feature', 'perm_bal_acc']].values)))
    T['rf_imp'] = md(pd.DataFrame(rows))
    lo = pd.read_csv('results/rf_loco.csv'); T['rf_loco'] = md(lo[lo.block.isin(['all', 'pfam_ab'])], 2)
    T['rf_loco_summary'] = md(lo.groupby('block').agg(median_of_class_median_maxprob=('median_maxprob', 'median'), max_class_median_maxprob=('median_maxprob', 'max'), AO_maxprob=('AO_maxprob', 'median'), AO_top_mode=('AO_top', lambda x: x.mode()[0])).reset_index(), 2)
os.makedirs('report', exist_ok=True); json.dump(T, open('report/tables.json', 'w'), indent=1); print(sorted(T))
