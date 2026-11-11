"""Build per-genome feature blocks for the comparison set from atlas intermediates."""
import duckdb, pandas as pd, numpy as np, os
os.makedirs('data/feat', exist_ok=True)
d = duckdb.connect()
U = pd.read_parquet('data/ao_set.parquet').set_index('LOCUSTAG')
d.register('U', U.reset_index()[['LOCUSTAG']])
def q(sql): return d.sql(sql).df()
def pivot(df, idx, col, val):
    return df.pivot_table(index=idx, columns=col, values=val, aggfunc='sum', fill_value=0).reindex(U.index).fillna(0)

# ---- (a) structure scalars
g = q("select g.* from 'data/A2/genome_gene_stats.parquet' g join U on g.genome=U.LOCUSTAG").set_index('genome').reindex(U.index)
s = pd.DataFrame(index=U.index)
mb = U.total_length_bp / 1e6
s['log10_genome_mb'] = np.log10(mb); s['gc_pct'] = U.gc_pct; s['masked_pct'] = U.masked_pct
s['log10_n_genes'] = np.log10(g.n_protein_coding)
s['genes_per_mb'] = g.n_protein_coding / mb
s['coding_fraction'] = g.mean_protein_len * 3 * g.n_proteins / U.total_length_bp
for c in ['median_gene_span','median_cds_len','median_protein_len','mean_protein_len','q10_protein_len','q90_protein_len','frac_short_lt100','frac_long_ge1000',
          'introns_per_transcript','frac_intronless','median_intron_len','q90_intron_len','mean_intron_gc','frac_gc_ag',
          'median_intergenic','q10_intergenic','q90_intergenic','mean_intergenic','partial_fraction']:
    s[c] = g[c]
s['intron_gc_minus_genome_gc'] = g.mean_intron_gc - U.gc_pct
s['n_trna'] = g.n_trna_genes; s['trna_per_mb'] = g.n_trna_genes / mb
s.to_parquet('data/feat/struct.parquet')
nprot = g.n_proteins

# ---- (b) codon, aa
cod = pivot(q("select species_prefix g, codon, frequency f from 'data/raw/codon_freq.parquet' join U on species_prefix=LOCUSTAG"), 'g', 'codon', 'f')
cod = cod.drop(columns=[c for c in cod.columns if c in ('TAA','TAG','TGA') or len(c)!=3 or set(c)-set('ACGT')])
cod.to_parquet('data/feat/codon.parquet')
aa = pivot(q("select species_prefix g, amino_acid a, frequency f from 'data/raw/aa_freq.parquet' join U on species_prefix=LOCUSTAG"), 'g', 'a', 'f')
aa = aa[[c for c in aa.columns if c in list('ACDEFGHIKLMNPQRSTVWY')]]
aa.to_parquet('data/feat/aa.parquet')
# GC3
gc3 = cod[[c for c in cod.columns if c[2] in 'GC']].sum(1) / cod.sum(1)
# ---- (c) tRNA
tr = q("select genome g, amino_acid||'_'||codon k, n from 'data/A5/trna_codon_counts.parquet' join U on genome=LOCUSTAG where amino_acid not in ('Undet','Sup','Pseudo')")
trm = pivot(tr, 'g', 'k', 'n'); trm = trm.loc[:, (trm > 0).mean() >= 0.02]
trm.to_parquet('data/feat/trna.parquet')
# ---- (d) pfam
pf = q("select genome g, pfam_acc||'|'||pfam_id k, n_proteins n from 'data/A7/pfam_genome_matrix.parquet' join U on genome=LOCUSTAG")
pfm = pivot(pf, 'g', 'k', 'n').astype('float32')
pfm.to_parquet('data/feat/pfam.parquet'); print('pfam', pfm.shape)
# ---- (e) cazy, merops, tf
cz = pivot(q("select genome g, cazy_family k, n_proteins n from 'data/A5/cazy_class_counts.parquet' join U on genome=LOCUSTAG"), 'g', 'k', 'n'); cz.to_parquet('data/feat/cazy.parquet')
czs = pivot(q("select genome g, substrate k, n_proteins n from 'data/A5/cazy_substrate_counts.parquet' join U on genome=LOCUSTAG where substrate is not null and substrate not in ('-','')"), 'g', 'k', 'n'); czs.to_parquet('data/feat/cazy_substrate.parquet')
czc = pivot(q("select genome g, cazy_class k, n_proteins n from 'data/A5/cazy_class_counts.parquet' join U on genome=LOCUSTAG"), 'g', 'k', 'n')
me = pivot(q("select genome g, merops_family k, n_proteins n from 'data/A5/merops_family_counts.parquet' join U on genome=LOCUSTAG where merops_family is not null"), 'g', 'k', 'n'); me.to_parquet('data/feat/merops.parquet')
tf = pivot(q("select genome g, tf_family k, n_proteins n from 'data/A5/tf_counts.parquet' join U on genome=LOCUSTAG"), 'g', 'k', 'n'); tf.to_parquet('data/feat/tf.parquet')
# ---- (f) cell-biology scalars
sec = q("select s.* from 'data/A5/secretome.parquet' s join U on genome=LOCUSTAG").set_index('genome').reindex(U.index)
wp = pivot(q("select genome g, localization k, n from 'data/A5/wolfpsort_counts.parquet' join U on genome=LOCUSTAG"), 'g', 'k', 'n')
wp.to_parquet('data/feat/wolfpsort.parquet')
idp = q("select i.* from 'data/A5/idp_genome.parquet' i join U on genome=LOCUSTAG").set_index('genome').reindex(U.index)
p4 = q("select i.* from 'data/A4/genome_protein_summary.parquet' i join U on genome=LOCUSTAG").set_index('genome').reindex(U.index)
tfg = tf.sum(1)
f = pd.DataFrame(index=U.index)
f['frac_signalp'] = sec.n_signalp / nprot; f['frac_secreted_strict'] = sec.n_secreted_strict / nprot
f['frac_gpi'] = sec.n_gpi / nprot; f['frac_tm1plus'] = sec.n_tm1 / nprot; f['frac_tm2plus'] = sec.n_tm2 / nprot
f['median_idp_fraction'] = idp.median_idp_fraction; f['frac_idp_ge50'] = idp.frac_idp_ge50
for c in ['nucl','mito','cyto','extr','plas','pero','E.R.','cysk','golg']:
    if c in wp: f['wolf_'+c] = wp[c] / wp.sum(1)
f['frac_cazy'] = cz.sum(1) / nprot
for c in czc.columns: f['frac_cazy_'+c] = czc[c] / nprot
f['frac_merops'] = me.sum(1) / nprot; f['frac_tf'] = tfg / nprot
f['frac_zn2cys6'] = tf.get('Zn(II)2Cys6 zinc cluster', 0) / nprot
f['frac_with_pfam'] = p4.frac_with_pfam; f['frac_multidomain'] = p4.frac_multidomain; f['frac_with_swissprot'] = p4.frac_with_swissprot
f['frac_best_hit_yeast'] = p4.frac_best_hit_yeast; f['frac_best_hit_neucr_emeni'] = p4.frac_best_hit_neucr_emeni
f['gc3'] = gc3
# signaling/NLR trait sets: recomputed from the Pfam matrix for every genome (atlas A11 table covers species representatives only);
# value = summed proteins-with-domain over the set's Pfam models (double counts multi-domain proteins identically for all genomes)
ts = pd.read_csv('data/labels/pfam_trait_sets.tsv', sep='\t')
acc = pd.Series([c.split('|')[0] for c in pfm.columns], index=pfm.columns)
for sname, grp in ts.groupby('set'):
    cols = acc.index[acc.isin(grp.pfam_acc)]
    f['sig_'+sname+'_per1k'] = 1000 * pfm[cols].sum(1) / nprot
f.to_parquet('data/feat/cell.parquet')
# ---- (g) wgd / telomere
w = q("select w.* from 'data/A5/wgd_genome.parquet' w join U on genome=LOCUSTAG").set_index('genome').reindex(U.index)
ws = q("select w.* from 'data/raw/wgd.ks.summary.parquet' w join U on species_prefix=LOCUSTAG").set_index('species_prefix').reindex(U.index)
t = q("select t.* from 'data/raw/telomere_summary.parquet' t join U on t.LOCUSTAG=U.LOCUSTAG").set_index('LOCUSTAG').reindex(U.index)
gg = pd.DataFrame(index=U.index)
gg['paralog_nodes_per1k'] = 1000 * w.n_tree_nodes / g.n_protein_coding; gg['young_nodes_per1k'] = 1000 * w.n_young_nodes / g.n_protein_coding
gg['nearzero_nodes_per1k'] = 1000 * w.n_nearzero_nodes / g.n_protein_coding; gg['paralog_families_per1k'] = 1000 * ws.n_families / g.n_protein_coding
gg['identical_groups'] = ws.n_identical_groups
gg['telomere_tracts'] = t.telomere_tracts; gg['telomere_top_monomer'] = t.telomere_top_monomer; gg['telomere_scaffolds'] = t.telomere_scaffolds
gg.to_parquet('data/feat/wgd_telo.parquet')
for n in ['struct','cell','codon','aa','trna','cazy','cazy_substrate','merops','tf','wolfpsort','wgd_telo']:
    x = pd.read_parquet(f'data/feat/{n}.parquet'); print(n, x.shape, int(x.isna().any(axis=1).sum()))
