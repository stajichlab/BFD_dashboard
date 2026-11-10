import duckdb, os, pandas as pd, numpy as np
pd.set_option('display.width',400); pd.set_option('display.max_rows',500); pd.set_option('display.max_columns',60)
W=os.path.expanduser('~/bigdata/tmp_kx_scratch/')
con=duckdb.connect()
d=con.sql(f"select * from '{W}g.parquet' where grp is not null and is_representative and q_ann").df()
d['mb']=d.total_length_bp/1e6
for c,src in [('sec_k','n_secreted_strict'),('sp_k','n_signalp'),('cazy_k','n_cazy'),('mer_k','n_merops'),('tf_k','n_tf_proteins'),('tm_k','n_tm2')]:
    d[c]=1000*d[src]/d.n_proteins
d['gene_density']=d.n_protein_coding/d.mb
mets=['mb','n_protein_coding','gene_density','busco_complete','N50_bp','partial_fraction','introns_per_transcript','frac_intronless','median_intron_len','median_intergenic','median_gene_span','median_protein_len','median_cds_len_complete','masked_pct','gc_pct','n_secreted_strict','sec_k','n_cazy','cazy_k','n_merops','mer_k','n_tf_proteins','tf_k','tm_k','median_idp_fraction']
def summ(df,label):
    gen=df.groupby(['grp','GENUS'])[mets].median().reset_index()
    out=gen.groupby('grp')[mets].median()
    out.insert(0,'n_gen',gen.groupby('grp').size()); out.insert(1,'n_sp',df.groupby('grp').size())
    print('=====',label); print(out.round(2).to_string())
    return out,gen
o_all,gen_all=summ(d,'ALL species reps, genus medians')
dq=d[d.q_busco_rel & d.q_dup & d.q_n50]
o_q,gen_q=summ(dq,'Q-filtered (busco_rel,dup,n50) genus medians')
o_all.to_csv(W+'grp_genus_medians_all.tsv',sep='\t'); o_q.to_csv(W+'grp_genus_medians_q.tsv',sep='\t')
# Coemansia species-level distribution: all vs high quality (N50>=50kb)
c=d[d.grp=='Kx:Coemansia']
print('Coemansia all n',len(c)); print(c[mets].quantile([.1,.25,.5,.75,.9]).round(2).T.to_string())
ch=c[c.N50_bp>=50000]; print('Coemansia N50>=50kb n',len(ch)); print(ch[mets].quantile([.1,.5,.9]).round(2).T.to_string())
# Kickxellomycota covariation, species reps (all) & hq
k=d[d.grp.str.startswith('Kx:')].copy()
cv=['mb','n_protein_coding','introns_per_transcript','median_intron_len','median_intergenic','median_gene_span','median_protein_len','masked_pct','gc_pct','busco_complete','N50_bp','sec_k','cazy_k']
print('Kx species reps n',len(k)); print(k[cv].corr(method='spearman').round(2).to_string())
kg=k.groupby('GENUS')[cv].median(); print('Kx genus medians n',len(kg)); print(kg.round(2).to_string()); print(kg.corr(method='spearman').round(2).to_string())
kh=k[k.N50_bp>=50000]; print('Kx HQ N50>=50kb n',len(kh), kh.grp.value_counts().to_dict()); print(kh[cv].corr(method='spearman').round(2).to_string())
k.to_csv(W+'kx_species_reps.tsv',sep='\t',index=False)
gen_all.to_csv(W+'genus_medians_all.tsv',sep='\t',index=False)
