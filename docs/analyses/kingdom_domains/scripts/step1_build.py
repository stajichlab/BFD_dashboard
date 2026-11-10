import duckdb, numpy as np, pandas as pd, sys, os
A='/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data'
W='/bigdata/stajichlab/jstajich/kw_fableB'
c=duckdb.connect(); c.sql("set threads=8")
L=c.sql(f"""select l.*, e.obligate_mode, e.obligate_lineage, e.hs_group, e.septation, p.n_proteins
 from '{A}/labels/genome_labels.parquet' l left join '{A}/labels/genome_labels_ext.parquet' e using(LOCUSTAG)
 left join '{A}/A4/genome_protein_summary.parquet' p on p.genome=l.LOCUSTAG
 where is_representative and q_ann""").df()
L['q']=L.q_busco_rel & L.q_dup & L.q_n50
L=L[L.n_proteins.notna()].reset_index(drop=True)
print('reps',len(L),'q',L.q.sum())
c.register('L',L[['LOCUSTAG']])
feats=c.sql(f"""
 select 'pfam' ftype, pfam_acc feat, genome, n_proteins n from '{A}/A7/pfam_genome_matrix.parquet' where genome in (select LOCUSTAG from L)
 union all select 'cazy_fam', cazy_family, genome, n_proteins from '{A}/A5/cazy_class_counts.parquet' where genome in (select LOCUSTAG from L)
 union all select 'cazy_class', cazy_class, genome, sum(n_proteins) from '{A}/A5/cazy_class_counts.parquet' where genome in (select LOCUSTAG from L) group by all
 union all select 'cazy_sub', substrate, genome, n_proteins from '{A}/A5/cazy_substrate_counts.parquet' where genome in (select LOCUSTAG from L)
 union all select 'merops', merops_family, genome, n_proteins from '{A}/A5/merops_family_counts.parquet' where genome in (select LOCUSTAG from L)
 union all select 'tf', tf_family, genome, n_proteins from '{A}/A5/tf_counts.parquet' where genome in (select LOCUSTAG from L)
 union all select 'trait', set_name, genome, n_proteins from '{A}/A11/pfam_trait_counts.parquet' where genome in (select LOCUSTAG from L)
""").df()
feats=feats[feats.feat.notna()]
print(feats.groupby('ftype').feat.nunique())
feats['fid']=feats.ftype+'|'+feats.feat
fcat=pd.Categorical(feats.fid); gcat=pd.Categorical(feats.genome,categories=L.LOCUSTAG)
M=np.zeros((len(L),len(fcat.categories)),dtype=np.int32)
M[gcat.codes,fcat.codes]=feats.n.values.astype(np.int32)
F=pd.DataFrame({'fid':fcat.categories}); F['ftype']=F.fid.str.split('|').str[0]; F['feat']=F.fid.str.split('|',n=1).str[1]
cl=pd.read_csv('/bigdata/stajichlab/shared/lib/funannotate_db/Pfam-A.clans.tsv',sep='\t',header=None,names=['feat','clan','clan_id','pfam_id','desc'])
F=F.merge(cl,how='left',on='feat')
pid=c.sql(f"select pfam_acc feat, any_value(pfam_id) pid2 from '{A}/A7/pfam_genome_matrix.parquet' group by 1").df()
F=F.merge(pid,how='left',on='feat'); F['name']=F.pid2.fillna(F.pfam_id).fillna(F.feat)
print('pfam w/o desc',(F.ftype.eq('pfam')&F.desc.isna()).sum())
# scalar per-genome
S=c.sql(f"""select s.genome, n_secreted_strict, n_signalp, n_gpi, n_tm2, i.median_idp_fraction, i.frac_idp_ge50, p.median_len_model, p.frac_with_pfam, p.frac_multidomain
 from '{A}/A5/secretome.parquet' s left join '{A}/A5/idp_genome.parquet' i using(genome) left join '{A}/A4/genome_protein_summary.parquet' p using(genome)""").df()
L=L.merge(S,how='left',left_on='LOCUSTAG',right_on='genome')
np.save(f'{W}/M_genome.npy',M); F.to_parquet(f'{W}/features.parquet'); L.to_parquet(f'{W}/genomes.parquet')
print(M.shape, M.sum())
