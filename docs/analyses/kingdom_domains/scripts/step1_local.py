import duckdb, numpy as np, pandas as pd, sys, os
A=os.environ.get('KW_BASE', '.') + '/kw/data'
W=os.environ.get('KW_BASE', '.') + '/kw/work'
c=duckdb.connect(); c.sql("set threads=8")
L=c.sql(f"""select l.*, e.obligate_mode, e.obligate_lineage, e.hs_group, e.septation, p.n_proteins
 from '{A}/genome_labels.parquet' l left join '{A}/genome_labels_ext.parquet' e using(LOCUSTAG)
 left join '{A}/genome_protein_summary.parquet' p on p.genome=l.LOCUSTAG
 where is_representative and q_ann""").df()
L['q']=L.q_busco_rel & L.q_dup & L.q_n50
L=L[L.n_proteins.notna()].reset_index(drop=True)
print('reps',len(L),'q',L.q.sum())
c.register('L',L[['LOCUSTAG']])
feats=c.sql(f"""
 select 'pfam' ftype, pfam_acc feat, genome, n_proteins n from '{A}/pfam_genome_matrix.parquet' where genome in (select LOCUSTAG from L)
 union all select 'cazy_fam', cazy_family, genome, n_proteins from '{A}/cazy_class_counts.parquet' where genome in (select LOCUSTAG from L)
 union all select 'cazy_class', cazy_class, genome, sum(n_proteins) from '{A}/cazy_class_counts.parquet' where genome in (select LOCUSTAG from L) group by all
 union all select 'cazy_sub', substrate, genome, n_proteins from '{A}/cazy_substrate_counts.parquet' where genome in (select LOCUSTAG from L)
 union all select 'merops', merops_family, genome, n_proteins from '{A}/merops_family_counts.parquet' where genome in (select LOCUSTAG from L)
 union all select 'tf', tf_family, genome, n_proteins from '{A}/tf_counts.parquet' where genome in (select LOCUSTAG from L)
 union all select 'trait', set_name, genome, n_proteins from '{A}/pfam_trait_counts.parquet' where genome in (select LOCUSTAG from L)
""").df()
feats=feats[feats.feat.notna()]
print(feats.groupby('ftype').feat.nunique())
feats['fid']=feats.ftype+'|'+feats.feat
fcat=pd.Categorical(feats.fid); gcat=pd.Categorical(feats.genome,categories=L.LOCUSTAG)
M=np.zeros((len(L),len(fcat.categories)),dtype=np.int32)
M[gcat.codes,fcat.codes]=feats.n.values.astype(np.int32)
F=pd.DataFrame({'fid':fcat.categories}); F['ftype']=F.fid.str.split('|').str[0]; F['feat']=F.fid.str.split('|',n=1).str[1]
cl=pd.read_csv(os.environ.get('KW_BASE', '.') + '/kw/data/Pfam-A.clans.tsv',sep='\t',header=None,names=['feat','clan','clan_id','pfam_id','desc'])
F=F.merge(cl,how='left',on='feat')
pid=c.sql(f"select pfam_acc feat, any_value(pfam_id) pid2 from '{A}/pfam_genome_matrix.parquet' group by 1").df()
F=F.merge(pid,how='left',on='feat'); F['name']=F.pid2.fillna(F.pfam_id).fillna(F.feat)
print('pfam w/o desc',(F.ftype.eq('pfam')&F.desc.isna()).sum())
# scalar per-genome
S=c.sql(f"""select s.genome, n_secreted_strict, n_signalp, n_gpi, n_tm2, i.median_idp_fraction, i.frac_idp_ge50, p.median_len_model, p.frac_with_pfam, p.frac_multidomain
 from '{A}/secretome.parquet' s left join '{A}/idp_genome.parquet' i using(genome) left join '{A}/genome_protein_summary.parquet' p using(genome)""").df()
L=L.merge(S,how='left',left_on='LOCUSTAG',right_on='genome')
np.save(f'{W}/M_genome.npy',M); F.to_parquet(f'{W}/features.parquet'); L.to_parquet(f'{W}/genomes.parquet')
print(M.shape, M.sum())
# ---- genus level
L['ORDER_']=L.ORDER_.fillna('NA'); 
def genus_agg(mask,tag):
    Ls=L[mask]; idx=Ls.index.values
    g=pd.Categorical(Ls.GENUS); ng=len(g.categories)
    P=np.zeros((ng,M.shape[1]),dtype=np.float32); Md=np.zeros((ng,M.shape[1]),dtype=np.float32); n=np.zeros(ng,int)
    for k in range(ng):
        rows=idx[g.codes==k]; sub=M[rows]
        P[k]=(sub>0).mean(0); Md[k]=np.median(sub,axis=0); n[k]=len(rows)
    G=Ls.groupby('GENUS').agg(PHYLUM=('PHYLUM','first'),CLASS=('CLASS','first'),ORDER_=('ORDER_','first'),FAMILY=('FAMILY','first'),
        growth_form=('growth_form','first'),primary_lifestyle=('primary_lifestyle','first'),yeast_origin=('yeast_origin','first'),
        obligate_mode=('obligate_mode','first'),obligate_lineage=('obligate_lineage','first'),hs_group=('hs_group','first'),
        n_genomes=('LOCUSTAG','size'),n_proteins=('n_proteins','median'),busco=('busco_complete','median'),genome_mb=('total_length_bp',lambda x: x.median()/1e6),
        sec_frac=('n_secreted_strict',lambda x:(x/Ls.loc[x.index,'n_proteins']).median()),gpi_frac=('n_gpi',lambda x:(x/Ls.loc[x.index,'n_proteins']).median()),
        tm_frac=('n_tm2',lambda x:(x/Ls.loc[x.index,'n_proteins']).median()),idp=('median_idp_fraction','median'),frac_pfam=('frac_with_pfam','median'),
        med_len=('median_len_model','median')).reindex(g.categories).reset_index()
    np.save(f'{W}/P_{tag}.npy',P); np.save(f'{W}/Md_{tag}.npy',Md); G.to_parquet(f'{W}/genera_{tag}.parquet'); print(tag,P.shape)
genus_agg(L.q.values,'q')
genus_agg(np.ones(len(L),bool),'all')
