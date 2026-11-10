"""Define the comparison set (universe) for the Antarctolithica obscura analysis.
Units: species representatives (atlas is_representative) passing atlas quality filters,
capped at 5 species/genus; genus representative = best BUSCO then N50 within genus."""
import duckdb, pandas as pd, numpy as np
d = duckdb.connect()
L = d.sql("select * from 'data/labels/genome_labels.parquet'").df()
have = set(d.sql("select genome from 'data/A2/genome_gene_stats.parquet'").df().genome)
L = L[L.LOCUSTAG.isin(have)].copy()
ao = L[L.SPECIES.str.contains('Antarctolith', na=False)].copy()
A = L[(L.PHYLUM == 'Ascomycota') & L.is_representative & ~L.LOCUSTAG.isin(ao.LOCUSTAG)].copy()
strict = A.q_ann & A.q_busco_rel & A.q_dup & A.q_n50
nstrict = A[strict].groupby('CLASS').size()
rare = set(A.CLASS.dropna().unique()) - set(nstrict[nstrict >= 15].index)
relaxed = A.CLASS.isin(rare) & A.q_ann & A.q_dup & (A.busco_complete >= 80)
unpl = A.CLASS.isna() & A.q_ann & A.q_dup & (A.busco_complete >= 80)
A = A[(strict & A.CLASS.notna()) | relaxed | unpl].copy()
A = A[~A.SPECIES.str.contains('fusant', na=False)]
A = A.sort_values(['busco_complete', 'N50_bp'], ascending=False)
A['genus_key'] = A.GENUS.fillna(A.SPECIES)
A['rank_in_genus'] = A.groupby('genus_key').cumcount()
A = A[A.rank_in_genus < 5].copy()
A['role'] = np.where(A.CLASS.isna(), 'unplaced', 'train')
B = L[(L.PHYLUM == 'Basidiomycota') & L.is_representative & L.q_ann & L.q_busco_rel & L.q_dup & L.q_n50 & L.CLASS.notna()].copy()
B = B.sort_values(['busco_complete', 'N50_bp'], ascending=False)
B['genus_key'] = B.GENUS.fillna(B.SPECIES); B['rank_in_genus'] = B.groupby('genus_key').cumcount()
B = B[B.rank_in_genus == 0].groupby('CLASS', group_keys=False).head(15).copy(); B['role'] = 'outgroup'
ao['genus_key'] = 'Antarctolithica'; ao['rank_in_genus'] = np.arange(len(ao)); ao['role'] = 'AO'
ao = ao.sort_values('busco_complete', ascending=False)
U = pd.concat([A, B, ao])
U['is_genus_unit'] = U.rank_in_genus == 0
taph = {'Taphrinomycetes','Schizosaccharomycetes','Pneumocystomycetes','Saitoellomycetes','Neolectomycetes','Novakomycetes','Archaeorhizomycetes'}
sacc = {'Saccharomycetes','Pichiomycetes','Dipodascomycetes','Lipomycetes','Trigonopsidomycetes','Alloascoideomycetes','Sporopachydermiomycetes'}
def sub(r):
    if r.PHYLUM != 'Ascomycota': return 'Basidiomycota'
    if r.CLASS is None or pd.isna(r.CLASS): return 'unplaced'
    if r.CLASS in taph: return 'Taphrinomycotina'
    if r.CLASS in sacc: return 'Saccharomycotina'
    return 'Pezizomycotina'
U['subphylum'] = U.apply(sub, axis=1)
U['hq_mag'] = U.LOCUSTAG.isin(['FF1FC3B3','F46C92EF','F7B4B958','F7E47E53'])
U.to_parquet('data/ao_set.parquet', index=False)
U[['LOCUSTAG','ASMID','PHYLUM','CLASS','ORDER_','GENUS','SPECIES','role','subphylum','is_genus_unit']].to_csv('data/ao_set.csv', index=False)
print(U.role.value_counts()); print(U[U.role=='train'].groupby(['subphylum','CLASS']).agg(n=('LOCUSTAG','size'), genera=('genus_key','nunique')).to_string())
print(len(U))
