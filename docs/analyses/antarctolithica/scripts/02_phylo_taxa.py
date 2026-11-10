"""Pick ~130 taxa for a quick phylogenomic placement of Antarctolithica."""
import pandas as pd, numpy as np
U = pd.read_parquet('data/ao_set.parquet'); L = pd.read_parquet('data/labels/genome_labels.parquet')
rng = np.random.default_rng(1)
want_genera = ['Friedmanniomyces','Cryomyces','Hortaea','Baudoinia','Coniosporium','Knufia','Exophiala','Cladophialophora','Rachicladosporium','Recurvomyces','Aureobasidium','Zymoseptoria','Neurospora','Aspergillus','Saccharomyces','Candida','Yarrowia','Lipomyces','Tuber','Botrytis','Xylona','Endocarpon','Umbilicaria','Cladonia','Arthonia','Lichina','Peltula','Orbilia','Pyronema','Trichoderma','Verrucaria','Capronia','Phaeomoniella','Symbiotaphrina','Trigonopsis','Geoglossum','Coccidioides','Penicillium','Cenococcum','Rhizodiscina','Lichenothelia','Saxomyces','Extremus','Neodevriesia','Acidomyces','Teratosphaeria','Piedraia','Myriangium','Elsinoe','Dothidea']
n_per = {'Dothideomycetes':8,'Eurotiomycetes':6,'Lecanoromycetes':6,'Leotiomycetes':4,'Sordariomycetes':4,'Pezizomycetes':4,'Orbiliomycetes':2,'Arthoniomycetes':4,'Lichinomycetes':4,'Geoglossomycetes':2,'Candelariomycetes':2,'Coniocybomycetes':2,'Saccharomycetes':2,'Pichiomycetes':2,'Dipodascomycetes':2,'Lipomycetes':2,'Trigonopsidomycetes':2}
g = U[U.is_genus_unit]
pick = set(U[U.role=='AO'].LOCUSTAG) | set(U[U.subphylum=='Taphrinomycotina'].LOCUSTAG) | set(U[U.role=='unplaced'].LOCUSTAG)
pick |= set(g[g.GENUS.isin(want_genera) & (g.role=='train')].LOCUSTAG)
for c, n in n_per.items():
    s = g[(g.CLASS==c)&(g.role=='train')&(g.busco_complete>=90 if c not in('Saccharomycetes','Pichiomycetes') else True)]
    pick |= set(s.sample(min(n,len(s)), random_state=1).LOCUSTAG)
pick |= set(g[(g.role=='train')&g.CLASS.isin(['Sareomycetes','Xylonomycetes','Xylobotryomycetes','Alloascoideomycetes','Sporopachydermiomycetes'])].LOCUSTAG)
out = g[g.role=='outgroup'].groupby('CLASS').head(1).head(7)
pick |= set(out.LOCUSTAG)
muc = L[(L.PHYLUM=='Mucoromycota')&L.is_representative&(L.busco_complete>90)&(L.busco_duplicated<10)].sort_values('N50_bp',ascending=False).head(2)
pick |= set(muc.LOCUSTAG)
T = L[L.LOCUSTAG.isin(pick)][['LOCUSTAG','PHYLUM','CLASS','ORDER_','GENUS','SPECIES','STRAIN','busco_complete']].copy()
pep = pd.read_csv('phylo/locustag_to_pep.tsv', sep='\t', header=None, names=['LOCUSTAG','name','pep'])
T = T.merge(pep, on='LOCUSTAG', how='left')
import re
T['label'] = [re.sub(r'[^A-Za-z0-9]+','_', f"{(c or 'unplaced')[:12]}|{s}|{l}".replace('|','__')) for c,s,l in zip(T.CLASS,T.SPECIES,T.LOCUSTAG)]
T.to_csv('phylo/taxa.tsv', sep='\t', index=False)
print(len(T), T.pep.isna().sum()); print(T.CLASS.fillna('NA').value_counts().to_string())
T.pep.dropna().to_csv('phylo/pep_paths.txt', index=False, header=False)
