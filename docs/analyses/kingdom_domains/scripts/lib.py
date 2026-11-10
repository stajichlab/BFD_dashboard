import os  # KW_BASE: directory holding kw/data, kw/work and fable_B_figs
import numpy as np, pandas as pd
W=os.environ.get('KW_BASE', '.') + '/kw/work'
pd.set_option('display.width',250); pd.set_option('display.max_columns',40); pd.set_option('display.max_rows',500); pd.set_option('display.max_colwidth',60)
F=pd.read_parquet(f'{W}/features.parquet')
SACC={'Pichiomycetes','Saccharomycetes','Dipodascomycetes','Lipomycetes','Trigonopsidomycetes','Sporopachydermiomycetes','Alloascoideomycetes'}
TAPH={'Taphrinomycetes','Pneumocystomycetes','Schizosaccharomycetes','Saitoellomycetes','Neolectomycetes','Novakomycetes'}
USTI={'Ustilaginomycetes','Exobasidiomycetes','Malasseziomycetes','Moniliellomycetes'}
PUCC={'Pucciniomycetes','Microbotryomycetes','Cystobasidiomycetes','Agaricostilbomycetes','Atractiellomycetes','Mixiomycetes','Classiculomycetes','Spiculogloeomycetes','Tritirachiomycetes'}
AGAR={'Agaricomycetes','Tremellomycetes','Dacrymycetes','Wallemiomycetes','Geminibasidiomycetes','Peribolosporomycetes'}
MUCO={'Mucoromycota','Mortierellomycota','Glomeromycota','Calcarisporiellomycota'}
ZOOP={'Zoopagomycota','Kickxellomycota','Entomophthoromycota','Basidiobolomycota'}
ZOOS={'Chytridiomycota','Blastocladiomycota','Neocallimastigomycota','Monoblepharomycota','Olpidiomycota','Sanchytriomycota'}
def load(tag):
    G=pd.read_parquet(f'{W}/genera_{tag}.parquet').rename(columns={'index':'GENUS'}); P=np.load(f'{W}/P_{tag}.npy'); Md=np.load(f'{W}/Md_{tag}.npy')
    G['CLASS']=G.CLASS.fillna(G.PHYLUM.fillna('NA')+'_incertae'); G['PHYLUM']=G.PHYLUM.fillna('NA')
    def sub(r):
        if r.PHYLUM=='Ascomycota': return 'Saccharomycotina' if r.CLASS in SACC else ('Taphrinomycotina' if r.CLASS in TAPH else 'Pezizomycotina')
        if r.PHYLUM=='Basidiomycota': return 'Ustilaginomycotina' if r.CLASS in USTI else ('Pucciniomycotina' if r.CLASS in PUCC else ('Agaricomycotina' if r.CLASS in AGAR else 'Basidio_other'))
        if r.PHYLUM in MUCO: return 'Mucoromycota_sl'
        if r.PHYLUM in ZOOP: return 'Zoopagomycota_sl'
        if r.PHYLUM in ZOOS: return 'zoosporic'
        return r.PHYLUM
    G['SUBGROUP']=G.apply(sub,axis=1)
    G['DIKARYA']=np.where(G.PHYLUM.isin(['Ascomycota','Basidiomycota']),'Dikarya','non-Dikarya')
    return G,P,Md
