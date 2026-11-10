from lib import *
L=pd.read_parquet(f'{W}/genomes.parquet'); M=np.load(f'{W}/M_genome.npy'); fidx={f:i for i,f in enumerate(F.fid)}
nm=dict(zip(F.feat,F.name))
doms={'FAS (Fas_alpha_ACP)':'PF18325','FAS_I_H':'PF18314','Piwi':'PF02171','Dicer_dimer':'PF03368','Complex1_49kDa':'PF00346','MoaC':'PF01967','Urease_alpha':'PF00449','Glycogen_syn':'PF05693','Glucan_synthase':'PF02364','ICL':'PF00463','Ku':'PF02735','HisG':'PF01634','Catalase':'PF00199','ATP-sulfurylase':'PF01747'}
doms={k:v for k,v in doms.items() if 'pfam|'+v in fidx}
print({k:nm[v] for k,v in doms.items()})
def lin(r):
    if r.GENUS in('Malassezia','Pneumocystis','Schizosaccharomyces','Hanseniaspora','Wallemia','Ustilago','Cryptococcus','Puccinia','Saccharomyces','Candida','Batrachochytrium'): return 'g:'+r.GENUS
    if r.PHYLUM in('Microsporidia','Glomeromycota','Neocallimastigomycota','Kickxellomycota','Zoopagomycota','Entomophthoromycota','Mortierellomycota','Mucoromycota','Chytridiomycota','Blastocladiomycota'): return r.PHYLUM
    if r.FAMILY in('Saccharomycetaceae','Saccharomycodaceae','Ustilaginaceae','Erysiphaceae'): return 'f:'+str(r.FAMILY)
    if r.CLASS in('Saccharomycetes','Pichiomycetes','Dipodascomycetes','Lipomycetes','Tremellomycetes','Pucciniomycetes','Microbotryomycetes','Ustilaginomycetes','Exobasidiomycetes','Agaricomycetes','Sordariomycetes','Lecanoromycetes','Eurotiomycetes','Pezizomycetes','Leotiomycetes'): return r.CLASS
    return None
L['lin']=L.apply(lin,axis=1)
T=pd.DataFrame({k:pd.Series((M[:,fidx['pfam|'+v]]>0)).groupby(L.lin.values).mean() for k,v in doms.items()}); T.insert(0,'n_species',L.lin.value_counts()); T.insert(1,'busco',L.groupby('lin').busco_complete.median())
print(T.round(2).to_string()); T.round(3).to_csv(f'{W}/key_loss_species_level.tsv',sep='\t')
# erysiphaceae etc: obligate biotroph losses
A=pd.read_csv(f'{W}/core_loss.tsv',sep='\t'); print('\ncore domains lost in >=12 lineages:'); print('; '.join((A[A.n_lin_lost>=12].sort_values('n_lin_lost',ascending=False).name+'('+A[A.n_lin_lost>=12].sort_values('n_lin_lost',ascending=False).n_lin_lost.astype(str)+')').tolist()))
# Pezizomycetes ankyrin, Glomeromycota expansions
G,P,Md=load('q')
for fid,cl in [('pfam|PF12796','Pezizomycetes'),('pfam|PF30695','Glomeromycetes'),('pfam|PF07534','Glomeromycetes'),('pfam|PF08238','Glomeromycetes'),('pfam|PF00069','Glomeromycetes'),('pfam|PF07714','Glomeromycetes'),('pfam|PF03106','Mucoromycetes'),('pfam|PF20152','Agaricomycetes'),('pfam|PF20684','Lecanoromycetes'),('pfam|PF00109','Lecanoromycetes')]:
    m=(G.CLASS==cl).values; v=Md[m,fidx[fid]]; o=Md[~m,fidx[fid]]
    print(fid,F.name[fidx[fid]],cl,'n',m.sum(),'median',np.median(v),'IQR',np.percentile(v,25),np.percentile(v,75),'| elsewhere median',np.median(o), 'top genera:',', '.join(f'{g}:{int(c)}' for g,c in sorted(zip(G.GENUS[m],v),key=lambda t:-t[1])[:4]))
R=pd.read_parquet(f'{W}/restricted_all.parquet')
print('\nzero-outside (out_n_all==0) & in_frac_q>=0.8 by clade:'); z=R[(R.out_n_all==0)&(R.in_frac_q>=0.8)]; print(z.groupby(['rank','clade']).size().sort_values(ascending=False).head(25).to_string())
print('restricted_best: DUF share', pd.read_csv(f'{W}/restricted_best.tsv',sep='\t').name.str.startswith('DUF').mean())
