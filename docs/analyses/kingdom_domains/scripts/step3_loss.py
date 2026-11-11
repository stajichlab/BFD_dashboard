from lib import *
import scipy.cluster.hierarchy as H
Gq,Pq,Mq=load('q'); Ga,Pa,Ma=load('all')
pf=(F.ftype=='pfam').values; Fp=F[pf].reset_index(drop=True); Bq=(Pq[:,pf]>=0.5); Ba=(Pa[:,pf]>=0.5)
# ancestral: >=50% genera in >=2 of 3 early-diverging groups (q set), or >=50% in one early group AND >=50% in both Asco and Basidio
ed=[ (Gq.SUBGROUP==g).values for g in ['zoosporic','Mucoromycota_sl','Zoopagomycota_sl']]
edf=np.array([Bq[m].mean(0) for m in ed]); asc=Bq[(Gq.PHYLUM=='Ascomycota').values].mean(0); bas=Bq[(Gq.PHYLUM=='Basidiomycota').values].mean(0)
anc=((edf>=0.5).sum(0)>=2)|(((edf>=0.5).sum(0)>=1)&(asc>=0.5)&(bas>=0.5))
print('n pfam',len(Fp),'ancestral',anc.sum(), 'ed n', [m.sum() for m in ed])
dik=Bq[(Gq.DIKARYA=='Dikarya').values].mean(0)
dl=anc&(dik<=0.05)
print('ancestral & Dikarya<=5%:',dl.sum())
D=Fp[dl].assign(zoosp=edf[0][dl],muco=edf[1][dl],zoop=edf[2][dl],dik=dik[dl])
D.to_csv(f'{W}/dikarya_lost.tsv',sep='\t',index=False)
pass
print('of which in zoosporic only (flagellar-like, muco<0.1):',((D.muco<0.1)).sum(),' retained in Mucoromycota_sl>=0.5:',(D.muco>=0.5).sum())
# lineage units
Ga['LIN']=Ga.CLASS
for g in ['Malassezia','Pneumocystis','Schizosaccharomyces','Wallemia','Hanseniaspora']:
    pass
Ga.loc[Ga.PHYLUM.isin(['Kickxellomycota']),'LIN']='Kickxellomycota'
Ga.loc[Ga.PHYLUM.isin(['Blastocladiomycota']),'LIN']='Blastocladiomycota'
Ga.loc[Ga.PHYLUM.isin(['Mucoromycota']),'LIN']='Mucoromycota'
Ga.loc[Ga.GENUS=='Hanseniaspora','LIN']='g_Hanseniaspora'
vc=Ga.LIN.value_counts()
special={'Malasseziomycetes','Pneumocystomycetes','Schizosaccharomycetes','Wallemiomycetes','g_Hanseniaspora','Neocallimastigomycetes','Pucciniomycetes','Taphrinomycetes'}
lins=[l for l in vc.index if (vc[l]>=5 or l in special) and 'NA' not in l and l!='Ascomycota_incertae']
pez=Bq[(Gq.SUBGROUP=='Pezizomycotina').values].mean(0); aga=Bq[(Gq.CLASS=='Agaricomycetes').values].mean(0)
anc0=anc.copy(); anc=anc&(pez>=0.7)&(aga>=0.7); print('core (ancestral & >=70% Pezizomycotina & >=70% Agaricomycetes):',anc.sum())
A=Fp[anc].reset_index(drop=True); Ba_a=Ba[:,anc]
LF=pd.DataFrame({l:Ba_a[(Ga.LIN==l).values].mean(0) for l in lins}); 
info=pd.DataFrame({'n_genera':[ (Ga.LIN==l).sum() for l in lins],'busco':[Ga[Ga.LIN==l].busco.median() for l in lins],'n_prot':[Ga[Ga.LIN==l].n_proteins.median() for l in lins]},index=lins)
lost=(LF<=0.1)
info['n_lost']=lost.sum(0).values; info['frac_lost']=info.n_lost/len(A)
# unique losses
nl=lost.sum(1)
info['n_unique_lost']=[ (lost[l]&(nl==1)).sum() for l in lins]
print(info.sort_values('n_lost',ascending=False).to_string())
A['n_lin_lost']=nl.values; A['lost_in']=[','.join([l for l in lins if lost.loc[i,l]]) for i in range(len(A))]
A.to_csv(f'{W}/core_loss.tsv',sep='\t',index=False); LF.to_parquet(f'{W}/lineage_frac.parquet')
# overlap matrix: log2 odds of co-loss
focal=['Microsporidia_incertae','Saccharomycetes','Pichiomycetes','Dipodascomycetes','Lipomycetes','g_Hanseniaspora','Schizosaccharomycetes','Pneumocystomycetes','Malasseziomycetes','Wallemiomycetes','Kickxellomycota','Zoopagomycetes','Entomophthoromycetes','Neocallimastigomycetes','Glomeromycetes','Pucciniomycetes','Tremellomycetes','Microbotryomycetes','Ustilaginomycetes','Exobasidiomycetes','Lecanoromycetes','Chytridiomycetes','Mucoromycota']
focal=[f for f in focal if f in lins]
N=len(A); rows=[]
import math
for a in focal:
    r=[]
    for b in focal:
        if a==b: r.append(np.nan); continue
        n11=(lost[a]&lost[b]).sum(); exp=lost[a].sum()*lost[b].sum()/N
        r.append(round(math.log2((n11+0.5)/(exp+0.5)),1))
    rows.append(r)
print('\nlog2(obs/exp) co-loss overlap'); print(pd.DataFrame(rows,index=focal,columns=[f[:6] for f in focal]).to_string())
rows=[[ (lost[a]&lost[b]).sum() for b in focal] for a in focal]
print('\nshared losses n'); print(pd.DataFrame(rows,index=focal,columns=[f[:6] for f in focal]).to_string())
# convergent losses: domains lost in >=4 lineages
print('\nDomains by number of lineages lost:'); print(nl.value_counts().sort_index().to_string())
# class-weighted co-absence modules at genus level (q set), excluding Microsporidia
keep=(Gq.PHYLUM!='Microsporidia').values
X=(~Bq[keep][:,anc]).astype(float); g=Gq[keep]
w=1/g.groupby('CLASS').GENUS.transform('size').values; w=w/w.sum()
mu=(w[:,None]*X).sum(0); var=mu*(1-mu)
sel=(mu>=0.03)&(mu<=0.7); print('variable ancestral domains for module analysis',sel.sum())
Xs=X[:,sel]-mu[sel]; C=(Xs*w[:,None]).T@Xs/np.sqrt(np.outer(var[sel],var[sel]))
np.fill_diagonal(C,1); d=1-C; d=(d+d.T)/2; d[d<0]=0
Z=H.linkage(d[np.triu_indices(len(d),1)],'average'); cl=H.fcluster(Z,0.35,'distance')
As=A[sel].reset_index(drop=True); As['module']=cl; As['w_abs']=mu[sel]
LFs=LF[sel].reset_index(drop=True)
sizes=As.module.value_counts(); As.to_csv(f'{W}/loss_modules.tsv',sep='\t',index=False)
print('modules size>=4:',(sizes>=4).sum())
for m in sizes[sizes>=4].index[:45]:
    sub=As[As.module==m]; lf=LFs[As.module==m].mean(0)
    absent=[f'{l}({info.n_genera[l]})' for l in lins if lf[l]<=0.2]
    print(f'\n# module {m} n={len(sub)} mean_w_abs={sub.w_abs.mean():.2f}\n  domains: '+'; '.join((sub.name+' ['+sub.desc.str[:38]+']').head(40))+f'\n  absent(<=0.2 mean presence) in: '+', '.join(absent))
np.save(f'{W}/coabs_C.npy',C); np.save(f'{W}/coabs_leaves.npy',H.leaves_list(Z)); info.to_csv(f'{W}/lineage_loss_info.tsv',sep='\t')
LFs.assign(module=As.module.values).groupby('module').mean().to_parquet(f'{W}/module_lineage_presence.parquet')
