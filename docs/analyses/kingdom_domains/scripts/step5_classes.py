from lib import *
from scipy import stats
G,P,Md=load('q')
G['n_proteins']=G.n_proteins.astype(float)
def grp(r):
    if r.SUBGROUP=='Pezizomycotina': return r.CLASS if r.CLASS in ('Sordariomycetes','Dothideomycetes','Eurotiomycetes','Leotiomycetes','Lecanoromycetes','Pezizomycetes','Orbiliomycetes') else 'Pezizomycotina_other'
    if r.SUBGROUP=='Agaricomycotina': return r.CLASS if r.CLASS in('Agaricomycetes','Tremellomycetes') else 'Agaricomycotina_other'
    return r.SUBGROUP
G['grp']=G.apply(grp,axis=1)
order=['Sordariomycetes','Leotiomycetes','Dothideomycetes','Eurotiomycetes','Lecanoromycetes','Pezizomycetes','Orbiliomycetes','Pezizomycotina_other','Saccharomycotina','Taphrinomycotina','Agaricomycetes','Tremellomycetes','Agaricomycotina_other','Ustilaginomycotina','Pucciniomycotina','Mucoromycota_sl','Zoopagomycota_sl','zoosporic','Microsporidia']
G.to_parquet(f'{W}/genera_q_grp.parquet')
sc=['n_proteins','genome_mb','sec_frac','gpi_frac','tm_frac','idp','frac_pfam','med_len']
t=G.groupby('grp')[sc].median().reindex(order); t.insert(0,'n_genera',G.grp.value_counts().reindex(order))
print(t.round(4).to_string()); t.to_csv(f'{W}/scalar_by_group.tsv',sep='\t')
# eta2 class and lifestyle-within-order for scalars
G['life']=G.primary_lifestyle; vc=G.life.value_counts(); G.loc[~G.life.isin(vc[vc>=10].index)&G.life.notna(),'life']='other'
def r2(y,D):
    D1=np.c_[np.ones(len(D)),D]; res=y-D1@np.linalg.lstsq(D1,y,rcond=None)[0]; return 1-(res**2).sum()/((y-y.mean())**2).sum(),res
cls=G.CLASS.where(G.groupby('CLASS').GENUS.transform('size')>=5,'small')
ordr=(G.CLASS+'|'+G.ORDER_); ordr=ordr.where(ordr.map(ordr.value_counts())>=3, G.CLASS+'|small')
print('\nscalar: eta2_class, eta2_life_raw, partialR2 life|order, lifestyle effects within order')
for s in sc:
    ok=G[s].notna().values&G.life.notna().values; y=G[s].values[ok].astype(float)
    if s in('n_proteins','genome_mb'): y=np.log10(y)
    e_c,_=r2(y,pd.get_dummies(cls[ok]).values[:,1:].astype(float)); e_l,_=r2(y,pd.get_dummies(G.life[ok]).values[:,1:].astype(float))
    D=pd.get_dummies(ordr[ok]).values[:,1:].astype(float); _,yr=r2(y,D); Ld=pd.get_dummies(G.life[ok]); 
    Lr=np.column_stack([r2(Ld.values[:,j].astype(float),D)[1] for j in range(1,Ld.shape[1])]); pr,_=r2(yr,Lr)
    eff={n:yr[(G.life[ok]==n).values].mean() for n in Ld.columns}; es=sorted(eff.items(),key=lambda kv:kv[1])
    print(f'{s:12s} eta2_class={e_c:.2f} eta2_life={e_l:.2f} pR2_life|order={pr:.3f}  low: {es[0][0]} {es[0][1]:+.4f}, {es[1][0]} {es[1][1]:+.4f}  high: {es[-1][0]} {es[-1][1]:+.4f}, {es[-2][0]} {es[-2][1]:+.4f}')
# selected families per 1000 proteins, group medians
selp={'PF05729':'NACHT','PF00931':'NB-ARC','PF06985':'HET','PF24883':'NPHP3_N (NLR NTPase)','PF00069':'Pkinase','PF00172':'Zn_clus','PF04082':'Fungal_trans','PF00067':'p450','PF00109':'ketoacyl-synt (PKS/FAS)','PF00668':'Condensation (NRPS)','PF00501':'AMP-binding','PF03936':'Terpene_synth_C','PF07690':'MFS_1','PF00083':'Sugar_tr','PF00005':'ABC_tran','PF01185':'Hydrophobin','PF06766':'Hydrophobin_2','PF01476':'LysM','PF05730':'CFEM','PF07249':'Cerato-platanin','PF20684':'Fung_rhodopsin (PTH11-like 7TM)','PF00082':'Peptidase_S8','PF00026':'Asp protease','PF03443':'AA9 LPMO','PF00734':'CBM_1','PF00704':'Glyco_hydro_18 chitinase','PF03106':'WRKY','PF20152':'DUF6534','PF00078':'RVT_1','PF12796':'Ank_2'}
rows=[]
fidx={f:i for i,f in enumerate(F.fid)}
per=Md/G.n_proteins.values[:,None]*1000
def add(fid,label):
    i=fidx.get(fid)
    if i is None: return
    r=pd.Series(per[:,i]).groupby(G.grp.values).median().reindex(order); r.name=label; rows.append(r)
for k,v in selp.items(): add('pfam|'+k,f'{k} {v}')
for k in ['cazy_class|GH','cazy_class|AA','cazy_class|CE','cazy_class|PL','cazy_class|GT','cazy_class|CBM','cazy_sub|cellulose','cazy_sub|xylan','cazy_sub|pectin','cazy_sub|lignin','cazy_sub|chitin','merops|A01A','merops|S08A','merops|S53','merops|M36','merops|S09X','merops|S10','tf|Zn(II)2Cys6 zinc cluster','tf|C2H2 zinc finger','tf|Homeodomain','tf|Basic helix-loop-helix','tf|GATA zinc finger','tf|High mobility group box']: add(k,k)
ST=pd.DataFrame(rows); ST.to_csv(f'{W}/selected_by_group.tsv',sep='\t')
print('\ncopies per 1000 proteins, median of genus medians'); print(ST.round(2).to_string())
# ECM sanity: n_proteins and frac_pfam by lifestyle within Agaricomycetes
a=G[G.CLASS=='Agaricomycetes']
print('\nAgaricomycetes by lifestyle'); print(a.groupby('life').agg(n=('GENUS','size'),n_prot=('n_proteins','median'),genome_mb=('genome_mb','median'),frac_pfam=('frac_pfam','median'),sec=('sec_frac','median')).round(3).to_string())
# absolute counts for key PCWDE in Agaricomycetes by lifestyle (not normalized)
for fid in ['cazy_sub|cellulose','cazy_sub|pectin','cazy_fam|AA9','cazy_fam|GH7','pfam|PF00106','pfam|PF00067','pfam|PF13409','pfam|PF17921','pfam|PF00069']:
    i=fidx[fid]; print(fid, pd.Series(Md[(G.CLASS=='Agaricomycetes').values,i]).groupby(a.life.values).median().round(1).to_dict())
