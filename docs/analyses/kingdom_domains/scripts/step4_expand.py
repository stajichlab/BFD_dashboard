from lib import *
from scipy import stats
G,P,Md=load('q')
G['life']=G.primary_lifestyle
vc=G.life.value_counts(); G.loc[~G.life.isin(vc[vc>=10].index)&G.life.notna(),'life']='other'
print(G.life.value_counts(dropna=False).to_string())
prev=(Md>0).mean(0); mean=Md.mean(0)
sel=(prev>=0.25)&(mean>=1.0)
Fs=F[sel].reset_index(drop=True); X=np.log10(Md[:,sel]+1); x=np.log10(G.n_proteins.astype(float).values)
print('features tested by type:\n',Fs.ftype.value_counts().to_string())
# class-balanced WLS allometry
w=1/G.groupby('CLASS').GENUS.transform('size').values; w=w/w.sum()
xm=(w*x).sum(); Xm=(w[:,None]*X).sum(0); sxx=(w*(x-xm)**2).sum()
b=((w*(x-xm))[:,None]*(X-Xm)).sum(0)/sxx; a=Xm-b*xm
R=X-(a+np.outer(x,b))
r2=1-(w[:,None]*R**2).sum(0)/(w[:,None]*(X-Xm)**2).sum(0)
def dummies(s):
    d=pd.get_dummies(s).values.astype(float); return d
def r2_fe(Y,D):
    D1=np.c_[np.ones(len(D)),D]; beta=np.linalg.lstsq(D1,Y,rcond=None)[0]; res=Y-D1@beta
    return 1-(res**2).sum(0)/((Y-Y.mean(0))**2).sum(0), res
cls=G.CLASS.where(G.groupby('CLASS').GENUS.transform('size')>=5,'small')
eta_class,_=r2_fe(R,dummies(cls)[:,1:])
# lifestyle within class / within order (genera with lifestyle label)
m=G.life.notna().values
def partial(fe):
    Y=R[m]; D=dummies(fe[m])[:,1:]; Lf=dummies(G.life[m]); names=pd.get_dummies(G.life[m]).columns
    _,Yr=r2_fe(Y,D); _,Lr=r2_fe(Lf[:,1:],D)
    beta=np.linalg.lstsq(Lr,Yr,rcond=None)[0]; res=Yr-Lr@beta
    ss0=(Yr**2).sum(0); ss1=(res**2).sum(0); k=Lr.shape[1]; df2=len(Y)-D.shape[1]-1-k
    Fst=((ss0-ss1)/k)/(ss1/df2); p=stats.f.sf(Fst,k,df2)
    # per-lifestyle mean of FE-centered residual
    Lm=np.array([Yr[(G.life[m]==n).values].mean(0) for n in names])
    return 1-ss1/ss0, p, Lm, names
ordr=(G.CLASS+'|'+G.ORDER_); ordr=ordr.where(ordr.map(ordr.value_counts())>=3, G.CLASS+'|small')
pe_c,p_c,Lm_c,names=partial(cls); pe_o,p_o,Lm_o,_=partial(ordr)
eta_life,_=r2_fe(R[m],dummies(G.life[m])[:,1:])
def bh(p):
    o=np.argsort(p); q=np.empty_like(p); q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]; return q
# class means
big=[c for c,n in G.CLASS.value_counts().items() if n>=5]
CM=np.array([R[(G.CLASS==c).values].mean(0) for c in big])
T=Fs[['ftype','feat','name','desc']].copy()
T['prev']=prev[sel]; T['med_all']=np.median(Md[:,sel],axis=0); T['slope']=b; T['allom_r2']=r2
T['eta2_class']=eta_class; T['sd_class_means']=CM.std(0); T['top_class']=[big[i] for i in CM.argmax(0)]; T['top_resid']=CM.max(0); T['low_class']=[big[i] for i in CM.argmin(0)]; T['low_resid']=CM.min(0)
T['eta2_life_raw']=eta_life; T['pR2_life_class']=pe_c; T['q_life_class']=bh(p_c); T['pR2_life_order']=pe_o; T['q_life_order']=bh(p_o)
T['top_life']=[names[i] for i in Lm_o.argmax(0)]; T['top_life_eff']=Lm_o.max(0); T['low_life']=[names[i] for i in Lm_o.argmin(0)]; T['low_life_eff']=Lm_o.min(0)
T.to_csv(f'{W}/expansion_stats.tsv',sep='\t',index=False)
np.save(f'{W}/resid_q.npy',R); Fs.to_parquet(f'{W}/resid_features.parquet'); pd.DataFrame(CM.T,columns=big).assign(fid=Fs.fid).to_parquet(f'{W}/class_mean_resid.parquet')
pd.DataFrame(Lm_o.T,columns=names).assign(fid=Fs.fid).to_parquet(f'{W}/life_eff_order.parquet')
c=['ftype','feat','name','desc','med_all','slope','eta2_class','sd_class_means','top_class','top_resid','low_class','low_resid','pR2_life_order','top_life','top_life_eff']
T['desc']=T.desc.fillna('').str[:40]
print('\nsummary of slopes by type'); print(T.groupby('ftype').slope.describe().round(2).to_string())
print('\nmedian eta2_class, pR2_life_class, pR2_life_order by type'); print(T.groupby('ftype')[['eta2_class','eta2_life_raw','pR2_life_class','pR2_life_order']].median().round(3).to_string())
print('n features q_life_order<0.01:',(T.q_life_order<0.01).sum(),'of',len(T), ' with pR2_life_order>=0.05:',((T.q_life_order<0.01)&(T.pR2_life_order>=0.05)).sum())
print('\n== TOP 60 pfam by between-class SD of residual (med_all>=2)'); print(T[(T.ftype=='pfam')&(T.med_all>=2)].sort_values('sd_class_means',ascending=False).head(60)[c].round(2).to_string(index=False))
print('\n== TOP 50 pfam by lifestyle partial R2 within order'); print(T[(T.ftype=='pfam')].sort_values('pR2_life_order',ascending=False).head(50)[c+['low_life','low_life_eff','q_life_order']].round(3).to_string(index=False))
for ft in ['cazy_sub','cazy_class','cazy_fam','merops','tf','trait']:
    print(f'\n== {ft}: by sd_class_means'); print(T[T.ftype==ft].sort_values('sd_class_means',ascending=False).head(22)[c+['low_life','low_life_eff']].round(2).to_string(index=False))
    if ft in('cazy_fam','merops','cazy_sub'):
        print(f'\n== {ft}: by lifestyle pR2 within order'); print(T[T.ftype==ft].sort_values('pR2_life_order',ascending=False).head(15)[c+['low_life','low_life_eff']].round(2).to_string(index=False))
