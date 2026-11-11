from lib import *
Gq,Pq,Mq=load('q'); Ga,Pa,Ma=load('all')
pf=(F.ftype=='pfam').values
Fp=F[pf].reset_index(drop=True); Bq=(Pq[:,pf]>=0.5); Ba=(Pa[:,pf]>=0.5); Mq=Mq[:,pf]
rows=[]
for rank in ['DIKARYA','PHYLUM','SUBGROUP','CLASS','ORDER_']:
    for cl,nq in Gq[rank].value_counts().items():
        if nq<5 or cl in ('NA','non-Dikarya','Basidio_other') or 'incertae' in cl: continue
        if rank=='ORDER_':
            pass
        iq=(Gq[rank]==cl).values; ia=(Ga[rank]==cl).values
        if rank=='ORDER_':  # order names may repeat as 'NA'
            if cl=='NA': continue
        inf=Bq[iq].mean(0); infa=Ba[ia].mean(0); outn=Ba[~ia].sum(0); outf=outn/(~ia).sum(); outnq=Bq[~iq].sum(0)
        med=np.median(Mq[iq],axis=0)
        d=pd.DataFrame({'rank':rank,'clade':cl,'n_in_q':nq,'n_in_all':ia.sum(),'pfam':Fp.feat,'name':Fp.name,'desc':Fp.desc,'clan':Fp.clan_id,'in_frac_q':inf,'in_frac_all':infa,'out_n_all':outn,'out_frac_all':outf,'out_n_q':outnq,'med_copies':med})
        rows.append(d[(d.in_frac_q>=0.5)])
R=pd.concat(rows); R['J']=R.in_frac_q-R.out_frac_all
R.to_parquet(f'{W}/restricted_all.parquet')
# best clade per domain
strict=R[(R.in_frac_q>=0.8)&(R.out_frac_all<=0.02)].copy()
best=strict.sort_values(['J','n_in_q'],ascending=False).drop_duplicates('pfam')
best.to_csv(f'{W}/restricted_best.tsv',sep='\t',index=False)
print('n strict domain-clade pairs',len(strict),'unique domains',len(best))
print(best.groupby(['rank','clade']).agg(n=('pfam','size'),n_in=('n_in_q','first')).sort_values('n',ascending=False).head(60))
for (rk,cl),d in best.groupby(['rank','clade']):
    if len(d)<1: continue
    print(f'\n## {rk} {cl} (n_q={d.n_in_q.iloc[0]}, n_all={d.n_in_all.iloc[0]}) : {len(d)} domains')
    print(d.sort_values('J',ascending=False).head(25)[['pfam','name','desc','in_frac_q','in_frac_all','out_n_all','med_copies']].to_string(index=False))
