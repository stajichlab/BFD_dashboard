"""Markdown tables for the report -> work/tables.json"""
from lib import *
import json, re
T={}
def md(df,floatfmt=2):
    df=df.copy()
    for c in df.columns:
        if df[c].dtype.kind=='f': df[c]=df[c].map(lambda v:'' if pd.isna(v) else (f'{v:.{floatfmt}f}'))
    h='| '+' | '.join(map(str,df.columns))+' |\n|'+'---|'*len(df.columns)+'\n'
    return h+'\n'.join('| '+' | '.join(str(v).replace('|','/') for v in r)+' |' for r in df.values)
best=pd.read_csv(f'{W}/restricted_best.tsv',sep='\t'); R=pd.read_parquet(f'{W}/restricted_all.parquet')
z=R[(R.out_n_all==0)&(R.in_frac_q>=0.8)].groupby(['rank','clade']).size().rename('n_zero_outside')
c=best.groupby(['rank','clade']).agg(n_genera_q=('n_in_q','first'),n_genera_all=('n_in_all','first'),n_restricted=('pfam','size'),DUF_share=('name',lambda s:s.str.startswith('DUF').mean())).join(z).fillna({'n_zero_outside':0}).sort_values('n_restricted',ascending=False).reset_index()
c['n_zero_outside']=c.n_zero_outside.astype(int); c['rank']=c['rank'].replace({'ORDER_':'order','CLASS':'class','PHYLUM':'phylum','SUBGROUP':'subphylum-level group','DIKARYA':'subkingdom'})
T['restr_counts']=md(c.head(26))
want=[('DIKARYA','Dikarya',5),('PHYLUM','Ascomycota',5),('SUBGROUP','Pezizomycotina',7),('CLASS','Sordariomycetes',2),('ORDER_','Pleosporales',2),('ORDER_','Hypocreales',1),('SUBGROUP','Saccharomycotina',6),('PHYLUM','Basidiomycota',6),('CLASS','Agaricomycetes',4),('CLASS','Ustilaginomycetes',4),('ORDER_','Ustilaginales',2),('PHYLUM','Mucoromycota',6),('PHYLUM','Mortierellomycota',5),('PHYLUM','Glomeromycota',6),('SUBGROUP','zoosporic',6),('PHYLUM','Kickxellomycota',3),('PHYLUM','Microsporidia',7)]
rows=[]
for rk,cl,k in want:
    d=best[(best['rank']==rk)&(best.clade==cl)].sort_values(['out_n_all','in_frac_q'],ascending=[True,False]).head(k); rows.append(d)
d=pd.concat(rows)[['clade','n_in_q','pfam','name','desc','in_frac_q','in_frac_all','out_n_all','med_copies']]; d['desc']=d.desc.str[:60]
d.columns=['clade','genera (q)','Pfam','name','description','frac genera in clade (q)','frac in clade (all)','genera outside with domain (of all)','median copies']
T['restr_top']=md(d)
info=pd.read_csv(f'{W}/lineage_loss_info.tsv',sep='\t',index_col=0).sort_values('n_lost',ascending=False).reset_index().rename(columns={'index':'lineage'})
info['lineage']=info.lineage.str.replace('_incertae','').str.replace('g_','genus ')
info['pct_lost']=info.frac_lost*100
T['loss_lineage']=md(info[['lineage','n_genera','busco','n_prot','n_lost','pct_lost','n_unique_lost']].rename(columns={'busco':'median BUSCO %','n_prot':'median proteins','n_lost':'core domains absent','pct_lost':'% of 3,568','n_unique_lost':'absent only here'}).head(24),1)
exec(re.search(r"mods=\{.*?\}\n",open('make_figures.py').read(),re.S).group(0))
MP=pd.read_parquet(f'{W}/module_lineage_presence.parquet'); LM=pd.read_csv(f'{W}/loss_modules.tsv',sep='\t'); li=pd.read_csv(f'{W}/lineage_loss_info.tsv',sep='\t',index_col=0)
rows=[]
for m,lab in mods.items():
    sub=LM[LM.module==m]; ab=[f"{l.replace('_incertae','').replace('g_','')} ({li.n_genera[l]})" for l in MP.columns if MP.loc[m,l]<=0.2]
    rows.append([lab.rsplit(' [',1)[0],len(sub),', '.join(sub.name.head(7)),len(ab),', '.join(ab)])
T['modules']=md(pd.DataFrame(rows,columns=['module (my label)','n domains','example Pfam domains','n lineages lost','lineages with mean presence <= 0.2 (n genera)']))
K=pd.read_csv(f'{W}/key_loss_species_level.tsv',sep='\t',index_col=0)
keep=['Microsporidia','g:Pneumocystis','g:Malassezia','g:Wallemia','g:Schizosaccharomyces','g:Hanseniaspora','f:Saccharomycetaceae','Saccharomycetes','Pichiomycetes','Dipodascomycetes','Lipomycetes','Tremellomycetes','g:Puccinia','f:Erysiphaceae','Glomeromycota','Zoopagomycota','Kickxellomycota','Neocallimastigomycota','Blastocladiomycota','Chytridiomycota','Mucoromycota','Sordariomycetes','Agaricomycetes']
K=K.loc[keep].drop(columns=['FAS_I_H','Dicer_dimer','HisG']).reset_index().rename(columns={'index':'lineage','FAS (Fas_alpha_ACP)':'FAS'}); K['n_species']=K.n_species.astype(int)
T['species_loss']=md(K)
E=pd.read_csv(f'{W}/expansion_stats.tsv',sep='\t'); E['desc']=E.desc.fillna('').str[:38]
cols=['feat','name','desc','med_all','slope','eta2_class','sd_class_means','top_class','top_resid','low_class','low_resid','pR2_life_order']
ren={'feat':'id','med_all':'median copies','eta2_class':'eta2 class','sd_class_means':'SD of class means','top_class':'highest class','top_resid':'resid','low_class':'lowest class','low_resid':'resid ','pR2_life_order':'pR2 lifestyle/order'}
T['exp_class']=md(E[(E.ftype=='pfam')&(E.med_all>=2)].sort_values('sd_class_means',ascending=False).head(30)[cols].rename(columns=ren))
c2=['ftype','feat','name','med_all','eta2_class','pR2_life_class','pR2_life_order','q_life_order','top_life','top_life_eff','low_life','low_life_eff']
x=E.sort_values('pR2_life_order',ascending=False); x=pd.concat([x[x.ftype=='pfam'].head(16),x[x.ftype=='cazy_fam'].head(6),x[x.ftype=='cazy_sub'].head(4),x[x.ftype=='merops'].head(3)]); x['q_life_order']=x.q_life_order.map(lambda v:f'{v:.1e}')
T['exp_life']=md(x[c2].rename(columns={'med_all':'median copies','eta2_class':'eta2 class','pR2_life_class':'pR2 life/class','pR2_life_order':'pR2 life/order','q_life_order':'BH q','top_life':'highest lifestyle','top_life_eff':'eff','low_life':'lowest lifestyle','low_life_eff':'eff '}),3)
s=E.groupby('ftype').agg(n=('feat','size'),median_slope=('slope','median'),eta2_class=('eta2_class','median'),eta2_life_raw=('eta2_life_raw','median'),pR2_life_class=('pR2_life_class','median'),pR2_life_order=('pR2_life_order','median')).reset_index()
T['var_summary']=md(s,3)
gc=['Sordariomycetes','Eurotiomycetes','Lecanoromycetes','Pezizomycetes','Saccharomycotina','Agaricomycetes','Tremellomycetes','Ustilaginomycotina','Pucciniomycotina','Mucoromycota_sl','Zoopagomycota_sl','zoosporic','Microsporidia']
ST=pd.read_csv(f'{W}/selected_by_group.tsv',sep='\t',index_col=0)[gc].reset_index().rename(columns={'index':'family'}); ST.columns=[c.replace('mycetes','.').replace('mycotina','.').replace('mycota_sl',' s.l.') for c in ST.columns]
T['selected']=md(ST)
SC=pd.read_csv(f'{W}/scalar_by_group.tsv',sep='\t',index_col=0).reset_index(); 
for c_ in ['sec_frac','gpi_frac','tm_frac','idp','frac_pfam']: SC[c_]=SC[c_].map(lambda v:f'{v:.3f}')
SC['n_proteins']=SC.n_proteins.round().astype(int); SC['genome_mb']=SC.genome_mb.map(lambda v:f'{v:.1f}'); SC['med_len']=SC.med_len.round().astype(int)
T['scalar']=md(SC.rename(columns={'grp':'group','n_proteins':'proteins','genome_mb':'genome Mb','sec_frac':'secretome','gpi_frac':'GPI','tm_frac':'TM>=2','idp':'median disorder','frac_pfam':'with Pfam','med_len':'median aa'}))
D=pd.read_csv(f'{W}/dikarya_lost.tsv',sep='\t'); T['n_dl']=str(len(D))
json.dump(T,open(f'{W}/tables.json','w'))
print({k:len(v) for k,v in T.items()})
