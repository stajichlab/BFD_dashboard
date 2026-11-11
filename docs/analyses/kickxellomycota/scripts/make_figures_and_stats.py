#!/usr/bin/env python3
"""Figures + genus-level effect sizes for the Kickxellomycota deep dive.
Inputs (kx_data/): genus_medians_all.tsv, kx_species_reps.tsv (from atlas A1/A2/A5 intermediates, scripts k1/k3/k5),
md_*.tsv (group-level results of the MotherDuck queries listed in the report).
Run: python3 kx_scripts/make_figures_and_stats.py   (from this analysis directory; kx_ paths map to data/, figs/, scripts/)"""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
D='kx_data/'; F='kx_figs/'
PAL=dict(blue='#2a78d6',orange='#eb6834',aqua='#1baf7a',yellow='#eda100',magenta='#e87ba4',green='#008300',violet='#6250d6',red='#e34948',grey='#b9b8b2')
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':'#e6e5e0','grid.linewidth':0.6,'axes.axisbelow':True,'figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb','savefig.dpi':170,'savefig.bbox':'tight'})
rng=np.random.default_rng(1)
gm=pd.read_csv(D+'genus_medians_all.tsv',sep='\t'); kx=pd.read_csv(D+'kx_species_reps.tsv',sep='\t')

def cls(g):
    if g.startswith('Kx:'): return 'Kickxellomycota'
    if g.startswith('Sacch:') and g not in('Sacch:Lipomycetes',): return 'Saccharomycotina yeasts'
    if g in('Schizosaccharomyces','Pneumocystis','Taphrinomycetes','Malassezia','Ustilaginomycetes','Wallemia'): return 'Other reduced Dikarya'
    if g in('Microsporidia','Rozello+Aphelid'): return 'Microsporidia / Rozello+Aphelid'
    if g in('Mucorales','Umbelopsidales','Endogonales','Mortierellales','Glomeromycota'): return 'Mucoro-/Mortierello-/Glomeromycota'
    if g in('Zoopagales','Entomophthorales','Basidiobolus'): return 'Zoopago-/Entomophthoro-/Basidiobolomycota'
    if g in('Chytrid+Monobleph','Blastocladiomycota','Neocallimastigomycota'): return 'Zoosporic'
    return 'Other Dikarya (filamentous etc.)'
gm['cls']=gm.grp.map(cls)
order=[('Other Dikarya (filamentous etc.)',PAL['grey'],'o'),('Zoosporic',PAL['aqua'],'s'),('Mucoro-/Mortierello-/Glomeromycota',PAL['blue'],'o'),('Zoopago-/Entomophthoro-/Basidiobolomycota',PAL['violet'],'D'),('Saccharomycotina yeasts',PAL['yellow'],'^'),('Other reduced Dikarya',PAL['green'],'v'),('Microsporidia / Rozello+Aphelid',PAL['magenta'],'X'),('Kickxellomycota',PAL['red'],'o')]
# ---- Fig 1
fig,ax=plt.subplots(1,2,figsize=(12,5))
for a,(x,y,xl,yl) in zip(ax,[('mb','introns_per_transcript','Genome size (Mb, log)','Introns per transcript'),('median_intergenic','median_intron_len','Median intergenic distance (bp, log)','Median intron length (bp, log)')]):
    for c,col,m in order:
        s=gm[gm.cls==c]; big=c=='Kickxellomycota'
        a.scatter(s[x],s[y],s=46 if big else 14,c=col,marker=m,alpha=1 if big else .75,edgecolor='#fcfcfb' if not big else '#0b0b0b',linewidth=.6,label=f'{c} (n={len(s)})',zorder=3 if big else 2)
    for _,r in gm[gm.cls=='Kickxellomycota'].iterrows():
        a.annotate(r.GENUS,(r[x],r[y]),xytext=(4,3),textcoords='offset points',fontsize=7,color='#0b0b0b')
    a.set_xscale('log'); a.set_xlabel(xl); a.set_ylabel(yl)
    if y=='median_intron_len': a.set_yscale('log')
ax[0].legend(fontsize=7,frameon=False,loc='upper left'); ax[0].set_ylim(-0.2,8)
fig.suptitle('Genus medians of species representatives (one point per genus; atlas A2 gene-structure intermediates)',fontsize=10)
fig.savefig(F+'fig1_genome_size_vs_introns.png'); plt.close(fig)
# ---- Fig 2: within Kickxellomycota
gorder=['Coemansia','Kickxella','Linderina','Martensiomyces','Dipsacomyces','Mycoemilia','Spiromyces','Ramicandelaber','Dimargaris','Dispira','Tieghemiomyces','Smittium','Furculomyces','Capniomyces','Zancudomyces']
gcol={**{g:PAL['red'] for g in gorder[:5]},**{g:PAL['orange'] for g in gorder[5:7]},'Ramicandelaber':PAL['yellow'],**{g:PAL['violet'] for g in gorder[8:11]},**{g:PAL['blue'] for g in gorder[11:]}}
mets=[('mb','Genome size (Mb)',True),('n_protein_coding','Protein-coding genes',False),('introns_per_transcript','Introns / transcript',False),('median_intron_len','Median intron length (bp)',False),('median_intergenic','Median intergenic (bp)',True),('masked_pct','Repeat-masked %',False)]
fig,axs=plt.subplots(1,6,figsize=(17,5.2),sharey=True)
for a,(m,lab,lg) in zip(axs,mets):
    for i,g in enumerate(gorder):
        s=kx[kx.GENUS==g]; y=i+rng.uniform(-.28,.28,len(s)); hq=s.N50_bp>=50000
        a.scatter(s[m][~hq],y[~hq],s=9,facecolor='none',edgecolor=gcol[g],linewidth=.7)
        a.scatter(s[m][hq],y[hq],s=16,c=gcol[g],edgecolor='#0b0b0b',linewidth=.4)
    if lg: a.set_xscale('log')
    a.set_xlabel(lab); a.set_yticks(range(len(gorder))); a.set_yticklabels([f'{g} (n={int((kx.GENUS==g).sum())})' for g in gorder]); a.invert_yaxis() if a is axs[0] else None
fig.suptitle('Kickxellomycota species representatives (n=169). Filled = N50 >= 50 kb; open = more fragmented. Red: core Kickxellales; orange: early Kickxellales; yellow: Ramicandelaberales; violet: Dimargaritales; blue: Harpellales',fontsize=9)
fig.savefig(F+'fig2_kickxellomycota_by_genus.png'); plt.close(fig)
# ---- Fig 3 quality
c=kx[kx.GENUS=='Coemansia']
fig,axs=plt.subplots(1,3,figsize=(13,3.8))
for a,(y,lab) in zip(axs,[('busco_complete','BUSCO complete % (fungi_odb12)'),('n_protein_coding','Protein-coding genes'),('introns_per_transcript','Introns / transcript')]):
    a.scatter(c.N50_bp,c[y],s=14,c=PAL['red'],edgecolor='#fcfcfb',linewidth=.4); a.set_xscale('log'); a.set_xlabel('Contig N50 (bp, log)'); a.set_ylabel(lab)
    r,p=stats.spearmanr(c.N50_bp,c[y]); a.set_title(f'Spearman rho={r:.2f} (n={len(c)})',fontsize=9)
fig.suptitle('Coemansia (146 species representatives): which statistics track assembly contiguity',fontsize=10)
fig.savefig(F+'fig3_coemansia_quality.png'); plt.close(fig)
# ---- Fig 4 position deciles
dec=pd.read_csv(D+'md_intron_position_deciles.tsv',sep='\t'); dc=[f'd{i}' for i in range(1,11)]
sel=[('Kx:Coemansia (N50>=50kb)',PAL['red'],'-'),('Kx:Harpellales',PAL['blue'],'-'),('Sacch:Dipodascomycetes',PAL['yellow'],'-'),('Schizosaccharomyces',PAL['green'],'-'),('Mucorales',PAL['violet'],'--'),('Zoopagales',PAL['aqua'],'--')]
fig,axs=plt.subplots(1,2,figsize=(11,4),sharey=True)
for a,sg,t in zip(axs,[1,0],['Genes with exactly one intron','Genes with >= 2 introns']):
    for g,col,ls in sel:
        r=dec[(dec.grp==g)&(dec.single==sg)].iloc[0]; a.plot(range(1,11),r[dc].values.astype(float),color=col,ls=ls,lw=2,marker='o',ms=4,label=g)
    a.axhline(.1,color='#52514e',lw=.8,ls=':'); a.set_title(t,fontsize=9); a.set_xlabel("Decile of CDS length (1 = 5' end)"); a.set_xticks(range(1,11))
axs[0].set_ylabel('Fraction of introns'); axs[1].legend(frameon=False,fontsize=8)
fig.suptitle('Intron position along the CDS (complete gene models; positions recomputed from gene_CDS)',fontsize=10)
fig.savefig(F+'fig4_intron_position.png'); plt.close(fig)
# ---- Fig 5 ortholog matched
om=pd.read_csv(D+'md_ortholog_matched_introns.tsv',sep='\t').sort_values('ratio_sum')
fig,ax=plt.subplots(figsize=(7.5,6.2))
cols=[PAL['red'] if g.startswith('Kx:') else PAL['yellow'] if g.startswith('Sacch') else PAL['grey'] for g in om.grp]
ax.barh(range(len(om)),om.ratio_sum,color=cols,height=.72); ax.axvline(1,color='#52514e',lw=.8,ls=':')
ax.set_yticks(range(len(om))); ax.set_yticklabels([f'{g} (n={n})' for g,n in zip(om.grp,om.n_fam)]); ax.invert_yaxis()
for i,(v,z) in enumerate(zip(om.ratio_sum,om.f_zero)): ax.text(v+.02,i,f'{v:.2f}  ({z:.0%} intronless)',va='center',fontsize=7,color='#52514e')
ax.set_xlabel('Introns per ortholog relative to Mucorales (ratio of sums over shared SwissProt-anchored families)'); ax.set_xlim(0,1.75)
ax.set_title('Ortholog-matched intron retention (n = families shared with Mucorales)\nred: Kickxellomycota; yellow: Saccharomycotina',fontsize=9)
fig.savefig(F+'fig5_ortholog_matched_introns.png'); plt.close(fig)
# ---- Fig 6 heatmap
pf=pd.read_csv(D+'md_pfam_presence_curated.tsv',sep='\t'); gc=list(pf.columns[2:])
fig,ax=plt.subplots(figsize=(9.5,10))
im=ax.imshow(pf[gc].values,cmap='Blues',vmin=0,vmax=1,aspect='auto')
ax.set_xticks(range(len(gc))); ax.set_xticklabels(gc,rotation=55,ha='right'); ax.set_yticks(range(len(pf))); ax.set_yticklabels([f'{c}: {p}' for c,p in zip(pf.category,pf.pfam_id)],fontsize=8); ax.grid(False)
for i in range(len(pf)):
    for j in range(len(gc)):
        v=pf[gc].values[i,j]; ax.text(j,i,f'{v:.2f}'.rstrip('0').rstrip('.') if v not in(0,1) else str(int(v)),ha='center',va='center',fontsize=6,color='#fcfcfb' if v>.55 else '#0b0b0b')
ax.axvline(4.5,color='#0b0b0b',lw=1); ax.axvline(7.5,color='#0b0b0b',lw=1)
plt.colorbar(im,ax=ax,shrink=.4,label='Fraction of species representatives with >= 1 hit')
ax.set_title('Presence of selected Pfam domains (Coemansia restricted to BUSCO >= 85 % genomes, n=41)',fontsize=9)
fig.savefig(F+'fig6_pfam_presence_heatmap.png'); plt.close(fig)
# ---- Fig 7 + Fisher: convergence of losses among 2,471 core Pfams
N=2471; nk=41
ov=pd.DataFrame([('Harpellales',85,33),('Dimargaritales',46,19),('Zoopagales',26,6),('Sacch. (Saccharomycetes+Pichiomycetes)',57,4),('Schizosaccharomyces',124,18),('Malassezia',63,4),('Ustilaginomycetes',13,1),('Pneumocystis',330,13),('Microsporidia',1271,36)],columns=['other','n_other','both'])
rows=[]
for _,r in ov.iterrows():
    a=r.both; b=nk-a; c_=r.n_other-a; d=N-nk-c_
    OR,p=stats.fisher_exact([[a,b],[c_,d]],alternative='two-sided'); rows.append((r.other,r.n_other,a,round(nk*r.n_other/N,2),round(OR,1),p,round(a/(nk+r.n_other-a),3)))
ft=pd.DataFrame(rows,columns=['lineage','n_lost_in_lineage','n_shared_with_core_Kickxellales','expected_if_independent','odds_ratio','fisher_p','jaccard'])
# reference pairs among yeasts
for nm,(n1,n2,b) in {'Sacch vs Schizosaccharomyces':(57,124,23),'Sacch vs Malassezia':(57,63,10),'Sacch vs Harpellales':(57,85,16),'Sacch vs Zoopagales':(57,26,1)}.items():
    OR,p=stats.fisher_exact([[b,n1-b],[n2-b,N-n1-n2+b]]); ft.loc[len(ft)]=[nm+' (reference pair)',n2,b,round(n1*n2/N,2),round(OR,1),p,round(b/(n1+n2-b),3)]
ft.to_csv(D+'stats_loss_overlap_fisher.tsv',sep='\t',index=False); print(ft.to_string())
f9=ft.iloc[:9]
fig,ax=plt.subplots(figsize=(8,4)); y=np.arange(9)
ax.barh(y-.19,f9.n_shared_with_core_Kickxellales,height=.36,color=PAL['red'],label='observed shared losses'); ax.barh(y+.19,f9.expected_if_independent,height=.36,color=PAL['grey'],label='expected if independent')
ax.set_yticks(y); ax.set_yticklabels([f'{l} (lost {n})' for l,n in zip(f9.lineage,f9.n_lost_in_lineage)]); ax.invert_yaxis(); ax.legend(frameon=False); ax.set_xlabel('Core Pfam domains (of 2,471) lost in both core Kickxellales (41 lost) and the lineage')
for i,(o,OR) in enumerate(zip(f9.n_shared_with_core_Kickxellales,f9.odds_ratio)): ax.text(o+.4,i-.19,f'OR {OR}',va='center',fontsize=7,color='#52514e')
fig.savefig(F+'fig7_loss_overlap.png'); plt.close(fig)
# ---- genus-level effect sizes
def cliff(a,b):
    a=np.asarray(a); b=np.asarray(b); return (np.sign(a[:,None]-b[None,:])).sum()/(len(a)*len(b))
def boot(a,b,n=5000):
    a=np.asarray(a); b=np.asarray(b); d=[np.median(rng.choice(a,len(a)))-np.median(rng.choice(b,len(b))) for _ in range(n)]; return np.percentile(d,[2.5,97.5])
core=['Coemansia','Dipsacomyces','Kickxella','Linderina','Martensiomyces']
sets={'core Kickxellales (5 genera)':gm[gm.GENUS.isin(core)&gm.grp.str.startswith('Kx')],'Harpellales (4 genera)':gm[gm.grp=='Kx:Harpellales'],'Dimargaritales (3 genera)':gm[gm.grp=='Kx:Dimargaritales']}
refs={'Mucorales':gm[gm.grp=='Mucorales'],'Mortierellales':gm[gm.grp=='Mortierellales'],'Zoopagales':gm[gm.grp=='Zoopagales'],'Sacch:Saccharomycetes':gm[gm.grp=='Sacch:Saccharomycetes'],'Sacch:Dipodascomycetes':gm[gm.grp=='Sacch:Dipodascomycetes']}
out=[]
for sn,s in sets.items():
    for rn,r in refs.items():
        for m in ['mb','n_protein_coding','introns_per_transcript','median_intron_len','median_intergenic','median_protein_len','masked_pct','sec_k','cazy_k','tf_k']:
            a=s[m].dropna(); b=r[m].dropna(); lo,hi=boot(a,b); u=stats.mannwhitneyu(a,b,alternative='two-sided')
            out.append((sn,rn,m,len(a),len(b),round(a.median(),2),round(b.median(),2),round(a.median()-b.median(),2),round(lo,2),round(hi,2),round(np.log2(a.median()/b.median()),2),round(cliff(a,b),2),u.pvalue))
es=pd.DataFrame(out,columns=['focal','reference','metric','n_focal_genera','n_ref_genera','median_focal','median_ref','median_diff','boot_lo','boot_hi','log2_ratio','cliffs_delta','mannwhitney_p'])
es.to_csv(D+'stats_genus_effect_sizes.tsv',sep='\t',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_rows',400)
print(es[es.reference.isin(['Mucorales','Zoopagales','Sacch:Dipodascomycetes'])&es.metric.isin(['mb','n_protein_coding','introns_per_transcript','median_intron_len','median_intergenic','median_protein_len'])].to_string())
# within-Coemansia covariation (HQ only) and across Kx genera
hq=c[c.N50_bp>=50000]; v=['mb','n_protein_coding','introns_per_transcript','median_intron_len','median_intergenic','masked_pct','gc_pct']
print('Coemansia HQ n=',len(hq)); print(hq[v].corr(method='spearman').round(2).to_string())
hq[v].corr(method='spearman').round(2).to_csv(D+'stats_coemansia_hq_spearman.tsv',sep='\t')
