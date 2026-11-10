"""Figures for fable_B_kingdom_domains.md. Run: ./venv/bin/python make_figures.py  (needs work/ outputs of step1-5)."""
from lib import *
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
OUT=os.environ.get('KW_BASE', '.') + '/fable_B_figs'
import os; os.makedirs(OUT,exist_ok=True)
plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb','savefig.facecolor':'#fcfcfb','axes.edgecolor':'#8a8980','xtick.color':'#55544e','ytick.color':'#55544e','axes.titlesize':9,'axes.titleweight':'bold'})
SEQ=LinearSegmentedColormap.from_list('seq',['#f4f3ef','#cde2fb','#6da7ec','#256abf','#0d366b'])
DIV=LinearSegmentedColormap.from_list('div',['#184f95','#6da7ec','#f0efec','#ee8a80','#b02a2a'])
CAT={'Pezizomycotina':'#2a78d6','Saccharomycotina':'#eb6834','Agaricomycotina':'#1baf7a','other Basidiomycota':'#eda100','early-diverging':'#6250d6','Microsporidia':'#e34948','Taphrinomycotina':'#e87ba4'}
G=pd.read_parquet(f'{W}/genera_q_grp.parquet'); _,P,Md=load('q')
order=['Sordariomycetes','Leotiomycetes','Dothideomycetes','Eurotiomycetes','Lecanoromycetes','Pezizomycetes','Orbiliomycetes','Pezizomycotina_other','Saccharomycotina','Taphrinomycotina','Agaricomycetes','Tremellomycetes','Agaricomycotina_other','Ustilaginomycotina','Pucciniomycotina','Mucoromycota_sl','Zoopagomycota_sl','zoosporic','Microsporidia']
ng=G.grp.value_counts()
fidx={f:i for i,f in enumerate(F.fid)}
def heat(ax,M,rows,cols,cmap,norm=None,vmin=None,vmax=None,cbl=''):
    im=ax.imshow(M,aspect='auto',cmap=cmap,norm=norm,vmin=vmin,vmax=vmax,interpolation='nearest')
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols,rotation=55,ha='right'); ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows)
    ax.set_xticks(np.arange(-.5,len(cols)),minor=True); ax.set_yticks(np.arange(-.5,len(rows)),minor=True); ax.grid(which='minor',color='#fcfcfb',lw=0.8); ax.tick_params(which='both',length=0)
    for s in ax.spines.values(): s.set_visible(False)
    cb=plt.colorbar(im,ax=ax,fraction=0.025,pad=0.01); cb.set_label(cbl); cb.outline.set_visible(False); return im
# ---------- Fig 1 lineage-restricted
best=pd.read_csv(f'{W}/restricted_best.tsv',sep='\t')
want=[('DIKARYA','Dikarya',5),('PHYLUM','Ascomycota',5),('SUBGROUP','Pezizomycotina',7),('CLASS','Sordariomycetes',2),('ORDER_','Pleosporales',2),('SUBGROUP','Saccharomycotina',6),('PHYLUM','Basidiomycota',6),('CLASS','Agaricomycetes',4),('CLASS','Ustilaginomycetes',3),('PHYLUM','Mucoromycota',5),('PHYLUM','Mortierellomycota',5),('PHYLUM','Glomeromycota',5),('SUBGROUP','zoosporic',6),('PHYLUM','Kickxellomycota',3),('PHYLUM','Microsporidia',6)]
rows=[];labs=[];sep=[]
B=(P>=0.5)
cols=order
# finer columns: add phyla for early-diverging
G['col']=G.grp
for ph in ['Mucoromycota','Mortierellomycota','Glomeromycota','Kickxellomycota']: G.loc[G.PHYLUM==ph,'col']=ph
cols2=[c for c in order if c not in('Mucoromycota_sl','Zoopagomycota_sl')]; cols2=cols2[:-2]+['Mucoromycota','Mortierellomycota','Glomeromycota','Kickxellomycota','Zoopagomycota_sl','zoosporic','Microsporidia']
for rk,cl,k in want:
    d=best[(best['rank']==rk)&(best.clade==cl)].sort_values(['out_n_all','in_frac_q'],ascending=[True,False]).head(k)
    for _,r in d.iterrows():
        i=fidx['pfam|'+r.pfam]; rows.append([B[(G.col==c).values,i].mean() for c in cols2]); labs.append(f"{r['name']} ({r.pfam})")
    sep.append(len(rows))
fig,ax=plt.subplots(figsize=(9.5,13)); heat(ax,np.array(rows),labs,[f'{c} ({(G.col==c).sum()})' for c in cols2],SEQ,vmin=0,vmax=1,cbl='fraction of genera with domain')
for s in sep[:-1]: ax.axhline(s-.5,color='#55544e',lw=0.8)
ax.set_title('Fig 1. Lineage-restricted Pfam domains: fraction of genera carrying each domain, \nby clade (quality-filtered genera; n in brackets)',loc='left')
fig.tight_layout(); fig.savefig(f'{OUT}/fig1_restricted_heatmap.png',dpi=170); plt.close()
# ---------- Fig 2 loss modules x lineages
MP=pd.read_parquet(f'{W}/module_lineage_presence.parquet'); info=pd.read_csv(f'{W}/lineage_loss_info.tsv',sep='\t',index_col=0)
mods={555:'Complex I (NADH dehydrogenase) [17]',787:'Mito. respiration assembly: SDH, COX11, SURF1 [15]',779:'Mito. ATP synthase / heme synthesis [13]',965:'Beta-oxidation + glyoxylate cycle [17]',742:'Type I fatty-acid synthase (FAS) [8]',980:'His/Arg/BCAA/Met biosynthesis [22]',891:'Sulfate assimilation / siroheme [9]',521:'Molybdenum cofactor biosynthesis [7]',429:'Mo-enzymes (sulfite/nitrate reductase type) [5]',430:'Urease + accessory proteins [7]',512:'Glycogen metabolism [9]',993:'Alpha-mannosidase GH38 [5]',728:'Beta-1,3-glucan synthase FKS/Kre9 [8]',366:'RNAi: Dicer / Argonaute [6]',847:'Ku70/Ku80 NHEJ [8]',994:'Autophagy ATG9/ATG12/ATG2 + misc [16]',428:'PAS / hemerythrin / CSN5 [7]',392:'Mechanosensitive channel / Pex16 [5]',438:'DOCK-ELMO GEF / animal haem peroxidase [5]',336:'Splicing factors SF3b1 / cactin [7]',722:'DASH kinetochore complex [8]',914:'SUR7 / KAR9 / ATG11 membrane-polarity set [19]',902:'Pom152 nucleoporin / DAK [7]',966:'Mixed oxidative metabolism (ETF, Cu-oxidase, FBPase) [35]'}
lin=['Sordariomycetes','Dothideomycetes','Eurotiomycetes','Lecanoromycetes','Agaricomycetes','Tremellomycetes','Ustilaginomycetes','Exobasidiomycetes','Malasseziomycetes','Wallemiomycetes','Microbotryomycetes','Pucciniomycetes','Lipomycetes','Dipodascomycetes','Pichiomycetes','Saccharomycetes','g_Hanseniaspora','Taphrinomycetes','Schizosaccharomycetes','Pneumocystomycetes','Mucoromycota','Mortierellomycetes','Glomeromycetes','Kickxellomycota','Zoopagomycetes','Entomophthoromycetes','Chytridiomycetes','Blastocladiomycota','Neocallimastigomycetes','Microsporidia_incertae']
M=MP.loc[list(mods.keys()),lin].values
fig,ax=plt.subplots(figsize=(11,6.8)); heat(ax,M,list(mods.values()),[f"{l.replace('_incertae','').replace('g_','')} ({info.n_genera[l]})" for l in lin],SEQ,vmin=0,vmax=1,cbl='mean fraction of genera retaining module domains')
ax.set_title('Fig 2. Co-absence loss modules (rows; [n domains])\nby lineage (columns; n genera). Pale = module lost.',loc='left')
fig.tight_layout(); fig.savefig(f'{OUT}/fig2_loss_modules.png',dpi=170); plt.close()
# ---------- Fig 3 co-absence correlation matrix + loss vs proteome
C=np.load(f'{W}/coabs_C.npy'); lv=np.load(f'{W}/coabs_leaves.npy')
fig,axs=plt.subplots(1,2,figsize=(12,5.6),gridspec_kw={'width_ratios':[1.05,1]})
im=axs[0].imshow(C[np.ix_(lv,lv)],cmap=DIV,norm=TwoSlopeNorm(0,-0.5,1),interpolation='nearest',rasterized=True); axs[0].set_xticks([]); axs[0].set_yticks([])
for s in axs[0].spines.values(): s.set_visible(False)
cb=plt.colorbar(im,ax=axs[0],fraction=0.04,pad=0.02); cb.set_label('class-weighted phi correlation of absence'); cb.outline.set_visible(False)
axs[0].set_title(f'a. Co-absence among {len(C)} variably-lost core domains\n(average-linkage order; Microsporidia excluded)',loc='left')
ax=axs[1]; x=info.n_prot; y=info.frac_lost*100
ax.scatter(x,y,s=np.clip(np.sqrt(info.n_genera)*14,18,200),color='#2a78d6',edgecolor='#fcfcfb',lw=1,zorder=3)
for l in info.index:
    if info.frac_lost[l]>0.09 or l in('Glomeromycetes','Agaricomycetes','Entomophthoromycetes','Saccharomycetes','Kickxellomycota'): ax.annotate(l.replace('_incertae','').replace('g_',''),(x[l],y[l]),xytext=(4,3),textcoords='offset points',fontsize=7,color='#33332f')
ax.set_xscale('log'); ax.set_yscale('symlog',linthresh=1); ax.set_xlabel('median proteome size (proteins, log)'); ax.set_ylabel('% of 3,568 core domains absent from lineage (symlog)'); ax.grid(color='#e6e5df',lw=0.6,zorder=0)
ax.set_title('b. Core-domain loss versus proteome size, per lineage\n(point area ~ number of genera)',loc='left')
fig.suptitle('Fig 3. Structure of repeated loss',x=0.01,ha='left',fontweight='bold'); fig.tight_layout(); fig.savefig(f'{OUT}/fig3_coabsence_and_loss_vs_size.png',dpi=170); plt.close()
# ---------- Fig 4 allometry scatters
R=np.load(f'{W}/resid_q.npy'); RF=pd.read_parquet(f'{W}/resid_features.parquet'); ridx={f:i for i,f in enumerate(RF.fid)}; ES=pd.read_csv(f'{W}/expansion_stats.tsv',sep='\t'); ES['fid']=ES.ftype+'|'+ES.feat; ES=ES.set_index('fid')
def sg(r):
    if r.PHYLUM=='Microsporidia': return 'Microsporidia'
    if r.SUBGROUP in('Pezizomycotina','Saccharomycotina','Agaricomycotina','Taphrinomycotina'): return r.SUBGROUP
    if r.PHYLUM=='Basidiomycota': return 'other Basidiomycota'
    return 'early-diverging'
G['sg']=G.apply(sg,axis=1); x=np.log10(G.n_proteins.astype(float).values)
sel=[('pfam|PF06985','HET'),('pfam|PF05729','NACHT'),('pfam|PF00172','Zn2Cys6 (Zn_clus)'),('pfam|PF00067','Cytochrome P450'),('pfam|PF00109','Ketoacyl-synt (PKS/FAS)'),('pfam|PF07690','MFS_1 transporter'),('pfam|PF00069','Protein kinase'),('pfam|PF03443','AA9 LPMO'),('pfam|PF20684','Fung_rhodopsin 7TM')]
fig,axs=plt.subplots(3,3,figsize=(11,10),sharex=True)
for ax,(fid,lab) in zip(axs.flat,sel):
    y=np.log10(Md[:,fidx[fid]]+1)
    for k,c in CAT.items():
        m=(G.sg==k).values; ax.scatter(x[m],y[m],s=9,color=c,alpha=0.75,lw=0,label=f'{k} ({m.sum()})',rasterized=True)
    r=ES.loc[fid]; resid=R[:,ridx[fid]]; a=np.mean(y-resid-r.slope*x); xx=np.array([x.min(),x.max()]); ax.plot(xx,a+r.slope*xx,color='#33332f',lw=1.5)
    ax.set_title(f'{lab}\nslope {r.slope:.2f}; eta2(class) {r.eta2_class:.2f}; pR2(lifestyle|order) {r.pR2_life_order:.2f}',loc='left',fontsize=8); ax.grid(color='#e6e5df',lw=0.5)
for ax in axs[-1]: ax.set_xlabel('log10 proteome size (genus median)')
for ax in axs[:,0]: ax.set_ylabel('log10 (proteins with domain + 1)')
h,l=axs[0,0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=7,frameon=False,markerscale=2.2,fontsize=8)
fig.suptitle('Fig 4. Allometry of selected families: genus-median copy number versus proteome size (class-balanced WLS fit, n = 1,173 genera)',x=0.01,ha='left',fontweight='bold')
fig.tight_layout(rect=(0,0.03,1,0.97)); fig.savefig(f'{OUT}/fig4_allometry_scatter.png',dpi=170); plt.close()
# ---------- Fig 5 class mean residual heatmap
CM=pd.read_parquet(f'{W}/class_mean_resid.parquet').set_index('fid')
feats=[('pfam|PF06985','HET'),('pfam|PF05729','NACHT'),('pfam|PF24883','NPHP3_N NLR NTPase'),('pfam|PF00931','NB-ARC'),('pfam|PF12796','Ank_2'),('pfam|PF20684','Fung_rhodopsin 7TM'),('pfam|PF05730','CFEM'),('pfam|PF00172','Zn2Cys6'),('pfam|PF04082','Fungal_trans'),('tf|C2H2 zinc finger','C2H2 TF'),('tf|Homeodomain','Homeodomain TF'),('tf|High mobility group box','HMG box TF'),('pfam|PF00069','Protein kinase'),('trait|two_comp','Two-component (HisKA/RR)'),('trait|g_alpha','G-alpha'),('pfam|PF00067','P450'),('pfam|PF00109','Ketoacyl-synt (PKS)'),('pfam|PF00668','Condensation (NRPS)'),('pfam|PF11807','UstYa (RiPP)'),('pfam|PF07690','MFS_1'),('pfam|PF00083','Sugar_tr'),('pfam|PF00005','ABC_tran'),('cazy_sub|cellulose','CAZy: cellulose'),('cazy_sub|xylan','CAZy: xylan'),('cazy_sub|pectin','CAZy: pectin'),('cazy_sub|lignin','CAZy: lignin'),('cazy_sub|chitin','CAZy: chitin'),('cazy_fam|AA9','AA9 LPMO'),('pfam|PF00734','CBM1'),('cazy_class|GT','CAZy: GT'),('merops|S53','MEROPS S53 sedolisin'),('merops|S08A','MEROPS S08A subtilisin'),('merops|A01A','MEROPS A01A pepsin'),('merops|S09X','MEROPS S09X'),('pfam|PF00264','Tyrosinase'),('trait|te','TE-derived domains'),('pfam|PF07534','TLDc'),('pfam|PF08238','Sel1 repeat')]
cc=['Sordariomycetes','Leotiomycetes','Dothideomycetes','Eurotiomycetes','Lecanoromycetes','Arthoniomycetes','Lichinomycetes','Geoglossomycetes','Pezizomycetes','Orbiliomycetes','Saccharomycetes','Pichiomycetes','Dipodascomycetes','Lipomycetes','Agaricomycetes','Tremellomycetes','Ustilaginomycetes','Exobasidiomycetes','Microbotryomycetes','Mucoromycetes','Mortierellomycetes','Glomeromycetes','Chytridiomycetes','Microsporidia_incertae']
cc=[c for c in cc if c in CM.columns]; feats=[f for f in feats if f[0] in CM.index]
M=CM.loc[[f[0] for f in feats],cc].values
fig,ax=plt.subplots(figsize=(10.5,9.5)); heat(ax,M,[f[1] for f in feats],[f"{c.replace('_incertae','')} ({(G.CLASS==c).sum()})" for c in cc],DIV,norm=TwoSlopeNorm(0,-1.2,1.2),cbl='mean allometry residual (log10 units; +0.3 = 2x expected)')
ax.set_title('Fig 5. Class-mean residuals from the proteome-size allometry \n(red = expanded relative to proteome size, blue = contracted)',loc='left')
fig.tight_layout(); fig.savefig(f'{OUT}/fig5_class_residual_heatmap.png',dpi=170); plt.close()
# ---------- Fig 6 taxonomy vs lifestyle
tc={'pfam':'#b9b8ae','cazy_fam':'#2a78d6','cazy_sub':'#eb6834','cazy_class':'#eb6834','merops':'#1baf7a','tf':'#6250d6','trait':'#eda100'}
fig,ax=plt.subplots(figsize=(8.5,6.2))
for t,c in tc.items():
    if t=='cazy_class': continue
    d=ES[ES.ftype.isin(['cazy_sub','cazy_class'] if t=='cazy_sub' else [t])]
    ax.scatter(d.eta2_class,d.pR2_life_order,s=10 if t=='pfam' else 26,color=c,alpha=0.6 if t=='pfam' else 0.95,lw=0.6 if t!='pfam' else 0,edgecolor='#fcfcfb',label={'pfam':'Pfam domain','cazy_fam':'CAZy family','cazy_sub':'CAZy substrate/class','merops':'MEROPS family','tf':'TF family','trait':'signalling/NLR set'}[t]+f' ({len(d)})',zorder=2 if t=='pfam' else 3,rasterized=True)
lab={'cazy_fam|GH10':'GH10','cazy_sub|pectin':'pectin','cazy_sub|cellulose':'cellulose','cazy_fam|AA9':'AA9','pfam|PF00248':'Aldo/keto red.','pfam|PF20684':'Fung_rhodopsin','pfam|PF26118':'DUF8035','tf|C6H2-type zinc finger':'Fungal_trans_2','pfam|PF06985':'HET','pfam|PF05729':'NACHT','pfam|PF00172':'Zn2Cys6','pfam|PF00069':'Pkinase','merops|S33':'S33','pfam|PF13409':'GST_N_2','pfam|PF17921':'Integrase_H2C2','cazy_fam|GT62':'GT62','pfam|PF00067':'P450','pfam|PF00109':'PKS KS'}
for f,l in lab.items():
    if f in ES.index: ax.annotate(l,(ES.eta2_class[f],ES.pR2_life_order[f]),xytext=(4,3),textcoords='offset points',fontsize=7,color='#33332f')
ax.set_xlabel('variance of allometry residual explained by CLASS (eta squared)'); ax.set_ylabel('partial R2 of lifestyle within ORDER'); ax.grid(color='#e6e5df',lw=0.5); ax.legend(frameon=False,loc='upper left',fontsize=7.5)
ax.set_title('Fig 6. Taxonomy versus lifestyle: 2,831 families, genus level\n(1,173 genera; 1,054 with a lifestyle label)',loc='left')
fig.tight_layout(); fig.savefig(f'{OUT}/fig6_taxonomy_vs_lifestyle.png',dpi=170); plt.close()
# ---------- Fig 7 lifestyle effects heatmap
LE=pd.read_parquet(f'{W}/life_eff_order.parquet').set_index('fid')
top=ES[(ES.ftype!='cazy_class')].sort_values('pR2_life_order',ascending=False)
pick=list(top[top.ftype=='pfam'].index[:14])+list(top[top.ftype=='cazy_fam'].index[:8])+list(top[top.ftype=='cazy_sub'].index[:5])+list(top[top.ftype=='merops'].index[:4])+['pfam|PF17921','pfam|PF00069','pfam|PF00172','pfam|PF05729','pfam|PF00067','pfam|PF00082','cazy_sub|chitin','cazy_fam|GH55']
pick=list(dict.fromkeys(pick)); lc=['ectomycorrhizal','lichenized','animal_endosymbiont','animal_parasite','mycoparasite','plant_pathogen','wood_saprotroph','litter_saprotroph','soil_saprotroph','dung_saprotroph','nectar/tap_saprotroph','unspecified_saprotroph']
nl=G.primary_lifestyle.value_counts()
fig,ax=plt.subplots(figsize=(8.5,9.5)); heat(ax,LE.loc[pick,lc].values,[f"{ES.name[p] if isinstance(ES.name[p],str) else ES.feat[p]} [{ES.ftype[p]}]" for p in pick],[f'{l} ({nl.get(l,0)})' for l in lc],DIV,norm=TwoSlopeNorm(0,-0.25,0.25),cbl='mean within-order residual (log10 units)')
ax.set_title('Fig 7. Lifestyle effects within orders (order means removed; allometry residuals)',loc='left')
fig.tight_layout(); fig.savefig(f'{OUT}/fig7_lifestyle_effects.png',dpi=170); plt.close()
# ---------- Fig 8 functional-class shift by group (log2 ratio to kingdom median of group medians)
ST=pd.read_csv(f'{W}/selected_by_group.tsv',sep='\t',index_col=0); SC=pd.read_csv(f'{W}/scalar_by_group.tsv',sep='\t',index_col=0)
rows=['cazy_class|GH','cazy_class|AA','cazy_class|CE','cazy_class|PL','cazy_class|GT','cazy_sub|cellulose','cazy_sub|xylan','cazy_sub|pectin','cazy_sub|lignin','cazy_sub|chitin','merops|A01A','merops|S08A','merops|S53','merops|S09X','merops|S10','tf|Zn(II)2Cys6 zinc cluster','tf|C2H2 zinc finger','tf|Homeodomain','tf|Basic helix-loop-helix','tf|GATA zinc finger','tf|High mobility group box']
A=pd.concat([ST.loc[rows,order], SC.loc[order,['sec_frac','gpi_frac','tm_frac','idp','frac_pfam']].T.rename(index={'sec_frac':'secretome fraction','gpi_frac':'GPI-anchored fraction','tm_frac':'TM (>=2 helices) fraction','idp':'median disorder fraction','frac_pfam':'fraction of proteins with Pfam'})])
ref=A.median(axis=1); Lr=np.log2((A.add(0.05*ref,axis=0)).div(ref*1.05,axis=0))
fig,ax=plt.subplots(figsize=(10,7.6)); heat(ax,Lr.values,[r.replace('cazy_class|','CAZy ').replace('cazy_sub|','CAZy substrate: ').replace('merops|','MEROPS ').replace('tf|','TF: ') for r in A.index],[f'{c} ({ng[c]})' for c in order],DIV,norm=TwoSlopeNorm(0,-2.5,2.5),cbl='log2 (group median / median of groups)')
ax.axhline(len(rows)-.5,color='#55544e',lw=0.8)
ax.set_title('Fig 8. Functional-class shifts: copies per 1,000 proteins (upper block) \nand proteome fractions (lower block), relative to the across-group median',loc='left')
fig.tight_layout(); fig.savefig(f'{OUT}/fig8_functional_class_shifts.png',dpi=170); plt.close()
# ---------- Fig 9 Dikarya-lost ancestral complexes
cur=[('PF10189','Integrator INTS3'),('PF24436','Integrator INTS7'),('PF14750','Integrator INTS2'),('PF25458','Integrator INTS4'),('PF14631','Fanconi FANCD2'),('PF14675','Fanconi FANCI'),('PF18891','Fanconi FANCL'),('PF09404','KICSTOR C12orf66'),('PF30693','KICSTOR kaptin'),('PF15907','KICSTOR ITFG2'),('PF16672','Ragulator LAMTOR5'),('PF02965','Met synthase B12 activation'),('PF02310','B12-binding'),('PF01923','Cobalamin adenosyltransferase'),('PF10229','MMADHC (B12 trafficking)'),('PF18140','Propionyl-CoA carboxylase BT'),('PF00351','Biopterin aromatic aa hydroxylase'),('PF01242','PTPS (biopterin synthesis)'),('PF15785','SMG1 kinase (NMD)'),('PF10220','Smg8/Smg9'),('PF08163','DNA-PKcs'),('PF06209','NELF-B (COBRA1)'),('PF01080','Presenilin'),('PF03493','BK K+ channel'),('PF01167','Tubby'),('PF12580','Tripeptidyl peptidase II'),('PF23659','UFL1 (UFMylation)'),('PF03106','WRKY'),('PF12437','Glutamine synthetase III')]
cur=[c for c in cur if 'pfam|'+c[0] in fidx]
M=np.array([[B[(G.col==c).values,fidx['pfam|'+p]].mean() for c in cols2] for p,_ in cur])
fig,ax=plt.subplots(figsize=(9.5,7.6)); heat(ax,M,[f'{n} ({p})' for p,n in cur],[f'{c} ({(G.col==c).sum()})' for c in cols2],SEQ,vmin=0,vmax=1,cbl='fraction of genera with domain')
ax.set_title('Fig 9. Ancestral (opisthokont-type) domains retained in early-diverging fungi \nand absent from Dikarya (curated from 203)',loc='left')
fig.tight_layout(); fig.savefig(f'{OUT}/fig9_dikarya_lost.png',dpi=170); plt.close()
print('done')
