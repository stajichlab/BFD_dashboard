#!/usr/bin/env python3
"""Markdown table fragments for the report (written to kx_data/frag_*.md)."""
import pandas as pd
D='kx_data/'
def md(df):
    rows=['| '+' | '.join(map(str,df.columns))+' |','|'+'---|'*len(df.columns)]
    rows+=['| '+' | '.join('' if (isinstance(v,float) and v!=v) else (str(int(v)) if isinstance(v,float) and v==int(v) and abs(v)>=100 else str(v)) for v in r)+' |' for r in df.itertuples(index=False)]
    return '\n'.join(rows)+'\n'
g=pd.read_csv(D+'grp_genus_medians_all.tsv',sep='\t')
keep=['Kx:Coemansia','Kx:Kickxellales-core-other','Kx:Kickxellales-early','Kx:Ramicandelaber','Kx:Dimargaritales','Kx:Harpellales','Zoopagales','Entomophthorales','Basidiobolus','Mortierellales','Mucorales','Umbelopsidales','Glomeromycota','Chytrid+Monobleph','Blastocladiomycota','Rozello+Aphelid','Microsporidia','Sacch:Saccharomycetes','Sacch:Pichiomycetes','Sacch:Dipodascomycetes','Sacch:Lipomycetes','Schizosaccharomyces','Taphrinomycetes','Pneumocystis','Malassezia','Ustilaginomycetes','Wallemia','Pez:Sordariomycetes','Agaricomycetes']
g=g.set_index('grp').loc[keep].reset_index()
a=g[['grp','n_gen','n_sp','mb','n_protein_coding','gene_density','busco_complete','N50_bp','partial_fraction','introns_per_transcript','frac_intronless','median_intron_len','median_intergenic','median_protein_len','masked_pct','gc_pct']].copy()
a.columns=['group','genera','species','Mb','genes','genes/Mb','BUSCO %','N50 bp','partial frac','introns/tx','intronless frac','intron bp','intergenic bp','protein aa','masked %','GC %']
for c in ['Mb','BUSCO %','masked %','GC %']: a[c]=a[c].round(1)
for c in ['genes','genes/Mb','N50 bp','intron bp','intergenic bp','protein aa']: a[c]=a[c].round(0).astype(int)
for c in ['partial frac','introns/tx','intronless frac']: a[c]=a[c].round(2)
open(D+'frag_group_structure.md','w').write(md(a))
b=g[['grp','n_gen','n_protein_coding','n_secreted_strict','sec_k','n_cazy','cazy_k','n_merops','mer_k','n_tf_proteins','tf_k','tm_k','median_idp_fraction']].copy()
b.columns=['group','genera','genes','secreted','secreted/1000','CAZymes','CAZy/1000','MEROPS','MEROPS/1000','TFs','TF/1000','TM>=2/1000','median IDP frac']
for c in b.columns[2:]: b[c]=b[c].round(1 if '/' in c else (2 if 'IDP' in c else 0))
open(D+'frag_group_function.md','w').write(md(b))
e=pd.read_csv(D+'stats_genus_effect_sizes.tsv',sep='\t')
e=e[e.focal.str.startswith('core')&e.reference.isin(['Mucorales','Mortierellales','Zoopagales','Sacch:Dipodascomycetes','Sacch:Saccharomycetes'])&e.metric.isin(['mb','n_protein_coding','introns_per_transcript','median_intron_len','median_intergenic','median_protein_len','masked_pct'])].copy()
e['95% CI']=e.boot_lo.astype(str)+' to '+e.boot_hi.astype(str); e['mannwhitney_p']=e.mannwhitney_p.map(lambda p:f'{p:.2g}')
open(D+'frag_effect_sizes.md','w').write(md(e[['reference','metric','n_focal_genera','n_ref_genera','median_focal','median_ref','median_diff','95% CI','log2_ratio','cliffs_delta','mannwhitney_p']]))
f=pd.read_csv(D+'stats_loss_overlap_fisher.tsv',sep='\t'); f['fisher_p']=f.fisher_p.map(lambda p:f'{p:.2g}')
open(D+'frag_fisher.md','w').write(md(f))
o=pd.read_csv(D+'md_ortholog_matched_introns.tsv',sep='\t'); open(D+'frag_ortho.md','w').write(md(o))
k=pd.read_csv(D+'kx_species_reps.tsv',sep='\t')
kk=k.groupby(['ORDER_','GENUS']).agg(species=('LOCUSTAG','size'),hq=('N50_bp',lambda x:(x>=50000).sum()),lifestyle=('primary_lifestyle','first'),Mb=('mb','median'),genes=('n_protein_coding','median'),busco=('busco_complete','median'),busco_max=('busco_complete','max'),n50=('N50_bp','median'),ipt=('introns_per_transcript','median'),ilen=('median_intron_len','median'),ig=('median_intergenic','median'),masked=('masked_pct','median'),gc=('gc_pct','median')).reset_index().round(2)
kk.columns=['order','genus','species reps','N50>=50kb','FungalTraits lifestyle','Mb','genes','BUSCO med','BUSCO max','N50 med','introns/tx','intron bp','intergenic bp','masked %','GC %']
open(D+'frag_kx_genera.md','w').write(md(kk))
