import duckdb, pandas as pd, numpy as np
pd.set_option('display.width',320); pd.set_option('display.max_rows',500); pd.set_option('display.max_columns',60)
con=duckdb.connect()
T='/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables'
import glob,os
print(sorted(os.listdir(T))[:80])
for t in ['gene_introns','gene_transcripts','gene_intergenic_distances','gene_CDS','pfam']:
    p=[x for x in glob.glob(f'{T}/{t}*')]
    print(t,p[:3])
R='/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/'
print(con.sql(f"describe select * from '{R}A3/protein_core.parquet'").df())
con.sql("create view g as select * from '/tmp/kx_g_all.parquet'")
print(con.sql("select PHYLUM,CLASS,count(*) n,count(distinct GENUS) ng, sum(is_representative::int) reps, any_value(yeast_origin), round(median(total_length_bp)/1e6,1) mb, median(n_protein_coding) genes,round(median(introns_per_transcript),2) ipt from g where PHYLUM in ('Ascomycota','Basidiomycota') group by all order by 1,2").df().to_string())
print(con.sql("select GENUS, count(*) n,  round(median(total_length_bp)/1e6,1) mb, median(n_protein_coding) genes,round(median(introns_per_transcript),2) ipt, round(median(busco_complete),1) b from g where GENUS in ('Pneumocystis','Malassezia','Taphrina','Schizosaccharomyces','Saitoella','Neolecta','Protomyces','Wallemia','Encephalitozoon') group by all").df())
# Coemansia quality relationships
d=con.sql("select * from g where PHYLUM='Kickxellomycota' and is_representative").df()
print(d[['busco_complete','busco_single','busco_duplicated','busco_fragmented','busco_missing','N50_bp','contig_count','total_length_bp','n_protein_coding','partial_fraction','frac_no_start','frac_no_stop','introns_per_transcript','frac_intronless','median_intron_len','median_intergenic','median_protein_len','median_gene_span','masked_pct','gc_pct','frac_gt_ag','frac_gc_ag']].describe().T.to_string())
c=d[d.GENUS=='Coemansia']
cols=['busco_complete','busco_fragmented','busco_missing','N50_bp','contig_count','total_length_bp','n_protein_coding','partial_fraction','introns_per_transcript','median_intron_len','median_intergenic','median_protein_len','masked_pct','gc_pct']
print('Coemansia species reps n=',len(c)); print(c[cols].corr(method='spearman').round(2).to_string())
print(c.sort_values('N50_bp',ascending=False)[['LOCUSTAG','SPECIES','STRAIN','N50_bp','contig_count','total_length_bp','busco_complete','busco_fragmented','busco_missing','n_protein_coding','partial_fraction','introns_per_transcript','frac_intronless','median_intron_len','median_intergenic','masked_pct','gc_pct']].head(25).to_string())
print(c.SPECIES.str.contains(' sp').sum(), c.SPECIES.head(20).tolist())
print(d[d.GENUS!='Coemansia'].sort_values(['ORDER_','GENUS'])[['LOCUSTAG','SPECIES','N50_bp','contig_count','total_length_bp','busco_complete','busco_duplicated','busco_fragmented','busco_missing','n_protein_coding','partial_fraction','introns_per_transcript','frac_intronless','median_intron_len','median_intergenic','masked_pct','gc_pct']].to_string())
