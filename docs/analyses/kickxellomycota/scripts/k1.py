import duckdb
R='/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/'
con=duckdb.connect()
con.sql(f"""create view g as select l.*, e.hs_group, e.septation, e.obligate_mode, s.* exclude(genome), x.n_signalp, x.n_secreted_strict, x.n_tm2, c.n_cazy, m.n_merops, t.n_tf_proteins, i.median_idp_fraction
from '{R}labels/genome_labels.parquet' l left join '{R}labels/genome_labels_ext.parquet' e using(LOCUSTAG)
left join '{R}A2/genome_gene_stats.parquet' s on s.genome=l.LOCUSTAG
left join '{R}A5/secretome.parquet' x on x.genome=l.LOCUSTAG
left join '{R}A5/cazy_genome.parquet' c on c.genome=l.LOCUSTAG
left join '{R}A5/merops_genome.parquet' m on m.genome=l.LOCUSTAG
left join '{R}A5/tf_genome.parquet' t on t.genome=l.LOCUSTAG
left join '{R}A5/idp_genome.parquet' i on i.genome=l.LOCUSTAG""")
con.sql("copy (select * from g) to '/tmp/kx_g_all.parquet'")
import pandas as pd
pd.set_option('display.width',300); pd.set_option('display.max_rows',500); pd.set_option('display.max_columns',50)
print(con.sql("""select CLASS,ORDER_,FAMILY,GENUS,count(*) n,count(distinct SPECIES) nsp,sum(is_representative::int) reps, sum((is_representative and q_busco_rel and q_dup and q_n50)::int) reps_q,
 any_value(growth_form) gf, any_value(primary_lifestyle) pl, any_value(hs_group) hs,
 round(median(total_length_bp)/1e6,1) mb, median(n_protein_coding) genes, round(median(busco_complete),1) busco, round(median(busco_duplicated),1) dup, round(median(N50_bp)/1e3) n50kb,
 round(median(introns_per_transcript),2) ipt, median(median_intron_len) ilen, median(median_intergenic) ig, round(median(gc_pct),1) gc, round(median(masked_pct),1) masked
 from g where PHYLUM='Kickxellomycota' group by all order by 1,2,3,4""").df().to_string())
print(con.sql("select busco_thr_level, busco_thr_p25, count(*) from g where PHYLUM='Kickxellomycota' group by all").df())
print(con.sql("""select PHYLUM,CLASS,ORDER_,count(*) n,count(distinct GENUS) ng,sum(is_representative::int) reps, sum((is_representative and q_busco_rel and q_dup and q_n50)::int) reps_q,
 round(median(total_length_bp)/1e6,1) mb, median(n_protein_coding) genes, round(median(busco_complete),1) busco, round(median(busco_duplicated),1) dup, round(median(N50_bp)/1e3) n50kb,
 round(median(introns_per_transcript),2) ipt, median(median_intron_len) ilen, median(median_intergenic) ig
 from g where PHYLUM not in ('Ascomycota','Basidiomycota') group by all order by 1,2,3""").df().to_string())
