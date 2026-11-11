import duckdb, os
W=os.path.expanduser('~/bigdata/tmp_kx_scratch/')
R='/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/'
con=duckdb.connect(); con.sql("set threads=4")
con.sql(f"copy (select m.* from '{R}A7/pfam_genome_matrix.parquet' m where genome in (select LOCUSTAG from '{W}scanset.parquet')) to '{W}pf_scan.parquet'")
for f,o in [('A5/cazy_class_counts','cazy_scan'),('A5/merops_family_counts','merops_scan'),('A5/tf_counts','tf_scan'),('A11/pfam_trait_counts','trait_scan')]:
    con.sql(f"copy (select m.* from '{R}{f}.parquet' m where genome in (select LOCUSTAG from '{W}scanset.parquet')) to '{W}{o}.parquet'")
print(con.sql(f"select count(*), count(distinct genome), count(distinct pfam_acc) from '{W}pf_scan.parquet'").fetchall())
