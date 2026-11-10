import duckdb, os, time
W=os.path.expanduser('~/bigdata/tmp_kx_scratch/')
T='/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables'
R='/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/'
con=duckdb.connect()
con.sql("set threads=16; set memory_limit='56GB'; set temp_directory='/scratch/jstajich/kxtmp'; set preserve_insertion_order=false")
t0=time.time()
con.sql(f"create table s as select LOCUSTAG from '{W}scanset.parquet'")
con.sql(f"""create table tx as select left(transcript_id,8) genome, transcript_id, CDS_length, "end"-"start"+1 span,
   (not is_partial and has_start_codon='TRUE' and has_stop_codon='TRUE') q_model
   from '{T}/gene_transcripts.parquet' where left(transcript_id,8) in (select LOCUSTAG from s)""")
print('tx',con.sql("select count(*) from tx").fetchone(), time.time()-t0, flush=True)
con.sql(f"""create table it as select left(transcript_id,8) genome, transcript_id, "end"-"start"+1 ilen, upper(splice_5) s5, upper(splice_3) s3, GC_content gc, codon_position cpos, codon_frame phase
   from '{T}/gene_introns.parquet' where left(transcript_id,8) in (select LOCUSTAG from s)""")
print('it',con.sql("select count(*) from it").fetchone(), time.time()-t0, flush=True)
print(con.sql("select phase,count(*) from it group by 1 order by 1").fetchall())
print(con.sql("select min(cpos),max(cpos),min(ilen),max(ilen) from it").fetchall())
# per-transcript intron counts
con.sql("""create table txi as select tx.*, coalesce(i.n,0) n_int, coalesce(i.sumlen,0) sum_ilen from tx left join
  (select transcript_id, count(*) n, sum(ilen) sumlen from it group by 1) i using(transcript_id)""")
con.sql(f"""copy (select genome, count(*) n_tx, sum(q_model::int) n_tx_complete, sum(CDS_length) sum_cds, sum(span) sum_span, sum(n_int) n_int, sum(sum_ilen) sum_ilen,
   avg(n_int) ipt, avg(n_int) filter(where q_model) ipt_complete,
   avg((n_int=0)::int) f0, avg((n_int=1)::int) f1, avg((n_int=2)::int) f2, avg((n_int>=3)::int) f3p,
   avg((n_int=0)::int) filter(where q_model) f0_complete,
   1000.0*sum(n_int)/sum(CDS_length) int_per_kb_cds,
   1000.0*sum(n_int) filter(where q_model)/sum(CDS_length) filter(where q_model) int_per_kb_cds_complete,
   median(CDS_length) med_cds, median(CDS_length) filter(where q_model) med_cds_complete,
   median(span) filter(where q_model) med_span_complete
   from txi group by 1) to '{W}scan_tx.parquet'""")
con.sql(f"""copy (select i.genome, count(*) n_int,
   avg((s5='GT' and s3='AG')::int) f_gtag, avg((s5='GC' and s3='AG')::int) f_gcag, avg((s5='AT' and s3='AC')::int) f_atac,
   avg((not ((s5='GT' and s3='AG') or (s5='GC' and s3='AG') or (s5='AT' and s3='AC')))::int) f_other,
   avg((phase=0)::int) ph0, avg((phase=1)::int) ph1, avg((phase=2)::int) ph2,
   quantile_cont(ilen,0.1) il_q10, quantile_cont(ilen,0.25) il_q25, median(ilen) il_med, quantile_cont(ilen,0.75) il_q75, quantile_cont(ilen,0.9) il_q90, avg(ilen) il_mean,
   avg((ilen<50)::int) f_lt50, avg((ilen>500)::int) f_gt500, avg(gc) intron_gc,
   avg(cpos::double/nullif(t.CDS_length,0)) relpos_mean, median(cpos::double/nullif(t.CDS_length,0)) relpos_med,
   avg((cpos::double/nullif(t.CDS_length,0) < 0.2)::int) f_first20,
   avg((cpos::double/nullif(t.CDS_length,0) > 0.8)::int) f_last20,
   avg((cpos<=60)::int) f_cpos_le60,
   avg(cpos::double/nullif(t.CDS_length,0)) filter(where t.n_int=1 and t.q_model) relpos_mean_single,
   avg((cpos::double/nullif(t.CDS_length,0) < 0.2)::int) filter(where t.n_int=1 and t.q_model) f_first20_single,
   avg((cpos<=60)::int) filter(where t.n_int=1 and t.q_model) f_cpos_le60_single
   from it i join txi t using(transcript_id) group by 1) to '{W}scan_intron.parquet'""")
# length histogram by group (log bins)
con.sql(f"""copy (select i.genome, least(floor(ilen/10)*10,1000) bin, count(*) n from it i group by 1,2) to '{W}scan_ilen_hist.parquet'""")
# relative position decile hist
con.sql(f"""copy (select i.genome, t.n_int=1 single, least(floor(10*cpos::double/nullif(t.CDS_length,0)),9) dec, count(*) n from it i join txi t using(transcript_id) where t.q_model group by 1,2,3) to '{W}scan_relpos_hist.parquet'""")
print('aggs',time.time()-t0, flush=True)
# ortholog-matched: strict swissprot best hit per protein
con.sql(f"""create table pc as select protein_id, genome, sp_id, org_code, length, CDS_length, bitscore, q_model, n_domains
   from '{R}A3/protein_core.parquet' where genome in (select LOCUSTAG from s) and sp_id is not null and qcov_x10>=800 and scov_x10>=800 and pident_x10>=300""")
print('pc',con.sql("select count(*) from pc").fetchone(), time.time()-t0, flush=True)
con.sql(f"""copy (select p.genome, p.sp_id, count(*) n_copies, arg_max(t.n_int,p.bitscore) n_int, arg_max(p.CDS_length,p.bitscore) cds_len, arg_max(t.sum_ilen,p.bitscore) sum_ilen, arg_max(p.q_model,p.bitscore) q_model, max(p.bitscore) bitscore
   from pc p join txi t on t.transcript_id=p.protein_id group by 1,2) to '{W}scan_sp_genome.parquet'""")
print('done',time.time()-t0, flush=True)
