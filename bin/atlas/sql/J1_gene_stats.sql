-- J1: per-genome gene structure statistics (A2). genome = LOCUSTAG = prefix of gene_id / transcript_id.
-- @output A2/gene_info_stats.parquet
COPY (
  SELECT LOCUSTAG AS genome,
         count(*) AS n_genes,
         count(*) FILTER (WHERE gene_type = 'protein_coding') AS n_protein_coding,
         count(*) FILTER (WHERE gene_type = 'tRNA_gene') AS n_trna_genes,
         count(*) FILTER (WHERE gene_type NOT IN ('protein_coding','tRNA_gene')) AS n_other_genes,
         median("end" - start + 1) FILTER (WHERE gene_type = 'protein_coding') AS median_gene_span,
         quantile_cont("end" - start + 1, 0.9) FILTER (WHERE gene_type = 'protein_coding') AS q90_gene_span,
         count(DISTINCT chrom) AS n_contigs_with_genes
  FROM '{T}/gene_info.parquet' GROUP BY 1
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A2/transcript_stats.parquet
COPY (
  SELECT split_part(gene_id,'_',1) AS genome,
         count(*) AS n_transcripts,
         count(DISTINCT gene_id) AS n_genes_with_transcript,
         count(*) FILTER (WHERE is_partial) AS n_partial,
         count(*) FILTER (WHERE is_partial) * 1.0 / count(*) AS partial_fraction,
         count(*) FILTER (WHERE has_start_codon = 'FALSE') * 1.0 / count(*) AS frac_no_start,
         count(*) FILTER (WHERE has_stop_codon = 'FALSE') * 1.0 / count(*) AS frac_no_stop,
         count(*) FILTER (WHERE has_start_codon = 'NULL' OR has_start_codon IS NULL) * 1.0 / count(*) AS frac_start_null,
         count(*) FILTER (WHERE NOT is_partial AND has_start_codon = 'TRUE' AND has_stop_codon = 'TRUE') AS n_complete_models,
         median(CDS_length) AS median_cds_len,
         median(CDS_length) FILTER (WHERE NOT is_partial AND has_start_codon = 'TRUE' AND has_stop_codon = 'TRUE') AS median_cds_len_complete,
         median("end" - start + 1) AS median_transcript_span
  FROM '{T}/gene_transcripts.parquet' GROUP BY 1
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A2/protein_stats.parquet
COPY (
  SELECT split_part(gene_id,'_',1) AS genome,
         count(*) AS n_proteins,
         count(DISTINCT gene_id) AS n_genes_with_protein,
         avg(length) AS mean_protein_len,
         quantile_cont(length, 0.10) AS q10_protein_len,
         quantile_cont(length, 0.25) AS q25_protein_len,
         median(length) AS median_protein_len,
         quantile_cont(length, 0.75) AS q75_protein_len,
         quantile_cont(length, 0.90) AS q90_protein_len,
         count(*) FILTER (WHERE length < 100) * 1.0 / count(*) AS frac_short_lt100,
         count(*) FILTER (WHERE length >= 1000) * 1.0 / count(*) AS frac_long_ge1000
  FROM '{T}/gene_proteins.parquet' GROUP BY 1
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A2/intron_stats.parquet
COPY (
  WITH per_tx AS (
    SELECT transcript_id, count(*) AS n_introns
    FROM '{T}/gene_introns.parquet' GROUP BY 1
  ), tx AS (
    SELECT split_part(transcript_id,'_',1) AS genome, count(*) AS n_tx_with_intron, sum(n_introns) AS n_introns
    FROM per_tx GROUP BY 1
  ), il AS (
    SELECT split_part(transcript_id,'_',1) AS genome,
           median("end" - start + 1) AS median_intron_len,
           quantile_cont("end" - start + 1, 0.9) AS q90_intron_len,
           avg(GC_content) AS mean_intron_gc,
           count(*) FILTER (WHERE upper(splice_5) = 'GT' AND upper(splice_3) = 'AG') * 1.0 / count(*) AS frac_gt_ag,
           count(*) FILTER (WHERE upper(splice_5) = 'GC' AND upper(splice_3) = 'AG') * 1.0 / count(*) AS frac_gc_ag
    FROM '{T}/gene_introns.parquet' GROUP BY 1
  )
  SELECT tx.genome, tx.n_tx_with_intron, tx.n_introns, il.median_intron_len, il.q90_intron_len, il.mean_intron_gc, il.frac_gt_ag, il.frac_gc_ag
  FROM tx JOIN il USING (genome)
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A2/intergenic_stats.parquet
COPY (
  SELECT species_prefix AS genome,
         count(*) AS n_intergenic,
         median(distance) AS median_intergenic,
         quantile_cont(distance, 0.1) AS q10_intergenic,
         quantile_cont(distance, 0.9) AS q90_intergenic,
         avg(distance) AS mean_intergenic,
         count(*) FILTER (WHERE distance < 0) * 1.0 / count(*) AS frac_overlapping
  FROM '{T}/gene_intergenic_distances.parquet' GROUP BY 1
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A2/genome_gene_stats.parquet
COPY (
  SELECT g.*, t.* EXCLUDE (genome), p.* EXCLUDE (genome), i.* EXCLUDE (genome), d.* EXCLUDE (genome),
         i.n_introns * 1.0 / t.n_transcripts AS introns_per_transcript,
         1.0 - i.n_tx_with_intron * 1.0 / t.n_transcripts AS frac_intronless,
         (t.median_transcript_span - t.median_cds_len) / t.median_transcript_span AS median_noncoding_fraction
  FROM '{A}/A2/gene_info_stats.parquet' g
  LEFT JOIN '{A}/A2/transcript_stats.parquet' t USING (genome)
  LEFT JOIN '{A}/A2/protein_stats.parquet' p USING (genome)
  LEFT JOIN '{A}/A2/intron_stats.parquet' i USING (genome)
  LEFT JOIN '{A}/A2/intergenic_stats.parquet' d USING (genome)
  ORDER BY genome
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
