-- J3: protein-level core table (A3) and family-by-genome length tables for H2.
-- Needs A3/protein_pfam_arch.parquet (J2) and labels/genome_labels.parquet (build_labels.py).
-- @output A3/protein_core.parquet
COPY (
  SELECT p.protein_id,
         split_part(p.gene_id,'_',1) AS genome,
         p.gene_id,
         p.length,
         t.is_partial,
         (t.has_start_codon = 'TRUE') AS has_start,
         (t.has_stop_codon = 'TRUE') AS has_stop,
         (NOT t.is_partial AND t.has_start_codon = 'TRUE' AND t.has_stop_codon = 'TRUE') AS q_model,
         t.CDS_length,
         coalesce(a.n_domains, 0) AS n_domains,
         a.pfam_arch,
         s.sp_id, s.org_code, s.pident_x10, s.qcov_x10, s.scov_x10, s.q_len, s.s_len, s.bitscore
  FROM '{T}/gene_proteins.parquet' p
  LEFT JOIN '{T}/gene_transcripts.parquet' t ON p.transcript_id = t.transcript_id
  LEFT JOIN '{A}/A3/protein_pfam_arch.parquet' a ON p.protein_id = a.protein_id
  LEFT JOIN '{T}/swissprot_best.parquet' s ON p.protein_id = s.protein_id
  ORDER BY genome, p.protein_id
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 500000);
-- @output A4/fam_sp_genome.parquet
-- SwissProt-family proxy (P-A): one row per (sp_id, representative genome), strict gates.
COPY (
  SELECT c.sp_id, c.genome,
         any_value(c.org_code) AS org_code,
         any_value(c.s_len) AS s_len,
         count(*) AS n_prot_all,
         median(c.length) AS med_len_all,
         count(*) FILTER (WHERE c.q_model) AS n_prot_model,
         median(c.length) FILTER (WHERE c.q_model) AS med_len_model,
         median(c.n_domains) AS med_n_domains
  FROM '{A}/A3/protein_core.parquet' c
  JOIN '{A}/labels/genome_labels.parquet' l ON c.genome = l.LOCUSTAG
  WHERE l.is_representative AND c.sp_id IS NOT NULL
    AND c.qcov_x10 >= 800 AND c.scov_x10 >= 800 AND c.pident_x10 >= 400
  GROUP BY 1,2
  ORDER BY 1,2
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A4/fam_arch_genome.parquet
-- Pfam-architecture proxy (P-C): architectures with >= 2 domains, representative genomes.
COPY (
  SELECT c.pfam_arch, c.genome,
         count(*) AS n_prot_all,
         median(c.length) AS med_len_all,
         count(*) FILTER (WHERE c.q_model) AS n_prot_model,
         median(c.length) FILTER (WHERE c.q_model) AS med_len_model
  FROM '{A}/A3/protein_core.parquet' c
  JOIN '{A}/labels/genome_labels.parquet' l ON c.genome = l.LOCUSTAG
  WHERE l.is_representative AND c.n_domains >= 2
  GROUP BY 1,2
  ORDER BY 1,2
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A4/genome_protein_summary.parquet
-- Per-genome protein summary used by V07, V31, V20 denominators.
COPY (
  SELECT genome,
         count(*) AS n_proteins,
         count(*) FILTER (WHERE q_model) AS n_model,
         median(length) AS median_len_all,
         median(length) FILTER (WHERE q_model) AS median_len_model,
         quantile_cont(length, 0.25) FILTER (WHERE q_model) AS q25_len_model,
         quantile_cont(length, 0.75) FILTER (WHERE q_model) AS q75_len_model,
         avg(length) FILTER (WHERE q_model) AS mean_len_model,
         count(*) FILTER (WHERE n_domains >= 1) * 1.0 / count(*) AS frac_with_pfam,
         count(*) FILTER (WHERE n_domains >= 2) * 1.0 / count(*) AS frac_multidomain,
         median(n_domains) FILTER (WHERE n_domains >= 1) AS median_domains_if_any,
         count(*) FILTER (WHERE sp_id IS NOT NULL) * 1.0 / count(*) AS frac_with_swissprot,
         count(*) FILTER (WHERE qcov_x10 >= 800 AND scov_x10 >= 800 AND pident_x10 >= 400) AS n_strict_hits,
         count(*) FILTER (WHERE org_code = 'YEAST') * 1.0 / nullif(count(*) FILTER (WHERE sp_id IS NOT NULL), 0) AS frac_best_hit_yeast,
         count(*) FILTER (WHERE org_code IN ('NEUCR','EMENI')) * 1.0 / nullif(count(*) FILTER (WHERE sp_id IS NOT NULL), 0) AS frac_best_hit_neucr_emeni
  FROM '{A}/A3/protein_core.parquet'
  GROUP BY 1 ORDER BY 1
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
