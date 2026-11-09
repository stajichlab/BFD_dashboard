-- J2: Pfam-derived intermediates. Domain gate: domain_i_evalue <= 1e-5. Accession version stripped.
-- @output A7/pfam_genome_matrix.parquet
COPY (
  SELECT split_part(protein_id,'_',1) AS genome,
         split_part(pfam_acc,'.',1) AS pfam_acc,
         any_value(pfam_id) AS pfam_id,
         count(DISTINCT protein_id) AS n_proteins,
         count(*) AS n_domains
  FROM '{T}/pfam.parquet'
  WHERE domain_i_evalue <= 1e-5
  GROUP BY 1,2
  ORDER BY 1,2
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A3/protein_pfam_arch.parquet
COPY (
  SELECT protein_id,
         count(*) AS n_domains,
         count(DISTINCT split_part(pfam_acc,'.',1)) AS n_distinct_domains,
         string_agg(split_part(pfam_acc,'.',1), ';' ORDER BY env_from, env_to) AS pfam_arch,
         min(env_from) AS first_env_from,
         max(env_to) AS last_env_to
  FROM '{T}/pfam.parquet'
  WHERE domain_i_evalue <= 1e-5
  GROUP BY 1
  ORDER BY 1
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
