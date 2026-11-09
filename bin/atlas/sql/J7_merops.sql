-- J7: MEROPS best hit per protein, family mapping (merops_fam.tsv, MER id -> family code), per-genome family counts.
CREATE OR REPLACE TABLE fam AS SELECT merops_id, nullif(merops_family,'') AS merops_family FROM read_csv('{A}/labels/merops_family_map.tsv', header=true, delim='\t', columns={'merops_id':'VARCHAR','merops_family':'VARCHAR'});
-- @output A5/merops_best.parquet
COPY (
  WITH b AS (
    SELECT species_prefix AS genome, protein_id, arg_max(merops_id, bitscore) AS merops_id, max(bitscore) AS bitscore,
           arg_max(percent_identity, bitscore) AS pident, arg_max(evalue, bitscore) AS evalue
    FROM '{T}/merops.parquet' WHERE evalue <= 1e-10 GROUP BY 1,2)
  SELECT b.*, f.merops_family, (s.protein_id IS NOT NULL) AS is_secreted
  FROM b LEFT JOIN fam f USING (merops_id) LEFT JOIN '{A}/A5/secreted_proteins.parquet' s USING (protein_id)
  ORDER BY genome, protein_id
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/merops_family_counts.parquet
COPY (
  SELECT genome, coalesce(merops_family, 'unmapped') AS merops_family, count(*) AS n_proteins,
         count(*) FILTER (WHERE is_secreted) AS n_secreted
  FROM '{A}/A5/merops_best.parquet' GROUP BY 1,2 ORDER BY 1,2
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/merops_mapping_check.parquet
COPY (
  SELECT count(*) AS n_best_hits, count(*) FILTER (WHERE merops_family IS NULL) AS n_unmapped,
         count(DISTINCT merops_id) AS n_ids, count(DISTINCT merops_id) FILTER (WHERE merops_family IS NULL) AS n_ids_unmapped
  FROM '{A}/A5/merops_best.parquet'
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
