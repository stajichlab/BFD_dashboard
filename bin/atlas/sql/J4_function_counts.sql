-- J4: per-genome function and localization counts (A5). genome = species_prefix.
-- @output A5/secretome.parquet
-- SignalP proteins with TM and GPI flags; strict secretome = prob >= 0.9, <= 1 TM helix, no GPI anchor.
COPY (
  WITH sp AS (
    SELECT species_prefix AS genome, protein_id, max(probability) AS sp_prob
    FROM '{T}/signalp.signal_peptide.parquet' GROUP BY 1,2
  ), tm AS (
    SELECT protein_id, max(PredHel) AS n_helix FROM '{T}/tmhmm.parquet' GROUP BY 1
  ), gpi AS (
    SELECT DISTINCT protein_id FROM '{T}/predgpi.parquet' WHERE feature = 'GPI-anchor'
  ), tm_genomes AS (
    SELECT species_prefix AS genome, count(*) AS n_tm_rows, count(*) FILTER (WHERE PredHel >= 2) AS n_tm2,
           count(*) FILTER (WHERE PredHel >= 1) AS n_tm1
    FROM '{T}/tmhmm.parquet' GROUP BY 1
  ), gpi_genomes AS (
    SELECT species_prefix AS genome, count(*) FILTER (WHERE feature = 'GPI-anchor') AS n_gpi,
           count(*) FILTER (WHERE feature = 'Chain') AS n_predgpi_rows
    FROM '{T}/predgpi.parquet' GROUP BY 1
  ), j AS (
    SELECT sp.genome, sp.protein_id, sp.sp_prob, tm.n_helix, (gpi.protein_id IS NOT NULL) AS is_gpi
    FROM sp LEFT JOIN tm USING (protein_id) LEFT JOIN gpi USING (protein_id)
  )
  SELECT j.genome,
         count(*) AS n_signalp,
         count(*) FILTER (WHERE sp_prob >= 0.9) AS n_signalp_09,
         count(*) FILTER (WHERE sp_prob >= 0.9 AND coalesce(n_helix, 0) <= 1 AND NOT is_gpi) AS n_secreted_strict,
         count(*) FILTER (WHERE sp_prob >= 0.9 AND coalesce(n_helix, 0) <= 1) AS n_secreted_noTM,
         count(*) FILTER (WHERE is_gpi) AS n_signalp_gpi,
         count(*) FILTER (WHERE n_helix IS NULL) AS n_signalp_no_tmhmm_row,
         any_value(tg.n_tm_rows) AS n_tmhmm_rows, any_value(tg.n_tm2) AS n_tm2, any_value(tg.n_tm1) AS n_tm1,
         any_value(gg.n_gpi) AS n_gpi, any_value(gg.n_predgpi_rows) AS n_predgpi_rows
  FROM j LEFT JOIN tm_genomes tg USING (genome) LEFT JOIN gpi_genomes gg USING (genome)
  GROUP BY 1 ORDER BY 1
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/secreted_proteins.parquet
-- protein list of the strict secretome (for CAZy / MEROPS secreted joins).
COPY (
  WITH sp AS (
    SELECT species_prefix AS genome, protein_id, max(probability) AS sp_prob
    FROM '{T}/signalp.signal_peptide.parquet' GROUP BY 1,2
  ), tm AS (SELECT protein_id, max(PredHel) AS n_helix FROM '{T}/tmhmm.parquet' GROUP BY 1),
  gpi AS (SELECT DISTINCT protein_id FROM '{T}/predgpi.parquet' WHERE feature = 'GPI-anchor')
  SELECT sp.genome, sp.protein_id FROM sp LEFT JOIN tm USING (protein_id) LEFT JOIN gpi USING (protein_id)
  WHERE sp.sp_prob >= 0.9 AND coalesce(tm.n_helix, 0) <= 1 AND gpi.protein_id IS NULL
  ORDER BY 1,2
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/wolfpsort_counts.parquet
COPY (
  SELECT species_prefix AS genome, localization, count(*) AS n
  FROM '{T}/wolfpsort.parquet' GROUP BY 1,2 ORDER BY 1,2
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/cazy_class_counts.parquet
-- cazyme_fam like 'GH27(12-262)+CBM35(...)'; class = leading letters of each family token. toolcount >= 2 (D7).
COPY (
  WITH u AS (
    SELECT species_prefix AS genome, protein_id, unnest(string_split(cazyme_fam, '+')) AS fam_tok
    FROM '{T}/cazy.overview.parquet' WHERE toolcount >= 2 AND cazyme_fam IS NOT NULL AND cazyme_fam <> '-'
  )
  SELECT genome, regexp_extract(fam_tok, '^([A-Za-z]+)', 1) AS cazy_class,
         regexp_extract(fam_tok, '^([A-Za-z]+[0-9]+)', 1) AS cazy_family,
         count(DISTINCT protein_id) AS n_proteins
  FROM u GROUP BY 1,2,3 ORDER BY 1,2,3
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/cazy_substrate_counts.parquet
COPY (
  WITH u AS (
    SELECT species_prefix AS genome, protein_id, trim(unnest(string_split(substrate, '+'))) AS substrate
    FROM '{T}/cazy.overview.parquet' WHERE toolcount >= 2 AND substrate IS NOT NULL AND substrate <> '-'
  )
  SELECT genome, substrate, count(DISTINCT protein_id) AS n_proteins FROM u GROUP BY 1,2 ORDER BY 1,2
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/cazy_genome.parquet
COPY (
  WITH c AS (
    SELECT species_prefix AS genome, protein_id FROM '{T}/cazy.overview.parquet' WHERE toolcount >= 2 GROUP BY 1,2
  )
  SELECT c.genome, count(*) AS n_cazy,
         count(*) FILTER (WHERE s.protein_id IS NOT NULL) AS n_cazy_secreted
  FROM c LEFT JOIN '{A}/A5/secreted_proteins.parquet' s USING (protein_id) GROUP BY 1 ORDER BY 1
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/merops_genome.parquet
COPY (
  WITH m AS (
    SELECT species_prefix AS genome, protein_id, arg_max(merops_id, bitscore) AS merops_id
    FROM '{T}/merops.parquet' WHERE evalue <= 1e-10 GROUP BY 1,2
  )
  SELECT m.genome, count(*) AS n_merops,
         count(*) FILTER (WHERE s.protein_id IS NOT NULL) AS n_merops_secreted
  FROM m LEFT JOIN '{A}/A5/secreted_proteins.parquet' s USING (protein_id) GROUP BY 1 ORDER BY 1
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/tf_counts.parquet
COPY (
  SELECT species_prefix AS genome, tf_family, count(DISTINCT protein_id) AS n_proteins
  FROM '{T}/tf_inventory.parquet' WHERE domain_i_evalue <= 1e-5 GROUP BY 1,2 ORDER BY 1,2
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/tf_genome.parquet
COPY (
  SELECT species_prefix AS genome, count(DISTINCT protein_id) AS n_tf_proteins
  FROM '{T}/tf_inventory.parquet' WHERE domain_i_evalue <= 1e-5 GROUP BY 1 ORDER BY 1
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/idp_genome.parquet
COPY (
  SELECT species_prefix AS genome, count(*) AS n_idp_rows,
         count(*) FILTER (WHERE protein_id LIKE '%-T%') AS n_rows_T_suffix,
         median(IDP_fraction) AS median_idp_fraction,
         count(*) FILTER (WHERE IDP_fraction >= 0.5) * 1.0 / count(*) AS frac_idp_ge50
  FROM '{T}/idp_summary.parquet' GROUP BY 1 ORDER BY 1
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/trna_codon_counts.parquet
COPY (
  SELECT split_part(gene_id,'_',1) AS genome, amino_acid, codon, count(*) AS n
  FROM '{T}/gene_trnas.parquet' GROUP BY 1,2,3 ORDER BY 1,2,3
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A5/wgd_genome.parquet
COPY (
  SELECT species_prefix AS genome,
         count(*) FILTER (WHERE kind = 'tree') AS n_tree_nodes,
         count(*) FILTER (WHERE kind = 'tree' AND n_valid >= 1 AND dS_mean BETWEEN 0.05 AND 1.0) AS n_young_nodes,
         count(*) FILTER (WHERE kind = 'tree' AND n_valid >= 1 AND dS_mean < 0.05) AS n_nearzero_nodes,
         count(*) FILTER (WHERE kind = 'identical') AS n_identical_nodes
  FROM '{T}/wgd.ks.nodes.parquet' GROUP BY 1 ORDER BY 1
) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
