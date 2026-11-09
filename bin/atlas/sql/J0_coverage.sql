-- J0: genome coverage of every large table (count distinct genome key). One small output per table.
-- @output coverage/pfam.parquet
COPY (SELECT 'pfam' AS tbl, 'split_part(protein_id)' AS key, count(DISTINCT split_part(protein_id,'_',1)) AS n_genomes, count(*) AS n_rows FROM '{T}/pfam.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/swissprot_best.parquet
COPY (SELECT 'swissprot_best' AS tbl, 'species_prefix' AS key, count(DISTINCT species_prefix) AS n_genomes, count(*) AS n_rows FROM '{T}/swissprot_best.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/gene_info.parquet
COPY (SELECT 'gene_info' AS tbl, 'LOCUSTAG' AS key, count(DISTINCT LOCUSTAG) AS n_genomes, count(*) AS n_rows FROM '{T}/gene_info.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/gene_transcripts.parquet
COPY (SELECT 'gene_transcripts' AS tbl, 'split_part(gene_id)' AS key, count(DISTINCT split_part(gene_id,'_',1)) AS n_genomes, count(*) AS n_rows FROM '{T}/gene_transcripts.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/gene_proteins.parquet
COPY (SELECT 'gene_proteins' AS tbl, 'split_part(gene_id)' AS key, count(DISTINCT split_part(gene_id,'_',1)) AS n_genomes, count(*) AS n_rows FROM '{T}/gene_proteins.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/gene_introns.parquet
COPY (SELECT 'gene_introns' AS tbl, 'split_part(transcript_id)' AS key, count(DISTINCT split_part(transcript_id,'_',1)) AS n_genomes, count(*) AS n_rows FROM '{T}/gene_introns.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/gene_intergenic_distances.parquet
COPY (SELECT 'gene_intergenic_distances' AS tbl, 'species_prefix' AS key, count(DISTINCT species_prefix) AS n_genomes, count(*) AS n_rows FROM '{T}/gene_intergenic_distances.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/merops.parquet
COPY (SELECT 'merops' AS tbl, 'species_prefix' AS key, count(DISTINCT species_prefix) AS n_genomes, count(*) AS n_rows FROM '{T}/merops.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/tmhmm.parquet
COPY (SELECT 'tmhmm' AS tbl, 'species_prefix' AS key, count(DISTINCT species_prefix) AS n_genomes, count(*) AS n_rows FROM '{T}/tmhmm.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/wolfpsort.parquet
COPY (SELECT 'wolfpsort' AS tbl, 'species_prefix' AS key, count(DISTINCT species_prefix) AS n_genomes, count(*) AS n_rows FROM '{T}/wolfpsort.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/predgpi.parquet
COPY (SELECT 'predgpi' AS tbl, 'species_prefix' AS key, count(DISTINCT species_prefix) AS n_genomes, count(*) AS n_rows FROM '{T}/predgpi.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/idp_summary.parquet
COPY (SELECT 'idp_summary' AS tbl, 'species_prefix' AS key, count(DISTINCT species_prefix) AS n_genomes, count(*) AS n_rows, count(*) FILTER (WHERE protein_id LIKE '%-T%') AS n_rows_with_T_suffix FROM '{T}/idp_summary.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/idp.parquet
COPY (SELECT 'idp' AS tbl, 'species_prefix' AS key, count(DISTINCT species_prefix) AS n_genomes, count(*) AS n_rows FROM '{T}/idp.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/gene_CDS.parquet
COPY (SELECT 'gene_CDS' AS tbl, 'split_part(transcript_id)' AS key, count(DISTINCT split_part(transcript_id,'_',1)) AS n_genomes, count(*) AS n_rows FROM '{T}/gene_CDS.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/gene_exons.parquet
COPY (SELECT 'gene_exons' AS tbl, 'split_part(transcript_id)' AS key, count(DISTINCT split_part(transcript_id,'_',1)) AS n_genomes, count(*) AS n_rows FROM '{T}/gene_exons.parquet') TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output coverage/swissprot_hits.parquet
COPY (SELECT 'swissprot_hits' AS tbl, 'species_prefix' AS key, count(DISTINCT species_prefix) AS n_genomes, count(*) AS n_rows FROM read_parquet('{T}/swissprot_hits/*.parquet')) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
