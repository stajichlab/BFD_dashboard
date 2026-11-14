-- th0: build every input table of the trehalose survey (th1-th5) from the BFD parquet tables with DuckDB.
-- Replaces export_queries.sql (MotherDuck, run from a laptop) so that the survey can be rebuilt on the HPCC
-- (workflows/atlas/J13_trehalose.sh). Run with bin/atlas/run_duckdb_job.py: blocks are restartable, each block
-- writes one CSV into <--atlas>/data/. Placeholders: {T} = tables dir, {A} = output root, {OUT} = block output.
--
-- Differences from export_queries.sql (stated, because they can change numbers; see README section 0):
--  * reads parquet, not MotherDuck;
--  * the Pfam domain matrix is the atlas A7/pfam_genome_matrix.parquet (domain_i_evalue <= 1e-5, version stripped),
--    so `qcov` is gated at 1e-5 everywhere (export_queries.sql counted ungated hits in th_univcov and gated hits
--    in th_uvmat: Opus review finding);
--  * kw_genus_pfam (analysis B scratch table, no saved SQL) is rebuilt from that matrix over all annotated
--    genomes with a genus (verified against the MotherDuck copy in the J13 log: Coemansia / PF00982 must give
--    146 genomes, 131 present, 2.5753 mean copies);
--  * th_pfprot (share of proteins with any Pfam hit) had no producer in export_queries.sql; it is rebuilt here as
--    distinct proteins with >= 1 gated domain / proteins.

CREATE TEMP VIEW pm AS
  SELECT * FROM read_parquet('/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/A7/pfam_genome_matrix.parquet');  -- read-only; genome, pfam_acc, pfam_id, n_proteins, n_domains

CREATE TEMP TABLE reps AS
WITH g AS (
  SELECT s.LOCUSTAG AS genome, s.ASMID, s.PHYLUM, s.CLASS, s."ORDER" AS ord, s.FAMILY, s.GENUS, s.SPECIES,
         a.total_length_bp, a.N50_bp, a.contig_count, a.gc_pct, a.masked_pct,
         b.complete_pct, b.duplicated_pct, b.fragmented_pct,
         row_number() OVER (PARTITION BY s.SPECIES ORDER BY b.complete_pct DESC, a.N50_bp DESC, s.LOCUSTAG) AS rn
  FROM read_parquet('{T}/species.parquet') s
  JOIN read_parquet('{T}/asm_stats.parquet') a ON a.ASMID = s.ASMID
  JOIN read_parquet('{T}/busco_genome.parquet') b ON b.ASMID = s.ASMID
  WHERE a.annotation_status = 'annotated' AND s.SPECIES IS NOT NULL AND s.GENUS IS NOT NULL)
SELECT * EXCLUDE (rn) FROM g WHERE rn = 1;

CREATE TEMP TABLE allg AS                                  -- every annotated genome with a genus (kw_genus_pfam universe)
SELECT s.LOCUSTAG AS genome, s.GENUS
FROM read_parquet('{T}/species.parquet') s
JOIN read_parquet('{T}/asm_stats.parquet') a ON a.ASMID = s.ASMID
WHERE a.annotation_status = 'annotated' AND s.GENUS IS NOT NULL;

CREATE TEMP TABLE nprot AS
SELECT split_part(protein_id,'_',1) AS genome, count(*) AS n_prot, median(length) AS med_len,
       avg(length) AS mean_len, count(*) FILTER (WHERE length < 100) AS n_lt100,
       count(*) FILTER (WHERE length >= 1000) AS n_ge1000
FROM read_parquet('{T}/gene_proteins.parquet') GROUP BY 1;

CREATE TEMP TABLE kw AS                                    -- genus x Pfam, all annotated genomes
WITH ng AS (SELECT GENUS, count(*) AS n_genomes FROM allg GROUP BY 1),
     gp AS (SELECT a.GENUS, m.pfam_acc, count(*) AS n_present, sum(m.n_proteins) AS copies
            FROM pm m JOIN allg a ON a.genome = m.genome GROUP BY 1,2)
SELECT gp.GENUS, gp.pfam_acc, ng.n_genomes, gp.n_present AS n_genomes_present,
       gp.n_present::DOUBLE / ng.n_genomes AS presence_frac, gp.copies::DOUBLE / ng.n_genomes AS mean_copies
FROM gp JOIN ng USING (GENUS);

CREATE TEMP TABLE univ AS                                  -- same rule as export_queries.sql step 5
SELECT DISTINCT pfam_acc FROM kw
WHERE presence_frac >= 0.95 AND n_genomes >= 5
  AND pfam_acc IN (SELECT pfam_acc FROM kw GROUP BY 1 HAVING count(*) >= 1600);

-- @output data/th_check_kw.csv
-- verification rows for the rebuild of kw_genus_pfam and th_univ (MotherDuck copy: 146 / 131 / 2.5753 and 1,409 families)
COPY (
  SELECT 'Coemansia_PF00982' AS check_name, n_genomes::DOUBLE AS v1, n_genomes_present::DOUBLE AS v2, mean_copies AS v3
  FROM kw WHERE GENUS = 'Coemansia' AND pfam_acc = 'PF00982'
  UNION ALL SELECT 'n_universal_families', count(*), NULL, NULL FROM univ
  UNION ALL SELECT 'n_genera', count(DISTINCT GENUS), NULL, NULL FROM kw
  UNION ALL SELECT 'n_reps', count(*), NULL, NULL FROM reps
  UNION ALL SELECT 'n_rep_genera', count(DISTINCT GENUS), NULL, NULL FROM reps
) TO '{OUT}' (FORMAT CSV, HEADER);

-- @output data/th_genomes.csv
COPY (SELECT r.*, n.* EXCLUDE (genome) FROM reps r LEFT JOIN nprot n USING (genome)) TO '{OUT}' (FORMAT CSV, HEADER);

-- @output data/th_dom.csv
COPY (
  SELECT m.genome, m.pfam_acc AS acc, any_value(m.pfam_id) AS pfam_id, m.n_proteins AS n_prot, m.n_domains AS n_dom
  FROM pm m JOIN reps r ON r.genome = m.genome
  WHERE m.pfam_acc IN (
   'PF00982','PF02358','PF01204','PF07492','PF03632','PF03636','PF03633',
   'PF05693','PF00343','PF00128','PF06202','PF01704','PF02878','PF02879','PF02880','PF00342','PF00349',
   'PF07479','PF01210','PF01232','PF08125','PF00230','PF00702',
   'PF00011','PF04119','PF01185','PF06766','PF00199','PF00080','PF00012','PF11765','PF05141','PF03907')
  GROUP BY m.genome, m.pfam_acc, m.n_proteins, m.n_domains
) TO '{OUT}' (FORMAT CSV, HEADER);

-- @output data/th_sp.csv
COPY (
  SELECT split_part(b.protein_id,'_',1) AS genome,
    count(*) FILTER (WHERE lower(b.protein_name) LIKE '%trehalose-phosphate synthase%') AS n_tps,
    count(*) FILTER (WHERE lower(b.protein_name) LIKE '%trehalose-phosphatase%'
                          OR lower(b.protein_name) LIKE '%trehalose-6-phosphate phosphatase%') AS n_tpp,
    count(*) FILTER (WHERE lower(b.protein_name) LIKE '%neutral trehalase%') AS n_nth,
    count(*) FILTER (WHERE lower(b.protein_name) LIKE '%acid trehalase%') AS n_ath,
    count(*) FILTER (WHERE lower(b.protein_name) LIKE '%glycogen synthase%') AS n_gsy,
    count(*) FILTER (WHERE lower(b.protein_name) LIKE '%glycogen phosphorylase%') AS n_gph,
    count(*) FILTER (WHERE lower(b.protein_name) LIKE '%glycerol-3-phosphate dehydrogenase%') AS n_gpd,
    count(*) FILTER (WHERE lower(b.protein_name) LIKE '%mannitol%dehydrogenase%') AS n_mfn
  FROM read_parquet('{T}/swissprot_best.parquet') b JOIN reps r ON r.genome = split_part(b.protein_id,'_',1)
  WHERE lower(b.protein_name) LIKE '%trehal%' OR lower(b.protein_name) LIKE '%glycogen%'
     OR lower(b.protein_name) LIKE '%glycerol-3-phosphate%' OR lower(b.protein_name) LIKE '%mannitol%'
  GROUP BY 1
) TO '{OUT}' (FORMAT CSV, HEADER);

-- @output data/th_pfprot.csv
COPY (
  SELECT split_part(protein_id,'_',1) AS genome, count(DISTINCT protein_id) AS n_prot_hit
  FROM read_parquet('{T}/pfam.parquet')
  WHERE domain_i_evalue <= 1e-5 AND split_part(protein_id,'_',1) IN (SELECT genome FROM reps)
  GROUP BY 1
) TO '{OUT}' (FORMAT CSV, HEADER);

-- @output data/th_univ.csv
COPY (SELECT pfam_acc FROM univ ORDER BY 1) TO '{OUT}' (FORMAT CSV, HEADER);

-- @output data/th_univcov.csv
COPY (
  SELECT r.genome, r.GENUS, r.CLASS, count(m.pfam_acc) AS n_univ_found
  FROM reps r LEFT JOIN (SELECT pm.genome, pm.pfam_acc FROM pm JOIN univ USING (pfam_acc)) m ON m.genome = r.genome
  GROUP BY 1,2,3
) TO '{OUT}' (FORMAT CSV, HEADER);

-- @output data/th_uvmiss.csv
-- genome x universal-family absences (gated), for genomes with n_univ_found >= 0.98 x number of universal families
COPY (
  WITH cov AS (SELECT r.genome, count(m.pfam_acc) AS n_found
               FROM reps r LEFT JOIN (SELECT pm.genome, pm.pfam_acc FROM pm JOIN univ USING (pfam_acc)) m ON m.genome = r.genome
               GROUP BY 1),
       g AS (SELECT genome FROM cov WHERE n_found >= 0.98 * (SELECT count(*) FROM univ))
  SELECT r.genome, u.pfam_acc, 1 AS missing
  FROM g JOIN reps r USING (genome) CROSS JOIN univ u
  WHERE NOT EXISTS (SELECT 1 FROM pm WHERE pm.genome = r.genome AND pm.pfam_acc = u.pfam_acc)
) TO '{OUT}' (FORMAT CSV, HEADER);

-- @output data/th_allobin.csv
-- allometry null: mean copies per family within deciles of genus proteome size (median n_prot over species representatives)
COPY (
  WITH gprot AS (SELECT GENUS, median(n_prot) AS n_prot FROM reps JOIN nprot USING (genome) GROUP BY 1),
       j AS (SELECT k.pfam_acc, k.GENUS, k.mean_copies, p.n_prot, ntile(10) OVER (ORDER BY p.n_prot) AS bin
             FROM kw k JOIN gprot p ON p.GENUS = k.GENUS WHERE p.n_prot > 0)
  SELECT pfam_acc, bin, count(*) AS n_bin_genus, avg(mean_copies) AS mean_copies, avg(n_prot) AS mean_prot
  FROM j GROUP BY 1,2
) TO '{OUT}' (FORMAT CSV, HEADER);
