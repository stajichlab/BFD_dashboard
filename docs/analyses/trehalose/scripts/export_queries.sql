-- MotherDuck queries behind docs/analyses/trehalose/ (run with duckdb CLI v1.5.6 +
-- motherduck extension; the repo's brew duckdb v2.0.0 is not yet supported by
-- MotherDuck). Scratch tables follow the analysis-prefix convention in
-- docs/FUNCTIONAL_ANALYSES.md; prefix `th_` = trehalose survey.
--
--   duckdb -c "ATTACH 'md:bfd' AS bfd (TYPE motherduck);
--              ATTACH 'md:bfd_scratch' AS scr (TYPE motherduck);" -f export_queries.sql
--
-- bfd is read-only throughout; everything written lives in bfd_scratch.th_*.

-- 1. species representatives, rebuilt from bfd to the atlas rule (highest
--    complete_pct, tie N50_bp, tie LOCUSTAG; annotated; genus+species known)
CREATE OR REPLACE TABLE scr.th_reps AS
WITH g AS (
  SELECT s.LOCUSTAG genome, s.ASMID, s.PHYLUM, s.CLASS, s."ORDER" ord, s.FAMILY, s.GENUS, s.SPECIES,
         a.total_length_bp, a.N50_bp, a.contig_count, a.gc_pct, a.masked_pct,
         b.complete_pct, b.duplicated_pct, b.fragmented_pct,
         row_number() OVER (PARTITION BY s.SPECIES ORDER BY b.complete_pct DESC, a.N50_bp DESC, s.LOCUSTAG) rn
  FROM bfd.species s
  JOIN bfd.asm_stats a ON a.ASMID = s.ASMID
  JOIN bfd.busco_genome b ON b.ASMID = s.ASMID
  WHERE a.annotation_status = 'annotated' AND s.SPECIES IS NOT NULL AND s.GENUS IS NOT NULL)
SELECT * FROM g WHERE rn = 1;                       -- 7,186 genomes, 1,669 genera

-- 2. module copy numbers: one scan of pfam (368 M rows) for 33 accessions.
--    Domain gate domain_i_evalue <= 1e-5 (atlas P-C rule); a protein counts once.
CREATE OR REPLACE TABLE scr.th_dom AS
SELECT split_part(p.protein_id,'_',1) genome, split_part(p.pfam_acc,'.',1) acc,
       any_value(p.pfam_id) pfam_id, count(DISTINCT p.protein_id) n_prot, count(*) n_dom
FROM bfd.pfam p JOIN scr.th_reps r ON r.genome = split_part(p.protein_id,'_',1)
WHERE split_part(p.pfam_acc,'.',1) IN (
 'PF00982','PF02358','PF01204','PF07492','PF03632','PF03636','PF03633',
 'PF05693','PF00343','PF00128','PF06202','PF01704','PF02878','PF02879','PF02880','PF00342','PF00349',
 'PF07479','PF01210','PF01232','PF08125','PF00230','PF00702',
 'PF00011','PF04119','PF01185','PF06766','PF00199','PF00080','PF00012','PF11765','PF05141','PF03907')
  AND p.domain_i_evalue <= 1e-5
GROUP BY 1,2;

-- 3. proteome size and protein-length descriptors (one scan of gene_proteins)
CREATE OR REPLACE TABLE scr.th_nprot AS
SELECT split_part(gp.protein_id,'_',1) genome, count(*) n_prot, median(gp.length) med_len,
       avg(gp.length) mean_len, count(*) FILTER (WHERE gp.length < 100) n_lt100,
       count(*) FILTER (WHERE gp.length >= 1000) n_ge1000
FROM bfd.gene_proteins gp JOIN scr.th_reps r ON r.genome = split_part(gp.protein_id,'_',1)
GROUP BY 1;

-- 4. second, name-based occurrence check that does not use Pfam at all:
--    diamond best SwissProt hit with a trehalose/glycogen/enzyme name
CREATE OR REPLACE TABLE scr.th_sp AS
SELECT split_part(b.protein_id,'_',1) genome,
  count(*) FILTER (WHERE lower(b.protein_name) LIKE '%trehalose-phosphate synthase%') n_tps,
  count(*) FILTER (WHERE lower(b.protein_name) LIKE '%trehalose-phosphatase%'
                        OR lower(b.protein_name) LIKE '%trehalose-6-phosphate phosphatase%') n_tpp,
  count(*) FILTER (WHERE lower(b.protein_name) LIKE '%neutral trehalase%') n_nth,
  count(*) FILTER (WHERE lower(b.protein_name) LIKE '%acid trehalase%') n_ath,
  count(*) FILTER (WHERE lower(b.protein_name) LIKE '%glycogen synthase%') n_gsy,
  count(*) FILTER (WHERE lower(b.protein_name) LIKE '%glycogen phosphorylase%') n_gph,
  count(*) FILTER (WHERE lower(b.protein_name) LIKE '%glycerol-3-phosphate dehydrogenase%') n_gpd,
  count(*) FILTER (WHERE lower(b.protein_name) LIKE '%mannitol%dehydrogenase%') n_mfn
FROM bfd.swissprot_best b JOIN scr.th_reps r ON r.genome = split_part(b.protein_id,'_',1)
WHERE lower(b.protein_name) LIKE '%trehal%' OR lower(b.protein_name) LIKE '%glycogen%'
   OR lower(b.protein_name) LIKE '%glycerol-3-phosphate%' OR lower(b.protein_name) LIKE '%mannitol%'
GROUP BY 1;

-- 5. annotation-completeness control: coverage of the 1,409 near-universal
--    families (present in >=95 % of genera of the kingdom-wide kw_genus_pfam)
CREATE OR REPLACE TABLE scr.th_univ AS
SELECT DISTINCT pfam_acc FROM scr.kw_genus_pfam
WHERE presence_frac >= 0.95 AND n_genomes >= 5
  AND pfam_acc IN (SELECT pfam_acc FROM (SELECT pfam_acc, count(*) n_gen, avg(presence_frac) mp
                                         FROM scr.kw_genus_pfam GROUP BY 1) WHERE n_gen >= 1600);
CREATE OR REPLACE TABLE scr.th_univcov AS
WITH h AS (SELECT split_part(p.protein_id,'_',1) g, split_part(p.pfam_acc,'.',1) acc
           FROM bfd.pfam p JOIN scr.th_univ u ON u.pfam_acc = split_part(p.pfam_acc,'.',1)
           WHERE split_part(p.protein_id,'_',1) IN (SELECT genome FROM scr.th_reps)
           GROUP BY 1,2)
SELECT r.genome, r.GENUS, r.CLASS, count(h.acc) n_univ_found
FROM scr.th_reps r LEFT JOIN h ON h.g = r.genome GROUP BY 1,2,3;

-- 6. loose (ungated) hit check: do the "absent" genomes carry sub-1e-5 records?
CREATE OR REPLACE TABLE scr.th_loose AS
WITH s AS (SELECT split_part(p.protein_id,'_',1) g, split_part(p.pfam_acc,'.',1) acc, p.domain_i_evalue ie
           FROM bfd.pfam p JOIN scr.th_reps r ON r.genome = split_part(p.protein_id,'_',1)
           WHERE split_part(p.pfam_acc,'.',1) IN ('PF00982','PF02358','PF18572','PF21141'))
SELECT g, acc, count(*) n_hits, min(ie) best_ie, count(*) FILTER (WHERE ie<=1e-5) n_gate
FROM s GROUP BY 1,2;

-- 7. allometry null: mean copies per family within deciles of proteome size,
--    from the kingdom-wide genus x Pfam table of analysis B
CREATE OR REPLACE TABLE scr.th_gprot AS
SELECT * FROM read_csv('<LOCAL>/genus_prot.csv', header=true);   -- GENUS, n_prot from step 3
CREATE OR REPLACE TABLE scr.th_allobin AS
WITH j AS (SELECT k.pfam_acc, k.GENUS, k.mean_copies, p.n_prot, ntile(10) OVER (ORDER BY p.n_prot) bin
           FROM scr.kw_genus_pfam k JOIN scr.th_gprot p ON p.GENUS = k.GENUS WHERE p.n_prot > 0)
SELECT pfam_acc, bin, count(*) n_bin_genus, avg(mean_copies) mean_copies, avg(n_prot) mean_prot
FROM j GROUP BY 1,2;

-- 8. genome x near-universal-family absence matrix for the co-absence null
CREATE OR REPLACE TABLE scr.th_uvmat AS
WITH g AS (SELECT genome FROM scr.th_univcov WHERE n_univ_found >= 0.98*1409),
     h AS (SELECT DISTINCT split_part(p.protein_id,'_',1) genome, split_part(p.pfam_acc,'.',1) acc
           FROM bfd.pfam p JOIN scr.th_univ u ON u.pfam_acc = split_part(p.pfam_acc,'.',1)
           WHERE p.domain_i_evalue <= 1e-5
             AND split_part(p.protein_id,'_',1) IN (SELECT genome FROM g))
SELECT r.genome, r.GENUS, r.CLASS, u.pfam_acc, CASE WHEN h.acc IS NULL THEN 1 ELSE 0 END missing
FROM scr.th_reps r CROSS JOIN scr.th_univ u
LEFT JOIN h ON h.genome = r.genome AND h.acc = u.pfam_acc
WHERE r.genome IN (SELECT genome FROM g);

-- exports consumed by th1..th4
COPY scr.th_univcov TO 'th_univcov.csv' (HEADER, DELIMITER ',');
COPY scr.th_sp      TO 'th_sp.csv'      (HEADER, DELIMITER ',');
COPY scr.th_allobin TO 'th_allobin.csv' (HEADER, DELIMITER ',');
COPY (SELECT genome, pfam_acc, missing FROM scr.th_uvmat WHERE missing=1) TO 'th_uvmiss.csv' (HEADER, DELIMITER ',');
COPY (SELECT r.*, n.* EXCLUDE (genome) FROM scr.th_reps r LEFT JOIN scr.th_nprot n USING (genome))
  TO 'th_genomes.csv' (HEADER, DELIMITER ',');
COPY scr.th_dom TO 'th_dom.csv' (HEADER, DELIMITER ',');

-- Reference note: bfd.pfam keeps weak hits (8.4 % of 368 M rows have
-- domain_i_evalue > 1e-5, max 1.8), but full_seq_e_value never exceeds 0.026, so
-- the upstream hmmsearch ran near Pfam's reporting threshold. A family can be
-- missing from the table because the protein is absent, because the gene model
-- is missing, or because it falls below the report threshold; all three need
-- tblastn against the assembly to separate.
