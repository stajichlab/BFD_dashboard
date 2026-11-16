-- 04_score.sql
-- Purpose : per-genome meiotic toolkit completeness score. Definition FIXED BEFORE any group value of the score was viewed (2026-11-12).
--           One point for each of seven informative genes detected by ANY of its specific evidence:
--             Hop2   : WHD_TBPIP PF07106 (<= 1e-3) or best-hit name 'homologous-pairing protein 2' / 'TBP-1-interacting'
--             Mnd1   : Mnd1 PF03962 (<= 1e-3) or best-hit name 'meiotic nuclear division protein 1'
--             Msh4   : best-hit name 'MutS protein homolog 4' / 'mismatch repair protein MSH4'  (no specific Pfam domain)
--             Msh5   : MutS5_N PF27035 (<= 1e-3) or best-hit name 'MutS protein homolog 5' / 'mismatch repair protein MSH5'
--             Spo11  : TOP6A-Spo11_Toprim PF21180 or TP6A_N PF04406 (<= 1e-3) or best-hit name Spo11 / rec12 / spo-11
--             Pch2   : TRIP13 lid PF23242 or TRIP13_N PF23563 (<= 1e-3) or best-hit name 'pachytene checkpoint' / 'TRIP13'
--             Mer3   : the Mer3/Hfm1 rule of 01_rules.sql (DEAD + Helicase_C + Sec63 and name) or best-hit name 'helicase MER3' / 'helicase HFM1'
--           Excluded as uninformative (not detectable by Pfam outside yeasts even in the outgroups): Mei4, Mei5, Mer2, Rec114, Red1, Spo16, Mek1.
--           Dmc1 is NOT in the score: it is classified by the RecA tree (../../family_trees/RECA_RESULT.md) for only 10 of the 107 Coemansia strains.
-- Name evidence: bfd_scratch.kx2_ku_sp best SwissProt hit with pident >= 25, scov >= 40, -log10(e) >= 10.
-- Inputs  : bfd_scratch.kx2_mt_pfam, bfd_scratch.kx2_ku_sp, bfd_scratch.kx2_mt_matrix, bfd_scratch.kx2_ku_genomes
-- Outputs : bfd_scratch.kx2_mt_score (one row per genome: seven flags, score 0-7, BUSCO, N50); a group summary query
-- Reports : README 'toolkit completeness score' table; also the Spearman correlation of score with BUSCO inside Coemansia (assembly effect).
CREATE OR REPLACE TABLE bfd_scratch.kx2_mt_score AS
WITH nm AS (
  SELECT s.locustag,
    bool_or(regexp_matches(s.protein_name,'(?i)homologous-pairing protein 2|TBP-1-interacting')) AS hop2_n,
    bool_or(regexp_matches(s.protein_name,'(?i)meiotic nuclear division protein 1')) AS mnd1_n,
    bool_or(regexp_matches(s.protein_name,'(?i)mutS protein homolog 4|mismatch repair protein MSH4')) AS msh4_n,
    bool_or(regexp_matches(s.protein_name,'(?i)mutS protein homolog 5|mismatch repair protein MSH5')) AS msh5_n,
    bool_or(regexp_matches(s.protein_name,'(?i)recombination protein SPO11|spo-11|recombination protein rec12|SPO11-')) AS spo11_n,
    bool_or(regexp_matches(s.protein_name,'(?i)pachytene checkpoint|trip13')) AS pch2_n,
    bool_or(regexp_matches(s.protein_name,'(?i)helicase MER3|helicase HFM1')) AS mer3_n
  FROM bfd_scratch.kx2_ku_sp s WHERE s.pident>=25 AND s.scov>=40 AND s.nlog10e>=10 GROUP BY s.locustag),
pf AS (
  SELECT locustag,
    bool_or(acc='PF07106') AS hop2_p, bool_or(acc='PF03962') AS mnd1_p, bool_or(acc='PF27035') AS msh5_p,
    bool_or(acc IN ('PF21180','PF04406')) AS spo11_p, bool_or(acc IN ('PF23242','PF23563')) AS pch2_p
  FROM bfd_scratch.kx2_mt_pfam GROUP BY locustag),
mr AS (SELECT locustag, bool_or(n_present>0) AS mer3_r FROM bfd_scratch.kx2_mt_matrix WHERE gene='Mer3/Hfm1' GROUP BY locustag)
SELECT g.locustag, g.genus, g.grp, g.busco_complete, g.n50, g.n_proteins,
  coalesce(pf.hop2_p,false) OR coalesce(nm.hop2_n,false) AS hop2,
  coalesce(pf.mnd1_p,false) OR coalesce(nm.mnd1_n,false) AS mnd1,
  coalesce(nm.msh4_n,false) AS msh4,
  coalesce(pf.msh5_p,false) OR coalesce(nm.msh5_n,false) AS msh5,
  coalesce(pf.spo11_p,false) OR coalesce(nm.spo11_n,false) AS spo11,
  coalesce(pf.pch2_p,false) OR coalesce(nm.pch2_n,false) AS pch2,
  coalesce(mr.mer3_r,false) OR coalesce(nm.mer3_n,false) AS mer3
FROM bfd_scratch.kx2_ku_genomes g
LEFT JOIN nm USING (locustag) LEFT JOIN pf USING (locustag) LEFT JOIN mr USING (locustag);

-- group summary: per gene detection (genomes with the gene / genomes) and the score (median, share of genomes with score >= 6)
SELECT grp, count(*) AS n,
  count(*) FILTER (WHERE hop2) AS hop2, count(*) FILTER (WHERE mnd1) AS mnd1, count(*) FILTER (WHERE msh4) AS msh4,
  count(*) FILTER (WHERE msh5) AS msh5, count(*) FILTER (WHERE spo11) AS spo11, count(*) FILTER (WHERE pch2) AS pch2, count(*) FILTER (WHERE mer3) AS mer3,
  median(hop2::INT + mnd1::INT + msh4::INT + msh5::INT + spo11::INT + pch2::INT + mer3::INT) AS median_score,
  count(*) FILTER (WHERE (hop2::INT + mnd1::INT + msh4::INT + msh5::INT + spo11::INT + pch2::INT + mer3::INT) >= 6) AS n_score_ge6
FROM bfd_scratch.kx2_mt_score GROUP BY grp ORDER BY grp;

-- assembly effect inside Coemansia: score and per-gene detection by BUSCO bin (separates real absence from annotation-in-fragments)
SELECT CASE WHEN busco_complete >= 90 THEN '1 BUSCO >= 90' WHEN busco_complete >= 80 THEN '2 BUSCO 80-90' WHEN busco_complete >= 70 THEN '3 BUSCO 70-80' ELSE '4 BUSCO < 70' END AS bin,
  count(*) AS n, round(median(hop2::INT + mnd1::INT + msh4::INT + msh5::INT + spo11::INT + pch2::INT + mer3::INT), 1) AS median_score,
  count(*) FILTER (WHERE hop2) AS hop2, count(*) FILTER (WHERE mnd1) AS mnd1, count(*) FILTER (WHERE msh4) AS msh4, count(*) FILTER (WHERE msh5) AS msh5,
  count(*) FILTER (WHERE spo11) AS spo11, count(*) FILTER (WHERE pch2) AS pch2, count(*) FILTER (WHERE mer3) AS mer3
FROM bfd_scratch.kx2_mt_score WHERE genus='Coemansia' GROUP BY 1 ORDER BY 1;
-- correlation of the score with assembly quality inside Coemansia
SELECT corr(hop2::INT + mnd1::INT + msh4::INT + msh5::INT + spo11::INT + pch2::INT + mer3::INT, busco_complete) AS pearson_score_busco,
       corr(hop2::INT + mnd1::INT + msh4::INT + msh5::INT + spo11::INT + pch2::INT + mer3::INT, ln(n50)) AS pearson_score_lnN50, count(*) AS n
FROM bfd_scratch.kx2_mt_score WHERE genus='Coemansia';
