-- Purpose: SwissProt-accession family shares and mean lengths per group (genus-weighted: each genus weight 1/ngenera, genomes equal within genus).
-- Inputs: kx2_h1_genomes, kx2_h1_prot, bfd.swissprot_best. Output: kx2_h1_fam_gen.
-- Filters: q genomes; groups core_Kickxellales, Mucorales, yeast_other, Microsporidia, Taphrinales;
--          hit neg_log10_evalue>=10 AND pident>=30%; unassigned proteins -> family 'NONE'.
-- Caveat: genomes with no rows in kx2_h1_prot are not represented, so Microsporidia shares sum to 0.913, not 1.
CREATE OR REPLACE TABLE kx2_h1_fam_gen AS
WITH GG AS (SELECT LOCUSTAG, grp, GENUS FROM kx2_h1_genomes WHERE q AND grp IN ('core_Kickxellales','Mucorales','yeast_other','Microsporidia','Taphrinales')),
W AS (SELECT LOCUSTAG, grp, GENUS, 1.0/count(*) OVER (PARTITION BY grp,GENUS) / count(distinct GENUS) OVER (PARTITION BY grp) AS w FROM GG),
PR AS (SELECT p.lt AS LOCUSTAG, p.protein_id, p.length, coalesce(s.accession,'NONE') AS fam FROM kx2_h1_prot p JOIN W ON W.LOCUSTAG=p.lt LEFT JOIN bfd.swissprot_best s ON s.protein_id=p.protein_id AND s.species_prefix=p.lt AND s.neg_log10_evalue>=10 AND s.pident_x10>=300),
GF AS (SELECT LOCUSTAG, fam, count(*) AS n, sum(length) AS sl FROM PR GROUP BY 1,2),
GT AS (SELECT LOCUSTAG, sum(n) AS tot FROM GF GROUP BY 1)
SELECT W.grp, GF.fam, sum(W.w*GF.n::double/GT.tot) AS fshare, sum(W.w*GF.sl/GT.tot)/sum(W.w*GF.n::double/GT.tot) AS fL
FROM GF JOIN GT USING(LOCUSTAG) JOIN W USING(LOCUSTAG) GROUP BY 1,2;
SELECT grp, count(*) fams, round(sum(fshare),3) ssum, round(sum(fshare*fL),1) meanL, round(max(fshare) FILTER (WHERE fam='NONE'),3) none_share FROM kx2_h1_fam_gen GROUP BY 1;
