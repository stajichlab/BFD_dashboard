-- 04_trna_tai.sql
-- Purpose : tRNA gene count, number of distinct tRNA species, genome-level tAI-like match, and
--           correlation between codon usage (RSCU) and tRNA supply within amino-acid families.
-- Inputs  : bfd.gene_trnas (gene_id, amino_acid, codon), bfd_scratch.kx2_h3_set, kx2_h3_codmap, bfd.codon_freq
-- Output  : bfd_scratch.kx2_h3_trna
-- Assumptions:
--   gene_trnas.codon is the COGNATE codon (Watson-Crick partner of the anticodon): verified by amino-acid/codon
--   consistency (e.g. Ile rows carry ATT/ATC/ATA). NULL codon rows, SeC, iMet, Sup are dropped.
--   n_antic = number of distinct cognate codons with >=1 tRNA gene (distinct anticodon species).
--   Wobble weights s (dos Reis, Savva & Wernisch 2004, NAR 32:5036; values and eukaryote A34->I34 rule:
--   citation to verify): cognate codon ending C (G34) reads C (0) and U (0.41); ending T (A34, treated as I34)
--   reads U (0), C (0.28), A (0.9999); ending A (U34) reads A (0) and G (0.68); ending G (C34) reads G (0).
--   Wobble is allowed only inside the same 2-nt box AND same amino acid.
--   W_codon = sum_j (1-s_ij) * tGCN_j ; w_i = W_i / max W (per genome); w_i=0 replaced by mean of non-zero w
--   (dos Reis convention). tai = exp( sum_c f_c ln w_c / sum_c f_c ) over sense codons excl. Met/Trp,
--   frequency-weighted geometric mean = pooled genome-level (not per gene) tAI.
--   rscu_supply_r = Pearson r across codons (multi-codon families) of RSCU (p*k) vs relative supply (q*k, q=w/sum w in aa).
-- Report  : tRNA gene counts, anticodons, tAI, supply-match r (Kx core 0.47 vs yeast 0.61).
CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_trna AS
WITH t AS (
  SELECT split_part(gene_id,'_',1) AS g, codon FROM bfd.gene_trnas
  WHERE codon IS NOT NULL AND amino_acid NOT IN ('SeC','iMet','Sup') AND split_part(gene_id,'_',1) IN (SELECT g FROM bfd_scratch.kx2_h3_set)),
tc AS (SELECT g, codon AS tcod, count(*) n FROM t GROUP BY 1,2),
w AS (SELECT * FROM (VALUES ('C','C',0.0),('C','T',0.41),('T','T',0.0),('T','C',0.28),('T','A',0.9999),('A','A',0.0),('A','G',0.68),('G','G',0.0)) v(tn3, cn3, s)),
dec AS (
  SELECT tc.g, m2.codon, sum((1-w.s)*tc.n) AS W
  FROM tc JOIN bfd_scratch.kx2_h3_codmap m1 ON m1.codon=tc.tcod
  JOIN w ON w.tn3=m1.n3
  JOIN bfd_scratch.kx2_h3_codmap m2 ON m2.box=m1.box AND m2.n3=w.cn3 AND m2.aa=m1.aa
  GROUP BY 1,2),
allc AS (SELECT s.g, m.codon, m.aa, m.n3, coalesce(d.W,0) AS W FROM bfd_scratch.kx2_h3_set s CROSS JOIN bfd_scratch.kx2_h3_codmap m LEFT JOIN dec d ON d.g=s.g AND d.codon=m.codon WHERE m.aa<>'*'),
nrm AS (SELECT *, W/nullif(max(W) OVER (PARTITION BY g),0) AS w0 FROM allc),
nrm2 AS (SELECT *, CASE WHEN w0>0 THEN w0 ELSE avg(w0) FILTER (WHERE w0>0) OVER (PARTITION BY g) END AS wi FROM nrm),
cf AS (SELECT n.g, n.codon, n.aa, n.wi, c.frequency fr,
   c.frequency/sum(c.frequency) OVER (PARTITION BY n.g, n.aa) AS p,
   wi/sum(wi) OVER (PARTITION BY n.g, n.aa) AS q,
   count(*) OVER (PARTITION BY n.g, n.aa) AS k
   FROM nrm2 n JOIN bfd.codon_freq c ON c.species_prefix=n.g AND c.codon=n.codon WHERE n.aa NOT IN ('M','W')),
stat AS (
  SELECT g, exp(sum(fr*ln(wi))/sum(fr)) AS tai, corr(p*k, q*k) AS rscu_supply_r
  FROM cf WHERE wi>0 GROUP BY g),
nt AS (SELECT g, count(*) n_trna, count(distinct codon) n_anticodons FROM t GROUP BY g)
SELECT s.g, coalesce(nt.n_trna,0) n_trna, coalesce(nt.n_anticodons,0) n_antic, stat.tai, stat.rscu_supply_r
FROM bfd_scratch.kx2_h3_set s LEFT JOIN nt USING(g) LEFT JOIN stat USING (g);
