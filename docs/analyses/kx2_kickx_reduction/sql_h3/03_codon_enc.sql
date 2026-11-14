-- 03_codon_enc.sql
-- Purpose : genome-wide codon metrics: ENC (Wright), GC3s, GC12, RSCU skew, expected ENC at GC3.
-- Inputs  : bfd_scratch.kx2_h3_set, bfd_scratch.kx2_h3_codmap, bfd.codon_freq (species_prefix, codon, frequency)
-- Output  : bfd_scratch.kx2_h3_cod
-- Formulas (genome-wide, pooled over all genes; NOT per gene):
--   p_c|aa = f_c / f_aa (within-amino-acid codon share); stops, Met, Trp excluded.
--   F_aa  ~= sum_c p^2   (Wright's homozygosity F=(n*sum p^2 -1)/(n-1); n unavailable, n is ~1e5-1e6 so
--            the correction is negligible and F ~= sum p^2).
--   Fbar_k = aa-frequency-weighted mean of F_aa over the degeneracy class k in {2,3,4,6}
--            (frequency weighting after Sun et al. 2013, citation to verify).
--   ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6     (Wright 1990, Gene 87:23-29; citation to verify).
--         Classes: 9 two-fold aa, Ile (3-fold), 5 four-fold aa, 3 six-fold aa (Leu, Ser, Arg treated as one family).
--   GC3s = share of sense codons (excl. Met/Trp/stops) with G/C at position 3.
--   ENC_exp = 2 + s + 29/(s^2+(1-s)^2), s=GC3s (Wright 1990). ENC departure = (ENC_exp-ENC)/ENC_exp (computed in 05).
--   rscu_skew = codon-frequency-weighted mean |RSCU-1|, RSCU = p * nf (multi-codon families).
--   gc_cds, gc12 = GC over all codons / over positions 1-2.
-- Limits  : codon_freq values are rounded to 4 decimals. Genome-level => confounded by GC/mutation bias.
-- Report  : GC3 (Kx core 67.2), ENC, RSCU skew.
CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_cod AS
WITH f AS (
  SELECT s.g, m.codon, m.aa, m.nf, m.n3, c.frequency AS fr
  FROM bfd_scratch.kx2_h3_set s JOIN bfd.codon_freq c ON c.species_prefix=s.g JOIN bfd_scratch.kx2_h3_codmap m USING (codon)
  WHERE m.aa NOT IN ('*','M','W')),
fa AS (SELECT *, sum(fr) OVER (PARTITION BY g, aa) AS faa, sum(fr) OVER (PARTITION BY g) AS ftot FROM f),
hom AS (
  SELECT g, aa, nf, max(faa) AS faa, sum(power(fr/faa,2)) AS F, max(ftot) AS ftot
  FROM fa WHERE faa>0 GROUP BY g, aa, nf),
cls AS (
  SELECT g, nf, sum(F*faa)/sum(faa) AS Fbar FROM hom GROUP BY g, nf),
enc AS (
  SELECT g, 2 + 9/max(Fbar) FILTER (WHERE nf=2) + 1/max(Fbar) FILTER (WHERE nf=3) + 5/max(Fbar) FILTER (WHERE nf=4) + 3/max(Fbar) FILTER (WHERE nf=6) AS enc,
         max(Fbar) FILTER (WHERE nf=2) F2, max(Fbar) FILTER (WHERE nf=4) F4
  FROM cls GROUP BY g),
gc3 AS (
  SELECT g, sum(fr) FILTER (WHERE n3 IN ('G','C'))/sum(fr) AS gc3s FROM fa GROUP BY g),
rs AS (
  SELECT g, sum(abs(fr/faa*nf - 1)*fr)/sum(fr) AS rscu_skew FROM fa WHERE faa>0 GROUP BY g),
gcc AS (
  SELECT c.species_prefix AS g,
    sum(c.frequency*((substr(codon,1,1) IN ('G','C'))::int + (substr(codon,2,1) IN ('G','C'))::int + (substr(codon,3,1) IN ('G','C'))::int))/3/sum(c.frequency) AS gc_cds,
    sum(c.frequency*((substr(codon,1,1) IN ('G','C'))::int + (substr(codon,2,1) IN ('G','C'))::int))/2/sum(c.frequency) AS gc12
  FROM bfd.codon_freq c WHERE c.species_prefix IN (SELECT g FROM bfd_scratch.kx2_h3_set) GROUP BY 1)
SELECT enc.g, enc.enc, enc.F2, enc.F4, gc3.gc3s, gcc.gc_cds, gcc.gc12, rs.rscu_skew,
  2 + gc3.gc3s + 29/(power(gc3.gc3s,2)+power(1-gc3.gc3s,2)) AS enc_exp
FROM enc JOIN gc3 USING (g) JOIN rs USING (g) JOIN gcc USING (g);
