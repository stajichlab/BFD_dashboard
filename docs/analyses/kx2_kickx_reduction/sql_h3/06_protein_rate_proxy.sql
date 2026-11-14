-- 06_protein_rate_proxy.sql
-- Purpose : relative protein-evolution-rate proxy from SwissProt best-hit identity at fixed anchor accessions.
-- Inputs  : bfd.swissprot_best, bfd_scratch.kx2_h3_set
-- Outputs : bfd_scratch.kx2_h3_anchor (anchor accessions), bfd_scratch.kx2_h3_prot (per-genome scores)
-- Definitions / thresholds:
--   hit filter: qcov_x10>=800 and scov_x10>=800 (>=80% both sides; controls for length/partial models).
--   best per (genome, accession) = max pident.
--   Anchor sets (source organism of the SwissProt accession, org_code):
--     SCHPO  : S. pombe accessions.
--     OUTGRP : animal/plant/amoeba outgroup accessions with org_code in
--              HUMAN, MOUSE, DROME, CAEEL, ARATH, DANRE, XENLA, BOVIN, RAT, DICDI.
--   An accession is an anchor if it passes the hit filter in >=30% of kx2_h3_set genomes.
--   (report: 707 SCHPO and 106 OUTGRP anchors; the list is exactly the output of the first statement and is
--   stored in bfd_scratch.kx2_h3_anchor. The 106 were NOT hand picked.)
--   delta = pident - median pident of that accession across all kx2_h3_set genomes (reference is the
--   pooled set, dominated by Pezizomycotina).
--   f_ge60/f_ge70 = fraction of retained anchors at >=60/>=70% identity.
--   cov = n anchors retained / total anchors (707 or 106): coverage collapses in yeasts/Microsporidia
--   (best hit moves to a different accession) so delta is survivorship-biased there.
-- Caveats : SCHPO-based delta is biased by distance to S. pombe (not used for conclusions); OUTGRP delta
--           (d_outgrp) is the one reported. Identity is only a proxy for rate; GC content also affects it.
-- Report  : "outgroup identity delta" column, Cliff's d vs Mucorales/Sacch, GC3-matched dout_res.
CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_anchor AS
WITH h AS (
  SELECT b.species_prefix g, b.accession, b.org_code, b.pident_x10/10.0 AS pid, b.bitscore, b.q_len, b.s_len
  FROM bfd.swissprot_best b JOIN bfd_scratch.kx2_h3_set s ON s.g=b.species_prefix
  WHERE b.qcov_x10>=800 AND b.scov_x10>=800 AND b.org_code IN ('SCHPO','HUMAN','MOUSE','DROME','CAEEL','ARATH','DANRE','XENLA','BOVIN','RAT','DICDI')),
best AS (SELECT g, accession, any_value(org_code) org_code, max(pid) pid FROM h GROUP BY g, accession),
pres AS (SELECT accession, org_code, count(*) n, median(pid) med_pid FROM best GROUP BY 1,2)
SELECT accession, org_code, n, med_pid, CASE WHEN org_code='SCHPO' THEN 'SCHPO' ELSE 'OUTGRP' END AS aset FROM pres
WHERE n >= 0.30*(SELECT count(*) FROM bfd_scratch.kx2_h3_set);

CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_prot AS
WITH h AS (
  SELECT b.species_prefix g, b.accession, max(b.pident_x10/10.0) pid
  FROM bfd.swissprot_best b JOIN bfd_scratch.kx2_h3_set s ON s.g=b.species_prefix
  JOIN bfd_scratch.kx2_h3_anchor a ON a.accession=b.accession
  WHERE b.qcov_x10>=800 AND b.scov_x10>=800
  GROUP BY 1,2),
d AS (SELECT h.g, a.aset, h.pid, h.pid-a.med_pid AS delta, a.accession FROM h JOIN bfd_scratch.kx2_h3_anchor a USING (accession))
SELECT g, aset, count(*) n_anchor_hit, median(pid) med_pid, median(delta) med_delta, avg((pid>=60)::int) f_ge60, avg((pid>=70)::int) f_ge70
FROM d GROUP BY 1,2;
