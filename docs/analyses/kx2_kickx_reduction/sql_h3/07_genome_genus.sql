-- 07_genome_genus.sql
-- Purpose : join all genome-level metrics; collapse to genus medians (all genomes, and BUSCO>=85 variant).
-- Inputs  : kx2_h3_set, kx2_h3_cod, kx2_h3_trna, kx2_h3_aa, kx2_h3_idp, kx2_h3_prot
-- Outputs : bfd_scratch.kx2_h3_genome, kx2_h3_genus (all reps), kx2_h3_genus_f (genomes with busco_complete>=85)
-- Notes   : enc_dep = (enc_exp-enc)/enc_exp. Genus value = median over the genus's species reps.
--           The BUSCO>=85 variant filters genomes BEFORE taking genus medians.
--           d_schpo/d_outgrp/cov_* come from kx2_h3_prot (707 SCHPO / 106 OUTGRP anchors).
-- Report  : the per-group medians table is built from these in 08.
CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_genome AS
SELECT s.*, c.enc, c.enc_exp, (c.enc_exp-c.enc)/c.enc_exp AS enc_dep, c.gc3s, c.gc_cds, c.gc12, c.rscu_skew, c.F2, c.F4,
 t.n_trna, t.n_antic, t.tai, t.rscu_supply_r, t.n_trna/s.mb AS trna_per_mb,
 a.fymink, a.garp, a.asn, a.lys, a.ile, a.cys, a.gln, a.lowcomp,
 i.idp_res_frac, i.idp_med_prot, i.n_prot, i.med_len,
 p1.med_delta AS d_schpo, p1.med_pid AS pid_schpo, p1.f_ge60 AS f60_schpo, p1.n_anchor_hit/707.0 AS cov_schpo,
 p2.med_delta AS d_outgrp, p2.n_anchor_hit/106.0 AS cov_outgrp
FROM bfd_scratch.kx2_h3_set s
LEFT JOIN bfd_scratch.kx2_h3_cod c USING (g) LEFT JOIN bfd_scratch.kx2_h3_trna t USING (g) LEFT JOIN bfd_scratch.kx2_h3_aa a USING (g) LEFT JOIN bfd_scratch.kx2_h3_idp i USING (g)
LEFT JOIN (SELECT * FROM bfd_scratch.kx2_h3_prot WHERE aset='SCHPO') p1 USING (g)
LEFT JOIN (SELECT * FROM bfd_scratch.kx2_h3_prot WHERE aset='OUTGRP') p2 USING (g);

CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_genus AS
SELECT hgrp, GENUS, count(*) n_sp, bool_or(busco_complete>=85) any85, median(busco_complete) busco,
 median(mb) mb, median(gc_pct) gc, median(gc3s) gc3, median(gc12) gc12, median(enc) enc, median(enc_dep) enc_dep, median(rscu_skew) rscu_skew,
 median(n_trna) n_trna, median(n_antic) n_antic, median(trna_per_mb) trna_per_mb, median(tai) tai, median(rscu_supply_r) rscu_supply_r,
 median(fymink) fymink, median(garp) garp, median(asn) asn, median(lys) lys, median(cys) cys, median(gln) gln, median(lowcomp) lowcomp, median(idp_res_frac) idp,
 median(d_schpo) d_schpo, median(f60_schpo) f60_schpo, median(cov_schpo) cov_schpo, median(d_outgrp) d_outgrp
FROM bfd_scratch.kx2_h3_genome GROUP BY 1,2;

CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_genus_f AS
SELECT hgrp, GENUS, count(*) n_sp, median(busco_complete) busco,
 median(mb) mb, median(gc_pct) gc, median(gc3s) gc3, median(gc12) gc12, median(enc) enc, median(enc_dep) enc_dep, median(rscu_skew) rscu_skew,
 median(n_trna) n_trna, median(n_antic) n_antic, median(trna_per_mb) trna_per_mb, median(tai) tai, median(rscu_supply_r) rscu_supply_r,
 median(fymink) fymink, median(garp) garp, median(asn) asn, median(lys) lys, median(cys) cys, median(gln) gln, median(lowcomp) lowcomp, median(idp_res_frac) idp,
 median(d_schpo) d_schpo, median(f60_schpo) f60_schpo, median(cov_schpo) cov_schpo, median(d_outgrp) d_outgrp
FROM bfd_scratch.kx2_h3_genome WHERE busco_complete>=85 GROUP BY 1,2;
