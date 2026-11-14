-- 05_aa_idp.sql
-- Purpose : amino-acid composition indicators and intrinsic-disorder fraction per genome.
-- Inputs  : bfd.aa_freq (species_prefix, amino_acid one-letter, frequency), bfd.idp_summary, bfd_scratch.kx2_h3_set
-- Outputs : bfd_scratch.kx2_h3_aa, bfd_scratch.kx2_h3_idp
-- Notes   : aa frequencies renormalised after dropping '*' and 'X'.
--           fymink = F+Y+M+I+N+K (AT-codon-rich aa); garp = G+A+R+P; lowcomp = S,T,P,Q,N,G.
--           idp_res_frac = sum(IDP_residues)/sum(length) (residue-weighted, per genome).
-- Report  : FYMINK %, Asn/Lys/Cys/Gln %, IDP % in the group table; lowcomp in residual table.
CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_aa AS
WITH a AS (SELECT species_prefix g, amino_acid aa, frequency fr FROM bfd.aa_freq WHERE species_prefix IN (SELECT g FROM bfd_scratch.kx2_h3_set) AND amino_acid NOT IN ('*','X')),
n AS (SELECT g, aa, fr/sum(fr) OVER (PARTITION BY g) f FROM a)
SELECT g,
 sum(f) FILTER (WHERE aa IN ('F','Y','M','I','N','K')) AS fymink,
 sum(f) FILTER (WHERE aa IN ('G','A','R','P')) AS garp,
 sum(f) FILTER (WHERE aa='N') AS asn, sum(f) FILTER (WHERE aa='K') AS lys, sum(f) FILTER (WHERE aa='I') AS ile,
 sum(f) FILTER (WHERE aa='C') AS cys, sum(f) FILTER (WHERE aa='Q') AS gln,
 sum(f) FILTER (WHERE aa IN ('N','Q','S','T')) AS polar_ns, sum(f) FILTER (WHERE aa IN ('S','T','P','Q','N','G')) AS lowcomp
FROM n GROUP BY g;

CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_idp AS
SELECT species_prefix g, sum(IDP_residues)::double/sum(length) AS idp_res_frac, median(IDP_fraction) AS idp_med_prot, count(*) n_prot, median(length) med_len
FROM bfd.idp_summary WHERE species_prefix IN (SELECT g FROM bfd_scratch.kx2_h3_set) GROUP BY 1;
