-- 15_verify.sql
-- Purpose : reproducibility check. Re-run the calls (05), presence matrix (06) and co-loss summary (13) from their SQL into
--           temporary bfd_scratch.kx2_ku_chk_* tables (same inputs: rules, kx2_ku_pfam, kx2_ku_sp, kx2_ku_genomes, kx2_ku_keymatrix,
--           kx2_ku_clade_tips), then compare to the reported tables with EXCEPT in both directions, then DROP the chk tables.
-- Inputs  : as above (inputs 01-04, 09, 11 are NOT regenerated here; see README for what was not re-run)
-- Output  : bfd_scratch.kx2_ku_chk_calls, kx2_ku_chk_matrix, kx2_ku_chk_coloss (dropped at the end); result rows = diff counts
-- Expected: all diff counts = 0.
-- Statement A (chk_calls) is the body of 05_calls.sql with the target renamed; B (chk_matrix) is 06 reading chk_calls;
-- C (chk_coloss) is 13 reading chk_matrix.
-- (The statements are not repeated here to avoid divergence: run 05/06/13 with the table names substituted by
--  sed 's/kx2_ku_calls/kx2_ku_chk_calls/g; s/kx2_ku_matrix/kx2_ku_chk_matrix/g; s/kx2_ku_coloss/kx2_ku_chk_coloss/g'
--  with the 06 and 13 reads re-pointed to the chk tables; the comparison queries are below.)
SELECT 'calls only in reported' AS what, count(*) AS n FROM (SELECT gene, protein_id, pfam_ok, name_ok, present FROM bfd_scratch.kx2_ku_calls EXCEPT SELECT gene, protein_id, pfam_ok, name_ok, present FROM bfd_scratch.kx2_ku_chk_calls)
UNION ALL SELECT 'calls only in chk', count(*) FROM (SELECT gene, protein_id, pfam_ok, name_ok, present FROM bfd_scratch.kx2_ku_chk_calls EXCEPT SELECT gene, protein_id, pfam_ok, name_ok, present FROM bfd_scratch.kx2_ku_calls)
UNION ALL SELECT 'matrix only in reported', count(*) FROM (SELECT gene, locustag, present, n_present FROM bfd_scratch.kx2_ku_matrix EXCEPT SELECT gene, locustag, present, n_present FROM bfd_scratch.kx2_ku_chk_matrix)
UNION ALL SELECT 'matrix only in chk', count(*) FROM (SELECT gene, locustag, present, n_present FROM bfd_scratch.kx2_ku_chk_matrix EXCEPT SELECT gene, locustag, present, n_present FROM bfd_scratch.kx2_ku_matrix)
UNION ALL SELECT 'coloss only in reported', count(*) FROM (SELECT gene, fA, fB, fC, fD, fCoem, class FROM bfd_scratch.kx2_ku_coloss EXCEPT SELECT gene, fA, fB, fC, fD, fCoem, class FROM bfd_scratch.kx2_ku_chk_coloss)
UNION ALL SELECT 'coloss only in chk', count(*) FROM (SELECT gene, fA, fB, fC, fD, fCoem, class FROM bfd_scratch.kx2_ku_chk_coloss EXCEPT SELECT gene, fA, fB, fC, fD, fCoem, class FROM bfd_scratch.kx2_ku_coloss);
-- then:
-- DROP TABLE bfd_scratch.kx2_ku_chk_calls; DROP TABLE bfd_scratch.kx2_ku_chk_matrix; DROP TABLE bfd_scratch.kx2_ku_chk_coloss;
