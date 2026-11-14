-- 05_calls.sql
-- Purpose : apply every rule in kx2_ku_rules to every protein of the 169 genomes -> protein-level calls.
-- Inputs  : bfd_scratch.kx2_ku_rules (04), kx2_ku_pfam (02), kx2_ku_sp (03)
-- Output  : bfd_scratch.kx2_ku_calls (gene, locustag, protein_id, mode, pfam_ok, name_ok, present, best-hit name/organism/pident/scov)
--           Only proteins with pfam_ok OR name_ok are kept.
-- Thresholds: Pfam rule: all req_all present, >=1 req_any present (if any), no excl, at the rule's ie_thr;
--           name rule: best hit pident>=25, scov>=40, -log10(e)>=10, name_re matches 'gene_name | protein_name', name_excl_re does not.
--           present = pfam_ok (mode pfam) | pfam_ok AND name_ok (both) | pfam_ok OR name_ok (either) | name_ok (name).
-- Reports : protein-level basis for every presence call.
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_calls AS
WITH pa AS (
  SELECT protein_id, locustag,
         list(acc) FILTER (WHERE min_ie <= 1e-5) AS strict_accs,
         list(acc) AS loose_accs
  FROM bfd_scratch.kx2_ku_pfam GROUP BY 1,2
),
univ AS (
  SELECT protein_id, locustag FROM pa
  UNION
  SELECT protein_id, locustag FROM bfd_scratch.kx2_ku_sp
),
feat AS (
  SELECT u.protein_id, u.locustag,
         coalesce(pa.strict_accs, []::VARCHAR[]) AS strict_accs, coalesce(pa.loose_accs, []::VARCHAR[]) AS loose_accs,
         sp.gene_name, sp.protein_name, sp.organism, sp.accession, sp.pident, sp.scov, sp.nlog10e
  FROM univ u LEFT JOIN pa USING (protein_id, locustag)
  LEFT JOIN bfd_scratch.kx2_ku_sp sp USING (protein_id)
),
ev AS (
  SELECT r.gene, r.category, r.mode, r.conf, r.is_ctrl, f.protein_id, f.locustag,
    (len(r.req_all) + len(r.req_any) > 0)
      AND list_has_all(CASE WHEN r.ie_thr <= 1e-5 THEN f.strict_accs ELSE f.loose_accs END, r.req_all)
      AND (len(r.req_any) = 0 OR list_has_any(CASE WHEN r.ie_thr <= 1e-5 THEN f.strict_accs ELSE f.loose_accs END, r.req_any))
      AND NOT list_has_any(CASE WHEN r.ie_thr <= 1e-5 THEN f.strict_accs ELSE f.loose_accs END, r.excl) AS pfam_ok,
    (r.name_re IS NOT NULL AND f.accession IS NOT NULL AND f.pident >= 25 AND f.scov >= 40 AND f.nlog10e >= 10
      AND regexp_matches(coalesce(f.gene_name,'') || ' | ' || f.protein_name, r.name_re)
      AND (r.name_excl_re IS NULL OR NOT regexp_matches(coalesce(f.gene_name,'') || ' | ' || f.protein_name, r.name_excl_re))) AS name_ok,
    f.gene_name AS sp_gene, f.protein_name AS sp_protein, f.organism AS sp_organism, f.accession AS sp_accession, f.pident, f.scov
  FROM bfd_scratch.kx2_ku_rules r CROSS JOIN feat f
)
SELECT gene, category, mode, conf, is_ctrl, protein_id, locustag, pfam_ok, name_ok,
       CASE mode WHEN 'pfam' THEN pfam_ok WHEN 'both' THEN (pfam_ok AND name_ok)
                 WHEN 'either' THEN (pfam_ok OR name_ok) ELSE name_ok END AS present,
       sp_gene, sp_protein, sp_organism, sp_accession, pident, scov
FROM ev WHERE pfam_ok OR name_ok;
