-- 02_calls_matrix.sql
-- Purpose : apply every rule in kx2_mt_rules to every protein of the 169 genomes (protein-level calls), then count proteins per
--           gene and genome (presence matrix). Adapted from ../../sql_coloss/05_calls.sql and 06_matrix.sql (same logic).
-- Inputs  : bfd_scratch.kx2_mt_rules (01), bfd_scratch.kx2_mt_pfam (Pfam hits, domain_i_evalue <= 1e-3, extracted from bfd.pfam for the
--           accession list in 01), bfd_scratch.kx2_ku_sp (best SwissProt hit per protein for the 169 genomes: gene_name, protein_name,
--           pident, scov, nlog10e; read-only reuse), bfd_scratch.kx2_ku_genomes (169 genomes with group labels; read-only reuse).
-- Outputs : bfd_scratch.kx2_mt_calls (protein level), bfd_scratch.kx2_mt_matrix (gene x genome: n_present proteins).
-- Thresholds: Pfam rule: every req_all accession present, >= 1 req_any present (if any), no excl accession, at the rule's ie_thr
--           (strict 1e-5 uses hits with i-evalue <= 1e-5; loose 1e-3 uses all extracted hits); name rule: SwissProt best hit pident >= 25,
--           scov >= 40, -log10(e) >= 10, name_re matches 'gene_name | protein_name' and name_excl_re does not.
--           present = pfam_ok (mode pfam) | pfam_ok AND name_ok (both) | pfam_ok OR name_ok (either) | name_ok (name).
-- Reports : the annotation-level presence calls used in the README tables.
CREATE OR REPLACE TABLE bfd_scratch.kx2_mt_calls AS
WITH pa AS (
  SELECT protein_id, locustag,
         list(acc) FILTER (WHERE min_ie <= 1e-5) AS strict_accs,
         list(acc) AS loose_accs
  FROM bfd_scratch.kx2_mt_pfam GROUP BY 1, 2
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
  SELECT r.gene, r.grp AS fn, r.mode, r.conf, r.is_ctrl, r.style, f.protein_id, f.locustag,
    (len(r.req_all) + len(r.req_any) > 0)
      AND list_has_all(CASE WHEN r.ie_thr <= 1e-5 THEN f.strict_accs ELSE f.loose_accs END, r.req_all)
      AND (len(r.req_any) = 0 OR list_has_any(CASE WHEN r.ie_thr <= 1e-5 THEN f.strict_accs ELSE f.loose_accs END, r.req_any))
      AND NOT list_has_any(CASE WHEN r.ie_thr <= 1e-5 THEN f.strict_accs ELSE f.loose_accs END, r.excl) AS pfam_ok,
    (r.name_re IS NOT NULL AND f.accession IS NOT NULL AND f.pident >= 25 AND f.scov >= 40 AND f.nlog10e >= 10
      AND regexp_matches(coalesce(f.gene_name,'') || ' | ' || f.protein_name, r.name_re)
      AND (r.name_excl_re IS NULL OR NOT regexp_matches(coalesce(f.gene_name,'') || ' | ' || f.protein_name, r.name_excl_re))) AS name_ok,
    f.gene_name AS sp_gene, f.protein_name AS sp_protein, f.accession AS sp_accession, f.pident, f.scov
  FROM bfd_scratch.kx2_mt_rules r CROSS JOIN feat f
)
SELECT gene, fn, mode, conf, is_ctrl, style, protein_id, locustag, pfam_ok, name_ok,
       CASE mode WHEN 'pfam' THEN pfam_ok WHEN 'both' THEN (pfam_ok AND name_ok)
                 WHEN 'either' THEN (pfam_ok OR name_ok) ELSE name_ok END AS present,
       sp_gene, sp_protein, sp_accession, pident, scov
FROM ev WHERE pfam_ok OR name_ok;

-- gene x genome matrix: number of proteins that satisfy the rule (0 = not seen in the annotation)
CREATE OR REPLACE TABLE bfd_scratch.kx2_mt_matrix AS
SELECT r.gene, r.grp AS fn, r.conf, r.is_ctrl, r.style, g.locustag, g.genus, g.grp AS grp, g.ord, g.phylum, g.busco_complete, g.n50, g.n_proteins,
       coalesce(c.n, 0) AS n_present
FROM bfd_scratch.kx2_mt_rules r CROSS JOIN bfd_scratch.kx2_ku_genomes g
LEFT JOIN (SELECT gene, locustag, count(*) FILTER (WHERE present) AS n FROM bfd_scratch.kx2_mt_calls GROUP BY 1, 2) c
  ON c.gene = r.gene AND c.locustag = g.locustag;
