-- 11b_dollo_keydomain.sql
-- Purpose : repeat the Dollo placement (11) using the relaxed single-key-domain presence calls (09), to see which loss
--           placements are robust to rule stringency.
-- Inputs  : bfd_scratch.kx2_ku_keymatrix (09), kx2_ku_clade_tips and kx2_ku_node_label (11), kx2_ku_genomes
-- Output  : bfd_scratch.kx2_ku_dollo_key (gene, node, n_tips, label, genera, kind)
-- Thresholds: i-evalue <= 1e-3 key domain; same maximal-absent-clade logic as 11 (parsimony on first-pass tree).
-- Reports : relaxed Ku loss count and placement.
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_dollo_key AS
WITH tipstate AS (
  SELECT gene, locustag AS tip, key_present AS present FROM bfd_scratch.kx2_ku_keymatrix WHERE set_name='backbone'
),
nodestate AS (
  SELECT s.gene, c.node, count(*) AS n_tips, sum(s.present::INT) AS n_pres
  FROM bfd_scratch.kx2_ku_clade_tips c JOIN tipstate s ON s.tip=c.tip GROUP BY 1,2
),
loss AS (
  SELECT n.gene, n.node, n.n_tips FROM nodestate n JOIN bfd_scratch.kx2_ku_tree_edges e ON e.node=n.node
  JOIN nodestate p ON p.gene=n.gene AND p.node=e.parent WHERE n.n_pres=0 AND p.n_pres>0
)
SELECT l.gene, l.node, l.n_tips, nl.label,
       (SELECT string_agg(g.genus, ', ' ORDER BY g.genus) FROM bfd_scratch.kx2_ku_clade_tips c JOIN bfd_scratch.kx2_ku_genomes g ON g.locustag=c.tip WHERE c.node=l.node) AS genera,
       CASE WHEN l.n_tips >= 2 THEN 'internal' ELSE 'terminal' END AS kind
FROM loss l LEFT JOIN bfd_scratch.kx2_ku_node_label nl USING (node);
