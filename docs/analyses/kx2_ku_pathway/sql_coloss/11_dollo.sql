-- 11_dollo.sql
-- Purpose : Dollo parsimony on the first-pass backbone tree: each gene is assumed present at the root (if present in
--           any tip) and gained once; a loss is placed on the branch above every MAXIMAL clade whose tips are all absent.
--           Minimum number of independent losses = number of such clades. THIS IS PARSIMONY ON A FIRST-PASS FASTTREE
--           TOPOLOGY (63 genus representatives), NOT A LIKELIHOOD RECONSTRUCTION; Linderina and Basidiobolus placement
--           is unstable, and the tree is rooted at its printed trifurcation (node 0).
-- Inputs  : bfd_scratch.kx2_ku_tree_edges (10_tree_clades.py -> 10_tree_clades.sql; the run used the same VALUES on one
--           line), kx2_ku_matrix (set_name='backbone'), kx2_ku_genomes
-- Output  : bfd_scratch.kx2_ku_node_label (hand labels of key internal nodes)
--           bfd_scratch.kx2_ku_clade_tips (node, tip)  descendant tips of every node
--           bfd_scratch.kx2_ku_dollo (gene, node, n_tips, label, genera, kind, parent_n_present)
--           kind = 'internal' (>=2 tips) or 'terminal' (single tip; also the signature of a possible annotation false absence)
--           bfd_scratch.kx2_ku_dollo_sum (gene, n_loss_total, n_loss_internal, n_loss_terminal, n_tips_present, n_tips)
-- Thresholds: none (uses 'present' from 06).
-- Reports : number of independent Ku losses and the branches of every other gene's losses.
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_node_label AS
SELECT * FROM (VALUES
 (10,'Zoopagales'),(26,'Entomophthorales'),(25,'Entomophthorales + Kickxellomycota'),(31,'Kickxellomycota stem'),
 (32,'Dimargaritales'),(37,'Ramicandelaber + (Kickxellales + Harpellales)'),(39,'Kickxellales + Harpellales stem'),
 (40,'Harpellales'),(47,'Kickxellales s.l. (incl. Mycoemilia, Spiromyces)'),(48,'Kickxellales core (Linderina + Coemansia clade)'),
 (49,'Coemansia/Kickxella/Dipsacomyces/Martensiomyces'),(57,'Mycoemilia + Spiromyces'),(60,'Mucoromycota s.l.'),
 (7,'Basidiobolus + Zoopagales + Entomophthorales + Kickxellomycota'),(9,'Zoopagales + Entomophthorales + Kickxellomycota'),
 (3,'Blastocladiomycota'),(1,'all but chytrid outgroups')
) t(node, label);

CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_clade_tips AS
WITH RECURSIVE up(node, tip) AS (
  SELECT node, tip FROM bfd_scratch.kx2_ku_tree_edges WHERE tip IS NOT NULL
  UNION ALL
  SELECT e.parent, u.tip FROM up u JOIN bfd_scratch.kx2_ku_tree_edges e ON e.node=u.node WHERE e.parent >= 0
)
SELECT node, tip FROM up;

CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_dollo AS
WITH tipstate AS (
  SELECT gene, locustag AS tip, present FROM bfd_scratch.kx2_ku_matrix WHERE set_name='backbone'
),
nodestate AS (
  SELECT s.gene, c.node, count(*) AS n_tips, sum(s.present::INT) AS n_pres
  FROM bfd_scratch.kx2_ku_clade_tips c JOIN tipstate s ON s.tip=c.tip GROUP BY 1,2
),
loss AS (
  SELECT n.gene, n.node, n.n_tips, p.n_pres AS parent_n_present
  FROM nodestate n JOIN bfd_scratch.kx2_ku_tree_edges e ON e.node=n.node
  JOIN nodestate p ON p.gene=n.gene AND p.node=e.parent
  WHERE n.n_pres = 0 AND p.n_pres > 0
)
SELECT l.gene, l.node, l.n_tips, nl.label,
       (SELECT string_agg(g.genus, ', ' ORDER BY g.genus) FROM bfd_scratch.kx2_ku_clade_tips c JOIN bfd_scratch.kx2_ku_genomes g ON g.locustag=c.tip WHERE c.node=l.node) AS genera,
       CASE WHEN l.n_tips >= 2 THEN 'internal' ELSE 'terminal' END AS kind, l.parent_n_present
FROM loss l LEFT JOIN bfd_scratch.kx2_ku_node_label nl USING (node);

CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_dollo_sum AS
SELECT r.gene, r.conf, r.category,
       count(d.node) AS n_loss_total, count(d.node) FILTER (WHERE d.kind='internal') AS n_loss_internal,
       count(d.node) FILTER (WHERE d.kind='terminal') AS n_loss_terminal,
       (SELECT sum(present::INT) FROM bfd_scratch.kx2_ku_matrix m WHERE m.gene=r.gene AND m.set_name='backbone') AS n_tips_present,
       (SELECT count(*) FROM bfd_scratch.kx2_ku_matrix m WHERE m.gene=r.gene AND m.set_name='backbone') AS n_tips
FROM bfd_scratch.kx2_ku_rules r LEFT JOIN bfd_scratch.kx2_ku_dollo d ON d.gene=r.gene
GROUP BY r.gene, r.conf, r.category;
