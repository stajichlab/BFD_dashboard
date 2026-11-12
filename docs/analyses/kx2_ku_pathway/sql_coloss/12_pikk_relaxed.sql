-- 12_pikk_relaxed.sql
-- Purpose : PIKK-family (PF00454 kinase domain) census per genome by SwissProt best-hit class, WITHOUT the coverage
--           filter used in 05. The strict ATM/Tel1 (TAN) and ATR/Mec1 (UME) Pfam rules call both kinases absent in all four
--           Harpellales, but their PIKK proteins hit ATM/rad3 over only 11-27% of the subject (truncated or fragmentary
--           gene models / diverged FAT region), so 'absent' may be 'truncated or diverged'. This table separates the cases.
-- Inputs  : bfd_scratch.kx2_ku_pfam (PF00454, <= 1e-3), kx2_ku_sp, kx2_ku_genomes
-- Output  : bfd_scratch.kx2_ku_pikk (locustag, pikk_class, n_prot, max_scov) ; bfd_scratch.kx2_ku_pikk_grp (group x class)
-- Thresholds: PF00454 i-evalue <= 1e-3; best-hit pident >= 25, -log10 e >= 10, NO subject-coverage minimum.
--           Classes (regex on 'gene_name | protein_name'): ATM/Tel1 (atm|tel1|mus-21|serine-protein kinase ATM),
--           ATR/Mec1 (rad3|mec1|ATR kinase), DNA-PKcs (PRKDC|DNA-dependent protein kinase), SMG1, TOR, TRRAP/Tra1, lipid-PI kinase (other).
-- Reports : relaxed ATM/ATR presence per group; DNA-PKcs name-based cross-check of the Pfam rule.
CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_pikk AS
SELECT f.locustag,
  CASE
    WHEN regexp_matches(coalesce(sp.gene_name,'')||' | '||coalesce(sp.protein_name,''), '(?i)\batm\b|tel1|mus-21|serine-protein kinase ATM') THEN 'ATM/Tel1'
    WHEN regexp_matches(coalesce(sp.gene_name,'')||' | '||coalesce(sp.protein_name,''), '(?i)\brad3\b|mec1|serine/threonine-protein kinase ATR|^atr ') THEN 'ATR/Mec1'
    WHEN regexp_matches(coalesce(sp.gene_name,'')||' | '||coalesce(sp.protein_name,''), '(?i)prkdc|DNA-dependent protein kinase') THEN 'DNA-PKcs'
    WHEN regexp_matches(coalesce(sp.gene_name,'')||' | '||coalesce(sp.protein_name,''), '(?i)smg-?1') THEN 'SMG1'
    WHEN regexp_matches(coalesce(sp.gene_name,'')||' | '||coalesce(sp.protein_name,''), '(?i)\btor[12]?\b|target of rapamycin|mTOR') THEN 'TOR'
    WHEN regexp_matches(coalesce(sp.gene_name,'')||' | '||coalesce(sp.protein_name,''), '(?i)tra1|trrap') THEN 'TRRAP/Tra1'
    ELSE 'other (PI3K/PI4K/unassigned)' END AS pikk_class,
  f.protein_id, sp.scov, sp.pident
FROM bfd_scratch.kx2_ku_pfam f LEFT JOIN bfd_scratch.kx2_ku_sp sp ON sp.protein_id=f.protein_id AND sp.pident>=25 AND sp.nlog10e>=10
WHERE f.acc='PF00454';

CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_pikk_grp AS
SELECT CASE WHEN g.genus='Coemansia' THEN 'Coemansia_strains' ELSE g.grp END AS grp2, k.pikk_class,
       count(DISTINCT g.locustag) AS n_genomes_with, count(*) AS n_prot,
       (SELECT count(*) FROM bfd_scratch.kx2_ku_genomes g2 WHERE (CASE WHEN g2.genus='Coemansia' THEN 'Coemansia_strains' ELSE g2.grp END) = (CASE WHEN g.genus='Coemansia' THEN 'Coemansia_strains' ELSE g.grp END) AND (g2.set_name='backbone' OR g2.genus='Coemansia')) AS n_genomes
FROM bfd_scratch.kx2_ku_pikk k JOIN bfd_scratch.kx2_ku_genomes g USING (locustag)
WHERE g.set_name='backbone' OR g.genus='Coemansia'
GROUP BY ALL;
