-- Purpose: per-protein table (length, partial flag, Pfam architecture, domain residues, disorder) for the chosen genomes.
-- Inputs: bfd_scratch.kx2_h1_genomes; bfd.gene_proteins, bfd.pfam, bfd.gene_transcripts, bfd.idp_summary
-- Output: bfd_scratch.kx2_h1_prot (~5.0M rows)
-- Filters: Pfam domain_i_evalue <= 1e-5 (accession version stripped). dom_res = union of env intervals (overlap merged via running max end).
--          gene_proteins ids carry a -T1 suffix; idp_summary ids do not (joined with regexp_replace).
--          IDP coverage is partial (2.56M of 5.0M proteins).
-- Report: all protein-level metrics downstream.
CREATE OR REPLACE TABLE kx2_h1_prot AS
WITH G AS (SELECT LOCUSTAG FROM kx2_h1_genomes),
P AS (SELECT p.protein_id, split_part(p.protein_id,'_',1) AS lt, p.gene_id, p.transcript_id, p.length FROM bfd.gene_proteins p JOIN G ON G.LOCUSTAG=split_part(p.protein_id,'_',1)),
D0 AS (SELECT f.protein_id, split_part(f.pfam_acc,'.',1) AS acc, f.env_from, f.env_to FROM bfd.pfam f JOIN G ON G.LOCUSTAG=split_part(f.protein_id,'_',1) WHERE f.domain_i_evalue<=1e-5),
D1 AS (SELECT *, max(env_to) OVER (PARTITION BY protein_id ORDER BY env_from, env_to ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS pm FROM D0),
DA AS (SELECT protein_id, count(*) n_dom, count(distinct acc) n_dist, string_agg(acc,';' ORDER BY env_from,env_to) arch,
  sum(greatest(0, env_to - greatest(env_from-1, coalesce(pm,0)))) cov FROM D1 GROUP BY 1)
SELECT P.*, t.is_partial, t.CDS_length, coalesce(DA.n_dom,0) n_dom, coalesce(DA.n_dist,0) n_dist, DA.arch, coalesce(DA.cov,0) dom_res, i.IDP_residues
FROM P LEFT JOIN DA USING(protein_id) LEFT JOIN bfd.gene_transcripts t ON t.transcript_id=P.transcript_id
LEFT JOIN bfd.idp_summary i ON i.protein_id=P.protein_id;
-- NOTE: the IDP_residues column from the join above is all NULL (id suffix mismatch); it is unused. The working column is idp_res, added next.
ALTER TABLE kx2_h1_prot ADD COLUMN idp_res BIGINT;
UPDATE kx2_h1_prot SET idp_res=i.IDP_residues FROM bfd.idp_summary i WHERE i.protein_id=regexp_replace(kx2_h1_prot.protein_id,'-T[0-9]+$','') AND i.species_prefix=kx2_h1_prot.lt;
SELECT count(*), count(idp_res) FROM kx2_h1_prot;
