-- 01_meta_nprot_gen.sql
-- Purpose: per-genome proteome size, lineage group assignment, and genus-level table.
-- Inputs : bfd.gene_proteins, bfd_scratch.kw_rep_genomes
-- Outputs: bfd_scratch.kx2_h4_nprot (genome, n_genes)
--          bfd_scratch.kx2_h4_meta  (one row per species-rep genome, with grp)
--          bfd_scratch.kx2_h4_gen   (genus x grp: n_sp, median n_genes, median BUSCO, median N50)
-- Thresholds/definitions:
--   genome key = LOCUSTAG = split_part(protein_id,'_',1); n_genes = count(distinct gene_id).
--   BUSCO>=85 variant is applied downstream (kx2_h4_gen.busco = genus median of complete_pct).
--   Groups: Kickxellales, Harpellales, Dimargaritales, Zoopagales, Mucorales, Microsporidia, Pneumocystis,
--   Taphrina, Hanseniaspora, Saccharomycetes (class, excl. Hanseniaspora), Saccharomycotina_other,
--   Pezizomycotina (listed classes), Kickxellomycota_other, Other_fungi.
-- Report numbers: group sizes/BUSCO/genes tables (e.g. Kickxellales 153 species, 7 genera, BUSCO 68.5,
--   6171 genes; only 3 core Kickxellales genera and 0 Harpellales/Hanseniaspora/Microsporidia pass BUSCO>=85).
-- Run with the query_rw MCP tool in database bfd_scratch (statements run one per call).

create or replace table kx2_h4_nprot as
select split_part(protein_id,'_',1) genome, count(distinct gene_id) n_genes
from bfd.gene_proteins
where split_part(protein_id,'_',1) in (select LOCUSTAG from kw_rep_genomes)
group by 1;

create or replace table kx2_h4_meta as
select r.LOCUSTAG genome, r.GENUS genus, r.CLASS, r.ORDER_, r.PHYLUM, r.complete_pct busco, r.N50_bp n50, r.total_length_bp, n.n_genes,
 case
  when r.ORDER_='Kickxellales' then 'Kickxellales'
  when r.ORDER_='Harpellales' then 'Harpellales'
  when r.ORDER_='Dimargaritales' then 'Dimargaritales'
  when r.ORDER_='Zoopagales' then 'Zoopagales'
  when r.ORDER_='Mucorales' then 'Mucorales'
  when r.PHYLUM='Microsporidia' then 'Microsporidia'
  when r.GENUS='Pneumocystis' then 'Pneumocystis'
  when r.GENUS='Taphrina' then 'Taphrina'
  when r.GENUS='Hanseniaspora' then 'Hanseniaspora'
  when r.CLASS='Saccharomycetes' then 'Saccharomycetes'
  when r.CLASS in ('Pichiomycetes','Dipodascomycetes','Lipomycetes','Trigonopsidomycetes','Alloascoideomycetes') then 'Saccharomycotina_other'
  when r.CLASS in ('Sordariomycetes','Eurotiomycetes','Dothideomycetes','Leotiomycetes','Pezizomycetes','Orbiliomycetes','Lecanoromycetes','Arthoniomycetes','Xylonomycetes','Lichinomycetes','Coniocybomycetes','Geoglossomycetes') then 'Pezizomycotina'
  when r.PHYLUM='Kickxellomycota' then 'Kickxellomycota_other'
  else 'Other_fungi' end grp
from kw_rep_genomes r join kx2_h4_nprot n on n.genome=r.LOCUSTAG;

create or replace table kx2_h4_gen as
select genus, grp, count(*) n_sp, median(n_genes) n_genes, median(busco) busco, median(n50) n50
from kx2_h4_meta group by genus, grp;

-- check queries (run, not stored):
-- select grp, count(*) n, count(distinct genus) g, round(median(busco),1) b, round(median(n_genes)) ng from kx2_h4_meta group by 1 order by 1;
-- select grp, count(*) g, sum((busco>=85)::int) g85, sum((busco>=75)::int) g75, sum((busco>=65)::int) g65, round(median(busco),1) b from kx2_h4_gen group by 1 order by 1;
