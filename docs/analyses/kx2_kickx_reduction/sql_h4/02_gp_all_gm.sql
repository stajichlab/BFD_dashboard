-- 02_gp_all_gm.sql
-- Purpose: genus median Pfam copy number (zero-padded over species reps), then group medians across genera,
--          for the unfiltered ('all') and BUSCO>=85 ('b85') genus sets.
-- Inputs : bfd_scratch.kw_pfam_matrix (n_proteins, gate domain_i_evalue<=1e-5), kx2_h4_meta, kx2_h4_gen
-- Outputs: bfd_scratch.kx2_h4_gp_all (genus, grp, pfam_acc, n_sp, n_pos, med_copies)
--          bfd_scratch.kx2_h4_gm     (flt, grp, pfam_acc, n_gen, prev, med)
-- Definitions:
--   med_copies = median over species of n_proteins, absent species counted as 0 (list_resize padding).
--   prev = share of genera in the group with genus median copies >= 1.
--   med  = median over genera of med_copies, absent genera counted as 0.
--   flt 'b85' keeps genera whose median BUSCO complete_pct >= 85.
-- Report numbers: basis for every Pfam table (loss split, lostset, gains).

create or replace table kx2_h4_gp_all as
with ng as (select genus, grp, count(*) n_sp, median(n_genes) n_genes, median(busco) busco from kx2_h4_meta group by 1,2),
pos as (select m.genus, m.grp, p.pfam_acc, list(p.n_proteins) pl
        from kx2_h4_meta m join kw_pfam_matrix p on p.genome=m.genome group by 1,2,3)
select pos.genus, pos.grp, pos.pfam_acc, ng.n_sp, len(pos.pl) n_pos,
 list_median(list_concat(pos.pl, list_resize([]::BIGINT[], ng.n_sp-len(pos.pl), 0))) med_copies
from pos join ng using(genus,grp);

create or replace table kx2_h4_gm as
with f as (select 'all' flt, genus, grp from kx2_h4_gen
           union all select 'b85', genus, grp from kx2_h4_gen where busco>=85),
ng as (select flt, grp, count(*) n_gen from f group by 1,2),
pos as (select f.flt, f.grp, p.pfam_acc, list(p.med_copies) pl, sum((p.med_copies>=1)::int) n_pres
        from f join kx2_h4_gp_all p on p.genus=f.genus and p.grp=f.grp group by 1,2,3)
select pos.flt, pos.grp, pos.pfam_acc, ng.n_gen, pos.n_pres/ng.n_gen::double prev,
 list_median(list_concat(pos.pl, list_resize([]::DOUBLE[], ng.n_gen-len(pos.pl), 0))) med
from pos join ng using(flt,grp);
