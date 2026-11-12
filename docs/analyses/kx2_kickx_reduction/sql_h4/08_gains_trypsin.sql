-- 08_gains_trypsin.sql
-- Purpose: Pfam families with copy number above both baselines in Kickxellales / Harpellales, targeted candidate
--          tables (AA9, laccase, tyrosinase, GH46, trypsin ...), and trypsin architecture check.
-- Inputs : kx2_h4_gm, kx2_h4_gp_all, kx2_h4_gen, kx2_h4_meta, kx2_h4_cazy_fam, bfd.pfam, kw_pfam_matrix
-- Outputs: bfd_scratch.kx2_h4_gains, kx2_h4_rep_gen, kx2_h4_trypsin_arch
-- Thresholds:
--   candidate: group median copies >= 2 (flt 'all' and 'b85'; Harpellales has no BUSCO>=85 genus)
--   gain: group median > 95th percentile (quantile_cont) of genus median copies in BOTH Mucorales and Pezizomycotina
--         genera (zero-padded), AND group median - max(Mucorales median, Pezizomycotina median) >= 2.
--   Trypsin architecture: best-BUSCO species rep per genus (rn=1) in Kickxellales, Harpellales, Mucorales; other Pfam
--   domains (domain_i_evalue<=1e-5) on proteins carrying PF00089; shown if on >=10% of trypsin proteins.
--   Caveats: copy counts are Pfam domain-bearing proteins, not orthologs; no phylogeny; annotation fragmentation can inflate.
-- Report numbers: top gained families (Trypsin 21.5 Kx / 72.5 Harp vs 2 Muc / 1 Pez; CAP 10; BAHD Transferase 17; GH46 6;
--   GH19 Harp 5.25), AA9/laccase/tyrosinase 'not gained', trypsin architecture (trypsin-only + PF13365 15%).

create or replace table kx2_h4_gains as
with ng as (select grp, count(*) n from kx2_h4_gen where grp in ('Mucorales','Pezizomycotina') group by 1),
base as (select p.grp, p.pfam_acc, list(p.med_copies) pl from kx2_h4_gp_all p join kx2_h4_gen g on g.genus=p.genus and g.grp=p.grp
         where p.grp in ('Mucorales','Pezizomycotina') group by 1,2),
bq as (select b.grp, b.pfam_acc, list_aggregate(list_concat(b.pl, list_resize([]::DOUBLE[], ng.n-len(b.pl),0)),'quantile_cont',0.95) q95,
        list_median(list_concat(b.pl, list_resize([]::DOUBLE[], ng.n-len(b.pl),0))) med from base b join ng using(grp)),
cand as (select flt, grp, pfam_acc, med from kx2_h4_gm where grp in ('Kickxellales','Harpellales') and med>=2),
j as (select c.flt, c.grp, c.pfam_acc, c.med kmed,
  coalesce(m.med,0) muc_med, coalesce(m.q95,0) muc_q95, coalesce(z.med,0) pez_med, coalesce(z.q95,0) pez_q95
 from cand c left join bq m on m.grp='Mucorales' and m.pfam_acc=c.pfam_acc left join bq z on z.grp='Pezizomycotina' and z.pfam_acc=c.pfam_acc)
select j.*, pn.pfam_id,
 j.kmed - greatest(muc_med,pez_med) delta,
 (j.kmed+0.5)/(greatest(muc_med,pez_med)+0.5) ratio
from j left join (select distinct pfam_acc, pfam_id from kw_pfam_matrix) pn using(pfam_acc)
where j.kmed > greatest(muc_q95, pez_q95) and j.kmed - greatest(muc_med,pez_med) >= 2;

-- Ranked gains listing (run, not stored):
-- select pfam_acc, pfam_id, max(case when grp='Kickxellales' and flt='all' then kmed end) kx_all, max(case when grp='Kickxellales' and flt='b85' then kmed end) kx_b85,
--  max(case when grp='Harpellales' then kmed end) harp, max(muc_med) muc, max(pez_med) pez, max(greatest(muc_q95,pez_q95)) q95max
-- from kx2_h4_gains group by 1,2
-- order by coalesce(max(case when grp='Kickxellales' and flt='all' then delta end),0)+coalesce(max(case when grp='Harpellales' then delta end),0) desc, 1 limit 60;

-- Targeted candidates (run, not stored): group medians for named families
-- with t as (select * from (values ('PF00089','Trypsin'),('PF00394','Cu-oxidase (laccase)'),('PF07731','Cu-oxidase_2'),('PF07732','Cu-oxidase_3'),('PF00264','Tyrosinase'),('PF01374','GH46'),('PF00188','CAP'),('PF02458','BAHD Transferase'),('PF02128','Peptidase_M36'),('PF00082','Peptidase_S8')) x(a,n))
-- select t.a, t.n, g.grp, g.flt, round(g.med,2) med, round(g.prev,2) prev, g.n_gen from t join kx2_h4_gm g on g.pfam_acc=t.a
-- where g.flt='all' and g.grp in ('Kickxellales','Harpellales','Dimargaritales','Zoopagales','Mucorales','Pezizomycotina','Saccharomycetes') order by t.a, g.grp;
-- CAZy families (positives-only medians; zero genera absent):
-- with g as (select grp, fam, count(*) n_gen, median(n) med_n from (select m.grp, m.genus, f.fam, median(f.n) n from kx2_h4_cazy_fam f join kx2_h4_meta m using(genome)
--   where f.fam in ('AA9','GH46','GH19','AA1','AA5','AA3','GH81','CBM96','GH18') group by 1,2,3) group by 1,2)
-- select fam, grp, med_n from g where grp in ('Kickxellales','Harpellales','Zoopagales','Mucorales','Pezizomycotina','Saccharomycetes') order by fam, grp;

create or replace table kx2_h4_rep_gen as
select genome, genus, grp, row_number() over (partition by genus, grp order by busco desc, n50 desc) rn
from kx2_h4_meta where grp in ('Kickxellales','Harpellales','Mucorales','Pezizomycotina');

create or replace table kx2_h4_trypsin_arch as
with tp as (select split_part(protein_id,'_',1) genome, protein_id from bfd.pfam
            where pfam_acc like 'PF00089%' and domain_i_evalue<=1e-5
              and split_part(protein_id,'_',1) in (select genome from kx2_h4_rep_gen where rn=1 and grp in ('Kickxellales','Harpellales','Mucorales'))),
o as (select p.protein_id, regexp_extract(p.pfam_acc,'^(PF[0-9]+)',1) acc from bfd.pfam p join tp using(protein_id) where p.domain_i_evalue<=1e-5)
select m.grp, o.acc, count(distinct o.protein_id) n_prot,
 (select count(distinct protein_id) from tp t2 join kx2_h4_meta m2 on m2.genome=t2.genome where m2.grp=m.grp) n_tryp_prot
from o join tp using(protein_id) join kx2_h4_meta m on m.genome=tp.genome group by 1,2;

-- Trypsin architecture listing (run, not stored):
-- select grp, acc, n_prot, n_tryp_prot, round(100.0*n_prot/n_tryp_prot) pct from kx2_h4_trypsin_arch where n_prot>=0.1*n_tryp_prot order by grp, n_prot desc;
