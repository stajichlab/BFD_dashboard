-- 07_lostset_overlap.sql
-- Purpose: per-family lost sets for each reduced group, dispensability-stratified overlap of Kickxellales losses with
--          each other reduced group, and functional-category enrichment of losses.
-- Inputs : kx2_h4_gm (flt='all'), kx2_h4_pfcat
-- Outputs: bfd_scratch.kx2_h4_lostset (grp, pfam_acc, mb_m, mb_p, med, prev, lost), bfd_scratch.kx2_h4_overlap
-- Thresholds:
--   core families = core in BOTH baselines (prevalence >= 0.5 and median >= 1 in Mucorales and in Pezizomycotina), 3,915 families.
--   lost = group median 0 AND group prevalence < 0.3 (unfiltered genus set, 'all').
--   Null: for each pair (Kickxellales A, partner B) stratify families by k_other = number of the OTHER 8 reduced
--   groups (excluding A and B) that lost the family, capped at 4; expected overlap = sum over strata nA*nB/N
--   (hypergeometric within stratum); enrich = observed/expected; jaccard_null from the same expectation.
--   Done entirely in SQL (no Python). The "other" groups include related lineages, so strata are not independent.
-- Report numbers: lost fractions (Kx 11.2%, Harp 15.7%, Sacch 12.6%, Hanseniaspora 31.5%, Pneumocystis 25.1%,
--   Microsporidia 70.3%); overlap enrichment (Saccharomycetes 0.93, Pneumocystis 1.03, Harpellales 1.87, ...);
--   category enrichment table.

create or replace table kx2_h4_lostset as
with core as (select a.pfam_acc, a.med mb_m, b.med mb_p from kx2_h4_gm a join kx2_h4_gm b on a.pfam_acc=b.pfam_acc and b.flt='all' and b.grp='Pezizomycotina'
              where a.flt='all' and a.grp='Mucorales' and a.prev>=0.5 and a.med>=1 and b.prev>=0.5 and b.med>=1),
g as (select distinct grp from kx2_h4_gm where flt='all' and grp in ('Kickxellales','Harpellales','Dimargaritales','Zoopagales','Saccharomycetes','Saccharomycotina_other','Hanseniaspora','Pneumocystis','Taphrina','Microsporidia'))
select g.grp, c.pfam_acc, c.mb_m, c.mb_p, coalesce(m.med,0) med, coalesce(m.prev,0) prev, (coalesce(m.prev,0)<0.3 and coalesce(m.med,0)=0) lost
from core c cross join g left join kx2_h4_gm m on m.flt='all' and m.grp=g.grp and m.pfam_acc=c.pfam_acc;

create or replace table kx2_h4_overlap as
with pairs as (select 'Kickxellales' a, b from (select distinct grp b from kx2_h4_lostset where grp<>'Kickxellales')),
f as (select p.a, p.b, la.pfam_acc, la.lost la, lb.lost lb,
  (select sum(o.lost::int) from kx2_h4_lostset o where o.pfam_acc=la.pfam_acc and o.grp not in (p.a,p.b)) k_other
 from pairs p join kx2_h4_lostset la on la.grp=p.a join kx2_h4_lostset lb on lb.grp=p.b and lb.pfam_acc=la.pfam_acc),
st as (select a,b, least(k_other,4) s, count(*) N, sum(la::int) nA, sum(lb::int) nB, sum((la and lb)::int) obs from f group by 1,2,3)
select a, b, sum(obs) obs_overlap, round(sum(nA*nB/N::double),1) exp_stratified, round(sum(obs)/nullif(sum(nA*nB/N::double),0),2) enrich,
 sum(nA) lostA, sum(nB) lostB, round(sum(obs)/(sum(nA)+sum(nB)-sum(obs)),3) jaccard_obs,
 round(sum(nA*nB/N::double)/(sum(nA)+sum(nB)-sum(nA*nB/N::double)),3) jaccard_null
from st group by 1,2 order by 1,2;

-- Lost fraction per group (run, not stored):
-- select grp, count(*) n, sum(lost::int) n_lost, round(100*avg(lost::int),1) pct from kx2_h4_lostset group by 1 order by 1;

-- Functional-category enrichment (run, not stored): enrichment = (share of lost families in category) / (share of core families in category)
-- with l as (select s.grp, coalesce(c.cat,'unclassified') cat, s.lost from kx2_h4_lostset s left join kx2_h4_pfcat c using(pfam_acc)),
-- r as (select grp, cat, count(*) n, sum(lost::int) nl from l group by 1,2),
-- t as (select grp, sum(n) TN, sum(nl) TNL from r group by 1)
-- select r.cat, max(case when r.grp='Kickxellales' then round(100.0*nl/n,1) end) kx_lossrate, max(case when r.grp='Kickxellales' then round((nl/TNL)/(n/TN),2) end) kx_enr,
--  max(case when r.grp='Harpellales' then round((nl/TNL)/(n/TN),2) end) harp_enr,
--  max(case when r.grp='Saccharomycetes' then round(100.0*nl/n,1) end) sac_lossrate, max(case when r.grp='Saccharomycetes' then round((nl/TNL)/(n/TN),2) end) sac_enr,
--  max(case when r.grp='Hanseniaspora' then round((nl/TNL)/(n/TN),2) end) han_enr, max(case when r.grp='Pneumocystis' then round((nl/TNL)/(n/TN),2) end) pneumo_enr,
--  max(case when r.grp='Microsporidia' then round((nl/TNL)/(n/TN),2) end) micro_enr, max(case when r.grp='Kickxellales' then n end) n_core
-- from r join t using(grp) group by 1 order by n_core desc;
