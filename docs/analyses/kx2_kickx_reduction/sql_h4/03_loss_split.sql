-- 03_loss_split.sql
-- Purpose: split lost gene copies into complete loss vs contraction relative to an unreduced baseline,
--          by baseline copy class; plus observed vs uniform-thinning expectation of complete-loss families.
-- Inputs : kx2_h4_gm
-- Outputs: bfd_scratch.kx2_h4_loss (flt, grp, bgrp, bclass, n_fam, base_copies, loss_copies, complete_copies,
--          contract_copies, n_complete, n_contract, n_kept, gain_copies)
-- Thresholds:
--   core baseline family: baseline group prevalence >= 0.5 and baseline median copies >= 1
--   baseline groups: Mucorales and Pezizomycotina (both run for every group; Mucorales is the sister baseline for
--   Kickxellomycota/Zoopagales, Pezizomycotina for yeasts/Pneumocystis/Taphrina)
--   bclass: 'single' b<=1; 'multi2-3' b<4; 'multi4+' b>=4 (b = baseline median copies)
--   complete loss: group median c = 0; contraction: 0 < c < b; loss copies = max(b-c,0)
--   BUSCO>=85 variant: flt='b85' (only groups with >=1 passing genus appear)
-- Report numbers: loss split table (Kickxellales 41% copies lost vs Pezizomycotina, 59% contraction / 41% complete;
--   46% / 81% / 19% vs Mucorales), single-copy family loss, obs/exp complete loss (Kx 0.89, 0.62 vs Muc).

create or replace table kx2_h4_loss as
with bl as (select flt, grp bgrp, pfam_acc, med b, prev bprev from kx2_h4_gm
            where grp in ('Mucorales','Pezizomycotina') and prev>=0.5 and med>=1),
tg as (select flt, grp, n_gen from kx2_h4_gm where grp not in ('Other_fungi') group by all),
x as (select b.flt, g.grp, b.bgrp, b.pfam_acc, b.b, coalesce(m.med,0) c,
  case when b.b<=1 then 'single' when b.b<4 then 'multi2-3' else 'multi4+' end bclass
  from bl b cross join (select distinct flt, grp from tg) g
  left join kx2_h4_gm m on m.flt=b.flt and m.grp=g.grp and m.pfam_acc=b.pfam_acc
  where g.flt=b.flt and g.grp<>b.bgrp)
select flt, grp, bgrp, bclass, count(*) n_fam, sum(b) base_copies, sum(greatest(b-c,0)) loss_copies,
 sum(case when c=0 then b else 0 end) complete_copies, sum(case when c>0 and c<b then b-c else 0 end) contract_copies,
 sum((c=0)::int) n_complete, sum((c>0 and c<b)::int) n_contract, sum((c>=b)::int) n_kept, sum(greatest(c-b,0)) gain_copies
from x group by all;

-- Reporting query 1 (loss split per group; produced the main loss table):
-- select flt, grp, bgrp, sum(n_fam) fams, sum(base_copies) bc,
--  round(100*sum(loss_copies)/sum(base_copies),1) pct_loss,
--  round(100*sum(complete_copies)/nullif(sum(loss_copies),0),1) pct_complete,
--  round(100*sum(contract_copies)/nullif(sum(loss_copies),0),1) pct_contract,
--  sum(n_complete) fam_lost, sum(n_contract) fam_contr
-- from kx2_h4_loss
-- where (grp,bgrp) in (('Kickxellales','Mucorales'),('Harpellales','Mucorales'),('Dimargaritales','Mucorales'),
--   ('Zoopagales','Mucorales'),('Microsporidia','Mucorales'),('Saccharomycetes','Pezizomycotina'),
--   ('Hanseniaspora','Pezizomycotina'),('Pneumocystis','Pezizomycotina'),('Taphrina','Pezizomycotina'),
--   ('Saccharomycotina_other','Pezizomycotina'),('Microsporidia','Pezizomycotina'),('Kickxellales','Pezizomycotina'))
-- group by all order by flt, bgrp, grp;

-- Reporting query 2 (by baseline copy class; single-copy family loss):
-- select flt, bgrp, grp, bclass, n_fam, round(100*n_complete/n_fam,1) pct_fam_lost,
--  round(100*loss_copies/base_copies,1) pct_copies_lost, round(100*complete_copies/nullif(loss_copies,0),0) pct_of_loss_complete,
--  round(100*gain_copies/base_copies,1) gain_pct
-- from kx2_h4_loss where flt='all' and grp in ('Kickxellales','Harpellales','Zoopagales','Saccharomycetes',
--  'Hanseniaspora','Pneumocystis','Taphrina','Microsporidia') order by bgrp, bclass, grp;

-- Reporting query 3 (observed vs uniform-thinning expectation; p = group copy loss fraction, expected P(lost)=p^round(b)):
-- with bl as (select grp bgrp, pfam_acc, med b from kx2_h4_gm where flt='all' and grp in ('Mucorales','Pezizomycotina') and prev>=0.5 and med>=1),
-- g as (select distinct grp from kx2_h4_gm where flt='all' and grp in ('Kickxellales','Harpellales','Dimargaritales','Zoopagales','Saccharomycetes','Saccharomycotina_other','Hanseniaspora','Pneumocystis','Taphrina','Microsporidia')),
-- x as (select g.grp, b.bgrp, b.pfam_acc, b.b, coalesce(m.med,0) c from bl b cross join g left join kx2_h4_gm m on m.flt='all' and m.grp=g.grp and m.pfam_acc=b.pfam_acc),
-- p as (select grp,bgrp, sum(greatest(b-c,0))/sum(b) p from x group by 1,2)
-- select x.bgrp, x.grp, round(p.p,3) p_copy_loss, sum((c=0)::int) obs_lost, round(sum(pow(p.p, round(b))),0) exp_lost_uniform,
--  round(sum((c=0)::int)/sum(pow(p.p, round(b))),2) obs_over_exp, sum((c=0 and b<=1)::int) obs_lost_single,
--  round(sum(case when b<=1 then p.p else 0 end),0) exp_single
-- from x join p using(grp,bgrp) group by x.bgrp,x.grp,p.p order by 1,2;
