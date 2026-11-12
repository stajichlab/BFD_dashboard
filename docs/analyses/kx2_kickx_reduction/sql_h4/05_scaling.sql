-- 05_scaling.sql
-- Purpose: kingdom-wide log-log scaling of functional-class counts vs proteome size (non-reduced genera),
--          residuals per genus and per group, class labels, and leverage of reduced groups. All in SQL (no Python).
-- Inputs : kx2_h4_cls, kx2_h4_meta, kx2_h4_gen
-- Outputs: kx2_h4_gcls (genus x class median count), kx2_h4_fit, kx2_h4_resid, kx2_h4_resid_grp, kx2_h4_scaling_summary
-- Method/thresholds:
--   y = log10(genus median count + 1), x = log10(genus median n_genes). Pept_C excluded (median 0).
--   Fit set: grp in (Pezizomycotina, Mucorales, Other_fungi) and genus median BUSCO >= 85 (1,308 genera). Reduced
--   groups (all Kickxellomycota, Zoopagales, yeasts, Microsporidia, Pneumocystis, Taphrina) are excluded from the fit.
--   OLS: regr_slope/regr_intercept. SMA (standardized major axis): slope = sign(r)*sd(y)/sd(x), through (xbar,ybar).
--   z = group median residual / sd of residuals in the fit set; fold = 10^median residual.
--   Label: over-contracted if z_ols<=-1 AND z_sma<=-1; expanded if both >=1; mixed_low/mixed_high if only one
--   fit crosses -1/+1; else on-curve.
--   frac_vs_muc / frac_vs_pez = (group median count / group median genes) / BUSCO>=85 baseline median fraction.
--   Residual groups are summarised with flt 'all' (no BUSCO filter) and 'b85'; scaling_summary uses 'all'.
-- Report numbers: scaling-residual / on-off-curve table; leverage line (fit proteome range 1.4k-15k, 5th pct 7.0k).

create or replace table kx2_h4_gcls as
select m.genus, m.grp, c.cls, median(c.n) n, any_value(g.n_genes) n_genes, any_value(g.busco) busco, any_value(g.n_sp) n_sp
from kx2_h4_cls c join kx2_h4_meta m using(genome) join kx2_h4_gen g on g.genus=m.genus and g.grp=m.grp
where c.cls<>'Pept_C' group by 1,2,3;

create or replace table kx2_h4_fit as
with d as (select cls, log10(n+1) y, log10(n_genes) x from kx2_h4_gcls
           where grp in ('Pezizomycotina','Mucorales','Other_fungi') and busco>=85)
select cls, count(*) n, regr_slope(y,x) ols_slope, regr_intercept(y,x) ols_int, regr_r2(y,x) r2,
 stddev(y)/stddev(x)*sign(corr(y,x)) sma_slope, avg(x) xbar, avg(y) ybar, stddev(x) sdx, stddev(y) sdy,
 stddev(y - 0) sd0 from d group by cls;

create or replace table kx2_h4_resid as
select g.genus, g.grp, g.cls, g.n, g.n_genes, g.busco, g.n_sp, log10(g.n+1) y, log10(g.n_genes) x,
 log10(g.n+1) - (f.ols_int + f.ols_slope*log10(g.n_genes)) res_ols,
 log10(g.n+1) - (f.ybar + f.sma_slope*(log10(g.n_genes)-f.xbar)) res_sma
from kx2_h4_gcls g join kx2_h4_fit f using(cls);

create or replace table kx2_h4_resid_grp as
with sd as (select cls, stddev(res_ols) sd_ols, stddev(res_sma) sd_sma from kx2_h4_resid
            where grp in ('Pezizomycotina','Mucorales','Other_fungi') and busco>=85 group by cls),
r as (select 'all' flt, * from kx2_h4_resid union all select 'b85', * from kx2_h4_resid where busco>=85)
select r.flt, r.grp, r.cls, count(*) n_gen, median(r.res_ols) res_ols, median(r.res_sma) res_sma,
 median(r.res_ols)/any_value(sd.sd_ols) z_ols, median(r.res_sma)/any_value(sd.sd_sma) z_sma,
 power(10,median(r.res_ols)) fold_ols, power(10,median(r.res_sma)) fold_sma, median(r.n) med_n, median(r.n_genes) med_genes
from r join sd using(cls) group by r.flt, r.grp, r.cls;

create or replace table kx2_h4_scaling_summary as
with b as (select cls, grp, median(n/n_genes) frac from kx2_h4_gcls where busco>=85 and grp in ('Mucorales','Pezizomycotina') group by 1,2),
r as (select * from kx2_h4_resid_grp where flt='all')
select r.grp, r.cls, r.n_gen, r.med_n, round(r.fold_ols,2) fold_ols, round(r.z_ols,1) z_ols, round(r.fold_sma,2) fold_sma, round(r.z_sma,1) z_sma,
 round((r.med_n/r.med_genes)/mu.frac,2) frac_vs_muc, round((r.med_n/r.med_genes)/pz.frac,2) frac_vs_pez,
 case when z_ols<=-1 and z_sma<=-1 then 'over-contracted' when z_ols>=1 and z_sma>=1 then 'expanded'
      when z_ols<=-1 or z_sma<=-1 then 'mixed_low' when z_ols>=1 or z_sma>=1 then 'mixed_high' else 'on-curve' end AS lab
from r join b mu on mu.cls=r.cls and mu.grp='Mucorales' join b pz on pz.cls=r.cls and pz.grp='Pezizomycotina'
where r.grp not in ('Other_fungi');

-- Leverage / extrapolation query (run, not stored):
-- with f as (select log10(n_genes) x from kx2_h4_gcls where cls='CAZy_GH' and grp in ('Pezizomycotina','Mucorales','Other_fungi') and busco>=85),
-- st as (select count(*) n, avg(x) xb, var_pop(x)*count(*) sxx, min(power(10,x)) mn, quantile_cont(power(10,x),0.05) q05, quantile_cont(power(10,x),0.5) q50 from f),
-- gg as (select grp, avg(log10(n_genes)) xg from kx2_h4_gcls where cls='CAZy_GH' and grp in ('Kickxellales','Harpellales','Zoopagales','Saccharomycetes','Hanseniaspora','Pneumocystis','Taphrina','Microsporidia') group by 1)
-- select gg.grp, round(power(10,xg)) genes, round(1.0/st.n + power(xg-st.xb,2)/st.sxx,3) h, round(st.mn) fit_min, round(st.q05) fit_q05, round(st.q50) fit_med, st.n from gg, st order by 1;
