-- Purpose: decompose mean-length difference vs Mucorales: within-family (shared accession families: avg share x delta length),
--          composition (avg length x delta share), lost (Mucorales-only), gained (group-only), unassigned ('NONE').
-- Inputs: kx2_h1_fam_gen. Output: kx2_h1_decomp.
-- Report: core dMean +7.5, within +13.1, composition +36.7, lost -88.3, gained +46, lost_long(>=1000 aa) -14.1; Microsporidia, yeasts, Taphrinales likewise.
CREATE OR REPLACE TABLE kx2_h1_decomp AS
WITH M AS (SELECT * FROM kx2_h1_fam_gen WHERE grp='Mucorales'), GR AS (SELECT DISTINCT grp FROM kx2_h1_fam_gen WHERE grp<>'Mucorales'),
C AS (SELECT * FROM kx2_h1_fam_gen WHERE grp<>'Mucorales'),
FAMS AS (SELECT GR.grp, f.fam FROM GR CROSS JOIN (SELECT fam FROM C UNION SELECT fam FROM M) f),
J AS (SELECT F.grp, F.fam, coalesce(C.fshare,0) sc, coalesce(M.fshare,0) sm, C.fL lc, M.fL lm FROM FAMS F LEFT JOIN C ON C.grp=F.grp AND C.fam=F.fam LEFT JOIN M ON M.fam=F.fam WHERE C.fam IS NOT NULL OR M.fam IS NOT NULL)
SELECT grp, fam, sc, sm, lc, lm,
 CASE WHEN fam='NONE' THEN 'unassigned' WHEN lc IS NULL THEN 'lost_in_group' WHEN lm IS NULL THEN 'gained_in_group' ELSE 'shared' END AS cls,
 CASE WHEN lc IS NOT NULL AND lm IS NOT NULL THEN (sc+sm)/2*(lc-lm) END AS within_c,
 CASE WHEN lc IS NOT NULL AND lm IS NOT NULL THEN (lc+lm)/2*(sc-sm) END AS comp_c,
 coalesce(sc*lc,0)-coalesce(sm*lm,0) AS tot FROM J;
SELECT grp, round(sum(tot),1) dMean, round(sum(within_c),1) within_shared, round(sum(comp_c),1) comp_shared, round(sum(tot) FILTER (WHERE cls='unassigned'),1) unassigned, round(sum(tot) FILTER (WHERE cls='lost_in_group'),1) lost_fams, round(sum(tot) FILTER (WHERE cls='gained_in_group'),1) gained_fams, round(sum(tot) FILTER (WHERE cls='lost_in_group' AND lm>=1000),1) lost_long, round(sum(sm) FILTER (WHERE cls='lost_in_group'),3) lost_share_M, round(sum(sm) FILTER (WHERE cls='lost_in_group' AND lm>=1000),4) lostlong_share_M, round(sum(sc) FILTER (WHERE cls='shared'),3) sh_c, round(sum(sm) FILTER (WHERE cls='shared'),3) sh_m FROM kx2_h1_decomp GROUP BY 1 ORDER BY 1;
