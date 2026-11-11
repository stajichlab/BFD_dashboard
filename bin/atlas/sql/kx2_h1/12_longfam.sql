-- Purpose: counts of long/multidomain families per genome (PKS, NRPS, multi-NBD ABC, >=1000/2000 aa, >=4 distinct domains).
-- Inputs: kx2_h1_prot. Output: kx2_h1_longfam.
-- PKS = PF00109+PF00698, no PF16073, length>=1500; NRPS = PF00501+PF00668; ABC = >=2 PF00005 copies (7-char accession string count).
-- Filter: q genomes plus all Harpellales/Dimargaritales. Report: long-family counts (PKS, NRPS, ABC, >=1000, >=2000, >=4 domains per 1000).
CREATE OR REPLACE TABLE kx2_h1_longfam AS
SELECT lt AS LOCUSTAG, count(*) n_prot,
 count(*) FILTER (WHERE arch LIKE '%PF00109%' AND arch LIKE '%PF00698%' AND arch NOT LIKE '%PF16073%' AND length>=1500) AS pks,
 count(*) FILTER (WHERE arch LIKE '%PF00501%' AND arch LIKE '%PF00668%') AS nrps,
 count(*) FILTER (WHERE (length(arch)-length(replace(arch,'PF00005','')))/7>=2) AS abc_ge2nbd,
 count(*) FILTER (WHERE length>=2000) AS n_ge2000, count(*) FILTER (WHERE length>=1000) AS n_ge1000,
 count(*) FILTER (WHERE n_dist>=4) AS n_ge4dist, count(*) FILTER (WHERE n_dist>=2) AS n_ge2dist
FROM kx2_h1_prot GROUP BY 1;
WITH x AS (SELECT g.grp, g.GENUS, l.* FROM kx2_h1_longfam l JOIN kx2_h1_genomes g USING(LOCUSTAG) WHERE g.q OR g.grp IN ('Harpellales','Dimargaritales')),
y AS (SELECT grp, GENUS, median(pks) pks, median(nrps) nrps, median(abc_ge2nbd) abc, median(n_ge2000) ge2000, median(n_ge1000) ge1000, median(n_ge4dist) ge4d, median(n_ge2dist) ge2d, median(1000.0*n_ge2dist/n_prot) ge2d_per1000, median(1000.0*n_ge4dist/n_prot) ge4d_per1000 FROM x GROUP BY 1,2)
SELECT grp, count(*) ngen, median(pks) pks, median(nrps) nrps, median(abc) abc, median(ge1000) ge1000, median(ge2000) ge2000, median(ge2d) n_multidom, median(ge4d) n_ge4dist, round(median(ge4d_per1000),1) ge4d_per1000 FROM y GROUP BY 1 ORDER BY 1;
