-- Purpose: genus medians per feature (replication unit = genus), clade group map, group summaries and Cliff's delta effect sizes.
-- Inputs: kx2_kc_wide
-- Outputs: bfd_scratch.kx2_kc_grp (PHYLUM, GENUS -> grp), bfd_scratch.kx2_kc_genus_long (genus x feature: median, frac genomes >0, n genomes),
--          bfd_scratch.kx2_kc_group_summary (group x feature: n genera, median of genus medians, share of genera with >=1 gene),
--          bfd_scratch.kx2_kc_cliff (KxHp vs each comparison group, Cliff's delta on genus medians)
-- Groups: KxHp = Kickxellales + Harpellales (Ku-lost clade); DimRam = Dimargaritales + Ramicandelaberales (Ku-retaining Kickxellomycota);
--         ZoopEnto = Zoopagomycota + Entomophthoromycota + Basidiobolomycota; Mucoro = Mucoromycota + Mortierellomycota + Glomeromycota;
--         EarlyOther = other early-diverging phyla; Dikarya = Ascomycota + Basidiomycota.
-- Thresholds: q_loose genomes only (BUSCO >= 60, >= 3000 genes, BUSCO dup <= 10). Cliff's delta = (P(a>b) - P(a<b)) over genus medians; |d|>=0.47 large (Romano convention).
-- Produces: per-group table of alternative-EJ / HR / telomerase features, effect sizes, proteome-size (n_genes) control.
CREATE OR REPLACE TABLE kx2_kc_grp AS
SELECT DISTINCT PHYLUM, GENUS,
 CASE WHEN ORD IN ('Kickxellales','Harpellales') THEN 'KxHp'
      WHEN ORD IN ('Dimargaritales','Ramicandelaberales') THEN 'DimRam'
      WHEN PHYLUM IN ('Zoopagomycota','Entomophthoromycota','Basidiobolomycota') THEN 'ZoopEnto'
      WHEN PHYLUM IN ('Mucoromycota','Mortierellomycota','Glomeromycota') THEN 'Mucoro'
      WHEN PHYLUM IN ('Ascomycota','Basidiomycota') THEN 'Dikarya'
      ELSE 'EarlyOther' END AS grp
FROM kx2_kc_wide WHERE q_loose;

CREATE OR REPLACE TABLE kx2_kc_genus_long AS
WITH u AS (
 UNPIVOT (SELECT * EXCLUDE (LOCUSTAG, ASMID, SUBPHYLUM, CLASS, ORD, FAMILY, SPECIES, STRAIN, complete_pct, duplicated_pct, N50_bp, contig_count, total_length_bp, t2t_scaffolds, n_prot, q_strict)
          FROM kx2_kc_wide WHERE q_loose)
 ON COLUMNS(* EXCLUDE (PHYLUM, GENUS, q_loose)) INTO NAME feature VALUE val)
SELECT PHYLUM, GENUS, feature, median(val) AS med, avg((val>0)::int) AS frac_pos, count(*) AS n_genomes FROM u GROUP BY ALL;

CREATE OR REPLACE TABLE kx2_kc_group_summary AS
SELECT g.grp, l.feature, count(*) AS n_genera, median(l.med) AS med_of_genus_med, avg((l.med>0)::int) AS frac_genera_pos, avg(l.med) AS mean_of_genus_med
FROM kx2_kc_genus_long l JOIN kx2_kc_grp g USING (PHYLUM, GENUS) GROUP BY ALL;

CREATE OR REPLACE TABLE kx2_kc_cliff AS
WITH a AS (SELECT l.feature, l.med FROM kx2_kc_genus_long l JOIN kx2_kc_grp g USING (PHYLUM, GENUS) WHERE g.grp='KxHp'),
 b AS (SELECT g.grp, l.feature, l.med FROM kx2_kc_genus_long l JOIN kx2_kc_grp g USING (PHYLUM, GENUS) WHERE g.grp<>'KxHp')
SELECT b.grp AS vs_group, a.feature, count(*) AS pairs,
 (sum((a.med>b.med)::int) - sum((a.med<b.med)::int))::double/count(*) AS cliff_delta
FROM a JOIN b USING (feature) GROUP BY ALL;
SELECT grp, count(*) FROM kx2_kc_grp GROUP BY 1 ORDER BY 1;
