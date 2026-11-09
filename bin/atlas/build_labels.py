#!/usr/bin/env python3
"""Build A1 genome_labels.parquet (one row per genome) for the BFD summary atlas.

Inputs (small tables only; runs on the login node):
  species, samples, asm_stats, busco_genome  (BFD tables)
  FungalTraits 1.2 CSV (latin-1 on disk; converted to UTF-8 here)
  growth_form_overrides.tsv (genus -> growth_form for genera FungalTraits lacks)

Outputs:
  <out>/labels/genome_labels.parquet
  <out>/labels/fungaltraits_1.2.utf8.csv
  <out>/labels/growth_form_overrides.tsv (written only if absent)

Rules (PLAN.md section 4):
  representative per species = highest busco complete_pct, tie N50_bp (D6)
  genus representative       = species representative with highest complete_pct, tie N50_bp
  q_busco_rel = complete_pct >= 25th percentile of its CLASS (classes with >= 20 genomes),
                else of its PHYLUM
  q_dup = duplicated_pct <= 10 ; q_n50 = N50_bp >= 50000
"""
import argparse
import os
import duckdb

SACCH = ['Saccharomycetes', 'Pichiomycetes', 'Dipodascomycetes', 'Lipomycetes', 'Trigonopsidomycetes',
         'Alloascoideomycetes', 'Sporopachydermiomycetes']
TAPHRINA = ['Schizosaccharomycetes', 'Taphrinomycetes', 'Pneumocystomycetes', 'Neolectomycetes']
PUCCINIO_YEAST = ['Microbotryomycetes', 'Cystobasidiomycetes', 'Agaricostilbomycetes', 'Mixiomycetes',
                  'Spiculogloeomycetes', 'Classiculomycetes', 'Tritirachiomycetes']
USTILAGO = ['Malasseziomycetes', 'Ustilaginomycetes', 'Exobasidiomycetes', 'Moniliellomycetes']
PEZIZO = ['Eurotiomycetes', 'Dothideomycetes', 'Sordariomycetes', 'Leotiomycetes', 'Lecanoromycetes',
          'Pezizomycetes', 'Orbiliomycetes', 'Arthoniomycetes', 'Xylonomycetes', 'Geoglossomycetes',
          'Lichinomycetes', 'Laboulbeniomycetes', 'Coniocybomycetes', 'Candelariomycetes']


def sql_list(xs):
    return ','.join("'%s'" % x for x in xs)


def write_default_overrides(path, con, species, ft_utf8):
    """Genera absent from FungalTraits: assign growth form from the class (source=class_rule)."""
    rows = con.execute(f"""
        SELECT s.GENUS, s.CLASS, count(*) AS n
        FROM '{species}' s
        ANTI JOIN read_csv('{ft_utf8}', header=true, all_varchar=true) f ON s.GENUS = f.GENUS
        WHERE s.GENUS IS NOT NULL GROUP BY 1,2 ORDER BY 3 DESC""").fetchall()
    with open(path, 'w') as o:
        o.write('GENUS\tgrowth_form\tsource\tn_genomes\tCLASS\n')
        for g, c, n in rows:
            if c in SACCH or c in TAPHRINA or c in PUCCINIO_YEAST or c == 'Tremellomycetes' or c == 'Malasseziomycetes':
                gf = 'yeast'
            elif c in ('Ustilaginomycetes', 'Exobasidiomycetes', 'Moniliellomycetes'):
                gf = 'dimorphic_yeast'
            elif c == 'Lecanoromycetes':
                gf = 'thallus_photosynthetic'
            elif c is None:
                gf = ''
            else:
                gf = 'filamentous_mycelium'
            o.write(f'{g}\t{gf}\tclass_rule\t{n}\t{c or ""}\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tables', default='/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables')
    ap.add_argument('--fungaltraits', default='/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/refs/lifestyle/FungalTraits_1.2.csv')
    ap.add_argument('--out', default='/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data')
    a = ap.parse_args()
    lab = os.path.join(a.out, 'labels')
    os.makedirs(lab, exist_ok=True)
    ft_utf8 = os.path.join(lab, 'fungaltraits_1.2.utf8.csv')
    if not os.path.exists(ft_utf8):
        with open(a.fungaltraits, encoding='latin-1') as i, open(ft_utf8, 'w', encoding='utf-8') as o:
            o.write(i.read())
    con = duckdb.connect()
    con.execute("PRAGMA threads=2; SET memory_limit='4GB'")
    T = a.tables
    species = f'{T}/species.parquet'
    ovr = os.path.join(lab, 'growth_form_overrides.tsv')
    if not os.path.exists(ovr):
        write_default_overrides(ovr, con, species, ft_utf8)

    con.execute(f"""
    CREATE TABLE ft AS
      SELECT GENUS,
             any_value(Growth_form_template) AS growth_form_ft,
             any_value(primary_lifestyle) AS primary_lifestyle,
             any_value(Secondary_lifestyle) AS secondary_lifestyle,
             any_value(Plant_pathogenic_capacity_template) AS plant_pathogenic_capacity,
             any_value(Decay_substrate_template) AS decay_substrate,
             any_value(Animal_biotrophic_capacity_template) AS animal_biotrophic_capacity
      FROM read_csv('{ft_utf8}', header=true, all_varchar=true) GROUP BY 1;
    CREATE TABLE ovr AS SELECT GENUS, nullif(growth_form,'') AS growth_form_ovr
      FROM read_csv('{ovr}', header=true, delim='\t', all_varchar=true);
    CREATE TABLE base AS
      SELECT s.LOCUSTAG, s.ASMID, s.PHYLUM, s.CLASS, s."ORDER" AS ORDER_, s.FAMILY, s.GENUS, s.SPECIES, s.STRAIN,
             s.NCBI_TAXONID, sa.TRANSL_TABLE AS transl_table,
             a.annotation_status, a.total_length_bp, a.contig_count, a.N50_bp, a.gc_pct, a.masked_pct,
             b.complete_pct AS busco_complete, b.single_pct AS busco_single, b.duplicated_pct AS busco_duplicated,
             b.fragmented_pct AS busco_fragmented, b.missing_pct AS busco_missing,
             coalesce(o.growth_form_ovr, f.growth_form_ft) AS growth_form,
             CASE WHEN o.growth_form_ovr IS NOT NULL THEN 'override'
                  WHEN f.growth_form_ft IS NOT NULL THEN 'fungaltraits' END AS growth_form_source,
             f.primary_lifestyle, f.secondary_lifestyle, f.plant_pathogenic_capacity, f.decay_substrate,
             f.animal_biotrophic_capacity,
             (f.GENUS IS NOT NULL) AS in_fungaltraits
      FROM '{species}' s
      LEFT JOIN '{T}/samples.parquet' sa ON s.LOCUSTAG = sa.LOCUSTAG
      LEFT JOIN '{T}/asm_stats.parquet' a ON s.ASMID = a.ASMID
      LEFT JOIN '{T}/busco_genome.parquet' b ON s.ASMID = b.ASMID
      LEFT JOIN ft f ON s.GENUS = f.GENUS
      LEFT JOIN ovr o ON s.GENUS = o.GENUS;
    """)
    # Taxonomy-only yeast definition (sensitivity) and yeast origin
    con.execute(f"""
    CREATE TABLE lab AS
    SELECT *,
      CASE WHEN CLASS IN ({sql_list(SACCH + TAPHRINA + PUCCINIO_YEAST)}) OR CLASS IN ('Tremellomycetes','Malasseziomycetes') THEN 'yeast'
           WHEN CLASS IN ('Ustilaginomycetes','Exobasidiomycetes','Moniliellomycetes') THEN 'dimorphic_yeast'
           WHEN CLASS IS NULL THEN NULL ELSE 'filamentous_mycelium' END AS growth_form_tax,
      CASE WHEN growth_form = 'yeast' AND CLASS IN ({sql_list(SACCH)}) THEN 'Saccharomycotina'
           WHEN growth_form = 'yeast' AND CLASS IN ({sql_list(TAPHRINA)}) THEN 'Taphrinomycotina'
           WHEN growth_form = 'yeast' AND CLASS = 'Tremellomycetes' THEN 'Tremellomycetes'
           WHEN growth_form = 'yeast' AND CLASS IN ({sql_list(PUCCINIO_YEAST)}) THEN 'Pucciniomycotina-yeast'
           WHEN growth_form = 'yeast' AND CLASS IN ({sql_list(USTILAGO)}) THEN 'Ustilaginomycotina-yeast'
           WHEN growth_form = 'yeast' THEN 'yeast-other'
           WHEN growth_form = 'dimorphic_yeast' AND CLASS IN ({sql_list(PEZIZO)}) THEN 'Pezizomycotina-dimorphic'
           WHEN growth_form = 'dimorphic_yeast' AND CLASS IN ({sql_list(USTILAGO)}) THEN 'Ustilaginomycotina-dimorphic'
           WHEN growth_form = 'dimorphic_yeast' THEN 'dimorphic-other'
           WHEN growth_form = 'filamentous_mycelium' THEN 'filamentous'
           ELSE 'other' END AS yeast_origin
    FROM base;
    """)
    # quality thresholds (class-relative), representatives
    con.execute("""
    CREATE TABLE cls_thr AS
      SELECT CLASS, quantile_cont(busco_complete, 0.25) AS p25, quantile_cont(busco_complete, 0.5) AS p50, count(*) AS n
      FROM lab WHERE annotation_status='annotated' AND CLASS IS NOT NULL GROUP BY 1 HAVING count(*) >= 20;
    CREATE TABLE phy_thr AS
      SELECT PHYLUM, quantile_cont(busco_complete, 0.25) AS p25, quantile_cont(busco_complete, 0.5) AS p50
      FROM lab WHERE annotation_status='annotated' AND PHYLUM IS NOT NULL GROUP BY 1;
    CREATE TABLE out AS
    SELECT l.*,
      coalesce(c.p25, p.p25) AS busco_thr_p25, coalesce(c.p50, p.p50) AS busco_thr_p50,
      CASE WHEN c.CLASS IS NOT NULL THEN 'class' WHEN p.PHYLUM IS NOT NULL THEN 'phylum' END AS busco_thr_level,
      (l.annotation_status = 'annotated') AS q_ann,
      (l.busco_complete >= coalesce(c.p25, p.p25)) AS q_busco_rel,
      (l.busco_complete >= coalesce(c.p50, p.p50)) AS q_busco_rel50,
      (l.busco_duplicated <= 10) AS q_dup,
      (l.N50_bp >= 50000) AS q_n50,
      row_number() OVER (PARTITION BY l.SPECIES ORDER BY l.busco_complete DESC NULLS LAST, l.N50_bp DESC NULLS LAST, l.LOCUSTAG) = 1
        AND l.annotation_status = 'annotated' AND l.SPECIES IS NOT NULL AS is_representative
    FROM lab l LEFT JOIN cls_thr c USING (CLASS) LEFT JOIN phy_thr p USING (PHYLUM);
    CREATE TABLE out2 AS
    SELECT o.*,
      (o.is_representative AND o.GENUS IS NOT NULL AND
       row_number() OVER (PARTITION BY o.GENUS, o.is_representative ORDER BY o.busco_complete DESC NULLS LAST, o.N50_bp DESC NULLS LAST, o.LOCUSTAG) = 1)
       AS is_genus_rep
    FROM out o;
    """)
    outp = os.path.join(lab, 'genome_labels.parquet')
    tmp = outp + '.tmp'
    con.execute(f"COPY (SELECT * FROM out2 ORDER BY LOCUSTAG) TO '{tmp}' (FORMAT PARQUET, COMPRESSION ZSTD)")
    os.replace(tmp, outp)
    for r in con.execute("""
        SELECT 'genomes', count(*) FROM out2 UNION ALL
        SELECT 'representatives', count(*) FILTER (WHERE is_representative) FROM out2 UNION ALL
        SELECT 'genus_reps', count(*) FILTER (WHERE is_genus_rep) FROM out2 UNION ALL
        SELECT 'growth_form:' || coalesce(growth_form,'NULL'), count(*) FROM out2 GROUP BY 1 UNION ALL
        SELECT 'rep growth_form:' || coalesce(growth_form,'NULL'), count(*) FROM out2 WHERE is_representative GROUP BY 1 UNION ALL
        SELECT 'yeast_origin:' || yeast_origin, count(*) FILTER (WHERE is_representative) FROM out2 GROUP BY 1
        """).fetchall():
        print(*r, sep='\t')
    print('wrote', outp)


if __name__ == '__main__':
    main()
