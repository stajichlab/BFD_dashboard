#!/usr/bin/env python3
"""H2b build: genome-level compaction traits, family-level consistency of shortening/lengthening,
and location of the length loss for the top families. Runs under SLURM (DuckDB). Writes atlas_data/A10/.
Restartable: every output is written to .tmp then renamed and skipped when present (use --force to redo).
"""
import argparse
import os
import sys
import time

import duckdb
import numpy as np
import pandas as pd

ROOT = '/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard'
A = f'{ROOT}/atlas_data'
T = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables'
OUT = f'{A}/A10'
GA = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation'
GATE = "qcov_x10 >= 800 AND scov_x10 >= 800 AND pident_x10 >= 400"
ORIGINS = ['Saccharomycotina', 'Taphrinomycotina', 'Tremellomycetes', 'Pucciniomycotina-yeast', 'Ustilaginomycotina-yeast']


def done(name):
    return os.path.exists(f'{OUT}/{name}')


def write(con, name, sql):
    if done(name) and not FORCE:
        print('skip', name, flush=True); return
    tmp = f'{OUT}/{name}.tmp'; t0 = time.time()
    con.execute(f"COPY ({sql}) TO '{tmp}' (FORMAT PARQUET, COMPRESSION ZSTD)")
    os.replace(tmp, f'{OUT}/{name}')
    n = con.execute(f"SELECT count(*) FROM '{OUT}/{name}'").fetchone()[0]
    print(f'done {name} rows={n} {time.time()-t0:.0f}s', flush=True)


def main():
    global FORCE
    ap = argparse.ArgumentParser(); ap.add_argument('--tmp', required=True); ap.add_argument('--threads', type=int, default=16)
    ap.add_argument('--memory', default='80GB'); ap.add_argument('--force', action='store_true'); a = ap.parse_args(); FORCE = a.force
    os.makedirs(OUT, exist_ok=True); os.makedirs(a.tmp, exist_ok=True)
    con = duckdb.connect(); con.execute(f"SET threads={a.threads}; SET memory_limit='{a.memory}'; SET temp_directory='{a.tmp}'; SET preserve_insertion_order=false")

    # ---- representative genomes with default filters (as render_views.Ctx.reps())
    con.execute(f"""
    CREATE TABLE reps AS
      SELECT LOCUSTAG AS genome, GENUS, CLASS, ORDER_, PHYLUM, SPECIES, growth_form, yeast_origin, busco_complete, busco_duplicated, N50_bp, total_length_bp, gc_pct,
             log10(N50_bp) AS log_n50, log10(total_length_bp) AS log_genome,
             CASE WHEN yeast_origin = 'Saccharomycotina' THEN 'Saccharomycotina'
                  WHEN growth_form = 'filamentous_mycelium' AND CLASS IN ('Eurotiomycetes','Dothideomycetes','Sordariomycetes','Leotiomycetes','Pezizomycetes','Orbiliomycetes','Xylonomycetes','Geoglossomycetes') THEN 'Pezizomycotina-filamentous'
                  WHEN growth_form = 'filamentous_mycelium' THEN 'filamentous-other'
                  ELSE yeast_origin END AS clade_group
      FROM '{A}/labels/genome_labels.parquet'
      WHERE is_representative AND q_ann AND coalesce(q_busco_rel, false) AND coalesce(q_dup, false) AND coalesce(q_n50, false);
    CREATE TABLE fil AS SELECT genome, GENUS FROM reps WHERE growth_form = 'filamentous_mycelium';
    """)
    print('reps', con.execute("SELECT count(*), count(DISTINCT GENUS) FROM reps").fetchone(), flush=True)

    # ---- 1. per-genome offset (P-A) + compaction traits + quality
    con.execute(f"""
    CREATE TABLE fam AS SELECT sp_id AS fam_key, genome, med_len_model AS len FROM '{A}/A4/fam_sp_genome.parquet' WHERE med_len_model IS NOT NULL;
    CREATE TABLE ref AS
      SELECT fam_key, median(gl) AS ref_len, count(*) AS n_fil_genera
      FROM (SELECT f.fam_key, l.GENUS, median(f.len) AS gl FROM fam f JOIN fil l USING (genome) GROUP BY 1,2) GROUP BY 1 HAVING count(*) >= 10;
    CREATE TABLE goff AS
      SELECT f.genome, median(log2(f.len / r.ref_len)) AS offset, count(*) AS n_fam
      FROM fam f JOIN ref r USING (fam_key) GROUP BY 1 HAVING count(*) >= 200;
    """)
    # RNA-seq training flag: existence of training/funannotate_train.pasa.gff3 for the genome's `out` folder
    m = pd.read_csv(f'{A}/phylo/locustag_to_pep.tsv', sep='\t')
    gen = con.execute("SELECT genome FROM goff").df()['genome']
    m = m[m['LOCUSTAG'].isin(gen)]
    flags = []
    for lt, out in zip(m['LOCUSTAG'], m['out']):
        d = f'{GA}/{out}/training'
        flags.append({'genome': lt, 'rnaseq_trained': os.path.exists(f'{d}/funannotate_train.pasa.gff3') or os.path.exists(f'{d}/kallisto.tsv')})
    con.register('rnaflag', pd.DataFrame(flags))
    write(con, 'genome_offset_traits.parquet', f"""
      SELECT g.*, r.* EXCLUDE (genome), a2.n_protein_coding, a2.n_protein_coding * 1e6 / r.total_length_bp AS genes_per_mb, a2.median_intergenic, a2.introns_per_transcript,
             a2.frac_intronless, a2.partial_fraction, a2.median_protein_len, rf.rnaseq_trained
      FROM goff g JOIN reps r USING (genome) JOIN '{A}/A2/genome_gene_stats.parquet' a2 USING (genome) LEFT JOIN rnaflag rf USING (genome)""")

    # ---- 2. family consistency tables for P-A, P-B (YEAST, NEUCR), P-C
    sources = {
        'P-A': f"SELECT sp_id::VARCHAR AS fam_key, genome, med_len_model AS len, med_n_domains AS ndom FROM '{A}/A4/fam_sp_genome.parquet' WHERE med_len_model IS NOT NULL",
        'P-B-YEAST': f"SELECT anchor_sp_id::VARCHAR AS fam_key, genome, med_len_model AS len, NULL::DOUBLE AS ndom FROM '{A}/A4/fam_anchor_genome.parquet' WHERE org_code='YEAST' AND med_len_model IS NOT NULL",
        'P-B-NEUCR': f"SELECT anchor_sp_id::VARCHAR AS fam_key, genome, med_len_model AS len, NULL::DOUBLE AS ndom FROM '{A}/A4/fam_anchor_genome.parquet' WHERE org_code='NEUCR' AND med_len_model IS NOT NULL",
        'P-C': f"SELECT pfam_arch AS fam_key, genome, med_len_model AS len, NULL::DOUBLE AS ndom FROM '{A}/A4/fam_arch_genome.parquet' WHERE med_len_model IS NOT NULL",
    }
    for name, src in sources.items():
        fn = f"fam_consistency_{name}.parquet"
        if done(fn) and not FORCE:
            print('skip', fn, flush=True); continue
        con.execute(f"""
        CREATE OR REPLACE TABLE gl AS
          SELECT s.fam_key, r.GENUS, r.clade_group, r.growth_form, r.yeast_origin, median(s.len) AS glen, median(s.ndom) AS gdom
          FROM ({src}) s JOIN reps r USING (genome) GROUP BY 1,2,3,4,5;
        CREATE OR REPLACE TABLE fstat AS
          SELECT fam_key,
            count(*) FILTER (WHERE growth_form='filamentous_mycelium') AS n_fil,
            median(glen) FILTER (WHERE growth_form='filamentous_mycelium') AS fil_median,
            quantile_cont(glen, 0.25) FILTER (WHERE growth_form='filamentous_mycelium') AS fil_p25,
            quantile_cont(glen, 0.10) FILTER (WHERE growth_form='filamentous_mycelium') AS fil_p10,
            quantile_cont(glen, 0.75) FILTER (WHERE growth_form='filamentous_mycelium') AS fil_p75,
            quantile_cont(glen, 0.90) FILTER (WHERE growth_form='filamentous_mycelium') AS fil_p90,
            median(gdom) FILTER (WHERE growth_form='filamentous_mycelium') AS fil_ndom,
            count(*) FILTER (WHERE growth_form='yeast') AS n_yeast,
            median(glen) FILTER (WHERE growth_form='yeast') AS yeast_median,
            median(gdom) FILTER (WHERE growth_form='yeast') AS yeast_ndom,
            count(*) FILTER (WHERE clade_group='Saccharomycotina') AS n_sacch,
            median(glen) FILTER (WHERE clade_group='Saccharomycotina') AS sacch_median,
            median(gdom) FILTER (WHERE clade_group='Saccharomycotina') AS sacch_ndom
          FROM gl GROUP BY 1 HAVING n_fil >= 10 AND n_yeast >= 10;
        CREATE OR REPLACE TABLE fcnt AS
          SELECT g.fam_key,
            count(*) FILTER (WHERE g.growth_form='yeast' AND g.glen < f.fil_p25) AS y_b25,
            count(*) FILTER (WHERE g.growth_form='yeast' AND g.glen < f.fil_p10) AS y_b10,
            count(*) FILTER (WHERE g.growth_form='yeast' AND g.glen > f.fil_p90) AS y_a90,
            count(*) FILTER (WHERE g.clade_group='Saccharomycotina' AND g.glen < f.fil_p25) AS s_b25,
            count(*) FILTER (WHERE g.clade_group='Saccharomycotina' AND g.glen < f.fil_p10) AS s_b10,
            count(*) FILTER (WHERE g.clade_group='Saccharomycotina' AND g.glen > f.fil_p90) AS s_a90,
            count(*) FILTER (WHERE g.clade_group='Saccharomycotina' AND g.glen > f.fil_p75) AS s_a75
          FROM gl g JOIN fstat f USING (fam_key) GROUP BY 1;
        CREATE OR REPLACE TABLE fcons AS
          SELECT f.*,
            c.y_b25 * 1.0 / nullif(f.n_yeast,0) AS yeast_frac_below_p25, c.y_b10 * 1.0 / nullif(f.n_yeast,0) AS yeast_frac_below_p10, c.y_a90 * 1.0 / nullif(f.n_yeast,0) AS yeast_frac_above_p90,
            c.s_b25 * 1.0 / nullif(f.n_sacch,0) AS sacch_frac_below_p25, c.s_b10 * 1.0 / nullif(f.n_sacch,0) AS sacch_frac_below_p10,
            c.s_a90 * 1.0 / nullif(f.n_sacch,0) AS sacch_frac_above_p90, c.s_a75 * 1.0 / nullif(f.n_sacch,0) AS sacch_frac_above_p75,
            f.fil_median - f.sacch_median AS loss_aa, (f.fil_median - f.sacch_median) / f.fil_median AS loss_frac,
            f.sacch_ndom - f.fil_ndom AS ndom_gain_sacch
          FROM fstat f JOIN fcnt c USING (fam_key);
        CREATE OR REPLACE TABLE orig AS
          SELECT fam_key, yeast_origin, count(*) AS n_gen, median(glen) AS omed FROM gl WHERE growth_form='yeast' AND yeast_origin IN ({','.join("'%s'" % o for o in ORIGINS)}) GROUP BY 1,2;
        CREATE OR REPLACE TABLE odir AS
          SELECT o.fam_key,
            count(*) FILTER (WHERE o.n_gen >= 3) AS n_origins_tested,
            count(*) FILTER (WHERE o.n_gen >= 3 AND o.omed < f.fil_median) AS n_origins_shorter,
            count(*) FILTER (WHERE o.n_gen >= 3 AND o.omed > f.fil_median) AS n_origins_longer,
            string_agg(o.yeast_origin || ':' || round(o.omed - f.fil_median, 1)::VARCHAR, '; ' ORDER BY o.yeast_origin) FILTER (WHERE o.n_gen >= 3) AS origin_deltas
          FROM orig o JOIN fstat f USING (fam_key) GROUP BY 1;
        """)
        write(con, fn, "SELECT c.*, d.n_origins_tested, d.n_origins_shorter, d.n_origins_longer, d.origin_deltas FROM fcons c LEFT JOIN odir d USING (fam_key)")

    # ---- 3. location of the loss for top P-A shrinkers and gainers (alignment coordinates in swissprot_best)
    pa = con.execute(f"SELECT * FROM '{OUT}/fam_consistency_P-A.parquet'").df()
    shr = pa[(pa['n_sacch'] >= 10) & (pa['sacch_frac_below_p10'] >= 0.8) & ((pa['loss_aa'] >= 20) | (pa['loss_frac'] >= 0.05))]
    gai = pa[(pa['n_sacch'] >= 10) & (pa['sacch_frac_above_p90'] >= 0.8) & ((-pa['loss_aa'] >= 20) | (-pa['loss_frac'] >= 0.05))]
    keys = sorted(set(shr['fam_key']) | set(gai['fam_key']))
    print('top P-A shrinkers', len(shr), 'gainers', len(gai), flush=True)
    con.register('topkeys', pd.DataFrame({'fam_key': keys}))
    if not done('top_family_proteins.parquet') or FORCE:
        con.execute(f"""
        CREATE OR REPLACE TABLE topp AS
          SELECT s.protein_id, s.species_prefix AS genome, s.sp_id::VARCHAR AS fam_key, s.q_len, s.s_len, s.q_start, s.q_end, s.s_start, s.s_end,
                 (s.s_start - s.q_start) AS nterm_loss, ((s.s_len - s.s_end) - (s.q_len - s.q_end)) AS cterm_loss,
                 ((s.s_end - s.s_start) - (s.q_end - s.q_start)) AS internal_loss,
                 c.q_model, c.gene_id, c.n_domains, c.pfam_arch, r.GENUS, r.clade_group, r.growth_form
          FROM '{T}/swissprot_best.parquet' s JOIN topkeys k ON s.sp_id::VARCHAR = k.fam_key
          JOIN reps r ON s.species_prefix = r.genome
          JOIN '{A}/A3/protein_core.parquet' c ON s.protein_id = c.protein_id
          WHERE {GATE.replace('qcov_x10','s.qcov_x10').replace('scov_x10','s.scov_x10').replace('pident_x10','s.pident_x10')} AND c.q_model;
        CREATE OR REPLACE TABLE idp AS
          SELECT i.protein_id AS gene_id, i.IDP_residues, i.IDP_fraction FROM '{T}/idp_summary.parquet' i
          WHERE i.protein_id IN (SELECT gene_id FROM topp);
        """)
        write(con, 'top_family_proteins.parquet', "SELECT t.*, i.IDP_residues, i.IDP_fraction FROM topp t LEFT JOIN idp i USING (gene_id)")
    # per family x group medians of the loss components, domain presence and disorder
    write(con, 'top_family_location.parquet', f"""
      WITH g AS (
        SELECT fam_key, GENUS, clade_group, growth_form, median(nterm_loss) AS nt, median(cterm_loss) AS ct, median(internal_loss) AS it,
               median(q_len) AS qlen, median(n_domains) AS nd, median(IDP_residues) AS idp_res, median(IDP_fraction) AS idp_frac, count(*) AS n_prot
        FROM '{OUT}/top_family_proteins.parquet' GROUP BY 1,2,3,4)
      SELECT fam_key,
        count(*) FILTER (WHERE growth_form='filamentous_mycelium') AS n_fil_genera, count(*) FILTER (WHERE clade_group='Saccharomycotina') AS n_sacch_genera,
        median(nt) FILTER (WHERE clade_group='Saccharomycotina') - median(nt) FILTER (WHERE growth_form='filamentous_mycelium') AS d_nterm,
        median(ct) FILTER (WHERE clade_group='Saccharomycotina') - median(ct) FILTER (WHERE growth_form='filamentous_mycelium') AS d_cterm,
        median(it) FILTER (WHERE clade_group='Saccharomycotina') - median(it) FILTER (WHERE growth_form='filamentous_mycelium') AS d_internal,
        median(qlen) FILTER (WHERE clade_group='Saccharomycotina') AS sacch_len, median(qlen) FILTER (WHERE growth_form='filamentous_mycelium') AS fil_len,
        median(nd) FILTER (WHERE clade_group='Saccharomycotina') AS sacch_ndom, median(nd) FILTER (WHERE growth_form='filamentous_mycelium') AS fil_ndom,
        median(idp_res) FILTER (WHERE clade_group='Saccharomycotina') AS sacch_idp_res, median(idp_res) FILTER (WHERE growth_form='filamentous_mycelium') AS fil_idp_res,
        median(idp_frac) FILTER (WHERE clade_group='Saccharomycotina') AS sacch_idp_frac, median(idp_frac) FILTER (WHERE growth_form='filamentous_mycelium') AS fil_idp_frac
      FROM g GROUP BY 1""")
    # domain presence per family x group (which Pfam domains are lost or gained)
    write(con, 'top_family_domains.parquet', f"""
      WITH u AS (
        SELECT fam_key, GENUS, clade_group, growth_form, unnest(string_split(pfam_arch, ';')) AS dom
        FROM '{OUT}/top_family_proteins.parquet' WHERE pfam_arch IS NOT NULL),
      gp AS (SELECT fam_key, clade_group, growth_form, GENUS, dom, count(*) AS n FROM u GROUP BY ALL),
      gsize AS (SELECT fam_key, GENUS, clade_group, growth_form, count(DISTINCT protein_id) AS np FROM '{OUT}/top_family_proteins.parquet' GROUP BY ALL),
      pres AS (SELECT gp.fam_key, gp.dom, gp.clade_group, gp.growth_form, gp.GENUS, gp.n * 1.0 / gs.np >= 0.5 AS present FROM gp JOIN gsize gs USING (fam_key, GENUS, clade_group, growth_form)),
      ng AS (SELECT fam_key, count(DISTINCT GENUS) FILTER (WHERE growth_form='filamentous_mycelium') AS n_fil, count(DISTINCT GENUS) FILTER (WHERE clade_group='Saccharomycotina') AS n_sacch FROM gsize GROUP BY 1),
      cnt AS (SELECT fam_key, dom, count(DISTINCT GENUS) FILTER (WHERE growth_form='filamentous_mycelium' AND present) AS fil_present,
                     count(DISTINCT GENUS) FILTER (WHERE clade_group='Saccharomycotina' AND present) AS sacch_present FROM pres GROUP BY 1,2)
      SELECT c.fam_key, c.dom, c.fil_present * 1.0 / nullif(ng.n_fil,0) AS fil_frac, c.sacch_present * 1.0 / nullif(ng.n_sacch,0) AS sacch_frac
      FROM cnt c JOIN ng USING (fam_key)""")
    print('all done', flush=True)


if __name__ == '__main__':
    main()
