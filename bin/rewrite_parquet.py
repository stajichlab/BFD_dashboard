#!/usr/bin/env python3
"""Rewrite BFD parquet tables for HTTP range-request access (DuckDB-WASM).

Reads the snapshot tables in --tables and writes browser-friendly copies to --out:

  * big tables re-sorted by genome so a per-genome query touches few row groups
  * ZSTD compression, bounded row-group sizes, column statistics (DuckDB default)
  * small per-genome / per-domain summary tables for the common questions

Genome key everywhere: species_prefix (== species.LOCUSTAG). pfam has no genome
column, so it is derived from protein_id ('F2D1C788_000724-T1' -> 'F2D1C788').

Usage:
    bin/rewrite_parquet.py --tables ../Fungi_BFD_runs/tables --out dist_parquet
    bin/rewrite_parquet.py --tables ... --out ... --only pfam_by_genome genome

Run it inside a SLURM job. Spill files go to $SCRATCH (node-local).
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import duckdb

ZSTD_LEVEL = 3


def q(path):
    return f"read_parquet('{path}')"


def build_jobs(t):
    """Return {name: (sql, row_group_size)}. `t` is the input tables directory."""
    T = lambda n: q(f"{t}/{n}.parquet")
    return {
        # --- dimension: one row per genome, all the per-genome scalars ---------
        "genome": (f"""
            SELECT s.LOCUSTAG AS species_prefix, s.ASMID, s.PHYLUM, s.SUBPHYLUM, s.CLASS,
                   s.SUBCLASS, s."ORDER" AS ORDER_, s.FAMILY, s.GENUS, s.SPECIES, s.STRAIN,
                   s.NCBI_TAXONID, s.BUSCO_LINEAGE,
                   a.contig_count, a.total_length_bp, a.N50_bp, a.L50, a.gc_pct, a.masked_pct,
                   a.t2t_scaffolds, b.complete_pct AS busco_complete_pct,
                   b.duplicated_pct AS busco_duplicated_pct
            FROM {T('species')} s
            LEFT JOIN {T('asm_stats')} a USING (ASMID)
            LEFT JOIN {T('busco_genome')} b USING (ASMID)
            ORDER BY s.LOCUSTAG""", 100_000),

        # --- protein dimension: length per protein (from gene_proteins) --------
        # tmhmm lists only proteins with >=1 predicted TM helix, so it must not be
        # the source of protein length; it only supplies PredHel (0 when absent).
        "protein": (f"""
            WITH tm AS (SELECT protein_id, max(PredHel) AS PredHel FROM {T('tmhmm')} GROUP BY 1)
            SELECT split_part(p.protein_id, '_', 1) AS species_prefix, p.protein_id,
                   p.length AS length_aa, coalesce(tm.PredHel, 0) AS PredHel
            FROM {T('gene_proteins')} p LEFT JOIN tm USING (protein_id)
            ORDER BY species_prefix, p.protein_id""", 250_000),

        # --- pfam, two sort orders ---------------------------------------------
        # by genome: per-genome profile / join to protein
        "pfam_by_genome": (f"""
            SELECT split_part(protein_id, '_', 1) AS species_prefix, protein_id,
                   pfam_acc, pfam_id, domain_num, domain_num_of,
                   domain_i_evalue, domain_score, hmm_from, hmm_to, ali_from, ali_to
            FROM {T('pfam')}
            ORDER BY species_prefix, protein_id, domain_num""", 122_880),
        # by domain: "all proteins with PF02364" reads few row groups
        "pfam_by_domain": (f"""
            SELECT split_part(protein_id, '_', 1) AS species_prefix, protein_id,
                   pfam_acc, pfam_id, domain_num, domain_num_of,
                   domain_i_evalue, domain_score, hmm_from, hmm_to, ali_from, ali_to
            FROM {T('pfam')}
            ORDER BY pfam_acc, species_prefix, protein_id, domain_num""", 122_880),

        # --- intergenic distances, sorted by genome ----------------------------
        "intergenic_by_genome": (f"""
            SELECT species_prefix, left_gene, right_gene, distance
            FROM {T('gene_intergenic_distances')}
            ORDER BY species_prefix, left_gene""", 250_000),

        # --- summaries ----------------------------------------------------------
        "pfam_dict": (f"""
            SELECT pfam_acc, any_value(pfam_id) AS pfam_id, any_value(pfam_len) AS pfam_len,
                   count(DISTINCT protein_id) AS n_proteins,
                   count(DISTINCT split_part(protein_id, '_', 1)) AS n_genomes
            FROM {T('pfam')} GROUP BY pfam_acc ORDER BY pfam_acc""", 100_000),

        "genome_pfam_counts": (f"""
            SELECT split_part(protein_id, '_', 1) AS species_prefix, pfam_acc,
                   count(DISTINCT protein_id) AS n_proteins, count(*) AS n_domains
            FROM {T('pfam')} GROUP BY 1, 2 ORDER BY 1, 2""", 250_000),

        "genome_protein_len": (f"""
            SELECT split_part(protein_id, '_', 1) AS species_prefix, count(*) AS n_proteins,
                   avg(length) AS mean_len, median(length) AS median_len,
                   quantile_cont(length, 0.10) AS p10_len, quantile_cont(length, 0.90) AS p90_len
            FROM {T('gene_proteins')} GROUP BY 1 ORDER BY 1""", 100_000),

        "genome_intergenic": (f"""
            SELECT species_prefix, count(*) AS n_pairs,
                   avg(distance) AS mean_dist, median(distance) AS median_dist,
                   quantile_cont(distance, 0.10) AS p10_dist,
                   quantile_cont(distance, 0.90) AS p90_dist
            FROM {T('gene_intergenic_distances')} WHERE distance >= 0
            GROUP BY 1 ORDER BY 1""", 100_000),
    }


# wide per-genome table: needs the summaries above, so it runs last from --out
def wide_job(out):
    O = lambda n: q(f"{out}/{n}.parquet")
    return ("genome_features", f"""
        SELECT g.*, p.n_proteins, p.mean_len, p.median_len, p.p10_len, p.p90_len,
               i.n_pairs, i.mean_dist, i.median_dist, i.p10_dist, i.p90_dist,
               f.n_pfam_domains  -- NULL = genome has no Pfam rows (used for coverage)
        FROM {O('genome')} g
        LEFT JOIN {O('genome_protein_len')} p USING (species_prefix)
        LEFT JOIN {O('genome_intergenic')} i USING (species_prefix)
        LEFT JOIN (SELECT species_prefix, sum(n_domains) AS n_pfam_domains
                   FROM {O('genome_pfam_counts')} GROUP BY 1) f USING (species_prefix)
        ORDER BY species_prefix""", 100_000)


def write(con, name, sql, rg, out):
    dest = Path(out) / f"{name}.parquet"
    tmp = dest.with_suffix(".parquet.tmp")
    t0 = time.time()
    con.execute(f"""COPY ({sql}) TO '{tmp}'
                    (FORMAT parquet, COMPRESSION zstd, COMPRESSION_LEVEL {ZSTD_LEVEL},
                     ROW_GROUP_SIZE {rg})""")
    tmp.replace(dest)
    rows = con.execute(f"SELECT count(*) FROM {q(dest)}").fetchone()[0]
    ngroups = con.execute(f"SELECT count(DISTINCT row_group_id) FROM parquet_metadata('{dest}')").fetchone()[0]
    info = {"table": name, "rows": rows, "row_groups": ngroups,
            "size_mb": round(dest.stat().st_size / 1e6, 1), "seconds": round(time.time() - t0)}
    print(json.dumps(info), flush=True)
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", required=True, help="directory of input *.parquet snapshots")
    ap.add_argument("--out", required=True, help="output directory (created)")
    ap.add_argument("--only", nargs="*", help="build only these tables")
    ap.add_argument("--memory", default="12GB")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_ON_NODE", 4)))
    a = ap.parse_args()

    scratch = os.environ.get("SCRATCH")
    if not scratch:
        sys.exit("SCRATCH is not set. Run inside a SLURM job.")
    tables, out = Path(a.tables).resolve(), Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()
    con.execute(f"SET memory_limit='{a.memory}'; SET threads={a.threads};")
    con.execute(f"SET temp_directory='{scratch}/duckdb_tmp'; SET preserve_insertion_order=true;")

    jobs = build_jobs(tables)
    names = a.only or list(jobs) + ["genome_features"]
    done = []
    for name, (sql, rg) in jobs.items():
        if name in names:
            done.append(write(con, name, sql, rg, out))
    if "genome_features" in names:
        n, sql, rg = wide_job(out)
        done.append(write(con, n, sql, rg, out))
    (out / "_build.json").write_text(json.dumps(done, indent=1))


if __name__ == "__main__":
    main()
