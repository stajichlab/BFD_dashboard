#!/usr/bin/env python3
"""T5 step 2: GH65 candidate proteins and genome table for the genus representatives, from the LOCAL BFD.duckdb (read-only).

Usage (HPCC, SLURM job):  02_candidates.py <BENCH> <genus_representatives.tsv>
  genus_representatives.tsv  the atlas genus representatives (columns LOCUSTAG, GENUS, CLASS, busco_complete, pep_path), the same 1,654
                             genomes as the atlas genus tree, so that tips of the tree and rows of the table agree.
Database: lib/bfd_connect.connect(scratch=False): BFD.duckdb attached read-only as `bfd` (override with BFD_DB).
Writes $BENCH/cand/genomes.tsv   one row per genus representative: taxonomy, BUSCO, N50, GC, number of GH65 proteins
       $BENCH/cand/candidates.tsv  one row per protein with the GH65 central domain PF03632 (domain_i_evalue <= 1e-5): the Pfam domains it
                                   carries (N-terminal PF03636, central PF03632, C-terminal PF03633), length and gene coordinates.
GH65 domains: PF03632 Glyco_hydro_65m (central catalytic), PF03636 Glyco_hydro_65N, PF03633 Glyco_hydro_65C. Accessions were verified
in bfd.pfam (pfam_id strings are printed by the script so the log shows them).
"""
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'lib'))
from bfd_connect import connect  # noqa: E402

bench, reps_tsv = sys.argv[1], sys.argv[2]
EVAL = 1e-5
DOMS = ['PF03632', 'PF03636', 'PF03633']
reps = list(csv.DictReader(open(reps_tsv), delimiter='\t'))
loci = sorted({r['LOCUSTAG'] for r in reps})
print('genus representatives:', len(loci), flush=True)

con = connect(scratch=False, threads=8, memory='24GB', tmp=os.environ.get('SCRATCH'))
con.execute('CREATE TEMP TABLE reps(locustag VARCHAR, genus VARCHAR)')
con.executemany('INSERT INTO reps VALUES (?,?)', [(r['LOCUSTAG'], r['GENUS']) for r in reps])
print('Pfam ids:', con.execute(
    "SELECT DISTINCT split_part(pfam_acc,'.',1), pfam_id FROM pfam WHERE protein_id >= 'F0000039_' AND protein_id < 'F0000039`' "
    "AND split_part(pfam_acc,'.',1) IN ('PF03632','PF03636','PF03633')").fetchall(), flush=True)

rng = ' OR '.join(f"(protein_id >= '{l}_' AND protein_id < '{l}`')" for l in loci)
acc_in = ','.join(f"'{a}'" for a in DOMS)
con.execute(f"""
CREATE TEMP TABLE hits AS
SELECT protein_id, split_part(protein_id,'_',1) AS locustag, split_part(pfam_acc,'.',1) AS acc, domain_i_evalue, env_from, env_to
FROM pfam
WHERE ({rng}) AND domain_i_evalue <= {EVAL} AND split_part(pfam_acc,'.',1) IN ({acc_in})""")
print('domain hits:', con.execute('SELECT count(*) FROM hits').fetchone()[0], flush=True)

os.makedirs(f'{bench}/cand', exist_ok=True)
con.execute(f"""
COPY (
  WITH c AS (
    SELECT protein_id, locustag,
           string_agg(DISTINCT acc, ';' ORDER BY acc) AS domains,
           count(*) FILTER (WHERE acc='PF03632') AS n_central,
           min(domain_i_evalue) FILTER (WHERE acc='PF03632') AS best_evalue_central,
           min(env_from) AS first_from, max(env_to) AS last_to
    FROM hits GROUP BY 1, 2
    HAVING count(*) FILTER (WHERE acc='PF03632') >= 1)
  SELECT c.protein_id, c.locustag, r.genus, c.domains, c.n_central, c.best_evalue_central, c.first_from, c.last_to,
         p.length AS prot_len, p.transcript_id, i.chrom, i.start AS gstart, i."end" AS gend, i.strand
  FROM c JOIN reps r ON r.locustag = c.locustag
  LEFT JOIN gene_proteins p ON p.protein_id = c.protein_id
  LEFT JOIN gene_info i ON i.gene_id = p.gene_id
  ORDER BY c.locustag, c.protein_id
) TO '{bench}/cand/candidates.tsv' (FORMAT CSV, DELIMITER '\t', HEADER)""")

con.execute(f"""
COPY (
  SELECT r.locustag, r.genus, s.PHYLUM, s."CLASS" AS class, s."ORDER" AS ord, s.FAMILY, s.SPECIES,
         b.complete_pct AS busco_complete, a.N50_bp, a.TOTAL_LENGTH AS total_length_bp, a.GC_PERCENT AS gc_pct, a.contig_count,
         coalesce(g.n_gh65, 0) AS n_gh65, coalesce(g.n_full, 0) AS n_gh65_full_arch
  FROM reps r
  JOIN species s ON s.LOCUSTAG = r.locustag
  LEFT JOIN asm_stats a ON a.ASMID = s.ASMID
  LEFT JOIN busco_genome b ON b.ASMID = s.ASMID
  LEFT JOIN (SELECT locustag, count(*) AS n_gh65,
                    count(*) FILTER (WHERE domains LIKE '%PF03636%' AND domains LIKE '%PF03633%') AS n_full
             FROM (SELECT protein_id, locustag, string_agg(DISTINCT acc, ';') AS domains FROM hits GROUP BY 1, 2
                   HAVING count(*) FILTER (WHERE acc='PF03632') >= 1) GROUP BY 1) g ON g.locustag = r.locustag
  ORDER BY r.locustag
) TO '{bench}/cand/genomes.tsv' (FORMAT CSV, DELIMITER '\t', HEADER)""")
n, g = con.execute("SELECT count(*), count(DISTINCT locustag) FROM (SELECT protein_id, locustag FROM hits WHERE acc='PF03632' GROUP BY 1,2)").fetchone()
print(f'GH65 (PF03632) proteins: {n} in {g} of {len(loci)} genus representatives')
tabs = [r[0] for r in con.execute("SELECT table_name FROM duckdb_tables() WHERE database_name='bfd' AND table_name ILIKE '%cazy%'").fetchall()]
tabs += [r[0] for r in con.execute("SELECT view_name FROM duckdb_views() WHERE database_name='bfd' AND view_name ILIKE '%cazy%'").fetchall()]
print('cazy tables/views:', tabs)
