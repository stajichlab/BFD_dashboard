#!/usr/bin/env python3
"""Step 2: annotated candidate proteins for each family, from the BFD Pfam table (read-only), in the genomes used for the trees.

Usage: 02_candidates.py <BENCH> <backbone_taxa.tsv> <strains genomes.tsv> [tables_dir]

Genomes: the 63 backbone genera (backbone_taxa.tsv, column LOCUSTAG) plus the additional Coemansia strains in genomes.tsv that are
not in the backbone (column locustag). Family domains are Pfam accessions verified in bfd.swissprot_subject / bfd.pfam:
  lig  : PF01068 DNA_ligase_A_M          (ATP-dependent ligase catalytic core)
  reca : PF08423 Rad51                   (RecA/Rad51 ATPase domain)
  polx : PF14792 DNA_pol_B_palm or PF14791 DNA_pol_B_thumb  (Pol X palm / thumb)
A protein is a candidate if it has the domain at domain_i_evalue <= 1e-5. Output $BENCH/cand/candidates.tsv with the proteins, their
lengths (gene_proteins) and gene coordinates (gene_info) so that miniprot rescue models can later be checked for overlap.
pfam.parquet is ordered by protein_id, so per-genome prefix ranges let DuckDB skip row groups instead of scanning the whole table.
"""
import csv, os, sys
import duckdb

bench, bb_tsv, st_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
tables = sys.argv[4] if len(sys.argv) > 4 else '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables'
FAM = {'lig': ['PF01068'], 'reca': ['PF08423'], 'polx': ['PF14792', 'PF14791']}
EVAL = 1e-5

loci = [r['LOCUSTAG'] for r in csv.DictReader(open(bb_tsv), delimiter='\t')]
loci += [r['locustag'] for r in csv.DictReader(open(st_tsv), delimiter='\t')]
loci = sorted(set(loci))
print('genomes:', len(loci), flush=True)
rng = ' OR '.join(f"(protein_id >= '{l}_' AND protein_id < '{l}`')" for l in loci)
doms = sorted({a for v in FAM.values() for a in v})
acc_in = ','.join(f"'{a}'" for a in doms)

con = duckdb.connect()
con.execute("SET threads=8; SET memory_limit='24GB'; SET preserve_insertion_order=false")
con.execute(f"""
CREATE TEMP TABLE hits AS
SELECT protein_id, split_part(protein_id,'_',1) AS locustag, split_part(pfam_acc,'.',1) AS acc, pfam_id, domain_i_evalue, env_from, env_to
FROM '{tables}/pfam.parquet'
WHERE ({rng}) AND domain_i_evalue <= {EVAL} AND split_part(pfam_acc,'.',1) IN ({acc_in})""")
print('domain hits:', con.execute('SELECT count(*) FROM hits').fetchone()[0], flush=True)
rows = []
for fam, accs in FAM.items():
    a = ','.join(f"'{x}'" for x in accs)
    q = f"""
    SELECT '{fam}' AS family, h.protein_id, h.locustag, string_agg(DISTINCT h.pfam_id || ':' || h.acc, ';' ORDER BY h.pfam_id || ':' || h.acc) AS domains,
           count(*) AS n_domain_hits, min(h.domain_i_evalue) AS best_evalue, min(h.env_from) AS first_from, max(h.env_to) AS last_to
    FROM hits h WHERE h.acc IN ({a}) GROUP BY h.protein_id, h.locustag"""
    rows += con.execute(q).fetchall()
con.execute("CREATE TEMP TABLE cand(family VARCHAR, protein_id VARCHAR, locustag VARCHAR, domains VARCHAR, n_domain_hits BIGINT, best_evalue DOUBLE, first_from BIGINT, last_to BIGINT)")
con.executemany("INSERT INTO cand VALUES (?,?,?,?,?,?,?,?)", rows)
os.makedirs(f'{bench}/cand', exist_ok=True)
con.execute(f"""
COPY (
  SELECT c.family, c.protein_id, c.locustag, c.domains, c.n_domain_hits, c.best_evalue, c.first_from, c.last_to,
         p.length AS prot_len, i.chrom, i.start AS gstart, i."end" AS gend, i.strand
  FROM cand c
  LEFT JOIN '{tables}/gene_proteins.parquet' p ON p.protein_id = c.protein_id
  LEFT JOIN '{tables}/gene_info.parquet' i ON i.gene_id = p.gene_id
  ORDER BY c.family, c.locustag, c.protein_id
) TO '{bench}/cand/candidates.tsv' (FORMAT CSV, DELIMITER '\t', HEADER)""")
for fam in FAM:
    n, g = con.execute(f"SELECT count(*), count(DISTINCT locustag) FROM cand WHERE family='{fam}'").fetchone()
    print(f'{fam}: {n} candidate proteins in {g} genomes')
