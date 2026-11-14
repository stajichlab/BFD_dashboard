#!/usr/bin/env python3
"""Step 7a: Pfam domain architecture of the annotated ligase candidates (BRCT, Lig4-specific domain, LIG3 BRCT and zinc finger).

Usage: 07a_ligase_architecture.py <BENCH> [tables_dir]
Reads trees/lig_members.tsv (annotated members) and the BFD Pfam table (read-only). Writes trees/lig_architecture.tsv with, per annotated
ligase protein, the ordered list of the informative domains (Pfam id:accession@start-end, domain_i_evalue <= 1e-3). Pfam names/accessions
come from the data (bfd.swissprot_subject shows them on the reference ligases): PF00533 BRCT, PF16589 BRCT_2, PF11411 DNA_ligase_IV (Lig4),
PF16759 LIG3_BRCT and PF00645 zf-PARP (LIG3), PF04675/PF01068/PF04679 ligase N/M/C. Models rescued by miniprot have no Pfam result.
"""
import csv, sys
import duckdb

import os
bench = sys.argv[1]
tdir = f"{bench}/trees_{os.environ.get('LAYER', 'annotated')}"
tables = sys.argv[2] if len(sys.argv) > 2 else '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables'
DOMS = ['PF00533', 'PF16589', 'PF11411', 'PF16759', 'PF00645', 'PF04675', 'PF01068', 'PF04679']
ids = [r['id'] for r in csv.DictReader(open(f'{tdir}/lig_members.tsv'), delimiter='\t') if r['kind'] == 'annotated']
locs = sorted({i.split('_')[0] for i in ids})
rng = ' OR '.join(f"(protein_id >= '{l}_' AND protein_id < '{l}`')" for l in locs)
con = duckdb.connect()
con.execute("SET threads=4")
rows = con.execute(f"""SELECT protein_id, split_part(pfam_acc,'.',1) AS acc, pfam_id, env_from, env_to
    FROM '{tables}/pfam.parquet' WHERE ({rng}) AND domain_i_evalue <= 1e-3 AND split_part(pfam_acc,'.',1) IN ({','.join(repr(d) for d in DOMS)})""").fetchall()
want = set(ids)
arch = {}
for pid, acc, name, s, e in rows:
    if pid in want:
        arch.setdefault(pid, []).append((s, e, f'{name}:{acc}'))
with open(f'{tdir}/lig_architecture.tsv', 'w') as f:
    f.write('protein_id\tdomains\tn_brct\thas_lig4_domain\thas_lig3_markers\n')
    for pid in ids:
        d = sorted(arch.get(pid, []))
        txt = ';'.join(f'{x[2]}@{x[0]}-{x[1]}' for x in d)
        n_brct = sum(1 for x in d if x[2].endswith(('PF00533', 'PF16589')))
        f.write(f"{pid}\t{txt}\t{n_brct}\t{any(x[2].endswith('PF11411') for x in d)}\t{any(x[2].endswith(('PF16759', 'PF00645')) for x in d)}\n")
print(len(ids), 'annotated ligases;', len(arch), 'with at least one informative domain')
