#!/usr/bin/env python3
"""Step 3: write the annotated candidate proteins of each family as FASTA, read from the BFD proteomes (read-only).

Usage: 03_extract_proteins.py <BENCH>
Reads $BENCH/cand/candidates.tsv (step 2) and the locustag -> proteome map used by the atlas phylogeny jobs. Writes
$BENCH/cand/<family>_annotated.faa with headers >ANN|<protein_id>|<locustag>. A protein that is not found in its proteome is an error.
"""
import csv, gzip, os, sys
from collections import defaultdict

bench = sys.argv[1]
MAP = '/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/locustag_to_pep.tsv'
pep = {r['LOCUSTAG']: r['pep_path'] for r in csv.DictReader(open(MAP), delimiter='\t')}
want = defaultdict(lambda: defaultdict(set))          # locustag -> protein_id -> {families}
for r in csv.DictReader(open(f'{bench}/cand/candidates.tsv'), delimiter='\t'):
    want[r['locustag']][r['protein_id']].add(r['family'])

seqs = {}
for loc, ids in want.items():
    name, buf = None, []
    with gzip.open(pep[loc], 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name in ids:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
        if name in ids:
            seqs[name] = ''.join(buf)
    missing = [i for i in ids if i not in seqs]
    if missing:
        raise SystemExit(f'{loc}: {len(missing)} proteins not found in {pep[loc]}, e.g. {missing[0]}')

out = {fam: open(f'{bench}/cand/{fam}_annotated.faa', 'w') for fam in ('lig', 'reca', 'polx')}
n = defaultdict(int)
for loc, ids in want.items():
    for pid, fams in ids.items():
        for fam in fams:
            s = seqs[pid].rstrip('*')
            out[fam].write(f'>ANN|{pid}|{loc}\n{s}\n')
            n[fam] += 1
for f in out.values():
    f.close()
print(dict(n))
