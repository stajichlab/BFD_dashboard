#!/usr/bin/env python3
"""T5 step 3: write the GH65 candidate proteins as FASTA, read from the BFD proteomes (read-only).

Usage: 03_extract_proteins.py <BENCH>
Reads $BENCH/cand/candidates.tsv (step 2) and the atlas locustag -> proteome map. Writes $BENCH/cand/gh65_annotated.faa with headers
>ANN|<protein_id>|<locustag>. A protein that is not found in its proteome is an error (the run stops).
"""
import csv
import gzip
import sys
from collections import defaultdict

MAP = '/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/locustag_to_pep.tsv'


def read_wanted(path, ids):
    """Sequences of the wanted protein ids from one (optionally gzipped) proteome FASTA."""
    op = gzip.open if str(path).endswith('.gz') else open
    seqs, name, buf = {}, None, []
    with op(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name in ids:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name in ids:
        seqs[name] = ''.join(buf)
    return seqs


def main(bench):
    pep = {r['LOCUSTAG']: r['pep_path'] for r in csv.DictReader(open(MAP), delimiter='\t')}
    want = defaultdict(set)
    for r in csv.DictReader(open(f'{bench}/cand/candidates.tsv'), delimiter='\t'):
        want[r['locustag']].add(r['protein_id'])
    n = 0
    with open(f'{bench}/cand/gh65_annotated.faa', 'w') as out:
        for loc in sorted(want):
            seqs = read_wanted(pep[loc], want[loc])
            missing = [i for i in want[loc] if i not in seqs]
            if missing:
                raise SystemExit(f'{loc}: {len(missing)} proteins not found in {pep[loc]}, e.g. {missing[0]}')
            for pid in sorted(want[loc]):
                out.write(f">ANN|{pid}|{loc}\n{seqs[pid].rstrip('*')}\n")
                n += 1
    print(f'{n} proteins written from {len(want)} genomes')


if __name__ == '__main__':
    main(sys.argv[1])
