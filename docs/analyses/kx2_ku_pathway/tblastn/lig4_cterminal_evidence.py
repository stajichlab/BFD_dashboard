#!/usr/bin/env python3
"""Evidence only: do tblastn hits of the Lig4 queries reach the C-terminal third of the query (where Lig4's BRCT region lies)?

Usage: lig4_cterminal_evidence.py <BENCH folder with hits_round3/ and queries_round3.faa> <genomes_round3.tsv> <out.tsv>

Lig4 and LIG1 share the ATP-dependent ligase core in the N-terminal part of the protein; the C-terminal region (BRCT repeats and
linker) is specific to Lig4. A hit that covers only the N-terminal core cannot tell the two apart, so this reports, per genome, the
best fraction of the C-terminal third (residues from 60% of the query length to the end) of any Lig4 query that is covered by merged
HSPs with e-value <= 1e-5 on one scaffold, and the best N-terminal-half coverage for comparison. It does not say whether a hit is
Lig4: that needs the gene-family trees. Self-genome queries are excluded as in the search.
"""
import csv
import os
import sys
from collections import defaultdict

bench, genomes_tsv, out_tsv = sys.argv[1:4]
genomes = list(csv.DictReader(open(genomes_tsv), delimiter='\t'))
qlen = {}
name = None
for line in open(f'{bench}/queries_round3.faa'):
    if line.startswith('>'):
        name = line[1:].strip()
        qlen[name] = 0
    else:
        qlen[name] += len(line.strip())


def merged(iv):
    iv = sorted(iv)
    tot, s0, e0 = 0, None, None
    for s, e in iv:
        if e0 is None or s > e0 + 1:
            if e0 is not None:
                tot += e0 - s0 + 1
            s0, e0 = s, e
        else:
            e0 = max(e0, e)
    return tot + (e0 - s0 + 1 if e0 is not None else 0)


rows = []
for g in genomes:
    loc = g['locustag']
    per_q = defaultdict(lambda: defaultdict(list))   # query -> scaffold -> [(qs, qe)]
    for l in open(f'{bench}/hits_round3/{loc}.tsv'):
        f = l.rstrip().split('\t')
        if not f[0].startswith('Lig4|') or float(f[4]) > 1e-5:
            continue
        per_q[f[0]][f[1]].append((min(int(f[6]), int(f[7])), max(int(f[6]), int(f[7]))))
    best_c, best_n, best_q = 0.0, 0.0, ''
    for q, scafs in per_q.items():
        L = qlen[q]
        cstart, nend = int(0.6 * L), int(0.5 * L)
        for scaf, iv in scafs.items():
            c = merged([(max(s, cstart), e) for s, e in iv if e >= cstart]) / max(1, L - cstart)
            n = merged([(s, min(e, nend)) for s, e in iv if s <= nend]) / max(1, nend)
            if c > best_c:
                best_c, best_q = c, q
            best_n = max(best_n, n)
    rows.append({'locustag': loc, 'genus': g['genus'], 'grp': g['grp'], 'strain': g.get('strain', ''), 'cterm_cov_best': round(best_c, 3),
                 'nterm_half_cov_best': round(best_n, 3), 'cterm_query': best_q})
with open(out_tsv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)
summ = defaultdict(list)
for r in rows:
    summ[r['grp']].append(r['cterm_cov_best'])
for grp, v in sorted(summ.items()):
    print(f'{grp:32s} n={len(v):3d}  C-terminal third covered >= 0.25: {sum(x >= 0.25 for x in v):3d}   median {sorted(v)[len(v) // 2]:.2f}  max {max(v):.2f}')
print('wrote', out_tsv)
