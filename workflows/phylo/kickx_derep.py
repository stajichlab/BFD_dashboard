#!/usr/bin/env python3
"""Dereplicate the Coemansia strains of a kickx_tree taxa table by ANI (fastANI all-vs-all output).

Usage: kickx_derep.py <coemansia_taxa.tsv> <fastani.out> <outdir> [primary_threshold=99.0]
Pre-stated rule (workflows/phylo/README.md): two ingroup strains are redundant when fastANI ANI >= 99.0 % in either direction;
clusters are the connected components (single linkage) at that threshold. One representative per cluster: highest BUSCO
complete_pct, then largest N50. Outgroup rows are always kept. Sensitivity: the same clustering at 97 % and 95 %.
fastANI reports no row for pairs below about 80 % ANI or with too few mapped fragments; those pairs count as "not redundant".
Writes <outdir>/derep_clusters_<thr>.tsv (one row per ingroup strain), <outdir>/derep_taxa.tsv (the retained taxa, taxa-table format)
and <outdir>/derep_summary.txt.
"""
import csv
import os
import sys
from collections import defaultdict

taxa_tsv, ani_out, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
PRIMARY = float(sys.argv[4]) if len(sys.argv) > 4 else 99.0
THRESHOLDS = sorted({PRIMARY, 97.0, 95.0}, reverse=True)
os.makedirs(outdir, exist_ok=True)

rows = list(csv.DictReader(open(taxa_tsv), delimiter='\t'))
ingroup = [r for r in rows if r['role'] == 'ingroup']
by_asm = {r['ASMID']: r for r in ingroup}


def asm_of(path):
    base = os.path.basename(path)
    for suf in ('.fa.gz', '.fa', '.fasta', '.fna', '.gz'):
        if base.endswith(suf):
            base = base[: -len(suf)]
    return base


ani = {}
for line in open(ani_out):
    f = line.rstrip('\n').split('\t')
    if len(f) < 5:
        continue
    a, b = asm_of(f[0]), asm_of(f[1])
    if a == b or a not in by_asm or b not in by_asm:
        continue
    key = tuple(sorted((a, b)))
    ani[key] = max(ani.get(key, 0.0), float(f[2]))


def components(thr):
    parent = {a: a for a in by_asm}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for (a, b), v in ani.items():
        if v >= thr:
            parent[find(a)] = find(b)
    groups = defaultdict(list)
    for a in by_asm:
        groups[find(a)].append(a)
    return list(groups.values())


def rank(asm):
    r = by_asm[asm]
    return (float(r['complete_pct'] or 0), int(r['N50_bp'] or 0))


summary = [f'ingroup strains: {len(ingroup)}; pairs with a fastANI value: {len(ani)}']
retained = {}
for thr in THRESHOLDS:
    comps = components(thr)
    reps = {}
    with open(f'{outdir}/derep_clusters_{thr:g}.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['threshold', 'cluster', 'representative', 'ASMID', 'LOCUSTAG', 'STRAIN', 'complete_pct', 'N50_bp', 'cluster_size',
                    'min_ani_to_rep'])
        for i, comp in enumerate(sorted(comps, key=lambda c: (-len(c), c[0])), 1):
            rep = max(comp, key=rank)
            reps[i] = rep
            for a in sorted(comp, key=lambda x: (x != rep, x)):
                r = by_asm[a]
                key = tuple(sorted((a, rep)))
                w.writerow([thr, i, 'yes' if a == rep else 'no', a, r['LOCUSTAG'], r['STRAIN'], r['complete_pct'], r['N50_bp'],
                            len(comp), '' if a == rep else ani.get(key, '')])
    sizes = sorted((len(c) for c in comps), reverse=True)
    summary.append(f'threshold {thr:g}% : {len(comps)} clusters; sizes {sizes[:12]}{" ..." if len(sizes) > 12 else ""}; singletons {sizes.count(1)}')
    retained[thr] = set(reps.values())

prim = retained[PRIMARY]
with open(f'{outdir}/derep_taxa.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(list(rows[0].keys()))
    for r in rows:
        if r['role'] != 'ingroup' or r['ASMID'] in prim:
            w.writerow(list(r.values()))
summary.append(f'primary threshold {PRIMARY:g}%: {len(prim)} Coemansia representatives + {len(rows) - len(ingroup)} outgroup taxa retained')
vals = sorted(ani.values(), reverse=True)
summary.append('highest pairwise ANI values: ' + ', '.join(f'{v:.2f}' for v in vals[:10]))
summary.append(f'pairs >= 99: {sum(v >= 99 for v in vals)}; >= 97: {sum(v >= 97 for v in vals)}; >= 95: {sum(v >= 95 for v in vals)}; '
               f'>= 90: {sum(v >= 90 for v in vals)}')
open(f'{outdir}/derep_summary.txt', 'w').write('\n'.join(summary) + '\n')
print('\n'.join(summary))
