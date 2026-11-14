#!/usr/bin/env python3
"""Subset the nf_phyling per-marker alignments to a taxa list and rebuild the concatenated alignment.

Usage: kickx_subset_aln.py <results/protein dir of a kickx_tree run> <taxa.tsv (kickx taxa table or list of LOCUSTAGs)> <outdir> [min_taxa=8]
Reads  <run>/filter/fungi_odb12-filter/<marker>.aa.mfa (the alignments that nf_phyling concatenated; header ">LOCUSTAG protein_id")
and the marker order of <run>/buildtree/*.partition.
Writes <outdir>/genes/<marker>.fa  (alignment restricted to the chosen taxa; all-gap columns removed; markers with fewer than
min_taxa sequences dropped), <outdir>/concat.faa (taxa without a marker get gaps, as in phyling), <outdir>/concat.partition
(IQ-TREE/RAxML style, 1-based), <outdir>/taxa.txt and <outdir>/subset_summary.txt.
Nothing is re-aligned or re-trimmed: this only removes taxa and columns that became all gaps.
"""
import csv
import glob
import os
import sys

run, taxa_file, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
MIN_TAXA = int(sys.argv[4]) if len(sys.argv) > 4 else 8
os.makedirs(f'{outdir}/genes', exist_ok=True)

first = open(taxa_file).readline()
if first.startswith('role\t'):
    taxa = [r['LOCUSTAG'] for r in csv.DictReader(open(taxa_file), delimiter='\t')]
else:
    taxa = [l.split()[0] for l in open(taxa_file) if l.strip()]
tset = set(taxa)

part = glob.glob(f'{run}/buildtree/*.partition')[0]
markers = [l.split(',')[1].split('=')[0] for l in open(part) if l.strip()]


def read_fa(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


concat = {t: [] for t in taxa}
lines_part, pos = [], 0
dropped, kept_len, ntax = [], 0, []
for m in markers:
    seqs = read_fa(f'{run}/filter/fungi_odb12-filter/{m}')
    sub = {t: s for t, s in seqs.items() if t in tset}
    if len(sub) < MIN_TAXA:
        dropped.append((m, len(sub)))
        continue
    L = len(next(iter(sub.values())))
    keep = [i for i in range(L) if any(s[i] != '-' for s in sub.values())]
    if not keep:
        dropped.append((m, 0))
        continue
    with open(f'{outdir}/genes/{m[:-len(".aa.mfa")]}.fa', 'w') as f:
        for t in taxa:
            if t in sub:
                f.write(f'>{t}\n{"".join(sub[t][i] for i in keep)}\n')
    for t in taxa:
        concat[t].append(''.join(sub[t][i] for i in keep) if t in sub else '-' * len(keep))
    lines_part.append(f'LG, {m[:-len(".aa.mfa")]} = {pos + 1}-{pos + len(keep)}')
    pos += len(keep)
    kept_len += len(keep)
    ntax.append(len(sub))
with open(f'{outdir}/concat.faa', 'w') as f:
    for t in taxa:
        f.write(f'>{t}\n{"".join(concat[t])}\n')
open(f'{outdir}/concat.partition', 'w').write('\n'.join(lines_part) + '\n')
open(f'{outdir}/taxa.txt', 'w').write('\n'.join(taxa) + '\n')
gap = {t: sum(c == '-' for c in ''.join(concat[t])) / max(1, pos) for t in taxa}
worst = sorted(gap.items(), key=lambda x: -x[1])[:5]
summ = [f'taxa {len(taxa)}; markers kept {len(ntax)} of {len(markers)}; dropped {len(dropped)} (fewer than {MIN_TAXA} taxa); concatenated sites {pos}',
        f'taxa per kept marker: min {min(ntax)}, median {sorted(ntax)[len(ntax)//2]}, max {max(ntax)}',
        f'gap fraction per taxon: median {sorted(gap.values())[len(gap)//2]:.3f}; highest ' + ', '.join(f'{t} {g:.2f}' for t, g in worst)]
open(f'{outdir}/subset_summary.txt', 'w').write('\n'.join(summ) + '\n')
print('\n'.join(summ))
