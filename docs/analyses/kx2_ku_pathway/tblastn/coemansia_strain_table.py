#!/usr/bin/env python3
"""Per-Coemansia-strain table of tblastn calls with assembly quality (round 3).

Usage: coemansia_strain_table.py <genomes_round3.tsv> <round3_calls.tsv> <out.tsv>

One row per Coemansia strain (genus == 'Coemansia'): the call (strong / weak / none / not_tested) for each gene, the strain name,
BUSCO complete, N50 and expected scaffold count. 'none_in_fragmented' lists genes called none in a low-contiguity assembly
(N50 < 50 kb, the atlas Q-n50 threshold, or BUSCO complete < 70), where a missing region of the assembly could explain the absence
and the 'none' is therefore weaker evidence. 'strong_or_weak_Ku' marks strains with any Ku hit so they can be inspected.
"""
import csv
import sys

genomes_tsv, calls_tsv, out_tsv = sys.argv[1:4]
GENES = ['Ku', 'Lig4', 'Pol4', 'Sae2_CtIP', 'DNA-PKcs', 'ATR_Mec1', 'CTRL_Rad51']
N50_MIN, BUSCO_MIN = 50000, 70.0

meta = {}
for r in csv.DictReader(open(genomes_tsv), delimiter='\t'):
    if r['genus'] == 'Coemansia':
        meta[r['locustag']] = r
calls = {}
for r in csv.DictReader(open(calls_tsv), delimiter='\t'):
    if r['locustag'] in meta and r['gene'] in GENES:
        calls[(r['locustag'], r['gene'])] = r

cols = ['locustag', 'strain', 'busco_complete', 'n50_bp', 'contigs_expected', 'fragmented'] + GENES + ['none_in_fragmented', 'strong_or_weak_Ku']
with open(out_tsv, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(cols)
    for loc, m in sorted(meta.items(), key=lambda kv: (kv[1]['strain'])):
        n50 = int(m['n50_bp']) if m['n50_bp'] not in ('NA', '') else None
        busco = float(m['busco_complete']) if m['busco_complete'] not in ('NA', '') else None
        frag = (n50 is not None and n50 < N50_MIN) or (busco is not None and busco < BUSCO_MIN)
        row_calls = [calls.get((loc, g), {}).get('call', 'NOT_RUN') for g in GENES]
        none_frag = [g for g, c in zip(GENES, row_calls) if c == 'none' and frag and g != 'CTRL_Rad51']
        ku = row_calls[0]
        w.writerow([loc, m['strain'], m['busco_complete'], m['n50_bp'], m['contigs_expected'], 'yes' if frag else 'no'] + row_calls +
                   [','.join(none_frag), 'yes' if ku in ('strong', 'weak') else 'no'])
print('wrote', out_tsv, 'for', len(meta), 'Coemansia strains')
