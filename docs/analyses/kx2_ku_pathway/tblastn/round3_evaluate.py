#!/usr/bin/env python3
"""Evaluate tblastn round 3: positive controls per gene, changes against round 2, and the Coemansia strain summary.

Usage: round3_evaluate.py <results dir containing round2_calls.tsv and round3_calls.tsv> <genomes_round3.tsv>

Controls (fixed before the round 3 results were seen) are the target genomes in which the Q1 Pfam calls
(bfd_scratch.kx2_ku_calls, present = true) found the gene. A gene's negatives are treated as informative only if at least 75% of its
control genomes are called 'strong' (a control is judged with its own protein left out). With a single control, one strong call is
reported as 'one control only' and a miss as 'uninformative'.
"""
import csv
import sys
from collections import Counter, defaultdict

res, genomes_tsv = sys.argv[1], sys.argv[2]
CONTROLS = {
    'Ku': ['FCDD8027', 'F0AFF821', 'F1C537C3', 'F1C9E2AA', 'F7FE4923', 'F8E635D0'],
    'Pol4': ['FCDD8027', 'F8E635D0'],
    'Lig4': ['FCDD8027'],
    'Aprataxin': ['FCDD8027'],
    'DNA-PKcs': ['F8E635D0'],
    'XLF_Nej1': ['F8E635D0'],
    'Sae2_CtIP': ['FCDD8027', 'F3D937C5', 'F0AFF821', 'F1C537C3', 'F1C9E2AA', 'F45BD9C7', 'F7FE4923', 'F8E635D0', 'FCCDAB75', 'F4FD0401'],
    'ATR_Mec1': ['FCDD8027', 'FDAC0301', 'F3F252DC', 'F1C537C3', 'F1C9E2AA', 'F9F51274', 'F18CA942', 'FC345CF7', 'FEB493EF', 'F7FE4923',
                 'F8E635D0', 'F4FD0401'],
}


def read(path):
    return list(csv.DictReader(open(path), delimiter='\t'))


r3 = read(f'{res}/round3_calls.tsv')
r2 = read(f'{res}/round2_calls.tsv')
c3 = {(r['locustag'], r['gene']): r for r in r3}
c2 = {(r['locustag'], r['gene']): r for r in r2}
genomes = {r['locustag']: r for r in read(genomes_tsv)}

print('== 1. Positive controls (round 3)')
for gene, locs in CONTROLS.items():
    calls = [(genomes[l]['genus'], c3[(l, gene)]['call'], c3[(l, gene)].get('best_evalue', '')) for l in locs if (l, gene) in c3]
    n_strong = sum(1 for c in calls if c[1] == 'strong')
    n_tested = sum(1 for c in calls if c[1] != 'not_tested')
    frac = n_strong / n_tested if n_tested else 0
    verdict = ('controls pass' if frac >= 0.75 and n_tested >= 2 else 'one control only (strong)' if n_tested == 1 and n_strong == 1
               else 'UNINFORMATIVE (controls fail)')
    print(f'  {gene:10s} {n_strong}/{n_tested} strong -> {verdict}')
    print('      ' + '; '.join(f'{g}:{c}' for g, c, _ in calls))

print('\n== 2. Rad51 control (CTRL_Rad51) strong fraction')
for label, pred in [('all 128 genomes', lambda l: True), ('107 Coemansia strains', lambda l: genomes[l]['genus'] == 'Coemansia'),
                    ('the other 21 genomes', lambda l: genomes[l]['genus'] != 'Coemansia')]:
    sel = [c3[(l, 'CTRL_Rad51')]['call'] for l in genomes if pred(l) and (l, 'CTRL_Rad51') in c3]
    cnt = Counter(sel)
    print(f'  {label}: n={len(sel)} {dict(cnt)}')

print('\n== 3. Changes from round 2 to round 3 in the original 31 genomes (calls differ)')
orig = [l for l, g in genomes.items() if g['min_priority'] in ('1', '2')]
changes = []
for l in orig:
    for gene in ['Ku', 'DNA-PKcs', 'Lig4', 'Pol4', 'XLF_Nej1', 'Aprataxin', 'ATR_Mec1']:
        a, b = c2.get((l, gene)), c3.get((l, gene))
        if a and b and a['call'] != b['call']:
            changes.append((gene, genomes[l]['genus'], l, a['call'], b['call'], b.get('query', ''), b.get('best_evalue', ''), b.get('query_cov', '')))
for ch in sorted(changes):
    print('  %-10s %-16s %s  %s -> %s  (best query %s, e=%s, cov=%s)' % ch)
print(f'  total changed calls: {len(changes)}; by direction: ' +
      str(Counter((c[3], c[4]) for c in changes)))

print('\n== 4. Coemansia strains (107): counts per gene')
coe = [l for l, g in genomes.items() if g['genus'] == 'Coemansia']
for gene in ['Ku', 'Lig4', 'Pol4', 'Sae2_CtIP', 'DNA-PKcs', 'ATR_Mec1', 'XLF_Nej1', 'Aprataxin', 'RecA_DMC1', 'RecA_RAD55', 'RecA_RAD57']:
    cnt = Counter(c3[(l, gene)]['call'] for l in coe if (l, gene) in c3)
    print(f'  {gene:11s} {dict(cnt)}')
ku_hits = [(l, genomes[l]['strain'], c3[(l, 'Ku')]) for l in coe if c3[(l, 'Ku')]['call'] in ('strong', 'weak')]
print(f'\n  Coemansia strains with any Ku hit: {len(ku_hits)}')
for l, strain, r in ku_hits:
    print(f"    {l} {strain}: {r['call']} query={r['query']} e={r['best_evalue']} cov={r['query_cov']} scaffold={r['scaffold']}")
