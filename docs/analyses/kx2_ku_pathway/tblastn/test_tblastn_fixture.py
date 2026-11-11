#!/usr/bin/env python3
"""Fixture test for summarize_tblastn.py and coemansia_strain_table.py (no BLAST needed).

Builds a tiny hits folder by hand with known best hits and checks the calls:
  strong  : e-value 1e-30, 60% of the query covered
  weak    : e-value 1e-7 (passes 1e-5, fails 1e-10)
  weak    : e-value 1e-30 but only 10% of the query covered (fails the 30% coverage rule, passes 1e-5)
  none    : no hit
  not_tested : the gene's only query comes from the genome itself (leave-one-out)
  merged  : two HSPs of one query on one scaffold, 20% + 20% (non-overlapping) = 40% coverage -> strong
Run:  python3 test_tblastn_fixture.py   (from this directory; exits non-zero on failure)
"""
import csv
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp()
os.makedirs(f'{tmp}/hits_round3')
# queries: gene key | accession; K1 has a self-genome query for genome G3
open(f'{tmp}/queries_round3.faa', 'w').write(
    '>GeneA|Q1\nMKV\n>GeneB|Q2\nMKV\n>GeneC|Q3\nMKV\n>GeneD|Q4\nMKV\n>GeneE|G3_000001-T1\nMKV\n>GeneF|Q6\nMKV\n')
genomes = [('G1', 'Coemansia', 'Kickxellales_core', '3', 'ASM1', '10', 'S1', '91.0', '200000'),
           ('G2', 'Coemansia', 'Kickxellales_core', '3', 'ASM2', '10', 'S2', '65.0', '7000'),
           ('G3', 'Dimargaris', 'Dimargaritales', '1', 'ASM3', '10', 'NA', 'NA', 'NA')]
with open(f'{tmp}/genomes.tsv', 'w') as f:
    f.write('locustag\tgenus\tgrp\tmin_priority\tasmid\tcontigs_expected\tstrain\tbusco_complete\tn50_bp\n')
    for g in genomes:
        f.write('\t'.join(g) + '\n')

# hits columns: qseqid sseqid pident length evalue bitscore qstart qend qlen sstart send frames
def hit(q, ev, qs, qe, qlen=100, scaf='sc1'):
    return f'{q}\t{scaf}\t40.0\t50\t{ev}\t100\t{qs}\t{qe}\t{qlen}\t1\t150\t+1\n'

# G1: GeneA strong (60% cover, 1e-30); GeneB weak (1e-7, 60%); GeneC weak (1e-30 but 10% cover); GeneD none (no hit);
#     GeneF strong from two HSPs 20% + 20% ... plus another 20% = 60%? use 2 x 20 = 40% -> strong
open(f'{tmp}/hits_round3/G1.tsv', 'w').write(
    hit('GeneA|Q1', '1e-30', 1, 60) + hit('GeneB|Q2', '1e-07', 1, 60) + hit('GeneC|Q3', '1e-30', 1, 10) +
    hit('GeneF|Q6', '1e-20', 1, 20) + hit('GeneF|Q6', '1e-15', 51, 70))
open(f'{tmp}/hits_round3/G2.tsv', 'w').write(hit('GeneA|Q1', '1e-30', 1, 60))
open(f'{tmp}/hits_round3/G3.tsv', 'w').write('')

subprocess.check_call([sys.executable, f'{HERE}/summarize_tblastn.py', tmp, f'{tmp}/genomes.tsv', f'{tmp}/round3', '_round3'])
calls = {(r['locustag'], r['gene']): r['call'] for r in csv.DictReader(open(f'{tmp}/round3_calls.tsv'), delimiter='\t')}
expected = {('G1', 'GeneA'): 'strong', ('G1', 'GeneB'): 'weak', ('G1', 'GeneC'): 'weak', ('G1', 'GeneD'): 'none',
            ('G1', 'GeneE'): 'none', ('G1', 'GeneF'): 'strong', ('G3', 'GeneE'): 'not_tested', ('G3', 'GeneA'): 'none',
            ('G2', 'GeneA'): 'strong', ('G2', 'GeneB'): 'none'}
bad = {k: (calls.get(k), v) for k, v in expected.items() if calls.get(k) != v}
assert not bad, f'unexpected calls (got, expected): {bad}'

# strain table: Ku/Lig4/... keys are not in the fixture; rename GeneA->Ku, GeneB->Lig4 to exercise the flags
calls_tsv = open(f'{tmp}/round3_calls.tsv').read().replace('\tGeneA\t', '\tKu\t').replace('\tGeneB\t', '\tLig4\t')
open(f'{tmp}/renamed_calls.tsv', 'w').write(calls_tsv)
subprocess.check_call([sys.executable, f'{HERE}/coemansia_strain_table.py', f'{tmp}/genomes.tsv', f'{tmp}/renamed_calls.tsv', f'{tmp}/strains.tsv'])
rows = {r['locustag']: r for r in csv.DictReader(open(f'{tmp}/strains.tsv'), delimiter='\t')}
assert set(rows) == {'G1', 'G2'}, 'only Coemansia strains expected'
assert rows['G1']['fragmented'] == 'no' and rows['G2']['fragmented'] == 'yes'
assert rows['G1']['Ku'] == 'strong' and rows['G1']['strong_or_weak_Ku'] == 'yes'
assert rows['G2']['Lig4'] == 'none' and 'Lig4' in rows['G2']['none_in_fragmented'].split(',')
assert rows['G1']['none_in_fragmented'] == '', 'a contiguous assembly must not be flagged'
print('fixture tests passed')
