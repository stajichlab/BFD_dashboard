#!/usr/bin/env python3
"""Fixture test for miniprot_introns.py: a hand-made GFF and genome with known answers. Run: python3 test_miniprot_introns.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import miniprot_introns as mi

# genome: scaffold s1, 1300 nt of A with the splice dinucleotides below written in at known 1-based positions
g = ['A'] * 1300


def put(pos, s):                      # 1-based position
    for i, c in enumerate(s):
        g[pos - 1 + i] = c


put(201, 'GT'); put(259, 'AG')        # Q1 intron 1 [201,260]  plus strand  -> GT-AG
put(331, 'GC'); put(409, 'AG')        # Q1 intron 2 [331,410]  plus strand  -> GC-AG
put(701, 'CT'); put(779, 'AC')        # Q2 intron   [701,780]  minus strand -> donor revcomp('AC')=GT, acceptor revcomp('CT')=AG
put(1251, 'AT'); put(1289, 'AC')      # Q5 intron   [1251,1290] plus strand -> AT-AC (other)
genome = {'s1': ''.join(g)}

gff = '\n'.join([
    '##gff-version 3',
    's1\tminiprot\tmRNA\t101\t500\t900\t+\t.\tID=MP1;Rank=1;Identity=0.90;Positive=0.95;Target=Q1 1 100',
    's1\tminiprot\tCDS\t101\t200\t.\t+\t0\tParent=MP1;Rank=1',
    's1\tminiprot\tCDS\t261\t330\t.\t+\t0\tParent=MP1;Rank=1',
    's1\tminiprot\tCDS\t411\t500\t.\t+\t0\tParent=MP1;Rank=1',
    's1\tminiprot\tmRNA\t101\t500\t999\t+\t.\tID=MP1x;Rank=2;Identity=0.99;Positive=0.99;Target=Q1 1 100',
    's1\tminiprot\tCDS\t101\t150\t.\t+\t0\tParent=MP1x;Rank=2',
    's1\tminiprot\tCDS\t201\t250\t.\t+\t0\tParent=MP1x;Rank=2',
    's1\tminiprot\tCDS\t301\t350\t.\t+\t0\tParent=MP1x;Rank=2',
    's1\tminiprot\tmRNA\t601\t850\t400\t-\t.\tID=MP2;Rank=1;Identity=0.70;Positive=0.80;Target=Q2 1 60',
    's1\tminiprot\tCDS\t601\t700\t.\t-\t0\tParent=MP2;Rank=1',
    's1\tminiprot\tCDS\t781\t850\t.\t-\t0\tParent=MP2;Rank=1',
    's1\tminiprot\tmRNA\t901\t950\t200\t+\t.\tID=MP3;Rank=1;Identity=0.80;Positive=0.90;Target=Q3 1 16',
    's1\tminiprot\tCDS\t901\t950\t.\t+\t0\tParent=MP3;Rank=1',
    's1\tminiprot\tmRNA\t1001\t1200\t300\t+\t.\tID=MP4;Rank=1;Identity=0.85;Positive=0.9;Target=Q4 1 66',
    's1\tminiprot\tCDS\t1001\t1100\t.\t+\t0\tParent=MP4;Rank=1',
    's1\tminiprot\tCDS\t1103\t1200\t.\t+\t0\tParent=MP4;Rank=1',
    's1\tminiprot\tmRNA\t1201\t1250\t100\t+\t.\tID=MP5a;Rank=1;Identity=0.90;Positive=0.95;Target=Q5 1 16',
    's1\tminiprot\tCDS\t1201\t1250\t.\t+\t0\tParent=MP5a;Rank=1',
    's1\tminiprot\tmRNA\t1201\t1300\t300\t+\t.\tID=MP5b;Rank=1;Identity=0.80;Positive=0.90;Target=Q5 1 33',
    's1\tminiprot\tCDS\t1201\t1250\t.\t+\t0\tParent=MP5b;Rank=1',
    's1\tminiprot\tCDS\t1291\t1300\t.\t+\t0\tParent=MP5b;Rank=1',
])
qlens = {'Q1': 100, 'Q2': 80, 'Q3': 16, 'Q4': 66, 'Q5': 40}

alns = mi.parse_gff(gff.split('\n'), qlens)
assert len(alns) == 7, len(alns)
best = mi.best_per_query(alns)
assert set(best) == {'Q1', 'Q2', 'Q3', 'Q4', 'Q5'}
assert best['Q1']['id'] == 'MP1', 'Rank=2 alignment must be ignored even though its score is higher'
assert best['Q5']['id'] == 'MP5b', 'the higher-scoring of two Rank=1 alignments must win'

rows = {q: mi.alignment_row(best[q], genome) for q in best}
r = rows['Q1']
assert r['n_cds'] == 3 and r['n_introns'] == 2 and r['intron_lens'] == '60,80' and r['motifs'] == 'GT-AG,GC-AG', r
assert abs(r['qcov'] - 1.0) < 1e-9 and r['cds_nt'] == 100 + 70 + 90
r = rows['Q2']
assert r['strand'] == '-' and r['n_introns'] == 1 and r['intron_lens'] == '80' and r['motifs'] == 'GT-AG', r
assert abs(r['qcov'] - 0.75) < 1e-9
r = rows['Q3']
assert r['n_introns'] == 0 and r['n_short_gaps'] == 0 and r['motifs'] == ''
r = rows['Q4']
assert r['n_introns'] == 0 and r['n_short_gaps'] == 1, 'a 2 nt gap is a frameshift/short gap, not an intron'
r = rows['Q5']
assert r['n_introns'] == 1 and r['intron_lens'] == '40' and r['motifs'] == 'AT-AC'
# without a genome the motifs are NA but the counts are unchanged
assert mi.alignment_row(best['Q1'], None)['motifs'] == 'NA,NA'
# the intron length threshold is honoured
assert mi.alignment_row(best['Q5'], genome, min_intron=50)['n_introns'] == 0

# --- real miniprot 0.2 output format (lines copied from a run on a Coemansia assembly): no Rank attribute, Frameshift=true on
# some models, stop_codon lines, and the Target on CDS lines. Every alignment must count as primary; best score wins.
real = '\n'.join([
    '##gff-version 3',
    '##PAF\t103013at4751\t305\t12\t305\t-\tJAIXMV010000018.1\t445344\t113911\t114823\t456\t918\t0\tAS:i:734',
    'JAIXMV010000018.1\tminiprot\tmRNA\t113912\t114823\t734\t-\t.\tID=MP000001;Identity=0.4967;Positive=0.6373;Target=103013at4751 13 305',
    'JAIXMV010000018.1\tminiprot\tCDS\t113912\t114823\t0\t-\t0\tParent=MP000001;Target=103013at4751 13 305',
    'JAIXMV010000004.1\tminiprot\tmRNA\t265138\t311856\t685\t+\t.\tID=MP000002;Frameshift=true;Identity=0.3734;Positive=0.5328;Target=3299670at4751 99 521',
    'JAIXMV010000004.1\tminiprot\tCDS\t265138\t265590\t0\t+\t0\tParent=MP000002;Frameshift=true;Target=3299670at4751 99 214',
    'JAIXMV010000004.1\tminiprot\tCDS\t311018\t311856\t0\t+\t2\tParent=MP000002;Target=3299670at4751 215 521',
    'JAIXMV010000007.1\tminiprot\tmRNA\t409868\t412378\t1499\t+\t.\tID=MP000003;Identity=0.4057;Positive=0.5749;Target=407231at4751 1 814',
    'JAIXMV010000007.1\tminiprot\tCDS\t409868\t412375\t0\t+\t0\tParent=MP000003;Target=407231at4751 1 814',
    'JAIXMV010000007.1\tminiprot\tstop_codon\t412376\t412378\t0\t+\t0\tParent=MP000003',
    # a second, lower-scoring alignment of the same protein at another locus must lose to MP000003
    'JAIXMV010000009.1\tminiprot\tmRNA\t1000\t3000\t400\t+\t.\tID=MP000004;Identity=0.9;Positive=0.9;Target=407231at4751 1 814',
    'JAIXMV010000009.1\tminiprot\tCDS\t1000\t1500\t0\t+\t0\tParent=MP000004;Target=407231at4751 1 300',
    'JAIXMV010000009.1\tminiprot\tCDS\t2000\t3000\t0\t+\t0\tParent=MP000004;Target=407231at4751 301 814',
])
rq = {'103013at4751': 305, '3299670at4751': 521, '407231at4751': 814}
ra = mi.parse_gff(real.split('\n'), rq)
assert len(ra) == 4 and all(a['rank'] == 1 for a in ra), 'a missing Rank must count as primary'
rb = mi.best_per_query(ra)
assert set(rb) == set(rq) and rb['407231at4751']['id'] == 'MP000003'
r = mi.alignment_row(rb['3299670at4751'], None)
assert rb['3299670at4751']['frameshift'] is True
assert r['n_cds'] == 2 and r['n_introns'] == 1 and r['intron_lens'] == '45427', r
assert abs(r['qcov'] - (521 - 99 + 1) / 521) < 1e-9 and abs(r['identity'] - 0.3734) < 1e-9
assert mi.alignment_row(rb['103013at4751'], None)['n_introns'] == 0
assert mi.alignment_row(rb['407231at4751'], None)['n_cds'] == 1, 'stop_codon lines are not CDS blocks'
print('all miniprot_introns fixture tests passed')
