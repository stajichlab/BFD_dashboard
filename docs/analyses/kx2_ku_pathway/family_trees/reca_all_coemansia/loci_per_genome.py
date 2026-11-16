#!/usr/bin/env python3
"""Distinct Rad51/Dmc1-like tblastn loci per genome (a per-genome version of ../../tblastn/reca_loci.py).

Usage: loci_per_genome.py <genomes_round3.tsv> <hits_round3_dir> <out.tsv>
Rule (same as reca_loci.py, fixed there): hits of the Rad51 (CTRL_Rad51) or Dmc1 (RecA_DMC1) query sets with e-value <= 1e-25 and >= 150
aligned residues; HSPs on the same scaffold within 20 kb are merged into one locus. For each locus the best bit score per query set is
reported, so a Dmc1-like second locus (Dmc1 queries score higher than Rad51 queries) can be told from a duplicated Rad51 locus.
"""
import csv
import os
import sys
from collections import defaultdict

EV, MINLEN, GAP = 1e-25, 150, 20000
KEYS = ('CTRL_Rad51', 'RecA_DMC1')


def loci_from_hits(rows):
    """rows: iterable of BLAST outfmt-6 style lists (qseqid sseqid pid aln mm go qs qe ss se evalue bits -> index as in the hit files:
    r[0] query, r[1] scaffold, r[3] aln length, r[4] e-value, r[5] bits, r[9], r[10] subject start/end). Returns list of loci dicts."""
    per = defaultdict(list)
    for r in rows:
        key = r[0].split('|')[0]
        ev, aln = float(r[4]), int(r[3])
        if key not in KEYS or ev > EV or aln < MINLEN:
            continue
        s, e = sorted((int(r[9]), int(r[10])))
        per[r[1]].append((s, e, key, float(r[5])))
    loci = []
    for scaf, iv in per.items():
        iv.sort()
        cur = None
        for s, e, k, b in iv:
            if cur is None or s > cur['end'] + GAP:
                if cur:
                    loci.append(cur)
                cur = {'scaf': scaf, 'start': s, 'end': e, 'bits': {'CTRL_Rad51': 0.0, 'RecA_DMC1': 0.0}}
            else:
                cur['end'] = max(cur['end'], e)
            cur['bits'][k] = max(cur['bits'][k], b)
        if cur:
            loci.append(cur)
    return loci


if __name__ == '__main__':
    genomes = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
    hits_dir = sys.argv[2]
    with open(sys.argv[3], 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['locustag', 'genus', 'grp', 'n_loci', 'locus_bits'])
        for g in genomes:
            p = f'{hits_dir}/{g["locustag"]}.tsv'
            if not os.path.exists(p):
                w.writerow([g['locustag'], g['genus'], g['grp'], '', 'no_hits_file'])
                continue
            loci = loci_from_hits(csv.reader(open(p), delimiter='\t'))
            bits = ' | '.join(f"R51={l['bits']['CTRL_Rad51']:.0f}/DMC={l['bits']['RecA_DMC1']:.0f}" for l in sorted(loci, key=lambda x: -max(x['bits'].values())))
            w.writerow([g['locustag'], g['genus'], g['grp'], len(loci), bits])
