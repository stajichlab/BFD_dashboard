#!/usr/bin/env python3
"""Per-strain RecA-family table for the Coemansia strains: tree call next to the tblastn locus count and assembly quality.

Usage: strain_summary.py <reca_classification.tsv> <genomes.tsv> <loci_per_genome.tsv> <out.tsv>
  reca_classification.tsv  from ../07_classify.py (annotated layer)
  genomes.tsv              ../../tblastn/genomes_round3.tsv (locustag genus grp ... busco_complete n50_bp)
  loci_per_genome.tsv      from loci_per_genome.py (distinct Rad51/Dmc1-like tblastn loci)
Columns: n_reca = annotated RecA-fold proteins in the tree input for that genome (>= 200 aa); n_dmc1 / n_rad51 = proteins assigned
(UFBoot >= 90 and SH-aLRT >= 80) to those clades; n_dmc1_weak = proteins whose nearest clade is Dmc1 but unsupported or ambiguous;
dmc1_call = 'present' (assigned), 'uncertain' (weak only) or 'not_seen' (annotation only, so not a proof of absence).
"""
import csv
import sys
from collections import defaultdict


def summarize(cls_rows, genome_rows, loci_rows):
    """cls_rows: dict rows of reca_classification; genome_rows: dict rows of genomes.tsv; loci_rows: dict rows of loci table."""
    per = defaultdict(lambda: {'n_reca': 0, 'n_dmc1': 0, 'n_rad51': 0, 'n_dmc1_weak': 0, 'n_other_assigned': 0})
    for r in cls_rows:
        if r['call'] == 'dropped_by_trim':
            per[r['locustag']]['n_dropped'] = per[r['locustag']].get('n_dropped', 0) + 1
            continue
        d = per[r['locustag']]
        d['n_reca'] += 1
        if r['call'] == 'assigned':
            if r['clade'] == 'Dmc1':
                d['n_dmc1'] += 1
            elif r['clade'] == 'Rad51':
                d['n_rad51'] += 1
            else:
                d['n_other_assigned'] += 1
        elif 'Dmc1' in (r.get('mixed_labels') or '') or r['clade'] == 'Dmc1':
            d['n_dmc1_weak'] += 1
    loci = {r['locustag']: r for r in loci_rows}
    out = []
    for g in genome_rows:
        loc = g['locustag']
        d = per.get(loc, {'n_reca': 0, 'n_dmc1': 0, 'n_rad51': 0, 'n_dmc1_weak': 0, 'n_other_assigned': 0})
        call = 'present' if d['n_dmc1'] > 0 else ('uncertain' if d['n_dmc1_weak'] > 0 else 'not_seen')
        lr = loci.get(loc, {})
        out.append({'locustag': loc, 'genus': g['genus'], 'group': g['grp'], 'busco_complete': g.get('busco_complete', ''),
                    'n50_bp': g.get('n50_bp', ''), 'n_reca': d['n_reca'], 'n_dropped_by_trim': d.get('n_dropped', 0),
                    'n_rad51': d['n_rad51'], 'n_dmc1': d['n_dmc1'], 'n_dmc1_weak': d['n_dmc1_weak'], 'dmc1_call': call,
                    'tblastn_loci': lr.get('n_loci', ''), 'tblastn_locus_bits_R51_DMC': lr.get('locus_bits', '')})
    return out


if __name__ == '__main__':
    cls = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
    gen = list(csv.DictReader(open(sys.argv[2]), delimiter='\t'))
    loci = list(csv.DictReader(open(sys.argv[3]), delimiter='\t'))
    rows = summarize(cls, gen, loci)
    with open(sys.argv[4], 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print(f'{len(rows)} genomes written')
