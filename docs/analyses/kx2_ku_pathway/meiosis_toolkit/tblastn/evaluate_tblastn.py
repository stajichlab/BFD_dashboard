#!/usr/bin/env python3
"""Evaluate the meiotic-toolkit tblastn calls: positive-control pass per gene and per-group strong / weak / none counts.

Usage: evaluate_tblastn.py <calls.tsv from summarize_tblastn.py> <positives.tsv> <out prefix>
Rule fixed before any result was seen (README): a gene is INFORMATIVE only if tblastn recovers it (strong or weak) in at least 75% of the
genomes where the annotation already calls it present in the reference groups (Harpellales, Mucoromycota, Root, Zoopagales, Entomophthorales
+ Basidiobolus, Dimargaritales); the list of those genomes is positives.tsv (from bfd_scratch.kx2_mt_matrix). A strong recovery rate is
reported as well. The mitotic controls CTRL_Msh2, CTRL_Mlh1, CTRL_Rad21 and CTRL_Rad51 (genes that must be present everywhere) are
evaluated on all genomes of the reference groups. Genes with only hit evidence (Msh4, Hop1, Rec8, Mlh3: paralogs cross-hit) are reported
without a pass/fail.
Outputs: <prefix>_control_pass.tsv and <prefix>_group_counts.tsv.
"""
import csv, sys
from collections import defaultdict

calls_tsv, pos_tsv, prefix = sys.argv[1], sys.argv[2], sys.argv[3]
REF = {'Harpellales', 'Mucoromycota_s.l.', 'Root', 'Zoopagales', 'Entomophthorales_Basidiobolus', 'Dimargaritales'}
rows = list(csv.DictReader(open(calls_tsv), delimiter='\t'))
call = {(r['gene'], r['locustag']): r['call'] for r in rows}
grp = {r['locustag']: r['grp'] for r in rows}
pos = {r['key']: set(r['annotated_positive_reference_genomes'].split(',')) for r in csv.DictReader(open(pos_tsv), delimiter='\t')}
for ctrl in ('CTRL_Msh2', 'CTRL_Mlh1', 'CTRL_Rad21', 'CTRL_Rad51'):
    pos[ctrl] = {l for l, g in grp.items() if g in REF}
keys = sorted({r['gene'] for r in rows if r['gene'] != 'NA'})
out = []
for k in keys:
    if k in pos:
        tested = [l for l in pos[k] if call.get((k, l)) in ('strong', 'weak', 'none')]
        n = len(tested)
        rec = sum(1 for l in tested if call[(k, l)] in ('strong', 'weak'))
        st = sum(1 for l in tested if call[(k, l)] == 'strong')
        frac = rec / n if n else 0.0
        out.append((k, n, rec, st, round(100 * frac, 1), round(100 * st / n, 1) if n else 0.0,
                    'PASS (informative)' if n and frac >= 0.75 else 'FAIL (uninformative)'))
    else:
        out.append((k, 0, 0, 0, '', '', 'hit evidence only (paralogs cross-hit)'))
with open(prefix + '_control_pass.tsv', 'w') as f:
    f.write('key\tn_positive_genomes_tested\trecovered_strong_or_weak\trecovered_strong\tpct_recovered\tpct_strong\tstatus\n')
    for o in out:
        f.write('\t'.join(str(x) for x in o) + '\n')
tab = defaultdict(lambda: defaultdict(int))
for r in rows:
    if r['gene'] != 'NA':
        tab[(r['gene'], r['grp'])][r['call']] += 1
with open(prefix + '_group_counts.tsv', 'w') as f:
    f.write('key\tgroup\tn_genomes\tstrong\tweak\tnone\tnot_tested\n')
    for (k, g), c in sorted(tab.items()):
        f.write(f"{k}\t{g}\t{sum(c.values())}\t{c['strong']}\t{c['weak']}\t{c['none']}\t{c['not_tested']}\n")
print('wrote', prefix + '_control_pass.tsv', prefix + '_group_counts.tsv')
