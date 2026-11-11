#!/usr/bin/env python3
"""Summarize tblastn hits per (genome, gene) for the Ku-pathway absence check.

Usage: summarize_tblastn.py <BENCH folder with hits/> <genomes.tsv> <out prefix>

For each genome and gene key (a gene may have several query proteins) the best query is the one with the highest total
non-overlapping aligned query span. Calls (thresholds chosen before looking at the results, see tblastn/README.md):
  strong : best e-value <= 1e-10 and at least 30% of the query covered by merged HSPs on one scaffold
  weak   : best e-value <= 1e-5 but not strong
  none   : otherwise
Coverage is the union of query intervals over all HSPs of that query on its best-scoring scaffold, divided by query length, so a
gene split by introns (several HSPs) is counted once. Query coverage is relative to the query, not the hit.
"""
import csv, os, sys
from collections import defaultdict

bench, genomes_tsv, prefix = sys.argv[1], sys.argv[2], sys.argv[3]
genomes = list(csv.DictReader(open(genomes_tsv), delimiter='\t'))


def merged_span(iv):
    iv = sorted(iv)
    tot, cur_s, cur_e = 0, None, None
    for s, e in iv:
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                tot += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        tot += cur_e - cur_s + 1
    return tot


rows = []
for g in genomes:
    loc = g['locustag']
    path = f'{bench}/hits/{loc}.tsv'
    if not os.path.exists(path):
        rows.append({'locustag': loc, 'genus': g['genus'], 'grp': g['grp'], 'gene': 'NA', 'call': 'NOT_RUN'})
        continue
    by_q = defaultdict(list)
    for r in csv.reader(open(path), delimiter='\t'):
        q, scaf, pid, aln, ev, bits, qs, qe, qlen = r[:9]
        by_q[q].append((scaf, float(pid), int(aln), float(ev), float(bits), int(qs), int(qe), int(qlen)))
    per_gene = {}
    for q, hs in by_q.items():
        key, acc = q.split('|')
        best_scaf = max({h[0] for h in hs}, key=lambda s: sum(h[4] for h in hs if h[0] == s))
        sub = [h for h in hs if h[0] == best_scaf]
        qlen = sub[0][7]
        cov = merged_span([(min(h[5], h[6]), max(h[5], h[6])) for h in sub]) / qlen
        best_ev = min(h[3] for h in sub)
        best_bits = max(h[4] for h in sub)
        pid = max(sub, key=lambda h: h[4])[1]
        rec = {'query': acc, 'scaffold': best_scaf, 'best_evalue': best_ev, 'best_bits': best_bits, 'pident_best_hsp': pid,
               'query_cov': round(cov, 3), 'n_hsps': len(sub)}
        if key not in per_gene or (rec['query_cov'] * rec['best_bits']) > (per_gene[key]['query_cov'] * per_gene[key]['best_bits']):
            per_gene[key] = rec
    keys = sorted({q.split('|')[0] for q in by_q} | set())
    allkeys = ['Ku', 'DNA-PKcs', 'Lig4', 'Pol4', 'XLF_Nej1', 'Aprataxin', 'ATR_Mec1', 'Sae2_CtIP', 'CTRL_Rad51']
    for key in allkeys:
        rec = per_gene.get(key)
        if rec is None:
            rows.append({'locustag': loc, 'genus': g['genus'], 'grp': g['grp'], 'gene': key, 'call': 'none', 'query': '',
                         'scaffold': '', 'best_evalue': '', 'best_bits': '', 'pident_best_hsp': '', 'query_cov': 0, 'n_hsps': 0})
            continue
        if rec['best_evalue'] <= 1e-10 and rec['query_cov'] >= 0.30:
            call = 'strong'
        elif rec['best_evalue'] <= 1e-5:
            call = 'weak'
        else:
            call = 'none'
        rows.append({'locustag': loc, 'genus': g['genus'], 'grp': g['grp'], 'gene': key, 'call': call, **rec})

cols = ['locustag', 'genus', 'grp', 'gene', 'call', 'query', 'scaffold', 'best_evalue', 'best_bits', 'pident_best_hsp', 'query_cov', 'n_hsps']
with open(prefix + '_calls.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore')
    w.writeheader()
    w.writerows(rows)

# per group x gene table: strong / weak / none counts
tab = defaultdict(lambda: defaultdict(int))
for r in rows:
    if r['gene'] != 'NA':
        tab[(r['grp'], r['gene'])][r['call']] += 1
with open(prefix + '_by_group.tsv', 'w') as f:
    f.write('grp\tgene\tn_genomes\tstrong\tweak\tnone\n')
    for (grp, gene), c in sorted(tab.items()):
        f.write(f"{grp}\t{gene}\t{sum(c.values())}\t{c['strong']}\t{c['weak']}\t{c['none']}\n")
print('wrote', prefix + '_calls.tsv', 'and', prefix + '_by_group.tsv')
