import os, csv, glob, sys, statistics
tree = sys.argv[1] if len(sys.argv) > 1 else 'backbone'
B = os.path.realpath(os.path.expanduser('~/projects/BFD_dashboard/bench/kickx_tree'))
taxa = list(csv.DictReader(open(f'{B}/taxa/{tree}_taxa.tsv'), delimiter='\t'))
A = f'{B}/{tree}/results/protein/align/fungi_odb12-align'
F = f'{B}/{tree}/results/protein/filter/fungi_odb12-filter'

def marker_counts(d):
    c = {}
    for f in glob.glob(f'{d}/*.mfa'):
        for l in open(f):
            if l.startswith('>'):
                n = l[1:].split()[0]
                c[n] = c.get(n, 0) + 1
    return c

n_align = len(glob.glob(f'{A}/*.mfa'))
ca = marker_counts(A)
cf = marker_counts(F) if os.path.isdir(F) else {}
rows = []
for t in taxa:
    k = t['LOCUSTAG']
    rows.append((ca.get(k, 0), cf.get(k, 0), float(t['duplicated_pct']), int(t['n_proteins']), k, t['GENUS'], t['role'],
                 'low_busco' if t.get('low_busco') == 'true' else '', 'high_dup' if t.get('high_dup') == 'true' else ''))
rows.sort()
print(f'tree={tree} taxa={len(taxa)} align_markers_total={n_align} taxa_with_any_marker={sum(1 for r in rows if r[0] > 0)} filter_dir={"yes" if cf else "not yet"}')
print('align_mk filt  dup%  n_prot  locustag genus role flags   (lowest 14)')
for r in rows[:14]:
    print('%6d %4d %5.1f %7d  %s %-15s %-11s %s %s' % r)
a = [r[0] for r in rows]
print('align markers per taxon: min %d, 10th pct %d, median %d' % (min(a), sorted(a)[len(a) // 10], statistics.median(a)))
d = [r[2] for r in rows]
def rank(x):
    s = sorted(range(len(x)), key=lambda i: x[i]); r = [0] * len(x)
    for p, i in enumerate(s):
        r[i] = p
    return r
ra, rd = rank(a), rank(d)
m = len(a); ma = sum(ra) / m; md = sum(rd) / m
cov = sum((x - ma) * (y - md) for x, y in zip(ra, rd))
den = (sum((x - ma) ** 2 for x in ra) * sum((y - md) ** 2 for y in rd)) ** 0.5
print('Spearman(align markers, BUSCO duplicated%%) = %.2f' % (cov / den))
print('flagged taxa:')
for r in rows:
    if r[7] or r[8]:
        print('  %-15s %-8s align=%d dup=%.1f %s %s' % (r[5], r[6], r[0], r[2], r[7], r[8]))
