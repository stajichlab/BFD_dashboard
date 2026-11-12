import os, csv, glob, statistics
from collections import defaultdict
B = os.path.realpath(os.path.expanduser('~/projects/BFD_dashboard/bench/kx2_ku_tblastn'))
gen = None
for cand in glob.glob(f'{B}/**/genomes_round3.tsv', recursive=True):
    gen = cand
    break
assert gen, 'genomes_round3.tsv not found'
rows = list(csv.DictReader(open(gen), delimiter='\t'))
EV = 1e-25          # strict: only close relatives of Rad51/Dmc1 (>~45% identity over most of the protein), not distant paralogs
MINLEN = 150        # aligned residues
GAP = 20000         # merge HSPs on the same scaffold within 20 kb (introns)
KEYS = ('CTRL_Rad51', 'RecA_DMC1')
res = {}
for g in rows:
    loc = g['locustag']
    p = f'{B}/hits_round3/{loc}.tsv'
    if not os.path.exists(p):
        continue
    per_scaf = defaultdict(list)
    keys_at = defaultdict(set)
    for r in csv.reader(open(p), delimiter='\t'):
        q, scaf, pid, aln, ev, bits = r[0], r[1], float(r[2]), int(r[3]), float(r[4]), float(r[5])
        key = q.split('|')[0]
        if key not in KEYS or ev > EV or aln < MINLEN:
            continue
        s, e = sorted((int(r[9]), int(r[10])))
        per_scaf[scaf].append((s, e, key))
    loci = []
    for scaf, iv in per_scaf.items():
        iv.sort()
        cs, ce, ks = None, None, set()
        for s, e, k in iv:
            if ce is None or s > ce + GAP:
                if ce is not None:
                    loci.append((scaf, cs, ce, ks))
                cs, ce, ks = s, e, {k}
            else:
                ce = max(ce, e); ks.add(k)
        if ce is not None:
            loci.append((scaf, cs, ce, ks))
    res[loc] = (g, loci)
def grp(g):
    return g['grp']
by = defaultdict(list)
for loc, (g, loci) in res.items():
    by[grp(g)].append(len(loci))
print('Distinct Rad51/Dmc1-like loci per genome (tblastn round 3; e<=1e-25, >=150 aa, HSPs merged within 20 kb)')
for k in sorted(by):
    v = by[k]
    print(f'{k:32s} n={len(v):3d} median={statistics.median(v):.1f} dist={ {x: v.count(x) for x in sorted(set(v))} }')
print('genomes searched', len(res), 'of', len(rows))
# per-locus best bit score by query key, for genomes with >= 2 loci (and a few with 1), to see whether the second locus is Dmc1-like
bits = defaultdict(lambda: defaultdict(float))
for g in rows:
    loc = g['locustag']
    p = f'{B}/hits_round3/{loc}.tsv'
    if not os.path.exists(p):
        continue
    for r in csv.reader(open(p), delimiter='\t'):
        key = r[0].split('|')[0]
        ev, aln = float(r[4]), int(r[3])
        if key in KEYS and ev <= EV and aln >= MINLEN:
            s, e = sorted((int(r[9]), int(r[10])))
            for i, (scaf, cs, ce, ks) in enumerate(res[loc][1]):
                if r[1] == scaf and s >= cs - 1 and e <= ce + 1:
                    bits[(loc, i)][key] = max(bits[(loc, i)][key], float(r[5]))
print('\nloci in genomes with 2+ loci: best bit score per query key (R51=Rad51 queries, DMC=Dmc1 queries); locus is Dmc1-like if DMC >> R51')
for loc, (g, loci) in res.items():
    if len(loci) >= 2:
        parts = []
        for i, (scaf, cs, ce, ks) in enumerate(loci):
            b = bits[(loc, i)]
            parts.append(f"R51={b.get('CTRL_Rad51', 0):.0f}/DMC={b.get('RecA_DMC1', 0):.0f}")
        print(f"{g['grp'][:18]:18s} {g['genus'][:15]:15s} {loc} " + ' | '.join(parts))
