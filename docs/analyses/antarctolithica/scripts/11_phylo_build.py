"""Quick phylogenomic placement: BUSCO fungi_odb10 markers -> per-marker MAFFT -> trimAl -> concatenation -> FastTree (LG+G).
Single-copy rule: best-scoring hit per marker if score >= BUSCO cutoff; marker dropped for a taxon if a 2nd hit scores >= 85% of the best (paralog ambiguity)."""
import os, glob, subprocess, sys, pandas as pd, numpy as np
from concurrent.futures import ThreadPoolExecutor
os.chdir('phylo'); os.makedirs('aln', exist_ok=True); os.makedirs('mk', exist_ok=True)
BIN = '../phyenv/bin/'
T = pd.read_csv('taxa.tsv', sep='\t'); T['base'] = T.pep.str.replace(r'.*/', '', regex=True).str.replace('.proteins.fa.gz', '')
cut = pd.read_csv('db/scores_cutoff', sep='\t', header=None, index_col=0)[1]
best = {}; pres = {}
for _, r in T.iterrows():
    f = f'tbl/{r.base}.tbl'
    if not os.path.exists(f): print('missing', f); continue
    rows = [l.split()[:6] for l in open(f) if not l.startswith('#')]
    h = pd.DataFrame(rows, columns=['prot', 'a', 'marker', 'b', 'ev', 'score']); h['score'] = h.score.astype(float)
    h = h[h.score >= h.marker.map(cut)].sort_values('score', ascending=False)
    b = {}
    for m, g in h.groupby('marker'):
        if len(g) > 1 and g.score.iloc[1] >= 0.85 * g.score.iloc[0]: b[m] = None      # ambiguous / duplicated
        else: b[m] = g.prot.iloc[0]
    best[r.label] = b; pres[r.label] = {m: (0 if m not in b else (2 if b[m] is None else 1)) for m in cut.index}
P = pd.DataFrame(pres).T; P.to_csv('marker_presence.csv')        # 0 missing, 1 single, 2 ambiguous
occ = (P == 1).mean(); ao = [l for l in P.index if 'Antarctolithica' in l]
hq = [l for l in ao if any(t in l for t in ['FF1FC3B3', 'F46C92EF', 'F7B4B958', 'F7E47E53'])]
sel = occ[(occ >= 0.85) & ((P.loc[hq] == 1).sum() >= 3)].sort_values(ascending=False).index[:300]
print('taxa', len(P), 'markers passing', len(sel)); print('single-copy markers per AO MAG:', (P.loc[ao] == 1).sum(1).to_dict())
seqs = {}
for _, r in T.iterrows():
    if r.label not in best: continue
    want = {p: m for m, p in best[r.label].items() if p and m in set(sel)}; name = None; d = {}
    for line in open(f'fa/{r.base}.fa'):
        if line[0] == '>': name = line[1:].split()[0]; name = name if name in want else None; 
        elif name: d[name] = d.get(name, '') + line.strip().rstrip('*')
    for p, s in d.items(): seqs.setdefault(want[p], {})[r.label] = s.replace('*', 'X')
def run(m):
    if os.path.exists(f'aln/{m}.trim.fa') and os.path.getsize(f'aln/{m}.trim.fa') > 0: return
    with open(f'mk/{m}.fa', 'w') as o:
        for l, s in seqs[m].items(): o.write(f'>{l}\n{s}\n')
    subprocess.run(f'{BIN}mafft --auto --thread 1 --quiet mk/{m}.fa > aln/{m}.aln.fa && {BIN}trimal -in aln/{m}.aln.fa -out aln/{m}.trim.fa -gappyout', shell=True, check=True)
with ThreadPoolExecutor(12) as ex: list(ex.map(run, sel))
cat = {l: [] for l in best}; parts = []
for m in sel:
    a = {}; name = None
    for line in open(f'aln/{m}.trim.fa'):
        if line[0] == '>': name = line[1:].split()[0]; a[name] = ''
        else: a[name] += line.strip()
    L = len(next(iter(a.values()))); parts.append((m, L))
    for l in cat: cat[l].append(a.get(l, '-' * L))
with open('concat.fa', 'w') as o:
    for l, p in cat.items():
        s = ''.join(p); o.write(f'>{l}\n{s}\n')
tot = sum(L for _, L in parts); print('concatenated columns:', tot, 'markers:', len(parts))
gapf = {l: ''.join(p).count('-') / tot for l, p in cat.items()}; print('gap fraction AO:', {l[-8:]: round(gapf[l], 3) for l in ao})
pd.Series(gapf).to_csv('concat_gapfrac.csv')
env = dict(os.environ, OMP_NUM_THREADS='12')
subprocess.run(f'{BIN}FastTree -lg -gamma -quiet concat.fa > concat.fasttree.nwk 2> fasttree.log', shell=True, check=True, env=env)
print('tree done')
