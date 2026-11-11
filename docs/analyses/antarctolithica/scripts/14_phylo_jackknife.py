"""Robustness of the placement: marker jackknife (4 x 100 random markers) and a reduced-taxon run without the longest-branch yeasts; FastTree LG+G and one WAG+G."""
import os, glob, subprocess, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
os.chdir('phylo'); os.makedirs('jk', exist_ok=True); BIN = '../phyenv/bin/'
T = pd.read_csv('taxa.tsv', sep='\t').set_index('label'); mk = sorted(m[4:-8] for m in glob.glob('aln/*.trim.fa')); rng = np.random.default_rng(7)
def load(m):
    a = {}; n = None
    for line in open(f'aln/{m}.trim.fa'):
        if line[0] == '>': n = line[1:].split()[0]; a[n] = ''
        else: a[n] += line.strip()
    return a
A = {m: load(m) for m in mk}
def write(name, markers, taxa):
    with open(f'jk/{name}.fa', 'w') as o:
        for t in taxa: o.write(f'>{t}\n' + ''.join(A[m].get(t, '-' * len(next(iter(A[m].values())))) for m in markers) + '\n')
alltx = list(T.index); drop = set(T.index[T.CLASS.isin(['Pneumocystomycetes', 'Schizosaccharomycetes', 'Saccharomycetes', 'Pichiomycetes', 'Mucoromycetes'])])
jobs = []
for i in range(4): write(f'jk{i}', list(rng.choice(mk, 100, replace=False)), alltx); jobs.append((f'jk{i}', '-lg'))
write('reduced', list(rng.choice(mk, 150, replace=False)), [t for t in alltx if t not in drop]); jobs.append(('reduced', '-lg'))
write('wag', list(rng.choice(mk, 100, replace=False)), alltx); jobs.append(('wag', '-wag'))
def run(j):
    env = dict(os.environ, OMP_NUM_THREADS='2'); subprocess.run(f'{BIN}FastTree {j[1]} -gamma -quiet jk/{j[0]}.fa > jk/{j[0]}.nwk 2> jk/{j[0]}.log', shell=True, env=env)
with ThreadPoolExecutor(6) as ex: list(ex.map(run, jobs))
print('JK_DONE')
