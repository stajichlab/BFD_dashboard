import os, csv, re, statistics
B = os.path.realpath(os.path.expanduser('~/projects/BFD_dashboard/bench/kickx_tree'))
R = f'{B}/coemansia/results/protein'
taxa = {r['LOCUSTAG']: r for r in csv.DictReader(open(f'{B}/taxa/coemansia_taxa.tsv'), delimiter='\t')}

fa = [f for f in os.listdir(f'{R}/buildtree') if f.endswith('.fa')][0]
seq, cur = {}, None
for l in open(f'{R}/buildtree/{fa}'):
    l = l.rstrip()
    if l.startswith('>'):
        cur = l[1:].split()[0]; seq[cur] = []
    else:
        seq[cur].append(l)
gaps = {k: 100 * sum(c in '-?Xx' for c in ''.join(v)) / len(''.join(v)) for k, v in seq.items()}
L = len(''.join(next(iter(seq.values()))))
print(f'concat: {len(seq)} sequences of {len(taxa)} taxa, {L} sites; median gap {statistics.median(gaps.values()):.1f}%')
print('missing from concat:', [k for k in taxa if k not in seq])
co = {k: g for k, g in gaps.items() if taxa[k]['role'] == 'ingroup'}
print('Coemansia strains: gap pct median %.1f, 90th pct %.1f, max %.1f' % (statistics.median(co.values()), sorted(co.values())[int(.9 * len(co))], max(co.values())))
print('strains with >60%% gaps: %d' % sum(g > 60 for g in co.values()))
print('outgroup gap %:', {taxa[k]['GENUS']: round(gaps[k]) for k in taxa if taxa[k]['role'] == 'outgroup' and k in gaps})

def parse(nwk):
    nwk = nwk.strip(); pos = 0
    def node():
        nonlocal pos
        if nwk[pos] == '(':
            pos += 1; kids = []
            while True:
                kids.append(node())
                if nwk[pos] == ',':
                    pos += 1; continue
                if nwk[pos] == ')':
                    pos += 1; break
            m = re.match(r'([^:,;()]*)', nwk[pos:]); lab = m.group(1); pos += len(lab)
        else:
            m = re.match(r'([^:,;()]+)', nwk[pos:]); lab = m.group(1); pos += len(lab); kids = []
        bl = 0.0
        if pos < len(nwk) and nwk[pos] == ':':
            m = re.match(r':([^,;()]*)', nwk[pos:]); bl = float(m.group(1)); pos += len(m.group(0))
        return {'kids': kids, 'label': lab, 'bl': bl}
    return node()

tf = [f for f in os.listdir(f'{R}/buildtree/fungi_odb12/fasttree') if f.endswith('fasttree.support.treefile')][0]
root = parse(open(f'{R}/buildtree/fungi_odb12/fasttree/{tf}').read())

def tips(n):
    return frozenset([n['label']]) if not n['kids'] else frozenset().union(*[tips(k) for k in n['kids']])
allt = tips(root)
print(f'\nFastTree: {len(allt)} tips')
COE = frozenset(k for k, r in taxa.items() if r['role'] == 'ingroup')
OUT = frozenset(k for k, r in taxa.items() if r['role'] == 'outgroup')
found = None
def walk(n):
    global found
    s = tips(n)
    if s == COE or (allt - s) == COE:
        found = n
    for k in n['kids']:
        walk(k)
walk(root)
print('Coemansia (107) monophyletic:', 'yes, support=' + (found['label'] or '-') if found else 'NO')
if not found:
    best = None
    def w2(n):
        global best
        s = tips(n)
        for side in (s, allt - s):
            if COE <= side and (best is None or len(side) < len(best)):
                best = side
        for k in n['kids']:
            w2(k)
    w2(root)
    print('smallest clade holding all Coemansia also holds:', sorted(taxa[k]['GENUS'] for k in best - COE))
    # where does Kickxella sit? find its parent and siblings
    kk = [k for k, r in taxa.items() if r['GENUS'] == 'Kickxella'][0]
    def findparent(n, target, path=()):
        for k in n['kids']:
            if not k['kids'] and k['label'] == target:
                return n, path + (n,)
            r = findparent(k, target, path + (n,))
            if r:
                return r
        return None
    par, path = findparent(root, kk)
    print('Kickxella parent node: support=%s, size=%d strains+taxa' % (par['label'] or '-', len(tips(par))))
    for k in par['kids']:
        s = tips(k)
        print('   child:', ('Kickxella' if s == frozenset([kk]) else f'{len(s)} taxa, {sum(1 for t in s if t in COE)} Coemansia strains, node support={k["label"] or "-"}'))
    for up in reversed(path[:-1][-3:]):
        s = tips(up)
        print('   ancestor node: %d taxa (%d Coemansia, %d outgroup) support=%s' % (len(s), sum(1 for t in s if t in COE), sum(1 for t in s if t in OUT), up['label'] or '-'))
    # locate the node whose tip set equals best (rooted orientation) and run redundancy stats on it
    def find_exact(n, S):
        if tips(n) == S:
            return n
        for k in n['kids']:
            r = find_exact(k, S)
            if r:
                return r
    found = find_exact(root, best)
    if found:
        print('analysing redundancy on the clade holding all strains: %d taxa' % len(best))

# near-clonal strains: leaf pairs (cherries) and short sister clades
def leaves_with_depth(n, d=0.0, out=None):
    out = [] if out is None else out
    if not n['kids']:
        out.append((n['label'], d + n['bl']))
    else:
        for k in n['kids']:
            leaves_with_depth(k, d + n['bl'], out)
    return out
if found:
    subtree = found
    cher = []
    sup = []
    def collect(n):
        if n['kids'] and all(not k['kids'] for k in n['kids']) and len(n['kids']) == 2:
            cher.append(n['kids'][0]['bl'] + n['kids'][1]['bl'])
        if n['kids'] and n['label']:
            try:
                sup.append((len(tips(n)), float(n['label'])))
            except ValueError:
                pass
        for k in n['kids']:
            collect(k)
    collect(subtree)
    print('cherries among Coemansia strains: %d; patristic distance between pair: median %.4f, min %.5f' % (len(cher), statistics.median(cher), min(cher)))
    for thr in (0.001, 0.005, 0.02):
        print('  cherries with distance < %.3f subst/site: %d' % (thr, sum(c < thr for c in cher)))
    big = sorted([(n, s) for n, s in sup if n >= 5 and s >= 0.95], reverse=True)
    print('clades inside Coemansia with >=5 strains and support >=0.95:', len(big), '| sizes (largest 12):', [n for n, _ in big[:12]])
    low = sum(1 for n, s in sup if s < 0.7)
    print('internal nodes with support < 0.70:', low, 'of', len(sup))
    # tree depth
    depths = [d for _, d in leaves_with_depth(subtree)]
    print('root-to-tip within Coemansia clade: median %.3f, range %.3f-%.3f' % (statistics.median(depths), min(depths), max(depths)))
# outgroup arrangement
def nearest_outgroups():
    pass
og = [(taxa[k]['GENUS'], taxa[k]['ord']) for k in OUT]
print('outgroups:', sorted(g for g, _ in og))
