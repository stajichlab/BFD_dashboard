import os, csv, re, sys, statistics
tree = sys.argv[1] if len(sys.argv) > 1 else 'backbone'
B = os.path.realpath(os.path.expanduser('~/projects/BFD_dashboard/bench/kickx_tree'))
R = f'{B}/{tree}/results/protein'
taxa = {r['LOCUSTAG']: r for r in csv.DictReader(open(f'{B}/taxa/{tree}_taxa.tsv'), delimiter='\t')}

# ---- concat coverage
fa = [f for f in os.listdir(f'{R}/buildtree') if f.endswith('.fa')][0]
seq, cur = {}, None
for l in open(f'{R}/buildtree/{fa}'):
    l = l.rstrip()
    if l.startswith('>'):
        cur = l[1:].split()[0]; seq[cur] = []
    else:
        seq[cur].append(l)
gaps = {}
for k, v in seq.items():
    s = ''.join(v)
    gaps[k] = 100 * sum(c in '-?Xx' for c in s) / len(s)
L = len(''.join(next(iter(seq.values()))))
print(f'concat: {len(seq)} sequences of {len(taxa)} taxa, {L} sites; median gap {statistics.median(gaps.values()):.1f}%')
print('missing from concat:', [k for k in taxa if k not in seq])
worst = sorted(gaps.items(), key=lambda x: -x[1])[:6]
print('highest gap %:', [(taxa.get(k, {}).get('GENUS', k), round(g, 1)) for k, g in worst])

# ---- newick parsing (unrooted splits with edge support)
def parse(nwk):
    nwk = nwk.strip()
    pos = 0
    def node():
        nonlocal pos
        if nwk[pos] == '(':
            pos += 1
            kids = []
            while True:
                kids.append(node())
                if nwk[pos] == ',':
                    pos += 1; continue
                if nwk[pos] == ')':
                    pos += 1; break
            m = re.match(r'([^:,;()]*)', nwk[pos:]); lab = m.group(1); pos += len(lab)
        else:
            m = re.match(r'([^:,;()]+)', nwk[pos:]); lab = m.group(1); pos += len(lab); kids = []
        if pos < len(nwk) and nwk[pos] == ':':
            m = re.match(r':[^,;()]*', nwk[pos:]); pos += len(m.group(0))
        return {'kids': kids, 'label': lab}
    return node()

def clades(root):
    out = []
    def walk(n):
        if not n['kids']:
            return frozenset([n['label']])
        s = frozenset().union(*[walk(k) for k in n['kids']])
        out.append((s, n['label']))
        return s
    walk(root)
    return out

def analyse(path, name):
    root = parse(open(path).read())
    cl = clades(root)
    allt = frozenset().union(*[c for c, _ in cl][:1]) if cl else frozenset()
    tips = max((c for c, _ in cl), key=len)
    print(f'\n== {name}: {len(tips)} tips')
    gen = lambda ks: sorted(taxa[k]['GENUS'] for k in ks if k in taxa)
    def mono(S, label):
        S = frozenset(S)
        if not S:
            return
        for c, sup in cl:
            if c == S or (tips - c) == S:
                print(f'  {label:42s} monophyletic (n={len(S):2d}) support={sup or "-"}')
                return
        # find smallest split containing S
        best = None
        for c, sup in cl:
            for side in (c, tips - c):
                if S <= side and (best is None or len(side) < len(best[0])):
                    best = (side, sup)
        extra = [taxa[k]['GENUS'] for k in best[0] - S if k in taxa] if best else []
        print(f'  {label:42s} NOT monophyletic (n={len(S):2d}); smallest clade holding it also has: {", ".join(sorted(extra))[:150]}')
    sel = lambda f: {k for k, r in taxa.items() if f(r)}
    mono(sel(lambda r: r['PHYLUM'] == 'Kickxellomycota'), 'Kickxellomycota')
    for o in ['Kickxellales', 'Dimargaritales', 'Harpellales', 'Zoopagales', 'Entomophthorales', 'Mucorales', 'Mortierellales', 'Glomerales']:
        mono(sel(lambda r, o=o: r['ord'] == o), o)
    mono(sel(lambda r: r['PHYLUM'] == 'Mucoromycota'), 'Mucoromycota')
    mono(sel(lambda r: r['PHYLUM'] in ('Mucoromycota', 'Mortierellomycota', 'Glomeromycota')), 'Mucoromycota+Mortierellomycota+Glomeromycota')
    mono(sel(lambda r: r['PHYLUM'] in ('Kickxellomycota', 'Zoopagomycota')), 'Kickxellomycota+Zoopagomycota')
    mono(sel(lambda r: r['PHYLUM'] in ('Kickxellomycota', 'Zoopagomycota', 'Entomophthoromycota', 'Basidiobolomycota')),
         'Kickx+Zoopago+Entomophth+Basidiobolus')
    mono(sel(lambda r: r['role'] == 'root'), 'root taxa (chytrid etc.)')
    mono(sel(lambda r: r['role'] != 'root'), 'all non-root taxa')

for name, p in [('FastTree (support)', f'{R}/buildtree/fungi_odb12/fasttree/' + [f for f in os.listdir(f'{R}/buildtree/fungi_odb12/fasttree') if f.endswith('fasttree.support.treefile')][0]),
                ('ASTRAL consensus/final', f'{R}/tree/fungi_odb12/final_tree.nw')]:
    try:
        analyse(p, name)
    except Exception as e:
        print('FAILED', name, repr(e))
