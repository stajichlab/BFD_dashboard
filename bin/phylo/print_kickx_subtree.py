import os, csv, re
B = os.path.realpath(os.path.expanduser('~/projects/BFD_dashboard/bench/kickx_tree'))
R = f'{B}/backbone/results/protein'
taxa = {r['LOCUSTAG']: r for r in csv.DictReader(open(f'{B}/taxa/backbone_taxa.tsv'), delimiter='\t')}
ING = frozenset(k for k, r in taxa.items() if r['PHYLUM'] == 'Kickxellomycota')

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
        if pos < len(nwk) and nwk[pos] == ':':
            m = re.match(r':[^,;()]*', nwk[pos:]); pos += len(m.group(0))
        return {'kids': kids, 'label': lab}
    return node()

def tips(n):
    return frozenset([n['label']]) if not n['kids'] else frozenset().union(*[tips(k) for k in n['kids']])

def find(n, S):
    if tips(n) == S:
        return n
    for k in n['kids']:
        r = find(k, S)
        if r:
            return r
    return None

def show(n, depth=0):
    if not n['kids']:
        t = taxa[n['label']]
        flag = ' [' + ','.join(f for f, c in (('lowBUSCO', t.get('low_busco') == 'true'), ('dup', t.get('high_dup') == 'true')) if c) + ']'
        flag = '' if flag == ' []' else flag
        return '  ' * depth + f"{t['GENUS']} ({t['ord']}){flag}"
    sup = n['label'] or '-'
    out = ['  ' * depth + f'+ node support={sup} ({len(tips(n))} taxa)']
    for k in sorted(n['kids'], key=lambda x: -len(tips(x))):
        out.append(show(k, depth + 1))
    return '\n'.join(out)

for name, p in [('FastTree', R + '/buildtree/fungi_odb12/fasttree/' + [f for f in os.listdir(R + '/buildtree/fungi_odb12/fasttree') if f.endswith('fasttree.support.treefile')][0]),
                ('ASTRAL', R + '/tree/fungi_odb12/final_tree.nw')]:
    root = parse(open(p).read())
    n = find(root, ING)
    print(f'\n===== {name}: Kickxellomycota subtree' + ('' if n else ' (clade not found in this rooting)'))
    if n:
        print(show(n))
