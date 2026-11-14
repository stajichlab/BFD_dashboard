import os, csv, re, statistics
B = os.path.realpath(os.path.expanduser('~/projects/BFD_dashboard/bench/kickx_tree'))
R = f'{B}/coemansia/results/protein'
taxa = {r['LOCUSTAG']: r for r in csv.DictReader(open(f'{B}/taxa/coemansia_taxa.tsv'), delimiter='\t')}
COE = {k for k, r in taxa.items() if r['role'] == 'ingroup'}

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

cher, pure = [], []
def walk(n):
    s = tips(n)
    if n['kids']:
        if len(n['kids']) == 2 and all(not k['kids'] for k in n['kids']) and s <= COE:
            cher.append(n['kids'][0]['bl'] + n['kids'][1]['bl'])
        if s <= COE and len(s) >= 2:
            try:
                pure.append((len(s), float(n['label'])))
            except ValueError:
                pass
    for k in n['kids']:
        walk(k)
walk(root)
print('cherries (two Coemansia strains that are each other\'s closest relatives): %d' % len(cher))
print('  patristic distance between the pair: min %.5f, median %.4f, max %.3f (substitutions per site)' % (min(cher), statistics.median(cher), max(cher)))
for thr in (0.001, 0.005, 0.02, 0.05):
    print('  pairs closer than %.3f: %d' % (thr, sum(c < thr for c in cher)))
# maximal pure-Coemansia clades with support >= 0.95 (no outgroup inside), count and size
sup = sorted([(n, s) for n, s in pure if s >= 0.95], reverse=True)
maximal = []
for n, s in sup:
    maximal.append(n)
print('pure Coemansia clades with support >= 0.95 (nested counted): %d; largest sizes: %s' % (len(sup), [n for n, _ in sup[:10]]))
low = [(n, s) for n, s in pure if s < 0.7]
print('pure Coemansia internal nodes with support < 0.70: %d of %d' % (len(low), len(pure)))
