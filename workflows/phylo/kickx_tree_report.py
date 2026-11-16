#!/usr/bin/env python3
"""Rooting, clade tests, comparison with the first-pass tree and figures for the improved Kickxellomycota trees.

Usage: kickx_tree_report.py <mode: coemansia|backbone> <iq dir> <taxa.tsv> <firstpass.treefile> <outdir>
  <iq dir>   $BENCH/iq/<name> with pmsf_ml.treefile (optional pmsf_ml_sup.suptree), guide.treefile, astral.tre, astral_c50.tre (missing files are skipped)
  <taxa.tsv> the kickx_tree taxa table of the run (columns role, LOCUSTAG, ord, GENUS, SPECIES, STRAIN, ...)
Outgroup rules, fixed before the trees were viewed:
  coemansia  the 121-taxon set has no taxon outside Kickxellomycota; root on Dimargaritales (Dimargaris, Dispira, Tieghemiomyces), which the
             63-taxon backbone tree (outgroups: chytrids, Blastocladiomycota, Monoblepharidomycota, Neocallimastigomycota) places as the earliest
             Kickxellomycota branch. Sensitivity: Dimargaritales + Ramicandelaber, Harpellales, midpoint.
  backbone   root on the 'root'-role taxa (chytrid-grade fungi). Sensitivity: without Blastocladiomycota, Zoopagomycota + Entomophthoromycota
             side excluded is not possible, so midpoint is the second option.
Writes <outdir>/<mode>_report.txt, rooted Newick files with readable labels, <mode>_tree.pdf/.png (PMSF tree, support on branches) and
<mode>_clades.tsv (clade tests per tree).
"""
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kickx_tree_tools as T
from Bio import Phylo

mode, iqdir, taxa_tsv, firstpass, outdir = sys.argv[1:6]
os.makedirs(outdir, exist_ok=True)
rows = list(csv.DictReader(open(taxa_tsv), delimiter='\t'))
info = {r['LOCUSTAG']: r for r in rows}
genus = {k: (v['GENUS'] or v['SPECIES'].split()[0]) for k, v in info.items()}
order = {k: v['ord'] for k, v in info.items()}
phylum = {k: v['PHYLUM'] for k, v in info.items()}
lines = []
P = lines.append


def label(k):
    r = info[k]
    sp = r.get('SPECIES') or genus[k]
    strain = r.get('STRAIN') or ''
    return f"{sp}{' ' + strain if strain and strain not in sp else ''} [{k}]"


def support(c):
    """UFBoot/bootstrap/local PP of an internal clade as (primary, secondary). IQ-TREE writes '-alrt' + '-B' labels as
    'SH-aLRT/UFBoot' (e.g. '97.3/100'), so for those the primary value is the second number (UFBoot) and the secondary the SH-aLRT."""
    if c.confidence is not None:
        return c.confidence, None
    if c.name and '/' in c.name:
        a, b = c.name.split('/')[:2]
        try:
            return float(b), float(a)
        except ValueError:
            pass
    return None, None


def load_tree(name):
    p = os.path.join(iqdir, name)
    return T.load(p) if os.path.exists(p) and os.path.getsize(p) > 0 else None


# PMSF ML tree = pmsf_ml.treefile (the 'pmsfml' stage, a run WITHOUT -b). NEVER use pmsf.treefile: with -b, IQ-TREE 3.0.1 overwrites it
# with the tree of every bootstrap replicate. Support, when available: iqtree3 -sup pmsf_ml.treefile -t pmsf.boottrees
# --prefix pmsf_ml_sup (written by the caller) gives pmsf_ml_sup.suptree; the label says how many replicates it rests on.
_nb = 0
_pm = 'pmsf_ml.treefile'
if os.path.exists(os.path.join(iqdir, 'pmsf_ml_sup.suptree')):
    _pm = 'pmsf_ml_sup.suptree'
    if os.path.exists(os.path.join(iqdir, 'pmsf.boottrees')):
        _nb = sum(1 for _ in open(os.path.join(iqdir, 'pmsf.boottrees')) if _.strip())
PMSF_LABEL = 'PMSF LG+C60+F+G ML' + (f' ({_nb} bootstrap replicates)' if _nb else ' (no bootstrap support)')
trees = {PMSF_LABEL: load_tree(_pm), 'LG+F+G (UFBoot)': load_tree('guide.treefile'),
         'ASTRAL (gene trees)': load_tree('astral.tre'), 'ASTRAL (UFBoot<50 collapsed)': load_tree('astral_c50.tre')}
fp = T.load(firstpass) if os.path.exists(firstpass) else None
if fp is not None:
    trees['first-pass FastTree'] = fp
trees = {k: v for k, v in trees.items() if v is not None}
sets = {k: set(T.tips(t.root)) for k, t in trees.items()}
primary_name = next(iter(trees))
taxa_in = sets[primary_name]
P(f'mode {mode}; tree set {primary_name} has {len(taxa_in)} tips; trees compared: ' + '; '.join(f'{k} ({len(s)} tips)' for k, s in sets.items()))

if mode == 'coemansia':
    OUT = [k for k in taxa_in if genus[k] in ('Dimargaris', 'Dispira', 'Tieghemiomyces')]
    ALT = {'Dimargaritales (primary)': OUT,
           'Dimargaritales + Ramicandelaber': OUT + [k for k in taxa_in if genus[k] == 'Ramicandelaber'],
           'Harpellales': [k for k in taxa_in if order[k] == 'Harpellales'],
           'Mycoemilia + Spiromyces': [k for k in taxa_in if genus[k] in ('Mycoemilia', 'Spiromyces')]}
else:
    OUT = [k for k in taxa_in if info[k]['role'] == 'root']
    ALT = {'root-role taxa (primary)': OUT,
           'root-role taxa without Blastocladiomycota': [k for k in OUT if phylum[k] != 'Blastocladiomycota']}
P(f'primary outgroup: {len(OUT)} taxa: ' + ', '.join(sorted({genus[k] for k in OUT})))

KU_NEG = lambda k: order[k] in ('Kickxellales', 'Harpellales')


def root_tree(t, outgroup):
    """Root on the outgroup; `mono` = the outgroup is separated by one branch of the unrooted tree (independent of the file's root)."""
    mono = T.has_split(t, outgroup)
    if mono:
        r = T.root_on(t, outgroup)
        return r, True, []
    r = T.root_on(t, outgroup)                       # rooted on the MRCA branch, but not a true outgroup rooting
    return r, False, T.is_monophyletic(r, outgroup)[1]


def sister_of(t, tip):
    """Genus composition (with counts) of the sister group of `tip` in the rooted tree t, and its support."""
    path = t.get_path(tip)
    parent = t.root if len(path) < 2 else path[-2]
    tc = [c for c in parent.clades if c is not path[-1]]
    comp = defaultdict(int)
    for c in tc:
        for n in T.tips(c):
            comp[genus[n]] += 1
    return ', '.join(f'{g}{" x" + str(n) if n > 1 else ""}' for g, n in sorted(comp.items(), key=lambda x: -x[1])), support(parent)[0]


def clade_test(t, name):
    out = {}
    tips_in = set(T.tips(t.root))
    co = [k for k in tips_in if genus[k] == 'Coemansia']
    if co:
        mono, extra = T.is_monophyletic(t, co)
        out['Coemansia monophyletic'] = f'{mono}' + ('' if mono else ' (MRCA also holds: ' + ', '.join(sorted({genus[k] for k in extra})) + ')')
    for g in ('Kickxella', 'Linderina', 'Dipsacomyces', 'Martensiomyces', 'Ramicandelaber', 'Basidiobolus'):
        k = next((x for x in tips_in if genus[x] == g), None)
        if k:
            s, sup = sister_of(t, k)
            out[f'sister of {g}'] = f'{s}' + (f' (support {sup:g})' if sup is not None else '')
    kx = [k for k in tips_in if phylum[k] == 'Kickxellomycota']
    if mode == 'backbone' and len(kx) > 1:
        monok, extrak = T.is_monophyletic(t, kx)
        mk = t.common_ancestor(kx)
        supk = support(mk)[0]
        out['Kickxellomycota monophyletic'] = str(monok) + (f' (support {supk:g})' if supk is not None and monok else '') + (
            '' if monok else ' (MRCA also holds: ' + ', '.join(sorted({genus[k] for k in extrak})) + ')')
        if monok:
            out['sister of Kickxellomycota'] = sister_comp(t, mk)
    neg = [k for k in tips_in if KU_NEG(k)]
    if neg:
        mono, extra = T.is_monophyletic(t, neg)
        m = t.common_ancestor(neg)
        sup = support(m)[0]
        out['Kickxellales+Harpellales monophyletic'] = f'{mono}' + (f' (support {sup:g})' if sup is not None and mono else '') + \
            ('' if mono else ' (MRCA also holds: ' + ', '.join(sorted({genus[k] for k in extra})) + ')')
        ke = [k for k in tips_in if order[k] == 'Kickxellales']
        ha = [k for k in tips_in if order[k] == 'Harpellales']
        out['Kickxellales monophyletic'] = str(T.is_monophyletic(t, ke)[0]) if ke else ''
        out['Harpellales monophyletic'] = str(T.is_monophyletic(t, ha)[0]) if ha else ''
        if mono:
            out['stem sister of Kickxellales+Harpellales'] = sister_comp(t, m)
        rami = [k for k in tips_in if genus[k] == 'Ramicandelaber']
        if mono and rami:
            khr = neg + rami
            mono2, extra2 = T.is_monophyletic(t, khr)
            out['Kickxellales+Harpellales+Ramicandelaber monophyletic'] = str(mono2) + (
                '' if mono2 else ' (MRCA also holds: ' + ', '.join(sorted({genus[k] for k in extra2})) + ')')
            if mono2:
                out['sister of Kickxellales+Harpellales+Ramicandelaber'] = sister_comp(t, t.common_ancestor(khr))
    return out


def sister_comp(t, m):
    """Genus composition of the sister group of clade m (rooted tree)."""
    path = t.get_path(m)
    parent = t.root if len(path) < 2 else path[-2]
    comp = defaultdict(int)
    for c in [c for c in parent.clades if c is not m]:
        for n in T.tips(c):
            comp[genus[n]] += 1
    return ', '.join(f'{g}{" x" + str(n) if n > 1 else ""}' for g, n in sorted(comp.items(), key=lambda x: -x[1]))


rooted = {}
P('\n== Rooting of the primary tree ==')
for name, t in trees.items():
    r, mono, extra = root_tree(t, [k for k in OUT if k in sets[name]])
    rooted[name] = r
    if name == primary_name or not mono:
        P(f'{name}: outgroup monophyletic in the unrooted tree = {mono}' + ('' if mono else f'; MRCA also holds {sorted({genus[k] for k in extra})}'))

clades = {n: clade_test(t, n) for n, t in rooted.items()}
keys = []
for c in clades.values():
    for k in c:
        if k not in keys:
            keys.append(k)
P('\n== Clade tests per tree (rooted on the primary outgroup) ==')
with open(f'{outdir}/{mode}_clades.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['question'] + list(rooted))
    for k in keys:
        w.writerow([k] + [clades[n].get(k, '') for n in rooted])
        P(f'-- {k}')
        for n in rooted:
            P(f'   {n:36s} {clades[n].get(k, "(not in this tree)")}')

P('\n== Robinson-Foulds distance between trees (unrooted, shared tips; normalized = RF / (splits1 + splits2)) ==')
names = list(trees)
for i, a in enumerate(names):
    for b in names[i + 1:]:
        d, m, n = T.rf(trees[a], trees[b])
        P(f'{a:36s} vs {b:36s}  RF {d:4d} / {m:4d} = {d / m:.3f}   ({n} shared tips)')

P('\n== Support in the primary tree ==')
t0 = trees[primary_name]
vals = [support(c)[0] for c in t0.get_nonterminals() if c is not t0.root and support(c)[0] is not None]
if vals:
    P(f'{len(vals)} internal branches; support >= 95: {sum(v >= 95 for v in vals)}; >= 70: {sum(v >= 70 for v in vals)}; < 50: {sum(v < 50 for v in vals)}')

P('\n== Rooting sensitivity (primary tree, alternative outgroups) ==')
tp = trees[primary_name]
for alt, og in ALT.items():
    og = [k for k in og if k in taxa_in]
    if not og:
        continue
    r, mono, extra = root_tree(tp, og)
    ct = clade_test(r, alt)
    P(f'{alt}: outgroup monophyletic in the unrooted tree = {mono}' + ('' if mono else f' (MRCA also holds {sorted({genus[k] for k in extra})})'))
    for k in ('Coemansia monophyletic', 'Kickxellales+Harpellales monophyletic', 'Kickxellales monophyletic', 'stem sister of Kickxellales+Harpellales',
              'sister of Kickxellales+Harpellales+Ramicandelaber'):
        if k in ct:
            P(f'   {k}: {ct[k]}')
import copy
mt0 = copy.deepcopy(tp)
mt0.root_at_midpoint()                       # Bio only finds the midpoint split; its rerooting mislabels supports, so reroot ourselves
mt = T.root_on(tp, T.tips(mt0.root.clades[0]))
ct = clade_test(mt, 'midpoint')
P('midpoint root:')
for k in ('Coemansia monophyletic', 'Kickxellales+Harpellales monophyletic', 'Kickxellales monophyletic', 'stem sister of Kickxellales+Harpellales',
              'sister of Kickxellales+Harpellales+Ramicandelaber'):
    if k in ct:
        P(f'   {k}: {ct[k]}')
P('   root splits: ' + ' | '.join(', '.join(sorted({genus[n] for n in T.tips(c)})[:6]) + f' ({len(T.tips(c))} tips)' for c in mt.root.clades))

# ---- strain redundancy (coemansia): shortest tip pairs in the primary tree
if mode == 'coemansia':
    co = [k for k in taxa_in if genus[k] == 'Coemansia']
    best = []
    tc = copy.deepcopy(tp)
    for i, a in enumerate(co):
        for b in co[i + 1:]:
            best.append((tc.distance(a, b), a, b))
    best.sort()
    P(f'\n== Coemansia representatives in the tree: {len(co)}; closest pairs (substitutions per site) ==')
    for d, a, b in best[:5]:
        P(f'   {d:.4f}  {label(a)}  {label(b)}')

# ---- write labelled rooted trees and figure
for name, r in rooted.items():
    slug = ''.join(ch if ch.isalnum() else '_' for ch in name).strip('_')[:40]
    rr = copy.deepcopy(r)
    for tip in rr.get_terminals():
        tip.name = label(tip.name).replace(' ', '_').replace('[', '').replace(']', '').replace(':', '').replace(',', '')
    Phylo.write(rr, f'{outdir}/{mode}_{slug}.rooted.nwk', 'newick')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

colors = {'Kickxellales': '#c0392b', 'Harpellales': '#e67e22', 'Dimargaritales': '#2471a3', 'Ramicandelaberales': '#1e8449',
          'Zoopagales': '#7d3c98', 'Entomophthorales': '#117a65', 'Basidiobolales': '#117a65', 'Mucorales': '#7f8c8d'}
plot_name = primary_name
if not vals and 'LG+F+G (UFBoot)' in rooted:        # the PMSF tree has no support until its bootstrap run finishes
    plot_name = 'LG+F+G (UFBoot)'
pt = copy.deepcopy(rooted[plot_name])
pt.ladderize()
n = len(pt.get_terminals())
fig, ax = plt.subplots(figsize=(11, max(8, 0.19 * n)))
lc = {}
for tip in pt.get_terminals():
    lc[label(tip.name)] = colors.get(order[tip.name], '#566573')
sup_map = {id(c): support(c)[0] for c in pt.get_nonterminals()}
for c in pt.get_nonterminals():
    c.name = None                      # IQ-TREE 'SH-aLRT/UFBoot' labels would otherwise be drawn as node names
    c.confidence = None
bl = lambda c: (f'{sup_map[id(c)]:.0f}' if id(c) in sup_map and sup_map[id(c)] is not None and sup_map[id(c)] < 100 else None)
for tip in pt.get_terminals():
    tip.name = label(tip.name)
Phylo.draw(pt, axes=ax, do_show=False, label_colors=lc, branch_labels=bl, show_confidence=False)
xmax = max(pt.distance(t) for t in pt.get_terminals())
ax.set_xlim(-0.02 * xmax, xmax * 1.55)
ax.set_xlabel('substitutions per site')
from matplotlib.lines import Line2D
used = sorted({order[t] for t in taxa_in if order[t] in colors})
ax.legend([Line2D([0], [0], color=colors[o], lw=3) for o in used] + [Line2D([0], [0], color='#566573', lw=3)], used + ['other'],
          loc='upper right', fontsize=7, frameon=False)
ax.set_title(f'{plot_name}; rooted on {"Dimargaritales" if mode == "coemansia" else "chytrid-grade outgroup"}; numbers = support below 100; '
             f'scale = substitutions per site', fontsize=9)
for lbl in ax.texts:
    if lbl.get_text().isdigit():
        lbl.set_fontsize(5)
    else:
        lbl.set_fontsize(6)
fig.tight_layout()
fig.savefig(f'{outdir}/{mode}_tree.pdf')
fig.savefig(f'{outdir}/{mode}_tree.png', dpi=130)
open(f'{outdir}/{mode}_report.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
