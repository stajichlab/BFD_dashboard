#!/usr/bin/env python3
"""T5 step 6: gate, classification of the GH65 candidates, and the tree-based tests R1, R2, R3 and R5(a) (see README.md, "Kill rules").

Usage: 06_classify.py <BENCH>
Reads  $BENCH/trees/gh65.treefile, gh65_members.tsv, $BENCH/refs/refs.tsv, $BENCH/cand/genomes.tsv, $BENCH/cand/candidates.tsv
Writes $BENCH/results/gh65_classification.tsv   one row per annotated protein (call, group, support, nearest reference, in_tclade, architecture)
       $BENCH/results/gh65_loo.tsv              leave-one-out of the references (the gate)
       $BENCH/results/gh65_by_genus.tsv         per genus representative: proteins, in T-clade, quality, taxonomy
       $BENCH/results/t5_tree_results.json      the numbers behind R1, R2, R3 and R5(a), and the verdicts
       $BENCH/results/t5_tree_results.txt       the same as text

Everything about the gene tree is ROOTING-FREE (the tree has no outgroup of its own): "the smallest clade holding the fungal trehalases" is the
smallest side of any edge (split) that contains all of them, and "nearest non-fungal neighbours" are found by absorbing, at the boundary node, the
smallest neighbouring subtree first and moving on across the largest, so the answer does not depend on where IQ-TREE printed its root.
Tip sets are Python integers used as bit masks.
"""
import csv
import json
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'kx2_ku_pathway', 'family_trees'))
import fam_lib as F  # noqa: E402

GATE_GROUP = 'TREHALASE_FUNGI'
DIKARYA = {'Ascomycota', 'Basidiomycota'}
MIN_BUSCO, MIN_N50 = 85.0, 50000
R1_KILL_FRACTION = 0.30
MIN_UFB, MIN_SH = 90.0, 80.0


# ------------------------------------------------------------------ tree structure with bit masks
class GTree:
    """A parsed Newick tree with, for every node, the bit mask of the tips below it."""

    def __init__(self, text):
        self.root = F.parse_newick(text)
        self.tips = [n.name for n in F.leaves(self.root)]
        self.bit = {t: 1 << i for i, t in enumerate(self.tips)}
        self.all = (1 << len(self.tips)) - 1
        self.down = {}
        self._fill(self.root)

    def _fill(self, n):
        # iterative post-order (the trees have thousands of tips and Python recursion depth is limited)
        order, stack = [], [n]
        while stack:
            x = stack.pop()
            order.append(x)
            stack.extend(x.kids)
        for x in reversed(order):
            self.down[id(x)] = self.bit[x.name] if not x.kids else sum_or([self.down[id(k)] for k in x.kids])

    def nodes(self):
        out, stack = [], [self.root]
        while stack:
            x = stack.pop()
            out.append(x)
            stack.extend(x.kids)
        return out

    def mask_of(self, names):
        m = 0
        for t in names:
            m |= self.bit[t]
        return m

    def tips_in(self, mask):
        return [t for t in self.tips if mask & self.bit[t]]


def sum_or(ms):
    r = 0
    for m in ms:
        r |= m
    return r


def popcount(x):
    return x.bit_count()


def smallest_side_containing(tree, tmask):
    """The smallest side (as a bit mask) of any edge that contains every bit of tmask, with the (node, 'down'|'up') that defines it.
    'down' = the tips below the node, 'up' = all other tips. Rooting-free."""
    best = None
    for v in tree.nodes():
        if v.parent is None:
            continue
        for side, which in ((tree.down[id(v)], 'down'), (tree.all ^ tree.down[id(v)], 'up')):
            if side & tmask == tmask and (best is None or popcount(side) < popcount(best[0])):
                best = (side, v, which)
    return best


def edge_support(node):
    """(SH-aLRT, UFBoot) of the edge above `node` as parsed by fam_lib, or (None, None)."""
    return node.support if node.support else (None, None)


def supported(sup, min_ufb=MIN_UFB, min_sh=MIN_SH):
    sh, ufb = sup
    return ufb is not None and ufb >= min_ufb and (sh is None or sh >= min_sh)


def maximal_pure_sides(tree, pure_mask):
    """Maximal sides (over all edges, both directions) whose tips are all inside pure_mask. Returns [(mask, node, which)]."""
    cands = []
    for v in tree.nodes():
        if v.parent is None:
            continue
        for side, which in ((tree.down[id(v)], 'down'), (tree.all ^ tree.down[id(v)], 'up')):
            if side and side & ~pure_mask == 0:
                cands.append((side, v, which))
    cands.sort(key=lambda c: -popcount(c[0]))
    kept = []
    for c in cands:
        if not any(c[0] & k[0] == c[0] for k in kept):
            kept.append(c)
    return kept


def nearest_nonfungal(tree, side_node, which, fungal_mask):
    """Expand outward from the split side defined by (side_node, which) and return the first neighbouring subtree that holds a non-fungal tip.
    At the boundary node the neighbouring subtrees are taken smallest first: an all-fungal one is absorbed, the first one with a non-fungal tip
    stops the search; if only the largest remains it is crossed and the search continues at its far end.
    Returns dict(mask=bit mask of the non-fungal tips in the stopping subtree, support=(sh, ufb) of that edge, depth=steps) or None."""
    if which == 'down':
        u, frm = side_node.parent, side_node
    else:
        u, frm = side_node, side_node.parent
    depth = 0
    while u is not None and frm is not None:
        neigh = []
        for k in u.kids:
            if k is not frm:
                neigh.append((tree.down[id(k)], k, edge_support(k)))
        if u.parent is not None and u.parent is not frm:
            neigh.append((tree.all ^ tree.down[id(u)], u.parent, edge_support(u)))
        if not neigh:
            return None
        neigh.sort(key=lambda t: popcount(t[0]))
        for mask, nb, sup in neigh[:-1]:                 # every neighbour but the largest: stop at the first with a non-fungal tip, else absorb it
            nonfungal = mask & ~fungal_mask
            if nonfungal:
                return {'mask': nonfungal, 'support': sup, 'depth': depth, 'side_mask': mask}
        mask, nb, sup = neigh[-1]                         # the largest remaining neighbour
        if not nb.kids and mask & ~fungal_mask:           # a single non-fungal tip
            return {'mask': mask & ~fungal_mask, 'support': sup, 'depth': depth, 'side_mask': mask}
        if not nb.kids and not (mask & ~fungal_mask):     # a single fungal tip: nothing beyond it
            return None
        u, frm = nb, u                                    # cross it and continue at its far end
        depth += 1
    return None


# ------------------------------------------------------------------ the analysis
def load_tsv(path):
    return list(csv.DictReader(open(path), delimiter='\t'))


def analyse(tree, members, refs, genomes, cands, log=print):
    """members: gh65_members rows; refs: {accession: row of refs.tsv}; genomes: {locustag: row}; cands: {protein_id: row}.
    Returns (results dict, per-protein rows, loo rows, per-genus rows)."""
    byl = {m['label']: m for m in members}
    ref_group = {m['label']: m['clade_or_group'] for m in members if m['kind'] == 'ref'}
    groups = Counter(ref_group.values())
    ann = [m['label'] for m in members if m['kind'] == 'annotated' and m['label'] in tree.bit]
    res = {'n_tips': len(tree.tips), 'n_annotated': len(ann), 'reference_groups': dict(groups)}

    # ---- gate: leave-one-out and the support of the smallest clade holding the fungal trehalase references
    ntest, ncorrect, details = F.loo_check(tree.root, {k: v for k, v in ref_group.items() if k in tree.bit}, min_ufb=MIN_UFB, min_sh=MIN_SH)
    loo_rows = [dict(zip(('reference', 'group', 'call', 'assigned', 'correct'), d)) for d in details]
    t_refs = [k for k, v in ref_group.items() if v == GATE_GROUP and k in tree.bit]
    t_loo = [d for d in details if d[1] == GATE_GROUP]
    side = smallest_side_containing(tree, tree.mask_of(t_refs))
    ref_clade_sup = edge_support(side[1]) if side else (None, None)
    ref_clade_ok = side is not None and supported(ref_clade_sup)
    gate = bool(t_loo) and all(d[4] for d in t_loo) and len(t_loo) == len(t_refs) and ref_clade_ok
    res['gate'] = {'loo_tested': ntest, 'loo_correct': ncorrect, 'trehalase_refs': len(t_refs), 'trehalase_loo_correct': sum(1 for d in t_loo if d[4]),
                   'trehalase_ref_clade_support': ref_clade_sup, 'trehalase_ref_clade_tips': popcount(side[0]) if side else None, 'passed': gate}
    log(f"gate: LOO {ncorrect}/{ntest}; fungal trehalase references {sum(1 for d in t_loo if d[4])}/{len(t_refs)} re-assigned; their smallest clade support {ref_clade_sup}; passed={gate}")

    # ---- classification of every annotated protein
    refs_only = {k: v for k, v in ref_group.items() if k in tree.bit}
    cls = F.classify(tree.root, refs_only, ann, min_ufb=MIN_UFB, min_sh=MIN_SH)
    rows = []
    for lab in ann:
        c, m = cls[lab], byl[lab]
        cand = cands.get(m['id'], {})
        g = genomes.get(m['locustag'], {})
        near = ref_group.get(c['nearest_ref'], '')
        rows.append({'label': lab, 'protein_id': m['id'], 'locustag': m['locustag'], 'genus': g.get('genus', ''), 'phylum': g.get('PHYLUM', ''),
                     'class': g.get('class', ''), 'length': m['length'], 'domains': cand.get('domains', ''), 'call': c['call'], 'group': c['clade'],
                     'ufboot': c['ufboot'], 'sh': c['sh'], 'nearest_ref': byl[c['nearest_ref']]['id'] if c['nearest_ref'] in byl else '',
                     'nearest_group': near, 'terminal_bl': round(c['terminal_bl'], 4)})
    n = len(rows)
    to_t = [r for r in rows if r['call'] == 'assigned' and r['group'] == GATE_GROUP]
    to_other = [r for r in rows if r['call'] == 'assigned' and r['group'] != GATE_GROUP]
    unres = [r for r in rows if r['call'] != 'assigned']
    res['R1'] = {'proteins': n, 'assigned_trehalase': len(to_t), 'assigned_other': len(to_other), 'unresolved': len(unres),
                 'other_groups': dict(Counter(r['group'] for r in to_other)), 'fraction_other': len(to_other) / n if n else None,
                 'unresolved_nearest_group': dict(Counter(r['nearest_group'] for r in unres))}
    res['R1']['verdict'] = ('tree cannot decide (gate failed)' if not gate else
                            'KILLED: >= 30% of fungal GH65 proteins are assigned to non-trehalase groups' if res['R1']['fraction_other'] >= R1_KILL_FRACTION else
                            'survives: fewer than 30% are assigned to non-trehalase groups')

    # ---- R2: one origin of the fungal trehalase lineage
    fungal_tips = [r['label'] for r in rows] + t_refs
    fungal_mask = tree.mask_of(fungal_tips)
    t_tips = [r['label'] for r in to_t] + t_refs
    tmask = tree.mask_of(t_tips)
    side = smallest_side_containing(tree, tmask)
    smask, snode, swhich = side
    inside = tree.tips_in(smask)
    nonfungal_inside = [t for t in inside if not (fungal_mask & tree.bit[t])]
    inside_fungal_ann = [t for t in inside if t in byl and byl[t]['kind'] == 'annotated']
    pure = maximal_pure_sides(tree, fungal_mask)
    big_pure = [(m, v, w) for m, v, w in pure if popcount(m & tmask) >= 5 and supported(edge_support(v))]
    sup = edge_support(snode)
    res['R2'] = {'tclade_tips': len(inside), 'tclade_assigned_tips': len(t_tips), 'tclade_fungal_annotated': len(inside_fungal_ann),
                 'tclade_nonfungal_tips': len(nonfungal_inside), 'tclade_nonfungal_labels': [byl[t]['id'] for t in nonfungal_inside][:20],
                 'tclade_support': sup, 'n_supported_fungal_pure_clades_with_5_assigned': len(big_pure)}
    if not gate:
        res['R2']['verdict'] = 'tree cannot decide (gate failed)'
    elif nonfungal_inside or len(big_pure) >= 2:
        res['R2']['verdict'] = 'REJECTED: not one fungal-only clade (non-fungal tips inside, or two supported fungal clades)'
    elif supported(sup):
        res['R2']['verdict'] = 'supported: one fungal-only clade, root supported (UFBoot >= 90)'
    else:
        res['R2']['verdict'] = 'unresolved: fungal-only clade but the root is not supported (UFBoot < 90)'

    # ---- R3: Dikarya stem
    in_t = set(inside_fungal_ann)
    for r in rows:
        r['in_tclade'] = r['label'] in in_t
    gen_in = defaultdict(list)
    for r in rows:
        if r['in_tclade']:
            gen_in[r['locustag']].append(r)
    non_dik = []
    for loc, rs in gen_in.items():
        g = genomes.get(loc, {})
        if g.get('PHYLUM') in DIKARYA:
            continue
        q = quality(g)
        full = any(('PF03636' in r['domains'] and 'PF03633' in r['domains']) for r in rs)
        non_dik.append({'locustag': loc, 'genus': g.get('genus'), 'phylum': g.get('PHYLUM'), 'class': g.get('class'), 'quality': q, 'full_architecture': full,
                        'counts': q and full})
    res['R3'] = {'genera_in_tclade': len(gen_in), 'dikarya_genera': sum(1 for loc in gen_in if genomes.get(loc, {}).get('PHYLUM') in DIKARYA),
                 'non_dikarya_genera': non_dik, 'non_dikarya_counting': sum(1 for x in non_dik if x['counts'])}
    if not gate:
        res['R3']['verdict'] = 'tree cannot decide (gate failed)'
    elif not res['R2']['verdict'].startswith('supported'):
        res['R3']['verdict'] = 'not evaluated: R2 did not support a single fungal-only clade'
    elif res['R3']['non_dikarya_counting'] >= 1:
        res['R3']['verdict'] = 'origin OLDER than the Dikarya stem: Dikarya-stem part of T5 rejected (>= 1 non-Dikarya genus with quality genome and full architecture)'
    else:
        res['R3']['verdict'] = 'supported: only Dikarya genera (quality, full architecture) inside the clade'

    # ---- R5(a): nearest non-fungal neighbours of the fungal trehalase clade
    nn = nearest_nonfungal(tree, snode, swhich, fungal_mask)
    if nn is None:
        res['R5a'] = {'verdict': 'no non-fungal neighbour found'}
    else:
        labs = tree.tips_in(nn['mask'])
        types = Counter(tip_type(byl[t], refs) for t in labs)
        groups_nn = Counter(byl[t]['clade_or_group'] for t in labs)
        prok = all(k in ('Bacteria', 'Archaea') for k in types)
        res['R5a'] = {'neighbour_tips': len(labs), 'neighbour_types': dict(types), 'neighbour_groups': dict(groups_nn), 'support': nn['support'],
                      'steps': nn['depth'], 'neighbour_ids': [byl[t]['id'] for t in labs][:15]}
        if not gate:
            res['R5a']['verdict'] = 'tree cannot decide (gate failed)'
        elif prok and supported(nn['support']):
            res['R5a']['verdict'] = 'HGT-consistent: nearest non-fungal neighbours are only bacterial/archaeal, edge supported'
        elif prok:
            res['R5a']['verdict'] = 'prokaryote-only neighbours but the edge is not supported (unresolved)'
        else:
            res['R5a']['verdict'] = 'vertical-consistent: nearest non-fungal neighbours include eukaryotes'
    genus_rows = []
    for loc, g in genomes.items():
        mine = [r for r in rows if r['locustag'] == loc]
        genus_rows.append({'locustag': loc, 'genus': g.get('genus'), 'phylum': g.get('PHYLUM'), 'class': g.get('class'), 'order': g.get('ord'),
                           'busco_complete': g.get('busco_complete'), 'N50_bp': g.get('N50_bp'), 'quality': quality(g), 'n_gh65': len(mine),
                           'n_in_tclade': sum(1 for r in mine if r['in_tclade']), 'n_assigned_trehalase': sum(1 for r in mine if r['call'] == 'assigned' and r['group'] == GATE_GROUP),
                           'n_assigned_other': sum(1 for r in mine if r['call'] == 'assigned' and r['group'] != GATE_GROUP)})
    return res, rows, loo_rows, genus_rows


def quality(g):
    try:
        return float(g.get('busco_complete') or 0) >= MIN_BUSCO and float(g.get('N50_bp') or 0) >= MIN_N50
    except ValueError:
        return False


def tip_type(member, refs):
    """Bacteria / Archaea / Eukaryota for a reference or context tip, from refs.tsv (UniProt lineage)."""
    r = refs.get(member['id'])
    return r['domain'] if r else 'Fungi(BFD)'


def write_table(path, rows):
    if not rows:
        open(path, 'w').write('')
        return
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main(bench):
    out = f'{bench}/results'
    os.makedirs(out, exist_ok=True)
    tree = GTree(open(f'{bench}/trees/gh65.treefile').read())
    members = load_tsv(f'{bench}/trees/gh65_members.tsv')
    refs = {r['accession']: r for r in load_tsv(f'{bench}/refs/refs.tsv')}
    genomes = {r['locustag']: r for r in load_tsv(f'{bench}/cand/genomes.tsv')}
    cands = {r['protein_id']: r for r in load_tsv(f'{bench}/cand/candidates.tsv')}
    lines = []
    res, rows, loo_rows, genus_rows = analyse(tree, members, refs, genomes, cands, log=lambda s: (print(s), lines.append(s)))
    write_table(f'{out}/gh65_classification.tsv', rows)
    write_table(f'{out}/gh65_loo.tsv', loo_rows)
    write_table(f'{out}/gh65_by_genus.tsv', genus_rows)
    json.dump(res, open(f'{out}/t5_tree_results.json', 'w'), indent=1, default=str)
    for k in ('R1', 'R2', 'R3', 'R5a'):
        lines.append(f'{k}: {json.dumps(res[k], default=str)}')
        print(f'{k}: {res[k].get("verdict")}')
    open(f'{out}/t5_tree_results.txt', 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main(sys.argv[1])
