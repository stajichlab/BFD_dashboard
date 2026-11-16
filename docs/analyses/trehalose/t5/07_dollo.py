#!/usr/bin/env python3
"""T5 step 7: minimum number of losses of the fungal acid trehalase on the genus tree (Dollo parsimony), rule R4 in README.md.

Usage: 07_dollo.py <BENCH> <label=tree.nw> [<label=tree.nw> ...]
Reads  $BENCH/results/gh65_by_genus.tsv (step 6: per genus representative, in_tclade count, quality, class) and each genus tree (tip names are
       LOCUSTAGs of the genus representatives; tips without a row in the table are ignored and counted).
Writes $BENCH/results/dollo_summary.tsv (per tree x scoring: tips, present, absent, unknown, gain node size, losses, classes touched, permutation
       mean/sd and the share of permutations with as few losses or fewer), $BENCH/results/dollo_losses.tsv (every loss clade), dollo_verdict.txt.

Rooting: on the early-diverging zoosporic fungi (Chytridiomycota, Blastocladiomycota, Neocallimastigomycota, Monoblepharidomycota/-omycota). If they
are not one branch of the unrooted tree, the branch whose side overlaps them best (Jaccard) is used and the overlap is reported. The rerooting
function is the tested one in workflows/phylo/kickx_tree_tools.py.
Dollo: presence = the genus representative has a protein inside the fungal trehalase clade. SCORED scoring: a genus is absent only if its
representative passes the quality filter (BUSCO >= 85, N50 >= 50 kb); other genera are unknown. UNFILTERED scoring: every genus is scored.
The gain is placed at the MRCA of the present genera; a loss is a maximal clade inside that subtree in which every scored genus is absent (clades of
unknown genera only are not losses), so losses are never nested. The permutation null shuffles the present/absent labels over the scored genera
(same number present, 1,000 permutations, seed 20261113) and repeats the whole procedure.
"""
import csv
import os
import random
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'workflows', 'phylo'))
sys.setrecursionlimit(100000)

OUTGROUP_PHYLA = {'Chytridiomycota', 'Blastocladiomycota', 'Neocallimastigomycota', 'Monoblepharidomycota', 'Monoblepharomycota'}
N_PERM, SEED = 1000, 20261113


# ------------------------------------------------------------------ pure Dollo on an array tree
def to_arrays(root):
    """Array form of a rooted Bio.Phylo tree: names (None for internal nodes), parent, children, postorder. Iterative."""
    names, parent, children = [], [], []
    stack = [(root, -1)]
    while stack:
        c, p = stack.pop()
        i = len(names)
        names.append(c.name if not c.clades else None)
        parent.append(p)
        children.append([])
        if p >= 0:
            children[p].append(i)
        for k in c.clades:
            stack.append((k, i))
    order = []
    st = [0]
    while st:
        x = st.pop()
        order.append(x)
        st.extend(children[x])
    return names, parent, children, order[::-1]


def count_states(children, order, tip_state):
    """present and absent counts below every node; tip_state[i] in {1, 0, -1} for tips (ignored for internal nodes)."""
    n = len(children)
    p, a = [0] * n, [0] * n
    for x in order:
        if not children[x]:
            p[x] = 1 if tip_state[x] == 1 else 0
            a[x] = 1 if tip_state[x] == 0 else 0
        else:
            p[x] = sum(p[c] for c in children[x])
            a[x] = sum(a[c] for c in children[x])
    return p, a


def dollo_losses(children, order, tip_state, want_clades=False):
    """Minimum losses under Dollo: gain at the MRCA of the present tips; losses = maximal clades below it with no present and >= 1 absent tip.
    Returns (n_losses, gain_node, [loss nodes])."""
    p, a = count_states(children, order, tip_state)
    root = order[-1]
    total = p[root]
    if total == 0:
        return 0, None, []
    g = root
    while True:                                   # descend while one child still holds all present tips
        nxt = [c for c in children[g] if p[c] == total]
        if not nxt:
            break
        g = nxt[0]
    losses, stack = [], list(children[g])
    while stack:
        x = stack.pop()
        if p[x] == 0:
            if a[x] > 0:
                losses.append(x)
        else:
            stack.extend(children[x])
    return len(losses), g, losses if want_clades else []


def tips_below(children, x):
    out, st = [], [x]
    while st:
        y = st.pop()
        if not children[y]:
            out.append(y)
        st.extend(children[y])
    return out


def permutation_null(children, order, tip_state, n_perm=N_PERM, seed=SEED):
    """Losses when the present/absent labels are shuffled over the scored tips (unknown tips stay unknown)."""
    scored = [i for i, s in enumerate(tip_state) if s in (0, 1)]
    labels = [tip_state[i] for i in scored]
    rng = random.Random(seed)
    out = []
    state = list(tip_state)
    for _ in range(n_perm):
        rng.shuffle(labels)
        for i, v in zip(scored, labels):
            state[i] = v
        out.append(dollo_losses(children, order, state)[0])
    return out


# ------------------------------------------------------------------ rooting and I/O
def outgroup_side(tree, outgroup_tips):
    """Tip set of the branch whose side best matches the outgroup (Jaccard), and that overlap."""
    from kickx_tree_tools import tips
    allt = set(tips(tree.root))
    g = set(outgroup_tips) & allt
    best = (0.0, g)
    for c in tree.find_clades():
        if c is tree.root:
            continue
        s = set(tips(c))
        for side in (s, allt - s):
            j = len(side & g) / len(side | g) if (side | g) else 0
            if j > best[0]:
                best = (j, side)
    return best[1], best[0]


def load_genus_rows(bench):
    return {r['locustag']: r for r in csv.DictReader(open(f'{bench}/results/gh65_by_genus.tsv'), delimiter='\t')}


def run_tree(label, path, genus, bench):
    from Bio import Phylo
    import kickx_tree_tools as T
    tree = Phylo.read(path, 'newick')
    in_tree = {t for t in T.tips(tree.root)}
    missing = sorted(set(genus) - in_tree)
    ignored = sorted(in_tree - set(genus))
    out_tips = [t for t in in_tree if genus.get(t, {}).get('phylum') in OUTGROUP_PHYLA]
    side, jac = outgroup_side(tree, out_tips)
    rooted = T.root_on(tree, sorted(side))
    names, parent, children, order = to_arrays(rooted.root)
    rows, loss_rows = [], []
    idx = {nm: i for i, nm in enumerate(names) if nm}
    for scoring in ('scored', 'unfiltered'):
        state = []
        for nm in names:
            if nm is None or nm not in genus:
                state.append(-1)
                continue
            g = genus[nm]
            present = int(g['n_in_tclade']) > 0
            if present:
                state.append(1)
            elif scoring == 'unfiltered' or g['quality'] == 'True':
                state.append(0)
            else:
                state.append(-1)
        n_loss, gnode, lnodes = dollo_losses(children, order, state, want_clades=True)
        null = permutation_null(children, order, state)
        classes = set()
        for x in lnodes:
            tl = [names[i] for i in tips_below(children, x)]
            sc = [t for t in tl if t in genus and state[idx[t]] == 0]
            cl = sorted({genus[t]['class'] for t in sc})
            classes.update(cl)
            loss_rows.append({'tree': label, 'scoring': scoring, 'n_scored_absent': len(sc), 'n_tips': len(tl), 'classes': ';'.join(cl),
                              'example_genera': ';'.join(sorted({genus[t]['genus'] for t in sc})[:6])})
        gsize = len(tips_below(children, gnode)) if gnode is not None else 0
        rows.append({'tree': label, 'scoring': scoring, 'tips_in_table': len(genus) - len(missing), 'tips_ignored': len(ignored),
                     'outgroup_jaccard': round(jac, 3), 'present': sum(1 for s in state if s == 1), 'absent': sum(1 for s in state if s == 0),
                     'unknown': sum(1 for i, s in enumerate(state) if s == -1 and names[i] in genus), 'gain_clade_tips': gsize, 'losses': n_loss,
                     'classes_touched': len(classes), 'perm_mean': round(statistics.mean(null), 1), 'perm_sd': round(statistics.pstdev(null), 1),
                     'share_perm_le_observed': round(sum(1 for x in null if x <= n_loss) / len(null), 4)})
    return rows, loss_rows


def verdict(rows):
    primary = next((r for r in rows if r['scoring'] == 'scored'), None)
    if primary is None:
        return 'no result'
    n, k = primary['losses'], primary['classes_touched']
    if n <= 2:
        return f"R4 REJECTED: only {n} loss(es) (scored genera, tree {primary['tree']})"
    if k >= 2:
        return f"R4 supported: {n} losses touching {k} classes (scored genera, tree {primary['tree']})"
    return f"R4 inconclusive: {n} losses but all in {k} class (scored genera, tree {primary['tree']})"


def main(bench, specs):
    genus = load_genus_rows(bench)
    allrows, alllosses = [], []
    for spec in specs:
        label, path = spec.split('=', 1)
        r, l = run_tree(label, path, genus, bench)
        allrows += r
        alllosses += l
    out = f'{bench}/results'
    for name, rows in (('dollo_summary.tsv', allrows), ('dollo_losses.tsv', alllosses)):
        with open(f'{out}/{name}', 'w', newline='') as f:
            if rows:
                w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n')
                w.writeheader()
                w.writerows(rows)
    v = verdict(allrows)
    open(f'{out}/dollo_verdict.txt', 'w').write(v + '\n')
    for r in allrows:
        print(r)
    print(v)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
