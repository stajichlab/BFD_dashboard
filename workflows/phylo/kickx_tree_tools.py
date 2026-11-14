#!/usr/bin/env python3
"""Small tree utilities for the Kickxellomycota trees (Biopython only; run with /usr/bin/python3.12 on the HPCC).

Subcommands
  collapse <in.tre> <out.tre> <min_support>   collapse internal branches with support < min_support in every tree of a multi-tree file
                                              (IQ-TREE gene trees: support = UFBoot; branch lengths kept; tips unchanged)
Also a library: load(), root_on(), clades(), rf(), is_monophyletic(), mrca_support()  used by kickx_tree_report.py
"""
import io
import sys

from Bio import Phylo


def load(path, which=0):
    trees = list(Phylo.parse(path, 'newick'))
    return trees[which]


def tips(clade):
    return [t.name for t in clade.get_terminals()]


def root_on(tree, outgroup):
    """Root on the branch leading to the MRCA of `outgroup` (a list of tip names). Returns a rooted copy. The outgroup must be
    monophyletic in the unrooted tree for this to be a true outgroup rooting; the caller checks that with is_monophyletic."""
    import copy
    from Bio.Phylo.BaseTree import Clade
    t = copy.deepcopy(tree)
    allt = set(tips(t.root))
    g = {n for n in outgroup if n in allt}
    # The file's root is arbitrary (IQ-TREE and FastTree print a trifurcating root, and Bio.Phylo.root_with_outgroup mishandles
    # an outgroup that is already a child of it), so reroot by hand on the branch that separates the outgroup from the rest.
    target = None
    for c in t.find_clades():
        if c is t.root:
            continue
        s = set(tips(c))
        if s == g or s == allt - g:
            target = c
            break
    if target is None or len(g) == 0:                  # outgroup is not one branch of the tree: fall back to Bio's MRCA rooting
        t.root_with_outgroup(sorted(g)) if len(g) > 1 else t.root_with_outgroup(sorted(g)[0])
        return t
    chain = [t.root] + t.get_path(target)               # old root ... target
    k = len(chain) - 1
    info = [(n.branch_length, n.confidence, n.name) for n in chain]   # attributes of the edge above each node
    C, P = chain[k], chain[k - 1]
    P.clades = [x for x in P.clades if x is not C]
    for i in range(k - 1, 0, -1):                       # reverse the path: chain[i-1] becomes a child of chain[i]
        chain[i - 1].clades = [x for x in chain[i - 1].clades if x is not chain[i]]
        chain[i].clades.append(chain[i - 1])
        chain[i - 1].branch_length, chain[i - 1].confidence, chain[i - 1].name = info[i]
    total, conf, name = info[k]
    P.branch_length, P.confidence, P.name = None, None, None
    rest = P

    def suppress(node):                                 # remove degree-2 nodes (merge branch lengths, keep the lower support)
        for i, ch in enumerate(node.clades):
            while len(ch.clades) == 1:
                only = ch.clades[0]
                if ch.branch_length is not None or only.branch_length is not None:
                    only.branch_length = (ch.branch_length or 0.0) + (only.branch_length or 0.0)
                node.clades[i] = only
                ch = only
            suppress(ch)
    while len(rest.clades) == 1:
        rest = rest.clades[0]
    suppress(rest)
    half = None if total is None else total / 2.0
    C.branch_length = half
    rest.branch_length, rest.confidence, rest.name = half, conf, name
    t.root = Clade(clades=[C, rest])
    return t


def clade_sets(tree):
    """Set of frozensets of tip names, one per internal clade of the rooted tree (root excluded)."""
    allt = frozenset(tips(tree.root))
    out = set()
    for c in tree.get_nonterminals():
        s = frozenset(tips(c))
        if 1 < len(s) < len(allt):
            out.add(s)
    return out


def bipartitions(tree):
    """Unrooted splits as frozensets of the smaller side (ties: side containing the sorted-first tip)."""
    allt = frozenset(tips(tree.root))
    first = min(allt)
    out = set()
    for c in tree.get_nonterminals():
        s = frozenset(tips(c))
        if len(s) < 2 or len(s) > len(allt) - 2:
            continue
        if first in s:
            s = allt - s
        out.add(s)
    return out


def rf(t1, t2, taxa=None):
    """Robinson-Foulds distance (unrooted) restricted to `taxa` (default: shared tips); returns (rf, max_rf, n_taxa)."""
    shared = set(tips(t1.root)) & set(tips(t2.root))
    if taxa is not None:
        shared &= set(taxa)
    import copy
    a, b = copy.deepcopy(t1), copy.deepcopy(t2)
    for t in (a, b):
        for n in [x for x in tips(t.root) if x not in shared]:
            t.prune(n)
    ba, bb = bipartitions(a), bipartitions(b)
    return len(ba ^ bb), len(ba) + len(bb), len(shared)


def is_monophyletic(tree, group):
    """True if `group` forms a clade of the rooted tree; also returns the tips of the MRCA that are not in the group."""
    group = [g for g in group if g in set(tips(tree.root))]
    if len(group) < 2:
        return True, []
    m = tree.common_ancestor(group)
    extra = sorted(set(tips(m)) - set(group))
    return not extra, extra


def has_split(tree, group):
    """True if the unrooted tree has a branch separating `group` from all other tips (independent of where the file's root sits)."""
    allt = frozenset(tips(tree.root))
    g = frozenset(x for x in group if x in allt)
    if len(g) < 2 or len(g) > len(allt) - 2:
        return True
    first = min(allt)
    side = allt - g if first in g else g
    return side in bipartitions(tree)


def collapse_trees(in_path, out_path, min_support):
    n_total = n_collapsed = 0
    with open(out_path, 'w') as out:
        for t in Phylo.parse(in_path, 'newick'):
            for c in list(t.get_nonterminals()):
                if c is t.root:
                    continue
                n_total += 1
                if c.confidence is not None and c.confidence < min_support:
                    t.collapse(c)
                    n_collapsed += 1
            s = io.StringIO()
            Phylo.write(t, s, 'newick')
            out.write(s.getvalue().strip() + '\n')
    print(f'collapsed {n_collapsed} of {n_total} internal branches with support < {min_support}')


if __name__ == '__main__':
    if len(sys.argv) >= 5 and sys.argv[1] == 'collapse':
        collapse_trees(sys.argv[2], sys.argv[3], float(sys.argv[4]))
    else:
        print(__doc__)
        sys.exit(2)
