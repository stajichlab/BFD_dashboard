#!/usr/bin/env python3
"""10_tree_clades.py
Purpose : turn the first-pass backbone Newick into an edge table for Dollo parsimony (done in SQL, 11_dollo.sql).
Input   : path to backbone_fasttree_support.nwk (tips = LOCUSTAGs; internal labels = FastTree local support)
Output  : writes SQL to stdout that creates bfd_scratch.kx2_ku_tree_edges(node INT, parent INT, support DOUBLE,
          tip VARCHAR (LOCUSTAG, NULL for internal nodes)). Node 0 is the printed root.
Usage   : python3 -I 10_tree_clades.py /path/to/backbone_fasttree_support.nwk > 10_tree_clades.sql
Note    : the tree is a first-pass FastTree topology, used as printed (rooted at the printed trifurcation, which is
          resolved by the caller: see README for the re-rooting statement). No data processing beyond tree text parsing.
"""
import re
import sys


def parse(s):
    s = s.strip().rstrip(";")
    pos = 0
    nodes = []  # (node, parent, support, tip)

    def newnode(parent, support, tip):
        nodes.append([len(nodes), parent, support, tip])
        return len(nodes) - 1

    def read_label_len():
        nonlocal pos
        m = re.compile(r"([^,:()]*)(?::([0-9.eE+-]+))?").match(s, pos)
        pos = m.end()
        return m.group(1), m.group(2)

    def read_node(parent):
        nonlocal pos
        if s[pos] == "(":
            me = newnode(parent, None, None)
            pos += 1
            while True:
                read_node(me)
                if s[pos] == ",":
                    pos += 1
                    continue
                assert s[pos] == ")"
                pos += 1
                break
            lab, _ = read_label_len()
            if lab:
                nodes[me][2] = float(lab)
            return me
        lab, _ = read_label_len()
        return newnode(parent, None, lab)

    read_node(-1)
    return nodes


if __name__ == "__main__":
    with open(sys.argv[1]) as fh:
        nodes = parse(fh.read())
    rows = []
    for n, p, sup, tip in nodes:
        rows.append("(%d,%d,%s,%s)" % (n, p, "NULL" if sup is None else repr(sup), "NULL" if tip is None else "'%s'" % tip))
    print("CREATE OR REPLACE TABLE bfd_scratch.kx2_ku_tree_edges AS SELECT * FROM (VALUES")
    print(",\n".join(rows))
    print(") t(node, parent, support, tip);")
