#!/usr/bin/env python3
"""Fixture tests for kickx_derep.py, kickx_subset_aln.py and kickx_tree_tools.py (run on the HPCC: /usr/bin/python3.12 -m unittest test_kickx_phylo -v)."""
import os
import subprocess
import sys
import tempfile
import unittest
from io import StringIO

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PY = sys.executable
import kickx_tree_tools as T  # noqa: E402
from Bio import Phylo  # noqa: E402

HDR = 'role\tLOCUSTAG\tASMID\tPHYLUM\tcls\tord\tFAMILY\tGENUS\tSPECIES\tSTRAIN\tcomplete_pct\tduplicated_pct\tsingle_pct\tn_proteins\tN50_bp\tn50_ge_50kb\thigh_dup\n'


def row(role, loc, asm, comp, n50):
    return f'{role}\t{loc}\t{asm}\tK\tc\tKickxellales\tF\tCoemansia\tCoemansia sp.\t{asm}\t{comp}\t0\t0\t0\t{n50}\ttrue\tfalse\n'


class Derep(unittest.TestCase):
    def test_clusters_and_representative(self):
        with tempfile.TemporaryDirectory() as d:
            taxa = HDR + row('ingroup', 'L1', 'A1', 80, 100) + row('ingroup', 'L2', 'A2', 90, 50) + row('ingroup', 'L3', 'A3', 85, 10) + \
                row('ingroup', 'L4', 'A4', 70, 10) + row('outgroup', 'L5', 'A5', 90, 10)
            open(f'{d}/t.tsv', 'w').write(taxa)
            # A1-A2 at 99.5 (redundant), A2-A3 at 98 (not redundant at 99, redundant at 97), A4 distinct; both directions listed
            ani = ['/x/A1.fa.gz\t/x/A2.fa.gz\t99.5\t900\t1000', '/x/A2.fa.gz\t/x/A1.fa.gz\t99.4\t900\t1000',
                   '/x/A2.fa.gz\t/x/A3.fa.gz\t98.0\t900\t1000', '/x/A1.fa.gz\t/x/A1.fa.gz\t100\t1000\t1000']
            open(f'{d}/ani.out', 'w').write('\n'.join(ani) + '\n')
            subprocess.run([PY, f'{HERE}/kickx_derep.py', f'{d}/t.tsv', f'{d}/ani.out', d, '99.0'], check=True, capture_output=True)
            kept = [l.split('\t')[1] for l in open(f'{d}/derep_taxa.tsv').read().splitlines()[1:]]
            self.assertEqual(sorted(kept), ['L2', 'L3', 'L4', 'L5'])           # A1 dropped (A2 has higher BUSCO); outgroup kept
            s = open(f'{d}/derep_summary.txt').read()
            self.assertIn('threshold 99% : 3 clusters', s)
            self.assertIn('threshold 97% : 2 clusters', s)                    # A1-A2-A3 chain at 97


class Subset(unittest.TestCase):
    def test_subset_removes_taxa_and_gap_columns(self):
        with tempfile.TemporaryDirectory() as d:
            run = f'{d}/run'
            os.makedirs(f'{run}/filter/fungi_odb12-filter')
            os.makedirs(f'{run}/buildtree')
            open(f'{run}/buildtree/x.partition', 'w').write('PROT,m1.aa.mfa=1-4\nPROT,m2.aa.mfa=5-8\n')
            open(f'{run}/filter/fungi_odb12-filter/m1.aa.mfa', 'w').write('>A p1\nAC-D\n>B p2\nAC-D\n>C p3\nACGD\n')
            open(f'{run}/filter/fungi_odb12-filter/m2.aa.mfa', 'w').write('>A p1\nWWWW\n>C p3\nWWYW\n')
            open(f'{d}/taxa.txt', 'w').write('A\nB\n')
            subprocess.run([PY, f'{HERE}/kickx_subset_aln.py', run, f'{d}/taxa.txt', f'{d}/out', '1'], check=True, capture_output=True)
            conc = dict(l.split('\n')[:2] for l in open(f'{d}/out/concat.faa').read().split('>')[1:])
            self.assertEqual(conc['A'], 'ACDWWWW')      # column 3 of m1 is all gaps in {A,B}: removed
            self.assertEqual(conc['B'], 'ACD----')      # B lacks m2: gaps
            self.assertIn('LG, m1 = 1-3', open(f'{d}/out/concat.partition').read())


class Trees(unittest.TestCase):
    NWK = '((((a1,a2)90,b1)80,c1)60,(d1,d2)100,e1);'

    def tree(self, s=None):
        return Phylo.read(StringIO(s or self.NWK), 'newick')

    def test_monophyly(self):
        t = T.root_on(self.tree(), ['e1'])
        self.assertTrue(T.is_monophyletic(t, ['a1', 'a2'])[0])
        mono, extra = T.is_monophyletic(t, ['a1', 'c1'])
        self.assertFalse(mono)
        self.assertIn('b1', extra)

    def test_split_and_rooting_independent_of_file_root(self):
        t = self.tree()
        self.assertTrue(T.has_split(t, ['d1', 'd2']))
        self.assertTrue(T.has_split(t, ['d1', 'd2', 'e1']))     # the complement {a1,a2,b1,c1} is the clade as printed
        self.assertFalse(T.has_split(t, ['a1', 'c1']))
        r = T.root_on(t, ['d1', 'd2', 'e1'])
        self.assertTrue(T.is_monophyletic(r, ['d1', 'd2', 'e1'])[0])
        self.assertTrue(T.is_monophyletic(r, ['a1', 'a2', 'b1', 'c1'])[0])

    def test_rerooting_preserves_the_unrooted_tree(self):
        nwk = '((((a1:1,a2:1)90:1,b1:2)80:1,c1:3)60:1,(d1:1,d2:1)100:2,e1:4);'
        t = self.tree(nwk)
        for og in (['e1'], ['d1', 'd2'], ['a1', 'a2'], ['d1', 'd2', 'e1'], ['a1', 'a2', 'b1', 'c1']):
            r = T.root_on(t, og)
            self.assertEqual(len(r.root.clades), 2, og)
            self.assertEqual(T.rf(t, r)[0], 0, og)                                           # same unrooted splits
            self.assertAlmostEqual(r.total_branch_length(), t.total_branch_length(), msg=str(og))   # no length lost or added
            self.assertTrue(T.is_monophyletic(r, og)[0] or T.is_monophyletic(r, [x for x in T.tips(r.root) if x not in og])[0])
        r = T.root_on(t, ['d1', 'd2'])
        sup = {frozenset(T.tips(c)): c.confidence for c in r.get_nonterminals() if c.confidence is not None}
        self.assertEqual(sup[frozenset(['a1', 'a2'])], 90)                                   # supports travel with their branch

    def test_rf(self):
        t1 = self.tree()
        t2 = self.tree('((((a1,b1)90,a2)80,c1)60,(d1,d2)100,e1);')
        self.assertEqual(T.rf(t1, t1)[0], 0)
        d, m, n = T.rf(t1, t2)
        self.assertEqual(d, 2)                 # one split differs on each side
        self.assertEqual(n, 7)

    def test_collapse(self):
        with tempfile.TemporaryDirectory() as d:
            open(f'{d}/in.tre', 'w').write(self.NWK + '\n' + self.NWK + '\n')
            T.collapse_trees(f'{d}/in.tre', f'{d}/out.tre', 85)
            out = open(f'{d}/out.tre').read().strip().split('\n')
            self.assertEqual(len(out), 2)
            t = Phylo.read(StringIO(out[0]), 'newick')
            sup = [c.confidence for c in t.get_nonterminals() if c.confidence is not None]
            self.assertEqual(sorted(sup), [90, 100])          # 80 and 60 collapsed


if __name__ == '__main__':
    unittest.main()
