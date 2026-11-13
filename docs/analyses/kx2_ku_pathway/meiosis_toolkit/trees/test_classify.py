"""Fixture test for the tree-placement logic as used in 04_classify.py (hand-made tree with known answers; no real data).
Run: python3 -I test_classify.py
Tree: two clades of references (A: a1 a2 a3; B: b1 b2 b3) with a supported split; candidates: c_a nests in A, c_b nests in B,
c_mix sits between the clades (ambiguous), c_low nests in A behind a weakly supported node.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fam_lib as F  # noqa: E402

NWK = ('(((a1:0.1,a2:0.1,c_a:0.1)99/100:0.2,((a3:0.1,c_low:0.3)40/50:0.1,a4:0.1)95/99:0.1)97/100:0.3,'
       '((b1:0.1,b2:0.1,c_b:0.1)98/100:0.2,b3:0.1)96/100:0.3,c_mix:0.5);')
REF = {'a1': 'A', 'a2': 'A', 'a3': 'A', 'a4': 'A', 'b1': 'B', 'b2': 'B', 'b3': 'B'}


def test_assignment():
    root = F.parse_newick(NWK)
    res = F.classify(root, REF, ['c_a', 'c_b', 'c_mix', 'c_low'])
    assert res['c_a']['call'] == 'assigned' and res['c_a']['clade'] == 'A'
    assert res['c_b']['call'] == 'assigned' and res['c_b']['clade'] == 'B'
    assert res['c_mix']['call'] == 'ambiguous'
    # c_low sits in a weak node of single-clade references; the walk climbs to the supported A node
    assert res['c_low']['clade'] == 'A'


def test_leave_one_out_recovers_references():
    root = F.parse_newick(NWK)
    tested, correct, det = F.loo_check(root, REF)
    assert tested == 7
    assert correct >= 5            # references nested in supported single-clade nodes re-assign; b3 sits behind the mixed root and fails
    wrong = {d[0] for d in det if not d[4]}
    assert 'b3' in wrong or correct == 7


if __name__ == '__main__':
    test_assignment()
    test_leave_one_out_recovers_references()
    print('classification fixture tests passed')
