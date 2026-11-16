"""Fixture tests for loci_per_genome.py and strain_summary.py (run: python3 -m unittest test_reca_all)."""
import unittest

import loci_per_genome as L
import strain_summary as S


def row(q, scaf, aln, ev, bits, s, e):
    # columns as read by loci_from_hits: 0 query, 1 scaffold, 3 aln, 4 e-value, 5 bits, 9 start, 10 end
    return [q, scaf, '60', str(aln), str(ev), str(bits), '0', '0', '0', str(s), str(e)]


class Loci(unittest.TestCase):
    def test_merge_within_gap_and_split_beyond(self):
        rows = [row('CTRL_Rad51|a', 's1', 200, 1e-50, 300, 1000, 2000),
                row('RecA_DMC1|b', 's1', 200, 1e-40, 250, 5000, 6000),      # 3 kb away -> same locus
                row('CTRL_Rad51|a', 's1', 200, 1e-50, 290, 100000, 101000),  # far away -> second locus
                row('CTRL_Rad51|a', 's2', 100, 1e-50, 200, 1, 500),          # too short -> ignored
                row('CTRL_Rad51|a', 's2', 300, 1e-10, 90, 1, 900)]            # e-value too weak -> ignored
        loci = L.loci_from_hits(rows)
        self.assertEqual(len(loci), 2)
        first = [l for l in loci if l['start'] == 1000][0]
        self.assertEqual(first['bits'], {'CTRL_Rad51': 300.0, 'RecA_DMC1': 250.0})

    def test_unrelated_queries_ignored(self):
        self.assertEqual(L.loci_from_hits([row('Lig4|x', 's1', 300, 1e-80, 400, 1, 900)]), [])


class Summary(unittest.TestCase):
    def test_calls(self):
        cls = [{'locustag': 'A', 'call': 'assigned', 'clade': 'Rad51', 'mixed_labels': ''},
               {'locustag': 'A', 'call': 'assigned', 'clade': 'Dmc1', 'mixed_labels': ''},
               {'locustag': 'B', 'call': 'assigned', 'clade': 'Rad51', 'mixed_labels': ''},
               {'locustag': 'B', 'call': 'low_support', 'clade': 'Dmc1', 'mixed_labels': ''},
               {'locustag': 'C', 'call': 'ambiguous', 'clade': '', 'mixed_labels': 'Rad51,Dmc1'},
               {'locustag': 'D', 'call': 'dropped_by_trim', 'clade': '', 'mixed_labels': ''}]
        gen = [{'locustag': x, 'genus': 'Coemansia', 'grp': 'g', 'busco_complete': '80', 'n50_bp': '1'} for x in 'ABCDE']
        loci = [{'locustag': 'A', 'n_loci': '2', 'locus_bits': 'x'}]
        out = {r['locustag']: r for r in S.summarize(cls, gen, loci)}
        self.assertEqual(out['A']['dmc1_call'], 'present')
        self.assertEqual(out['B']['dmc1_call'], 'uncertain')
        self.assertEqual(out['C']['dmc1_call'], 'uncertain')
        self.assertEqual(out['D']['dmc1_call'], 'not_seen')
        self.assertEqual(out['D']['n_dropped_by_trim'], 1)
        self.assertEqual(out['E']['n_reca'], 0)
        self.assertEqual(out['A']['tblastn_loci'], '2')


if __name__ == '__main__':
    unittest.main()
