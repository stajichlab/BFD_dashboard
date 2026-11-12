"""Fixture tests for fam_lib (tiny made-up sequences and trees; no real data). Run: python3 -I tests/test_fam_lib.py"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import fam_lib as F  # noqa: E402


def test_translate_and_revcomp():
    assert F.translate('ATGAAATTTTAA') == 'MKF*'
    assert F.translate('ATGNNNTTT') == 'MXF'
    assert F.translate('ATGAAAT') == 'MK'                      # trailing partial codon dropped
    assert F.revcomp('AACG') == 'CGTT'


# protein MKFG*: ATG AAA TTT GGT TAA, split as exon1 ATGAA + intron + exon2 ATTTGGTTAA
EX1, INTRON, EX2 = 'ATGAA', 'GTAAGTTTTAG', 'ATTTGGTTAA'
FWD = 'CCC' + EX1 + INTRON + EX2 + 'CCC'


def test_model_protein_forward_spliced():
    genome = {'s1': FWD}
    # block 1 ATGAA leaves 2 bases of codon AAA pending; block 2 starts with the 1 base that completes it, so its phase is 1
    m = {'scaf': 's1', 'strand': '+', 'cds': [(4, 8, 0), (20, 29, 1)]}
    assert F.model_protein(m, genome) == ('MKFG', 0, True, 0)


def test_model_protein_reverse_spliced():
    L = len(FWD)
    genome = {'s1': F.revcomp(FWD)}
    # forward pieces (4,8) and (20,29) map to (L-8+1, L-4+1) = (25,29) and (L-29+1, L-20+1) = (4,13); the first block in transcript order
    # is the high-coordinate one (phase 0), the second has phase 1
    m = {'scaf': 's1', 'strand': '-', 'cds': [(L - 29 + 1, L - 20 + 1, 1), (L - 8 + 1, L - 4 + 1, 0)]}
    assert F.model_protein(m, genome) == ('MKFG', 0, True, 0)


def test_model_protein_frameshift_is_counted_not_translated_across():
    genome = {'s1': FWD}
    # block 2 given phase 2 although the pending 2 bases need a lead of 1: the frame changed there. The incomplete codon is dropped, the
    # remaining 8 bases of block 2 (TTGGTTAA) give L V with 2 bases left over.
    m = {'scaf': 's1', 'strand': '+', 'cds': [(4, 8, 0), (20, 29, 2)]}
    prot, stops, ok, breaks = F.model_protein(m, genome)
    assert (prot, stops, ok, breaks) == ('MLV', 0, True, 1), (prot, stops, ok, breaks)


def test_norm_strand():
    assert [F.norm_strand(x) for x in ('1', '-1', '+', '-', 1, -1)] == ['+', '-', '+', '-', '+', '-']


def test_dedup_against_annotated_with_gene_info_strand():
    # annotated gene (gene_info: strand '1') vs a miniprot model on '+': the real F0DE0BF4 case that was missed before the fix
    ann = {('L', 's', F.norm_strand('1')): [(351610, 356721)]}
    model = (351820, 356563)
    assert any(F.overlap_frac(model, a) >= 0.5 for a in ann.get(('L', 's', '+'), []))


def test_max_intron():
    assert F.max_intron({'cds': [(4, 8, 0), (20, 29, 2)]}) == 11
    assert F.max_intron({'cds': [(1, 10, 0), (21, 30, 0), (100, 120, 0)]}) == 69
    assert F.max_intron({'cds': [(1, 10, 0)]}) == 0


def test_internal_stop_is_counted():
    genome = {'s': 'ATGTAAGGTTGA'}
    m = {'scaf': 's', 'strand': '+', 'cds': [(1, 12, 0)]}
    prot, stops, _, breaks = F.model_protein(m, genome)
    assert prot == 'M*G' and stops == 1 and breaks == 0, (prot, stops)


def test_phase_skips_leading_bases():
    genome = {'s': 'GGATGAAATAA'}                              # first CDS piece starts 2 bases before the ATG, phase 2
    m = {'scaf': 's', 'strand': '+', 'cds': [(1, 11, 2)]}
    prot, stops, _, breaks = F.model_protein(m, genome)
    assert prot == 'MK' and stops == 0 and breaks == 0, prot


def test_parse_gff_and_loci():
    gff = ['##gff-version 3',
           's1\tminiprot\tmRNA\t4\t29\t100\t+\t.\tID=MP1;Rank=1;Identity=0.8;Positive=0.9;Target=P1 1 40',
           's1\tminiprot\tCDS\t4\t8\t100\t+\t0\tParent=MP1;Rank=1;Target=P1 1 2',
           's1\tminiprot\tCDS\t20\t29\t100\t+\t1\tParent=MP1;Rank=1;Target=P1 2 5',
           's1\tminiprot\tmRNA\t5\t28\t50\t+\t.\tID=MP2;Rank=2;Identity=0.5;Positive=0.6;Target=P2 1 40',
           's1\tminiprot\tmRNA\t5\t28\t50\t-\t.\tID=MP3;Rank=1;Identity=0.5;Positive=0.6;Target=P2 1 40',
           's1\tminiprot\tmRNA\t100\t200\t60\t+\t.\tID=MP4;Rank=1;Identity=0.5;Positive=0.6;Target=P2 1 40']
    models = F.parse_miniprot_gff(gff)
    assert len(models) == 4 and models[0]['cds'] == [(4, 8, 0), (20, 29, 1)] and models[0]['target'] == 'P1'
    kept = {m['id'] for m in F.greedy_loci(models)}
    # MP2 overlaps the better MP1 on the same strand (dropped); MP3 is the other strand; MP4 is elsewhere
    assert kept == {'MP1', 'MP3', 'MP4'}, kept


T1 = ('((A1:0.1,A2:0.1)100/100:0.1,((B1:0.1,B2:0.1)99/99:0.1,(c1:0.1,B3:0.1)95/96:0.1)90/92:0.1,c2:0.5);')


def test_newick_support():
    r = F.parse_newick(T1)
    assert sorted(l.name for l in F.leaves(r)) == ['A1', 'A2', 'B1', 'B2', 'B3', 'c1', 'c2']
    assert F.parse_support('87.3/95') == (87.3, 95.0)
    assert F.parse_support('') == (None, None)


def test_classify_assigned_ambiguous():
    r = F.parse_newick(T1)
    refs = {'A1': 'A', 'A2': 'A', 'B1': 'B', 'B2': 'B', 'B3': 'B'}
    res = F.classify(r, refs, ['c1', 'c2'])
    assert res['c1']['call'] == 'assigned' and res['c1']['clade'] == 'B', res['c1']
    assert res['c2']['call'] == 'ambiguous', res['c2']         # sister to a mixed A/B group
    assert res['c1']['nearest_ref'] == 'B3'


def test_classify_low_support():
    r = F.parse_newick('((A1:0.1,A2:0.1,c4:0.1)50/60:0.1,(B1:0.1,B2:0.1)100/100:0.1,(C1:0.1,C2:0.1)100/100:0.1);')
    refs = {'A1': 'A', 'A2': 'A', 'B1': 'B', 'B2': 'B', 'C1': 'C', 'C2': 'C'}
    res = F.classify(r, refs, ['c4'])['c4']
    assert res['call'] == 'low_support' and res['clade'] == 'A', res


def test_leave_one_out():
    r = F.parse_newick(T1)
    refs = {'A1': 'A', 'A2': 'A', 'B1': 'B', 'B2': 'B', 'B3': 'B'}
    tested, correct, det = F.loo_check(r, refs)
    assert tested == 5 and correct == 5, det


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for fn in fns:
        fn()
        print('ok', fn.__name__)
    print(len(fns), 'tests passed')
