"""Fixture tests for the T5 scripts (toy trees with known answers; no HPCC data, no network).
Run:  python3 -m pytest -q test_t5.py     (on the HPCC: /usr/bin/python3.12 -m pytest -q test_t5.py; needs Biopython for the Dollo tests)
"""
import csv
import importlib.util
import io
import math
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'kx2_ku_pathway', 'family_trees'))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', '..', 'workflows', 'phylo'))
sys.setrecursionlimit(100000)


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace('.py', '').replace('0', 'm0', 1) if name[0].isdigit() else name, os.path.join(HERE, name))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------ 01 references
def uni(entry, org, seq, ec='', lineage='Eukaryota, Fungi'):
    return {'Entry': entry, 'Organism': org, 'Sequence': seq, 'Length': str(len(seq)), 'EC number': ec, 'Taxonomic lineage': lineage,
            'Reviewed': 'reviewed', 'Protein names': 'x', 'Organism (ID)': '1'}


def test_group_of_and_duplicates_and_context():
    m = load('01_fetch_references.py')
    assert m.group_of('3.2.1.28', 'Fungi') == 'TREHALASE_FUNGI'
    assert m.group_of('3.2.1.107', 'Metazoa') == 'ATHL1_METAZOA'
    assert m.group_of('2.4.1.64', '') == 'TREHALOSE_PHOSPHORYLASE'
    assert m.group_of('3.2.1.-', '') == 'GH65_OTHER'
    s1 = 'A' * 100
    s2 = 'A' * 99 + 'C'                       # 99% identical, same species -> collapsed
    s3 = 'C' * 100                              # another species -> kept
    rev = [uni('P2', 'Saccharomyces cerevisiae (strain X)', s2, '3.2.1.28'), uni('P1', 'Saccharomyces cerevisiae (strain Y)', s1, '3.2.1.28'),
           uni('P3', 'Candida albicans', s3, '3.2.1.28')]
    kept = m.collapse_near_duplicates(rev)
    assert [r['Entry'] for r in kept] == ['P1', 'P3']      # the lower accession of the near-duplicate pair survives
    pool = [uni(f'B{i}', f'Genus{chr(97 + i % 3)} sp{i}', 'D' * 50, lineage='Bacteria, Firmicutes') for i in range(9)]
    pool += [uni('F1', 'Aspergillus niger', 'D' * 50, lineage='Eukaryota, Fungi'), uni('C1', 'Candida zzz', 'D' * 50, lineage='Bacteria')]
    out = m.select(rev, {'bact': (pool, 2)})
    ctx = [r for r, d, k, role in out if role == 'context']
    assert len(ctx) == 2                                     # cap
    assert len({m.genus_of(r['Organism']) for r in ctx}) == 2   # one per genus
    assert all(r['Entry'] != 'F1' for r in ctx)              # fungi never come from UniProt context
    assert all(m.genus_of(r['Organism']) != 'Candida' for r in ctx)   # a genus with a reference is not reused


# ------------------------------------------------------------------ 04 inputs
def test_build_inputs():
    m = load('04_build_inputs.py')
    with tempfile.TemporaryDirectory() as d:
        os.makedirs(f'{d}/refs')
        os.makedirs(f'{d}/cand')
        open(f'{d}/refs/refs.faa', 'w').write('>TREHALASE_FUNGI|P1|Sce\nAAAA\n>CONTEXT_BACTERIA|B1|Bsub\nCCCC\n')
        open(f'{d}/cand/gh65_annotated.faa', 'w').write('>ANN|F1_1|F1\n' + 'M' * 500 + '\n>ANN|F2_1|F2\n' + 'M' * 100 + '\n')
        nref, nctx, nann, nex = m.build(d, 400)
        assert (nref, nctx, nann, nex) == (1, 1, 1, 1)
        rows = list(csv.DictReader(open(f'{d}/trees/gh65_members.tsv'), delimiter='\t'))
        assert [r['label'] for r in rows] == ['REF_1', 'CTX_1', 'ANN_1']
        assert rows[0]['clade_or_group'] == 'TREHALASE_FUNGI' and rows[2]['id'] == 'F1_1'
        assert open(f'{d}/trees/excluded_short.tsv').read().count('F2_1') == 1


# ------------------------------------------------------------------ 06 classification and rules
MEMBERS = ([('REF_%d' % i, 'ref', '', 'Q%d' % i, g, 900, '') for i, g in
            [(1, 'TREHALASE_FUNGI'), (2, 'TREHALASE_FUNGI'), (3, 'TREHALASE_FUNGI'), (4, 'ATHL1_METAZOA'), (5, 'ATHL1_METAZOA'),
             (6, 'KOJIBIOSE_PHOSPHORYLASE'), (7, 'KOJIBIOSE_PHOSPHORYLASE')]] +
           [('CTX_1', 'context', '', 'C1', 'CONTEXT_BACTERIA', 700, '')] +
           [('ANN_%d' % i, 'annotated', 'L%d' % i, 'P%d' % i, '', 900, '') for i in range(1, 5)])
MEMBERS = [dict(zip(['label', 'kind', 'locustag', 'id', 'clade_or_group', 'length', 'info'], m)) for m in MEMBERS]
REFS = {'Q1': {'domain': 'Eukaryota'}, 'Q2': {'domain': 'Eukaryota'}, 'Q3': {'domain': 'Eukaryota'}, 'Q4': {'domain': 'Eukaryota'},
        'Q5': {'domain': 'Eukaryota'}, 'Q6': {'domain': 'Bacteria'}, 'Q7': {'domain': 'Bacteria'}, 'C1': {'domain': 'Bacteria'}}
FULL = 'PF03636;PF03632;PF03633'


def genomes(phyla):
    return {f'L{i + 1}': {'genus': f'G{i + 1}', 'PHYLUM': p, 'class': f'C{i + 1}', 'ord': 'o', 'busco_complete': '95', 'N50_bp': '200000'} for i, p in enumerate(phyla)}


def cands():
    return {f'P{i}': {'domains': FULL, 'transcript_id': f'L{i}_t'} for i in range(1, 5)}


T = '99/100'
TCLADE = f'((REF_1:1,REF_2:1){T}:1,(REF_3:1,ANN_1:1){T}:1,(ANN_2:1,ANN_3:1){T}:1){T}:1'


def run(newick, phyla=('Ascomycota', 'Basidiomycota', 'Ascomycota', 'Ascomycota')):
    m = load('06_classify.py')
    tree = m.GTree(newick)
    return m.analyse(tree, MEMBERS, REFS, genomes(phyla), cands(), log=lambda s: None)[0]


def test_vertical_scenario_all_rules():
    # fungal clade next to the animal glucosidases, bacteria further out; ANN_4 sits with a bacterium
    nw = f'(({TCLADE},(REF_4:1,REF_5:1){T}:1){T}:1,(REF_6:1,REF_7:1){T}:1,(CTX_1:1,ANN_4:1)50/50:1);'
    r = run(nw)
    assert r['gate']['passed'] and r['gate']['trehalase_loo_correct'] == 3
    assert r['R1']['assigned_trehalase'] == 3 and r['R1']['unresolved'] == 1 and r['R1']['fraction_other'] == 0
    assert r['R1']['verdict'].startswith('survives')
    assert r['R2']['tclade_tips'] == 6 and r['R2']['tclade_nonfungal_tips'] == 0 and r['R2']['verdict'].startswith('supported')
    assert r['R3']['non_dikarya_counting'] == 0 and r['R3']['verdict'].startswith('supported')
    assert r['R5a']['verdict'].startswith('vertical-consistent') and r['R5a']['neighbour_types'] == {'Eukaryota': 2}


def test_hgt_neighbours_and_older_origin():
    nw = f'(({TCLADE},(REF_6:1,REF_7:1){T}:1){T}:1,(REF_4:1,REF_5:1){T}:1,(CTX_1:1,ANN_4:1)50/50:1);'
    r = run(nw, phyla=('Ascomycota', 'Basidiomycota', 'Mortierellomycota', 'Ascomycota'))
    assert r['R5a']['verdict'].startswith('HGT-consistent') and r['R5a']['neighbour_types'] == {'Bacteria': 2}
    assert r['R3']['non_dikarya_counting'] == 1 and 'OLDER' in r['R3']['verdict']        # ANN_3 is in a Mortierellomycota genome, quality, full architecture


def test_nonfungal_tip_inside_rejects_single_origin():
    inner = f'((REF_1:1,REF_2:1){T}:1,(REF_3:1,ANN_1:1){T}:1,(ANN_2:1,ANN_3:1){T}:1,CTX_1:1){T}:1'
    nw = f'(({inner},(REF_4:1,REF_5:1){T}:1){T}:1,(REF_6:1,REF_7:1){T}:1,ANN_4:1);'
    r = run(nw)
    assert r['R2']['tclade_nonfungal_tips'] == 1 and r['R2']['verdict'].startswith('REJECTED')


def test_gate_fails_when_a_reference_is_misplaced():
    # REF_3 (fungal trehalase) nests with the bacterial references: leave-one-out cannot re-assign it
    nw = f'((((REF_1:1,REF_2:1){T}:1,(ANN_1:1,ANN_2:1,ANN_3:1){T}:1){T}:1,(REF_4:1,REF_5:1){T}:1){T}:1,((REF_6:1,REF_7:1,REF_3:1){T}:1,CTX_1:1,ANN_4:1){T}:1);'
    r = run(nw)
    assert not r['gate']['passed']
    assert r['R1']['verdict'].startswith('tree cannot decide') and r['R2']['verdict'].startswith('tree cannot decide')


def test_r1_killed_when_many_proteins_assign_elsewhere():
    # ANN_1..3 sit with the kojibiose references, ANN_4 with the animal references: 4 of 4 assigned to non-trehalase groups
    nw = f'((REF_1:1,REF_2:1,REF_3:1){T}:1,((REF_6:1,REF_7:1,ANN_1:1,ANN_2:1,ANN_3:1){T}:1,(REF_4:1,REF_5:1,ANN_4:1){T}:1){T}:1);'
    r = run(nw)
    assert r['gate']['passed']
    assert r['R1']['assigned_other'] == 4 and r['R1']['verdict'].startswith('KILLED')


def test_rooting_free_smallest_side():
    m = load('06_classify.py')
    a = m.GTree(f'((A:1,B:1){T}:1,(C:1,D:1){T}:1,E:1);')
    b = m.GTree(f'(((A:1,B:1){T}:1,E:1){T}:1,C:1,D:1);')       # same unrooted tree printed with another root
    for t in (a, b):
        mask, node, which = m.smallest_side_containing(t, t.mask_of(['A', 'B']))
        assert sorted(t.tips_in(mask)) == ['A', 'B']


# ------------------------------------------------------------------ 07 Dollo
def arrays(newick):
    from Bio import Phylo
    m = load('07_dollo.py')
    tree = Phylo.read(io.StringIO(newick), 'newick')
    names, parent, children, order = m.to_arrays(tree.root)
    return m, names, children, order


def states(names, spec):
    return [spec.get(n, -1) if n else -1 for n in names]


def test_dollo_losses_known_answer():
    m, names, children, order = arrays('(((A,B),(C,D)),((E,F),(G,H)));')
    st = states(names, {'A': 1, 'E': 1, 'B': 0, 'C': 0, 'D': 0, 'F': 0, 'G': 0, 'H': 0})
    n, g, losses = m.dollo_losses(children, order, st, want_clades=True)
    assert n == 4 and g == order[-1]                                   # gain at the root; B, (C,D), F, (G,H)
    st2 = states(names, {'A': 1, 'E': 1, 'B': 0, 'F': 0, 'G': 0, 'H': 0})   # C, D unknown: that clade is not a loss
    assert m.dollo_losses(children, order, st2)[0] == 3
    st3 = states(names, {'A': 1, 'B': 1, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0})   # present in one clade only: gain placed there, one loss clade outside is NOT counted
    n3, g3, _ = m.dollo_losses(children, order, st3, want_clades=True)
    assert n3 == 0 and len(m.tips_below(children, g3)) == 2
    assert m.dollo_losses(children, order, states(names, {'A': 0, 'B': 0}))[0] == 0   # nothing present: no gain, no losses


def test_permutation_null_is_deterministic_and_bounded():
    m, names, children, order = arrays('(((A,B),(C,D)),((E,F),(G,H)));')
    st = states(names, {'A': 1, 'B': 1, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0})
    a = m.permutation_null(children, order, st, n_perm=50, seed=1)
    b = m.permutation_null(children, order, st, n_perm=50, seed=1)
    assert a == b and len(a) == 50 and max(a) <= 6 and min(a) >= 0
    allpresent = states(names, {n: 1 for n in 'ABCDEFGH'})
    assert set(m.permutation_null(children, order, allpresent, n_perm=5)) == {0}


def test_verdict_rules():
    m = load('07_dollo.py')
    mk = lambda losses, classes: [{'scoring': 'scored', 'tree': 't', 'losses': losses, 'classes_touched': classes}]
    assert m.verdict(mk(2, 2)).startswith('R4 REJECTED')
    assert m.verdict(mk(5, 3)).startswith('R4 supported')
    assert m.verdict(mk(5, 1)).startswith('R4 inconclusive')


# ------------------------------------------------------------------ 08 composition
def test_zscore_percentile_and_verdict():
    m = load('08_hgt_gc.py')
    assert m.zscore(5, 3, 1) == 2 and math.isnan(m.zscore(1, 1, 0))
    assert m.percentile_rank(3, [1, 2, 3, 4]) == 0.75
    rows = [{'locustag': f'L{i}', 'z': (3.0 if i < 4 else 0.1), 'quality': True} for i in range(10)]      # 4 of 10 genomes beyond 2 sd
    s, v = m.summarize(rows, 'HGT-consistent: ...')
    assert s['quality_genomes']['share_beyond_2sd'] == 0.4 and v.startswith('HGT supported')
    assert m.summarize(rows, 'vertical-consistent: ...')[1].startswith('inconclusive')
    calm = [{'locustag': f'L{i}', 'z': 0.2, 'quality': True} for i in range(10)]
    assert m.summarize(calm, 'vertical-consistent: ...')[1].startswith('vertical origin supported')
    assert m.summarize(calm, 'HGT-consistent: ...')[1].startswith('inconclusive')
