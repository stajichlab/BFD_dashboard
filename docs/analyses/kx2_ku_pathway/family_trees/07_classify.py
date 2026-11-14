#!/usr/bin/env python3
"""Step 7: classify tree members by placement and tabulate clade presence per genome and group.

Usage: 07_classify.py <BENCH> <backbone_taxa.tsv> <strain genomes.tsv> [tblastn_round2_calls.tsv]
For each family (lig, reca, polx) reads trees/<family>.treefile and trees/<family>_members.tsv and writes, in trees/:
  <family>_classification.tsv  every non-reference member with its call (assigned / low_support / ambiguous), clade, support, nearest reference
  <family>_loo.tsv             leave-one-out test of the references (each reference classified as if unknown)
  <family>_by_genome.tsv       per genome and clade: members assigned (annotated / rescued) and a presence call
  <family>_by_group.tsv        per group and clade: genomes with the clade / genomes
Classification rule (fam_lib.classify): a member is assigned to a reference clade when the nearest ancestor containing references holds
only that clade's references and has UFBoot >= 90 and SH-aLRT >= 80; unsupported but single-clade -> 'low_support'; mixed -> 'ambiguous'.
'Long branch' = terminal branch longer than the 95th percentile of all terminal branches in that tree.
Presence per genome: in the 31 genomes searched with miniprot, 'present' = at least one assigned member (annotated or rescued),
'uncertain' = only low_support / ambiguous members, 'absent' = no member of the family tree placed near that clade at all.
The other genomes were not searched for missed loci, so their calls are 'annotation_only' (present / not seen in annotation).
"""
import csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fam_lib as F
from collections import defaultdict

bench, bb_tsv, st_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
tb_calls = sys.argv[4] if len(sys.argv) > 4 else None
tdir = f"{bench}/trees_{os.environ.get('LAYER', 'annotated')}"     # annotated (references + annotated BFD proteins) or rescue (adds miniprot models)
CORE = {'Coemansia', 'Kickxella', 'Linderina', 'Dipsacomyces', 'Martensiomyces'}
CLADES = {'lig': ['LIG1', 'LIG3', 'LIG4', 'LIG6'],
          'reca': ['Rad51', 'Dmc1', 'Rad55', 'Rad57', 'Rad51B', 'Rad51C', 'Rad51D', 'XRCC2', 'XRCC3', 'RadA'],
          'polx': ['Pol4', 'PolLambda', 'PolMu', 'TdT', 'PolBeta']}


def group_of(genus, ord_, phylum, role):
    if phylum == 'Kickxellomycota':
        if ord_ == 'Harpellales':
            return 'Harpellales'
        if ord_ == 'Dimargaritales':
            return 'Dimargaritales'
        if ord_ == 'Ramicandelaberales':
            return 'Ramicandelaber'
        return 'Kickxellales_core' if genus in CORE else 'Mycoemilia_Spiromyces'
    if ord_ == 'Zoopagales':
        return 'Zoopagales'
    if phylum in ('Entomophthoromycota', 'Basidiobolomycota'):
        return 'Entomophthorales_Basidiobolus'
    return 'Mucoromycota_s.l.' if role == 'mucoromycot' else 'Root'


info = {}
for r in csv.DictReader(open(bb_tsv), delimiter='\t'):
    info[r['LOCUSTAG']] = (r['GENUS'], group_of(r['GENUS'], r['ord'], r['PHYLUM'], r['role']))
searched = set()
for r in csv.DictReader(open(st_tsv), delimiter='\t'):
    searched.add(r['locustag'])
    info.setdefault(r['locustag'], (r['genus'], r['grp']))
genomes = sorted(info)


def pct(values, q):
    v = sorted(values)
    return v[min(len(v) - 1, int(q * len(v)))]


summary = []
for fam in ('lig', 'reca', 'polx'):
    if not os.path.exists(f'{tdir}/{fam}.treefile'):
        summary.append(f'{fam}: tree not finished, skipped')
        continue
    tree = F.parse_newick(open(f'{tdir}/{fam}.treefile').read())
    mem = {r['label']: r for r in csv.DictReader(open(f'{tdir}/{fam}_members.tsv'), delimiter='\t')}
    ref_clade = {l: r['clade_or_best_ref'] for l, r in mem.items() if r['kind'] == 'ref'}
    cands_all = [l for l, r in mem.items() if r['kind'] != 'ref']
    in_tree = {n.name for n in F.leaves(tree)}
    dropped = [l for l in cands_all if l not in in_tree]            # removed by the alignment trimming (mostly-gap sequences)
    cands = [l for l in cands_all if l in in_tree]
    res = F.classify(tree, ref_clade, cands)
    long_thr = pct([n.bl for n in F.leaves(tree)], 0.95)
    # ---- classification table
    with open(f'{tdir}/{fam}_classification.tsv', 'w') as f:
        f.write('label\tkind\tlocustag\tgenus\tgroup\tid\tlength\tcall\tclade\tn_refs\tufboot\tsh\tnearest_ref\tnearest_dist\tterminal_bl\tlong_branch\tmixed_labels\tinfo\n')
        for l in cands:
            m, r = mem[l], res[l]
            g, grp = info.get(m['locustag'], ('?', '?'))
            nref = mem[r['nearest_ref']]['id'] + '|' + mem[r['nearest_ref']]['clade_or_best_ref'] if r['nearest_ref'] else ''
            f.write('\t'.join(str(x) for x in [l, m['kind'], m['locustag'], g, grp, m['id'], m['length'], r['call'], r['clade'], r['n_refs'],
                    r['ufboot'], r['sh'], nref, '%.3f' % r['nearest_dist'], '%.3f' % r['terminal_bl'], r['terminal_bl'] > long_thr,
                    r['labels'], m['info']]) + '\n')
        for l in dropped:
            m = mem[l]
            g, grp = info.get(m['locustag'], ('?', '?'))
            f.write('\t'.join(str(x) for x in [l, m['kind'], m['locustag'], g, grp, m['id'], m['length'], 'dropped_by_trim', '', '', '', '', '', '', '', '', '', m['info']]) + '\n')
    # ---- leave-one-out
    tested, correct, det = F.loo_check(tree, ref_clade)
    with open(f'{tdir}/{fam}_loo.tsv', 'w') as f:
        f.write('reference\tclade\tcall\tassigned_clade\tcorrect\tacc\n')
        for r, lab, call, got, ok in det:
            f.write(f"{r}\t{lab}\t{call}\t{got}\t{ok}\t{mem[r]['id']}\n")
    summary.append(f'{fam}: tree {len(F.leaves(tree))} leaves ({len(dropped)} members dropped by trimming); leave-one-out {correct}/{tested} '
                   f'references re-assigned to their own clade; long-branch threshold {long_thr:.3f}')
    # ---- per genome
    by = defaultdict(lambda: defaultdict(lambda: {'annotated': 0, 'rescued': 0, 'weak': 0}))
    seen = defaultdict(int)
    for l in cands:
        m, r = mem[l], res[l]
        seen[m['locustag']] += 1
        if r['call'] == 'assigned':
            by[m['locustag']][r['clade']][m['kind']] += 1
        else:
            by[m['locustag']][r['clade'] or 'unplaced']['weak'] += 1
    with open(f'{tdir}/{fam}_by_genome.tsv', 'w') as f:
        f.write('locustag\tgenus\tgroup\tclade\tn_annotated\tn_rescued\tn_uncertain\tcall\tsearched_for_rescue\n')
        for loc in genomes:
            g, grp = info[loc]
            for cl in CLADES[fam]:
                d = by[loc].get(cl, {'annotated': 0, 'rescued': 0, 'weak': 0})
                if d['annotated'] + d['rescued'] > 0:
                    call = 'present'
                elif d['weak'] > 0:
                    call = 'uncertain'
                else:
                    call = 'absent' if loc in searched else 'not_seen_in_annotation'
                f.write(f"{loc}\t{g}\t{grp}\t{cl}\t{d['annotated']}\t{d['rescued']}\t{d['weak']}\t{call}\t{loc in searched}\n")
    # ---- per group
    grp_g = defaultdict(set)
    for loc in genomes:
        grp_g[info[loc][1]].add(loc)
    rows = list(csv.DictReader(open(f'{tdir}/{fam}_by_genome.tsv'), delimiter='\t'))
    with open(f'{tdir}/{fam}_by_group.tsv', 'w') as f:
        f.write('group\tclade\tn_genomes\tpresent\tpresent_annotated_only\tpresent_rescued\tuncertain\tabsent_or_not_seen\tn_searched_for_rescue\n')
        for grp in sorted(grp_g):
            for cl in CLADES[fam]:
                rr = [r for r in rows if r['group'] == grp and r['clade'] == cl]
                pres = [r for r in rr if r['call'] == 'present']
                resc = [r for r in pres if int(r['n_rescued']) > 0]
                ann = [r for r in pres if int(r['n_rescued']) == 0]
                unc = [r for r in rr if r['call'] == 'uncertain']
                ab = [r for r in rr if r['call'] in ('absent', 'not_seen_in_annotation')]
                ns = sum(1 for r in rr if r['searched_for_rescue'] == 'True')
                f.write(f'{grp}\t{cl}\t{len(rr)}\t{len(pres)}\t{len(ann)}\t{len(resc)}\t{len(unc)}\t{len(ab)}\t{ns}\n')

with open(f'{tdir}/classification_summary.txt', 'w') as f:
    f.write('\n'.join(summary) + '\n')
print('\n'.join(summary))

# ---- comparison with the earlier tblastn round 2 calls (Lig4)
if tb_calls and os.path.exists(tb_calls) and os.path.exists(f'{tdir}/lig_by_genome.tsv'):
    tb = {r['locustag']: r for r in csv.DictReader(open(tb_calls), delimiter='\t') if r['gene'] == 'Lig4'}
    bg = {(r['locustag'], r['clade']): r for r in csv.DictReader(open(f'{tdir}/lig_by_genome.tsv'), delimiter='\t')}
    with open(f'{tdir}/lig4_vs_tblastn.tsv', 'w') as f:
        f.write('locustag\tgenus\tgroup\ttblastn_round2_Lig4\ttree_LIG4\ttree_LIG1\ttree_LIG3\n')
        for loc in sorted(tb):
            g = info.get(loc, ('?', '?'))
            f.write(f"{loc}\t{g[0]}\t{g[1]}\t{tb[loc]['call']}\t{bg[(loc, 'LIG4')]['call']}\t{bg[(loc, 'LIG1')]['call']}\t{bg[(loc, 'LIG3')]['call']}\n")
    print('wrote lig4_vs_tblastn.tsv')
