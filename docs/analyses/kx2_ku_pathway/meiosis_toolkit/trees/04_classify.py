#!/usr/bin/env python3
"""Step 4: classify the annotated candidates of each family by tree placement and tabulate clade presence per genome and group.

Usage: 04_classify.py <BENCH> <backbone_taxa.tsv> <coemansia_strains.txt>
For each family (muts, mutl, kleisin, horma) reads $BENCH/trees/<family>.treefile and <family>_members.tsv and writes in $BENCH/trees/:
  <family>_classification.tsv  every annotated member: call (assigned / low_support / ambiguous), clade, support, nearest reference, long-branch flag
  <family>_loo.tsv             leave-one-out test of the references (each reference classified as if unknown)
  <family>_by_genome.tsv       per genome and clade: number of members assigned / low_support / ambiguous and a presence call
  <family>_by_group.tsv        per group and clade: genomes present / uncertain / not seen
  classification_summary.txt   tree size, trimmed width, leave-one-out score per clade
Classification (fam_lib.classify, unchanged from ../../family_trees/): a member is assigned to a reference clade when the nearest ancestor
containing references holds only that clade's references and has UFBoot >= 90 and SH-aLRT >= 80; unsupported but single-clade ->
'low_support'; mixed -> 'ambiguous'. 'Long branch' = terminal branch longer than the 95th percentile in that tree.
Presence per genome and clade: 'present' = >= 1 assigned member; 'uncertain' = only low_support / ambiguous members that name the clade;
'not_seen' = no member assigned or named for that clade. Annotated layer only: a gene missed by annotation shows as 'not_seen'.
GATE: a clade's absence calls are used only if every reference of that clade is re-assigned correctly in leave-one-out (clades with a
single reference cannot be tested). The gate is applied in the README tables, not hidden here.
"""
import csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fam_lib as F
from collections import defaultdict

bench, bb_tsv, strains_txt = sys.argv[1], sys.argv[2], sys.argv[3]
tdir = f'{bench}/trees'
CORE = {'Coemansia', 'Kickxella', 'Linderina', 'Dipsacomyces', 'Martensiomyces'}
CLADES = {'muts': ['Msh2', 'Msh3', 'Msh6', 'Msh4', 'Msh5', 'Msh1'], 'mutl': ['Mlh1', 'Pms1', 'Mlh3', 'Mlh2'],
          'kleisin': ['Rad21', 'Rec8'], 'horma': ['Hop1', 'Mad2', 'Rev7']}


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
for l in open(strains_txt):
    if l.strip():
        info.setdefault(l.strip(), ('Coemansia', 'Kickxellales_core'))
genomes = sorted(info)


def pct(values, q):
    v = sorted(values)
    return v[min(len(v) - 1, int(q * len(v)))]


summary = []
for fam, clades in CLADES.items():
    if not os.path.exists(f'{tdir}/{fam}.treefile'):
        summary.append(f'{fam}: tree not finished, skipped')
        continue
    tree = F.parse_newick(open(f'{tdir}/{fam}.treefile').read())
    mem = {r['label']: r for r in csv.DictReader(open(f'{tdir}/{fam}_members.tsv'), delimiter='\t')}
    ref_clade = {l: r['clade_or_best_ref'] for l, r in mem.items() if r['kind'] == 'ref'}
    cands_all = [l for l, r in mem.items() if r['kind'] != 'ref']
    in_tree = {n.name for n in F.leaves(tree)}
    dropped = [l for l in cands_all if l not in in_tree]
    cands = [l for l in cands_all if l in in_tree]
    res = F.classify(tree, ref_clade, cands)
    long_thr = pct([n.bl for n in F.leaves(tree)], 0.95)
    with open(f'{tdir}/{fam}_classification.tsv', 'w') as fh:
        fh.write('label\tlocustag\tid\tlength\tcall\tclade\tn_refs\tufboot\tsh_alrt\tnearest_ref\tnearest_dist\tterminal_bl\tlong_branch\tlabels\n')
        for l in cands:
            r, m = res[l], mem[l]
            fh.write('\t'.join(str(x) for x in [l, m['locustag'], m['id'], m['length'], r['call'], r['clade'], r['n_refs'], r['ufboot'], r['sh'],
                                                  r['nearest_ref'], r['nearest_dist'], r['terminal_bl'], r['terminal_bl'] > long_thr, r['labels']]) + '\n')
    tested, correct, det = F.loo_check(tree, ref_clade)
    with open(f'{tdir}/{fam}_loo.tsv', 'w') as fh:
        fh.write('reference\tclade\tcall\tassigned_clade\tcorrect\n')
        for r in det:
            fh.write('\t'.join(str(x) for x in r) + '\n')
    per_clade_ok = {}
    for c in clades:
        rows = [d for d in det if d[1] == c]
        per_clade_ok[c] = (sum(1 for d in rows if d[4]), len(rows))
    bygen = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))          # (locustag) -> clade -> [assigned, low, ambiguous-naming]
    for l in cands:
        loc = mem[l]['locustag']
        r = res[l]
        if r['call'] == 'assigned':
            bygen[loc][r['clade']][0] += 1
        elif r['call'] == 'low_support' and r['clade']:
            bygen[loc][r['clade']][1] += 1
        else:
            for c in (r['labels'] or '').split(','):
                if c:
                    bygen[loc][c][2] += 1
    with open(f'{tdir}/{fam}_by_genome.tsv', 'w') as fh:
        fh.write('locustag\tgenus\tgroup\tclade\tn_assigned\tn_low_support\tn_ambiguous\tcall\n')
        for loc in genomes:
            for c in clades:
                a, lo, am = bygen[loc][c] if loc in bygen else (0, 0, 0)
                call = 'present' if a else ('uncertain' if (lo or am) else 'not_seen')
                fh.write(f'{loc}\t{info[loc][0]}\t{info[loc][1]}\t{c}\t{a}\t{lo}\t{am}\t{call}\n')
    grp = defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0]))
    for loc in genomes:
        g = info[loc][1]
        for c in clades:
            a, lo, am = bygen[loc][c] if loc in bygen else (0, 0, 0)
            idx = 0 if a else (1 if (lo or am) else 2)
            grp[g][c][idx] += 1
            grp[g][c][3] += 1
    with open(f'{tdir}/{fam}_by_group.tsv', 'w') as fh:
        fh.write('group\tclade\tn_genomes\tpresent\tuncertain\tnot_seen\n')
        for g in sorted(grp):
            for c in clades:
                p, u, n, tot = grp[g][c]
                fh.write(f'{g}\t{c}\t{tot}\t{p}\t{u}\t{n}\n')
    summary.append(f"{fam}: tree {len(in_tree)} leaves ({len(dropped)} members dropped by trimming); leave-one-out {correct}/{tested} references re-assigned; "
                   f"per clade {', '.join(f'{c} {a}/{b}' for c, (a, b) in per_clade_ok.items())}; long-branch threshold {long_thr:.3f}")
open(f'{tdir}/classification_summary.txt', 'w').write('\n'.join(summary) + '\n')
print('\n'.join(summary))
