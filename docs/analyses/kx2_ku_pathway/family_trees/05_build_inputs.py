#!/usr/bin/env python3
"""Step 5: combine references, annotated candidates and miniprot-rescued models into one tree input per family.

Usage: 05_build_inputs.py <BENCH> [layer]     layer = annotated (default) or rescue
Writes $BENCH/trees_<layer>/<family>_input.faa (unique labels REF_n / ANN_n / MP_n, safe for tree software) and
$BENCH/trees_<layer>/<family>_members.tsv (label -> kind, locustag, protein id or model position, clade or best reference, length, domains).
Layer 'annotated' = references + annotated BFD candidates only (does not depend on the miniprot rescue). Layer 'rescue' = the same plus the
rescued miniprot models, tagged kind=rescued, which is only meaningful if the rescue passed validation. Genomes whose miniprot task did not
finish (no <locustag>.done) are written to $BENCH/trees_<layer>/rescue_failed_genomes.txt and are never silently treated as 'no model'.
Rules, fixed before any result was seen:
  - annotated candidates and rescued models must be >= MIN_LEN aa (ligase 300, RecA/Rad51 200, Pol X 200), as in step 4b
  - a miniprot model is a duplicate (dropped) if it overlaps an annotated candidate of the same family on the same scaffold and strand by
    >= 50% of the shorter locus (gene_info coordinates); otherwise it is 'rescued' and kept
  - references are always kept
Also reports how many annotated candidates' scaffold names match miniprot scaffold names (a naming check for the overlap test).
"""
import csv, glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fam_lib as F
from collections import defaultdict

bench = sys.argv[1]
layer = sys.argv[2] if len(sys.argv) > 2 else 'annotated'
assert layer in ('annotated', 'rescue'), layer
MIN_LEN = {'lig': 300, 'reca': 200, 'polx': 200}
tdir = f'{bench}/trees_{layer}'
os.makedirs(tdir, exist_ok=True)
if layer == 'rescue':
    searched = [l.split('\t')[0] for l in list(open(f'{bench}/scripts/strain_genomes.tsv'))[1:]]
    failed = [l for l in searched if not os.path.exists(f'{bench}/miniprot/{l}.done')]
    open(f'{tdir}/rescue_failed_genomes.txt', 'w').write('\n'.join(failed) + ('\n' if failed else ''))
    print(f'rescue failed in {len(failed)} of {len(searched)} genomes: {failed}')


def fasta_records(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name is not None:
        yield name, ''.join(buf)


cand = defaultdict(dict)                                   # family -> protein_id -> row
for r in csv.DictReader(open(f'{bench}/cand/candidates.tsv'), delimiter='\t'):
    cand[r['family']][r['protein_id']] = r

for fam in ('lig', 'reca', 'polx'):
    members, out = [], open(f'{tdir}/{fam}_input.faa', 'w')
    n = defaultdict(int)
    for name, seq in fasta_records(f'{bench}/refs/{fam}.faa'):
        clade, acc, gene, org = name.split('|')
        n['REF'] += 1
        lab = f'REF_{n["REF"]}'
        out.write(f'>{lab}\n{seq}\n')
        members.append((lab, 'ref', '', acc, clade, len(seq), gene + ' ' + org))
    ann_loci = defaultdict(list)                           # (locustag, scaffold, strand) -> [(start, end)]
    dropped_short = 0
    for name, seq in fasta_records(f'{bench}/cand/{fam}_annotated.faa'):
        _, pid, loc = name.split('|')
        if len(seq) < MIN_LEN[fam]:
            dropped_short += 1
            continue
        n['ANN'] += 1
        lab = f'ANN_{n["ANN"]}'
        out.write(f'>{lab}\n{seq}\n')
        r = cand[fam][pid]
        members.append((lab, 'annotated', loc, pid, '', len(seq), r['domains']))
        if r['chrom'] and r['gstart']:
            ann_loci[(loc, r['chrom'], F.norm_strand(r['strand']))].append((int(r['gstart']), int(r['gend'])))
    # rescued models
    dup = rescued = 0
    scaf_seen = set()
    for tsv in (sorted(glob.glob(f'{bench}/miniprot/*.{fam}.tsv')) if layer == 'rescue' else []):
        for r in csv.DictReader(open(tsv), delimiter='\t'):
            scaf_seen.add((r['locustag'], r['scaffold']))
    names_ok = sum(1 for (loc, ch, st) in ann_loci if (loc, ch) in scaf_seen)
    for fa in sorted(glob.glob(f'{bench}/miniprot/*.{fam}.faa')) if layer == 'rescue' else []:
        for name, seq in fasta_records(fa):
            _, loc, pos, clade = name.split('|')
            scaf, rng, strand = pos.rsplit(':', 2)[0], pos.rsplit(':', 2)[1], pos.rsplit(':', 2)[2]
            s, e = (int(x) for x in rng.split('-'))
            if any(F.overlap_frac((s, e), a) >= 0.5 for a in ann_loci.get((loc, scaf, strand), [])):
                dup += 1
                continue
            rescued += 1
            lab = f'MP_{rescued}'
            out.write(f'>{lab}\n{seq}\n')
            members.append((lab, 'rescued', loc, pos, clade, len(seq), 'miniprot'))
    out.close()
    with open(f'{tdir}/{fam}_members.tsv', 'w') as f:
        f.write('label\tkind\tlocustag\tid\tclade_or_best_ref\tlength\tinfo\n')
        for m in members:
            f.write('\t'.join(str(x) for x in m) + '\n')
    print(f'{fam}: {n["REF"]} refs, {n["ANN"]} annotated (+{dropped_short} dropped as short), {rescued} rescued models '
          f'({dup} miniprot models duplicated an annotated candidate); annotated loci with matching scaffold names: {names_ok}/{len(ann_loci)}')
