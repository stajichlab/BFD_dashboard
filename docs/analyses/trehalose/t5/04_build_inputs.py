#!/usr/bin/env python3
"""T5 step 4: combine references, context taxa and the BFD candidates into one FASTA for the tree.

Usage: 04_build_inputs.py <BENCH> [min_len]
Candidates shorter than min_len aa (default 400; the GH65 catalytic domain alone is about 400 aa, so shorter proteins are fragments that
cannot align over most of the family) are excluded and listed in $BENCH/trees/excluded_short.tsv. Writes $BENCH/trees/gh65_input.faa and
$BENCH/trees/gh65_members.tsv (label, kind ref/context/annotated, locustag, id, group, length, info). Sequence labels are short (REF_n,
CTX_n, ANN_n) because IQ-TREE and the Newick parser dislike '|' and long names; gh65_members.tsv maps them back.
"""
import csv
import os
import sys


def read_fasta(path):
    name, buf, out = None, [], []
    for line in open(path):
        if line.startswith('>'):
            if name is not None:
                out.append((name, ''.join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name is not None:
        out.append((name, ''.join(buf)))
    return out


def build(bench, min_len=400):
    os.makedirs(f'{bench}/trees', exist_ok=True)
    members, seqs, excl = [], [], []
    nref = nctx = nann = 0
    for name, s in read_fasta(f'{bench}/refs/refs.faa'):
        grp, acc, org = name.split('|', 2)
        if grp.startswith('CONTEXT_'):
            nctx += 1
            label, kind = f'CTX_{nctx}', 'context'
        else:
            nref += 1
            label, kind = f'REF_{nref}', 'ref'
        members.append((label, kind, '', acc, grp, len(s), org))
        seqs.append((label, s))
    for name, s in read_fasta(f'{bench}/cand/gh65_annotated.faa'):
        _, pid, loc = name.split('|')
        if len(s) < min_len:
            excl.append((pid, loc, len(s)))
            continue
        nann += 1
        label = f'ANN_{nann}'
        members.append((label, 'annotated', loc, pid, '', len(s), ''))
        seqs.append((label, s))
    with open(f'{bench}/trees/gh65_input.faa', 'w') as f:
        f.writelines(f'>{label}\n{s}\n' for label, s in seqs)
    with open(f'{bench}/trees/gh65_members.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['label', 'kind', 'locustag', 'id', 'clade_or_group', 'length', 'info'])
        w.writerows(members)
    with open(f'{bench}/trees/excluded_short.tsv', 'w') as f:
        f.write('protein_id\tlocustag\tlength\n')
        f.writelines(f'{a}\t{b}\t{c}\n' for a, b, c in excl)
    return nref, nctx, nann, len(excl)


if __name__ == '__main__':
    ml = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    nref, nctx, nann, nex = build(sys.argv[1], ml)
    print(f'references {nref}, context {nctx}, annotated {nann} (excluded shorter than {ml}: {nex}); total {nref + nctx + nann}')
