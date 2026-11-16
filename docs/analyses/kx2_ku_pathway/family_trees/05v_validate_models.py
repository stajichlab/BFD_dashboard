#!/usr/bin/env python3
"""Step 5v: validate the miniprot rescue against genes that are already annotated.

Usage: 05v_validate_models.py <BENCH>
For every annotated candidate protein (step 2) in a genome that was searched, find the miniprot locus that covers it (same scaffold and
strand, reciprocal overlap >= 90% of the shorter locus). Reports, per family:
  - how many annotated loci were recovered by a model at all
  - of those, how many translated cleanly (status 'kept': no frameshift, no internal stop, length ok)
  - for the clean ones, the identity of the model protein to the annotated protein (difflib ratio, exact-match count)
Pass criterion, fixed before the results were seen: >= 80% of recovered annotated loci translate cleanly AND the median identity of the
clean ones to the annotated protein is >= 0.95. If the criterion fails the rescued models are not used as a layer. Writes
$BENCH/validation_models.tsv (one row per annotated locus) and $BENCH/validation_summary.txt.
"""
import csv, glob, os, sys, difflib, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fam_lib as F
from collections import defaultdict

bench = sys.argv[1]


def fasta(path):
    d, name, buf = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if name:
                d[name] = ''.join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name:
        d[name] = ''.join(buf)
    return d


def model_seqs(fam):
    out = {}
    for fa in glob.glob(f'{bench}/miniprot/*.{fam}.faa'):
        for h, s in fasta(fa).items():
            _, loc, pos, _c = h.split('|')
            out[(loc, pos)] = s
    return out


searched = {l.split('\t')[0] for l in list(open(f'{bench}/scripts/strain_genomes.tsv'))[1:]}
done = {l for l in searched if os.path.exists(f'{bench}/miniprot/{l}.done')}
cand = defaultdict(list)
for r in csv.DictReader(open(f'{bench}/cand/candidates.tsv'), delimiter='\t'):
    if r['locustag'] in done and r['chrom'] and r['gstart']:
        cand[r['family']].append(r)
lines, summary = [], []
for fam in ('lig', 'reca', 'polx'):
    ann = fasta(f'{bench}/cand/{fam}_annotated.faa')
    ann = {k.split('|')[1]: v for k, v in ann.items()}
    seqs = model_seqs(fam)
    models = defaultdict(list)                                       # (loc, scaffold, strand) -> model rows
    for tsv in glob.glob(f'{bench}/miniprot/*.{fam}.tsv'):
        for r in csv.DictReader(open(tsv), delimiter='\t'):
            models[(r['locustag'], r['scaffold'], r['strand'])].append(r)
    n_ann = len(cand[fam])
    rec = clean = 0
    ids = []
    exact = 0
    for r in cand[fam]:
        key = (r['locustag'], r['chrom'], F.norm_strand(r['strand']))
        a = (int(r['gstart']), int(r['gend']))
        best = None
        for m in models.get(key, []):
            b = (int(m['start']), int(m['end']))
            if F.overlap_frac(a, b) >= 0.9 and F.overlap_frac(b, a) >= 0.9 and (best is None or float(m['score']) > float(best['score'])):
                best = m
        status, ident = 'not_recovered', ''
        if best:
            rec += 1
            status = best['status']
            if status == 'kept':
                clean += 1
                pos = f"{best['scaffold']}:{best['start']}-{best['end']}:{best['strand']}"
                s = seqs.get((r['locustag'], pos))
                p = ann.get(r['protein_id'])
                if s and p:
                    ratio = difflib.SequenceMatcher(None, s, p, autojunk=False).ratio()
                    ids.append(ratio)
                    ident = f'{ratio:.3f}'
                    exact += (s == p)
        lines.append(f"{fam}\t{r['locustag']}\t{r['protein_id']}\t{status}\t{ident}")
    med = statistics.median(ids) if ids else float('nan')
    ok = rec > 0 and clean / rec >= 0.8 and med >= 0.95
    summary.append(f'{fam}: {n_ann} annotated loci in {len(done)} finished genomes; recovered by a model: {rec} ({100*rec/max(n_ann,1):.0f}%); '
                   f'translated cleanly: {clean} ({100*clean/max(rec,1):.0f}% of recovered); identical to annotation: {exact}; '
                   f'median identity of clean models to annotation: {med:.3f}; criterion {"PASSED" if ok else "FAILED"}')
with open(f'{bench}/validation_models.tsv', 'w') as f:
    f.write('family\tlocustag\tprotein_id\tmodel_status\tidentity_to_annotation\n' + '\n'.join(lines) + '\n')
failed = sorted(searched - done)
summary.append(f'genomes with a finished rescue: {len(done)} of {len(searched)}; rescue failed in: {failed}')
open(f'{bench}/validation_summary.txt', 'w').write('\n'.join(summary) + '\n')
print('\n'.join(summary))
