#!/usr/bin/env python3
"""T5 step 8: composition test of the fungal GH65 genes against their host genomes (rule R5(b) in README.md), from the LOCAL BFD.duckdb.

Usage (HPCC, SLURM job):  08_hgt_gc.py <BENCH>
Reads  $BENCH/results/gh65_classification.tsv (step 6; proteins with in_tclade = True), $BENCH/cand/candidates.tsv (transcript ids),
       $BENCH/results/gh65_by_genus.tsv (quality flag), $BENCH/results/t5_tree_results.json (R5(a) verdict), and bfd.gene_CDS / bfd.gene_introns.
For every protein inside the fungal trehalase clade: the length-weighted GC of its coding sequence (gene_CDS.GC_content weighted by CDS segment
length) and its z-score against all transcripts of the same genome (mean and sd over the genome's transcripts); the percentile; the number of introns,
and the share of the genome's transcripts that have no intron. Whole-CDS GC is used because the database holds no sequence (GC3 would be more
sensitive; this is a limit). Primary scope: genomes passing the quality filter; all genomes are shown as a sensitivity row.
Rule: HGT supported only if R5(a) is HGT-consistent AND in more than 25% of genomes (five times the 5% expected by chance) at least one clade gene has
|z| > 2; vertical-consistent needs R5(a) eukaryotic AND at most 10% of genomes beyond |z| > 2; anything else is inconclusive.
Writes $BENCH/results/hgt_gc.tsv, hgt_summary.json, hgt_verdict.txt.
"""
import csv
import json
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'lib'))

Z_CUT, HGT_SHARE, VERT_SHARE = 2.0, 0.25, 0.10


def zscore(x, mean, sd):
    return (x - mean) / sd if sd and sd > 0 else float('nan')


def percentile_rank(x, values):
    """Share of values at or below x (values is a sorted list)."""
    import bisect
    return bisect.bisect_right(values, x) / len(values)


def summarize(rows, tree_verdict):
    """rows: per protein dicts with locustag, z, quality. Returns the summary dict and the verdict string."""
    def scope(sel):
        by = {}
        for r in sel:
            if not math.isnan(r['z']):
                by.setdefault(r['locustag'], []).append(r['z'])
        n = len(by)
        beyond = sum(1 for zs in by.values() if any(abs(z) > Z_CUT for z in zs))
        medz = statistics.median([statistics.median(zs) for zs in by.values()]) if by else float('nan')
        medabs = statistics.median([max(abs(z) for z in zs) for zs in by.values()]) if by else float('nan')
        return {'genomes': n, 'genomes_beyond_2sd': beyond, 'share_beyond_2sd': beyond / n if n else float('nan'), 'median_z': medz, 'median_max_abs_z': medabs}
    out = {'quality_genomes': scope([r for r in rows if r['quality']]), 'all_genomes': scope(rows)}
    share = out['quality_genomes']['share_beyond_2sd']
    if tree_verdict.startswith('HGT-consistent') and share > HGT_SHARE:
        v = 'HGT supported: bacterial/archaeal nearest neighbours and clade genes deviate from host GC in more than 25% of genomes'
    elif tree_verdict.startswith('vertical-consistent') and share <= VERT_SHARE:
        v = 'vertical origin supported: eukaryotic nearest neighbours and clade genes match host GC (<= 10% of genomes beyond 2 sd)'
    else:
        v = f'inconclusive: tree says "{tree_verdict[:60]}", GC share beyond 2 sd = {share:.2f}'
    return out, v


def main(bench):
    from bfd_connect import connect
    cls = [r for r in csv.DictReader(open(f'{bench}/results/gh65_classification.tsv'), delimiter='\t') if r['in_tclade'] == 'True']
    cands = {r['protein_id']: r for r in csv.DictReader(open(f'{bench}/cand/candidates.tsv'), delimiter='\t')}
    qual = {r['locustag']: r['quality'] == 'True' for r in csv.DictReader(open(f'{bench}/results/gh65_by_genus.tsv'), delimiter='\t')}
    tree = json.load(open(f'{bench}/results/t5_tree_results.json'))
    loci = sorted({r['locustag'] for r in cls})
    print('genomes with a protein in the fungal trehalase clade:', len(loci), flush=True)
    con = connect(scratch=False, threads=8, memory='48GB', tmp=os.environ.get('SCRATCH'))
    rng = ' OR '.join(f"(transcript_id >= '{l}_' AND transcript_id < '{l}`')" for l in loci)
    con.execute(f"""
    CREATE TEMP TABLE tx AS
    SELECT transcript_id, split_part(transcript_id,'_',1) AS locustag,
           sum(GC_content * ("end" - start + 1)) / sum("end" - start + 1) AS gc, sum("end" - start + 1) AS cds_len
    FROM gene_CDS WHERE {rng} GROUP BY 1, 2""")
    con.execute(f"""
    CREATE TEMP TABLE intr AS
    SELECT transcript_id, count(*) AS n_introns FROM gene_introns WHERE {rng.replace('transcript_id', 'transcript_id')} GROUP BY 1""")
    stats = {r[0]: r[1:] for r in con.execute("""
        SELECT t.locustag, avg(t.gc), stddev_pop(t.gc), count(*), avg((i.n_introns IS NULL)::int)
        FROM tx t LEFT JOIN intr i USING (transcript_id) GROUP BY 1""").fetchall()}
    gcs = {}
    for loc, gc in con.execute('SELECT locustag, gc FROM tx ORDER BY locustag, gc').fetchall():
        gcs.setdefault(loc, []).append(gc)
    txrows = {r[0]: r[1:] for r in con.execute('SELECT t.transcript_id, t.gc, coalesce(i.n_introns, 0) FROM tx t LEFT JOIN intr i USING (transcript_id)').fetchall()}
    rows = []
    for r in cls:
        tid = cands.get(r['protein_id'], {}).get('transcript_id', '')
        loc = r['locustag']
        if tid not in txrows or loc not in stats:
            continue
        gc, nint = txrows[tid]
        mean, sd, ntx, intronless = stats[loc]
        rows.append({'locustag': loc, 'genus': r['genus'], 'class': r['class'], 'protein_id': r['protein_id'], 'quality': qual.get(loc, False),
                     'cds_gc': round(gc, 4), 'genome_mean_gc': round(mean, 4), 'genome_sd_gc': round(sd, 4), 'z': zscore(gc, mean, sd),
                     'percentile': round(percentile_rank(gc, gcs[loc]), 3), 'n_introns': nint, 'genome_share_intronless': round(intronless, 3),
                     'genome_transcripts': ntx})
    with open(f'{bench}/results/hgt_gc.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows({**r, 'z': round(r['z'], 3)} for r in rows)
    summ, verdict = summarize(rows, tree.get('R5a', {}).get('verdict', ''))
    summ['n_proteins'] = len(rows)
    summ['clade_intronless_share'] = sum(1 for r in rows if r['n_introns'] == 0) / len(rows) if rows else None
    summ['genome_intronless_share_median'] = statistics.median([r['genome_share_intronless'] for r in rows]) if rows else None
    summ['verdict'] = verdict
    json.dump(summ, open(f'{bench}/results/hgt_summary.json', 'w'), indent=1, default=str)
    open(f'{bench}/results/hgt_verdict.txt', 'w').write(verdict + '\n')
    print(json.dumps(summ, indent=1, default=str))


if __name__ == '__main__':
    main(sys.argv[1])
