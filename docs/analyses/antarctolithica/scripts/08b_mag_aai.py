"""Strain relatedness of the annotated MAGs: DIAMOND blastp of each MAG proteome against the Mars-82.9 proteome.
Best hit per query, kept if the alignment covers >= 90 % of both query and subject; AAI = mean/median % identity of
those pairs; frac_identical = share of pairs at 100 % identity. Writes results/mag_aai.csv (read by 16_report_tables.py).
Needs `diamond` on PATH (HPCC: module load diamond). Temporary DIAMOND files go to $SCRATCH (or TMPDIR)."""
import glob, os, re, subprocess, tempfile, pandas as pd

REF = 'Mars-82.9'
COLS = ['q', 's', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore']


def strain(path):
    return re.search(r'obscura_(Mars-[\d.]+?)\.proteins', os.path.basename(path)).group(1)


def summarise(hits, mag, cov=0.9):
    """hits: DIAMOND tabular output (COLS). Best hit per query by bitscore, then the two-way coverage filter."""
    h = hits.sort_values('bitscore', ascending=False).drop_duplicates('q')
    h = h[(h.length >= cov * h.qlen) & (h.length >= cov * h.slen)]
    return dict(mag=mag, vs=REF, n_pairs=len(h), mean_AAI=h.pident.mean(), median_AAI=h.pident.median(),
                frac_identical=(h.pident >= 100).mean(), frac_ge99=(h.pident >= 99).mean(), frac_lt90=(h.pident < 90).mean())


if __name__ == '__main__':
    peps = {strain(p): p for p in sorted(glob.glob('mag/*/predict_results/*.proteins.fa.gz'))}
    ncpu = os.environ.get('NCPU', '4')
    tmp = tempfile.mkdtemp(dir=os.environ.get('SCRATCH') or os.environ.get('TMPDIR') or None)
    subprocess.run(['diamond', 'makedb', '--in', peps[REF], '-d', f'{tmp}/ref', '--quiet'], check=True)
    rows = []
    for m, p in peps.items():
        if m == REF: continue
        out = f'{tmp}/{m}.tsv'
        subprocess.run(['diamond', 'blastp', '-q', p, '-d', f'{tmp}/ref', '-o', out, '--threads', ncpu, '--quiet', '--more-sensitive',
                        '--max-target-seqs', '1', '--evalue', '1e-10', '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore'], check=True)
        rows.append(summarise(pd.read_csv(out, sep='\t', header=None, names=COLS), m))
    R = pd.DataFrame(rows); R.to_csv('results/mag_aai.csv', index=False); print(R.round(3).to_string(index=False))
