#!/usr/bin/env python3
"""Call introns from a miniprot GFF3 (no gene finder involved).

Library functions (tested in test_miniprot_introns.py) and a command line:

  miniprot_introns.py --gff aln.gff[.gz] --ref ref.faa --out table.tsv [--genome asm.fa[.gz]] [--min-intron 20]

For each reference protein the best-scoring alignment with Rank=1 is kept. Its CDS blocks are sorted by coordinate; a gap
between consecutive blocks of >= --min-intron nt is an intron, a shorter gap is counted as a frameshift/short gap and is not
an intron. If --genome is given the donor and acceptor dinucleotides of each intron are read from the assembly (reverse
complemented on the minus strand) and reported as GT-AG, GC-AG or other.
"""
import argparse
import gzip
import sys

COMP = str.maketrans('ACGTN', 'TGCAN')


def _open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def read_fasta(path):
    """Return {id: sequence} using the first token of the header as id; sequences are upper case."""
    seqs, cur, parts = {}, None, []
    with _open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if cur is not None:
                    seqs[cur] = ''.join(parts).upper()
                cur, parts = line[1:].split()[0], []
            elif cur is not None:
                parts.append(line)
    if cur is not None:
        seqs[cur] = ''.join(parts).upper()
    return seqs


def parse_attrs(s):
    d = {}
    for kv in s.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
    return d


def parse_gff(lines, qlens):
    """Yield one dict per miniprot mRNA with its CDS blocks.

    qlens: {reference protein id: length}. Alignments whose Target is not in qlens get qlen None.
    """
    cur = None
    out = {}
    for line in lines:
        if not line or line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        seqid, _src, typ, start, end, score, strand, _ph, attrs = f[:9]
        a = parse_attrs(attrs)
        if typ == 'mRNA':
            tgt = a.get('Target', '').split()
            q = tgt[0] if tgt else None
            try:
                sc = float(score)
            except ValueError:
                sc = 0.0
            rec = {'id': a.get('ID'), 'query': q, 'scaffold': seqid, 'strand': strand, 'start': int(start), 'end': int(end),
                   # miniprot 0.2 writes no Rank attribute (newer versions do): a missing Rank counts as primary
                   'score': sc, 'rank': int(a.get('Rank', '1') or 1), 'identity': float(a.get('Identity', 'nan')),
                   'positive': float(a.get('Positive', 'nan')), 'frameshift': a.get('Frameshift') == 'true',
                   'qstart': int(tgt[1]) if len(tgt) > 2 else None, 'qend': int(tgt[2]) if len(tgt) > 2 else None,
                   'cds': []}
            rec['qlen'] = qlens.get(q)
            out[rec['id']] = rec
        elif typ == 'CDS':
            parent = out.get(a.get('Parent'))
            if parent is not None:
                parent['cds'].append((int(start), int(end)))
    return list(out.values())


def best_per_query(alns):
    """Best-scoring alignment per reference protein (ties: higher identity). An alignment with an explicit Rank other than 1 is
    ignored; miniprot 0.2 writes no Rank, so there the best score alone decides."""
    best = {}
    for a in alns:
        if a['rank'] != 1 or a['query'] is None:
            continue
        key = (a['score'], a['identity'])
        if a['query'] not in best or key > (best[a['query']]['score'], best[a['query']]['identity']):
            best[a['query']] = a
    return best


def introns_of(aln, genome=None, min_intron=20):
    """Return (introns, short_gaps) for one alignment; introns is a list of dicts with length and motif."""
    cds = sorted(aln['cds'])
    introns, short = [], 0
    for (s1, e1), (s2, e2) in zip(cds, cds[1:]):
        gap_start, gap_end = e1 + 1, s2 - 1          # 1-based inclusive intron coordinates on the forward strand
        glen = gap_end - gap_start + 1
        if glen < min_intron:
            if glen >= 0:
                short += 1
            continue
        motif = None
        if genome is not None and aln['scaffold'] in genome:
            seq = genome[aln['scaffold']]
            first2 = seq[gap_start - 1:gap_start + 1]
            last2 = seq[gap_end - 2:gap_end]
            if aln['strand'] == '+':
                motif = first2 + '-' + last2
            else:
                motif = last2.translate(COMP)[::-1] + '-' + first2.translate(COMP)[::-1]
        introns.append({'len': glen, 'motif': motif})
    return introns, short


def alignment_row(aln, genome=None, min_intron=20):
    ints, short = introns_of(aln, genome, min_intron)
    cds_nt = sum(e - s + 1 for s, e in aln['cds'])
    qcov = None
    if aln['qlen'] and aln['qstart'] is not None:
        qcov = (aln['qend'] - aln['qstart'] + 1) / aln['qlen']
    return {'query': aln['query'], 'scaffold': aln['scaffold'], 'strand': aln['strand'], 'score': aln['score'],
            'identity': aln['identity'], 'positive': aln['positive'], 'qlen': aln['qlen'], 'qcov': qcov,
            'n_cds': len(aln['cds']), 'n_introns': len(ints), 'n_short_gaps': short, 'cds_nt': cds_nt,
            'intron_lens': ','.join(str(i['len']) for i in ints),
            'motifs': ','.join(i['motif'] or 'NA' for i in ints)}


COLS = ['query', 'scaffold', 'strand', 'score', 'identity', 'positive', 'qlen', 'qcov', 'n_cds', 'n_introns', 'n_short_gaps',
        'cds_nt', 'intron_lens', 'motifs']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gff', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--genome')
    ap.add_argument('--min-intron', type=int, default=20)
    a = ap.parse_args()
    qlens = {k: len(v) for k, v in read_fasta(a.ref).items()}
    genome = read_fasta(a.genome) if a.genome else None
    with _open(a.gff) as f:
        alns = parse_gff(f, qlens)
    best = best_per_query(alns)
    with open(a.out, 'w') as o:
        o.write('\t'.join(COLS) + '\n')
        for q in sorted(best):
            r = alignment_row(best[q], genome, a.min_intron)
            o.write('\t'.join('' if r[c] is None else (f'{r[c]:.4f}' if isinstance(r[c], float) else str(r[c])) for c in COLS) + '\n')
    print(f'{len(best)} proteins with a Rank=1 alignment from {len(alns)} alignments', file=sys.stderr)


if __name__ == '__main__':
    main()
