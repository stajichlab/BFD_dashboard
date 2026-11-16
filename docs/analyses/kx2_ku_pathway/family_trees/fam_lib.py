"""Shared helpers for the gene-family tree classification (translation of miniprot models, Newick handling, placement).

Everything here is plain Python 3 without third-party packages, and is covered by tests/test_fam_lib.py.
"""
import gzip
import re
from collections import defaultdict

# ---------------------------------------------------------------- sequences
_BASES = 'TCAG'
_AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate((x, y, z) for x in _BASES for y in _BASES for z in _BASES)}
_COMP = str.maketrans('ACGTUacgtuNnRYKMBDHVrykmbdhv', 'TGCAAtgcaaNnYRMKVHDByrmkvhdb')


def revcomp(s):
    return s.translate(_COMP)[::-1]


def translate(dna, table=None):
    """Standard genetic code (NCBI table 1). Codons with ambiguity codes give 'X'. Trailing bases that do not make a codon are dropped."""
    t = table or CODON
    s = dna.upper().replace('U', 'T')
    return ''.join(t.get(s[i:i + 3], 'X') for i in range(0, len(s) - len(s) % 3, 3))


def read_fasta(path):
    op = gzip.open if str(path).endswith('.gz') else open
    seqs, name, buf = {}, None, []
    with op(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(buf)
    return seqs


# ---------------------------------------------------------------- miniprot
def parse_miniprot_gff(lines):
    """Parse `miniprot --gff` output into a list of models: dict(id, scaf, strand, start, end, score, identity, positive, target, cds=[(start,end,phase)]).
    Coordinates are 1-based inclusive as in GFF3."""
    models, by_id = [], {}
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 9:
            continue
        scaf, _src, typ, s, e, score, strand, phase, attr = c[:9]
        a = dict(kv.split('=', 1) for kv in attr.split(';') if '=' in kv)
        if typ == 'mRNA':
            m = {'id': a['ID'], 'scaf': scaf, 'strand': strand, 'start': int(s), 'end': int(e),
                 'score': float(score) if score != '.' else 0.0, 'identity': float(a.get('Identity', 'nan')),
                 'positive': float(a.get('Positive', 'nan')), 'target': a.get('Target', '').split(' ')[0], 'cds': []}
            models.append(m)
            by_id[m['id']] = m
        elif typ == 'CDS' and a.get('Parent') in by_id:
            by_id[a['Parent']]['cds'].append((int(s), int(e), int(phase) if phase.isdigit() else 0))
    return models


def model_protein(model, genome):
    """Translate a miniprot model from the genome, following the GFF3 phase of every CDS block.
    Returns (protein, n_internal_stops, ok, n_frame_breaks).

    Blocks are taken in transcript order (descending coordinates on the reverse strand, each reverse complemented). The phase of a block is
    the number of its leading bases that complete the codon begun at the end of the previous block. If the pending bases from the previous
    block plus that lead do not make exactly one codon, the alignment changed frame there (miniprot reports this as a frameshift): the
    incomplete codon is dropped and the break is counted. The phase of the first block (normally 0) is the number of leading bases to skip.
    A final stop codon is removed and not counted as internal."""
    seq = genome.get(model['scaf'])
    if seq is None or not model['cds']:
        return '', 0, False, 0
    pieces = sorted(model['cds'])
    order = list(reversed(pieces)) if model['strand'] == '-' else pieces
    pending, prot, breaks = '', [], 0
    for i, (s, e, ph) in enumerate(order):
        b = seq[s - 1:e]
        if model['strand'] == '-':
            b = revcomp(b)
        lead, body = b[:ph], b[ph:]
        if i > 0:
            if len(pending) + len(lead) == 3 and (pending or lead):
                prot.append(translate(pending + lead))
            elif pending or lead:
                breaks += 1                       # frame changed between blocks
        k = len(body) - len(body) % 3
        prot.append(translate(body[:k]))
        pending = body[k:]
    p = ''.join(prot)
    if p.endswith('*'):
        p = p[:-1]
    return p, p.count('*'), True, breaks


def norm_strand(s):
    """gene_info stores the strand as '1' / '-1'; miniprot and GFF3 use '+' / '-'. Returns '+' or '-'."""
    return '-' if str(s).strip() in ('-', '-1', 'minus', 'reverse') else '+'


def max_intron(model):
    """Longest gap between consecutive CDS pieces of a model (0 for a single-exon model)."""
    p = sorted(model['cds'])
    return max((p[i + 1][0] - p[i][1] - 1 for i in range(len(p) - 1)), default=0)


def overlap_frac(a, b):
    """Overlap of two (start,end) intervals as a fraction of the shorter one."""
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    if hi < lo:
        return 0.0
    return (hi - lo + 1) / min(a[1] - a[0] + 1, b[1] - b[0] + 1)


def greedy_loci(models, min_overlap=0.5):
    """Keep the best-scoring model per locus: models on the same scaffold and strand that overlap by >= min_overlap of the shorter
    one are one locus. Returns the kept models, best score first."""
    kept = []
    for m in sorted(models, key=lambda x: -x['score']):
        if not any(k['scaf'] == m['scaf'] and k['strand'] == m['strand'] and
                   overlap_frac((k['start'], k['end']), (m['start'], m['end'])) >= min_overlap for k in kept):
            kept.append(m)
    return kept


# ---------------------------------------------------------------- newick
class Node:
    __slots__ = ('name', 'bl', 'kids', 'parent', 'support')

    def __init__(self):
        self.name, self.bl, self.kids, self.parent, self.support = '', 0.0, [], None, None


def parse_support(label):
    """IQ-TREE writes 'SH-aLRT/UFBoot' (e.g. '87.3/95'); a single number is taken as UFBoot. Returns (sh, ufb) or (None, None)."""
    if not label:
        return None, None
    m = re.fullmatch(r'([\d.]+)/([\d.]+)', label)
    if m:
        return float(m.group(1)), float(m.group(2))
    m = re.fullmatch(r'[\d.]+', label)
    if m:
        v = float(label)
        return None, (v * 100 if v <= 1.0 else v)
    return None, None


def parse_newick(text):
    text = text.strip().rstrip(';')
    pos = 0

    def node(parent):
        nonlocal pos
        n = Node()
        n.parent = parent
        if text[pos] == '(':
            pos += 1
            while True:
                n.kids.append(node(n))
                if text[pos] == ',':
                    pos += 1
                    continue
                if text[pos] == ')':
                    pos += 1
                    break
                raise ValueError('bad newick at %d' % pos)
            m = re.match(r"([^:,;()]*)", text[pos:])
            lab = m.group(1)
            pos += len(lab)
            n.support = parse_support(lab)
        else:
            m = re.match(r"([^:,;()]+)", text[pos:])
            n.name = m.group(1)
            pos += len(n.name)
        if pos < len(text) and text[pos] == ':':
            m = re.match(r':([^,;()]*)', text[pos:])
            n.bl = float(m.group(1)) if m.group(1) else 0.0
            pos += len(m.group(0))
        return n

    return node(None)


def leaves(n):
    if not n.kids:
        return [n]
    out = []
    for k in n.kids:
        out += leaves(k)
    return out


def depth_to_root(n):
    d = 0.0
    while n.parent is not None:
        d += n.bl
        n = n.parent
    return d


def path_distance(a, b):
    anc = {}
    d = 0.0
    x = a
    while x is not None:
        anc[id(x)] = d
        d += x.bl
        x = x.parent
    d = 0.0
    x = b
    while x is not None:
        if id(x) in anc:
            return anc[id(x)] + d
        d += x.bl
        x = x.parent
    raise ValueError('not in the same tree')


# ---------------------------------------------------------------- classification
def classify(root, ref_clade, candidates, min_ufb=90.0, min_sh=80.0):
    """Place candidate leaves by the references they nest with.

    ref_clade: {leaf name -> clade label} for the reference leaves. candidates: iterable of leaf names to classify (they are never
    treated as references, even if present in ref_clade, which is what the leave-one-out test uses).
    For each candidate, walk up the tree. The first ancestor that contains references decides: if all of its references carry one
    label L the candidate is assigned to L provided the node is supported (UFBoot >= min_ufb and SH-aLRT >= min_sh); if not
    supported, keep climbing while the references stay single-label; mixed labels before a supported node -> 'ambiguous'.
    Returns {candidate: dict(call, clade, n_refs, ufboot, sh, nearest_ref, nearest_dist, terminal_bl, labels)}.
    """
    cand = set(candidates)
    byname = {n.name: n for n in leaves(root)}
    refs = {k: v for k, v in ref_clade.items() if k in byname and k not in cand}
    out = {}
    for c in cand:
        node = byname[c]
        res = {'call': 'ambiguous', 'clade': '', 'n_refs': 0, 'ufboot': None, 'sh': None, 'labels': ''}
        x = node.parent
        seen_single = False
        while x is not None:
            labs = sorted({refs[l.name] for l in leaves(x) if l.name in refs})
            nref = sum(1 for l in leaves(x) if l.name in refs)
            if labs:
                if len(labs) > 1:
                    res['labels'] = ','.join(labs)
                    res['call'] = 'ambiguous' if not seen_single else res['call']
                    if not seen_single:
                        res['n_refs'] = nref
                    break
                sh, ufb = x.support if x.support else (None, None)
                supported = (ufb is None or ufb >= min_ufb) and (sh is None or sh >= min_sh)
                seen_single = True
                res.update(clade=labs[0], labels=labs[0], n_refs=nref, ufboot=ufb, sh=sh)
                if supported and x.parent is not None:
                    res['call'] = 'assigned'
                    break
                res['call'] = 'low_support'
            x = x.parent
        dists = sorted((path_distance(node, byname[r]), r) for r in refs)
        res['nearest_ref'] = dists[0][1] if dists else ''
        res['nearest_dist'] = dists[0][0] if dists else None
        res['terminal_bl'] = node.bl
        out[c] = res
    return out


def loo_check(root, ref_clade, **kw):
    """Leave-one-out: classify each reference as if it were a candidate. Returns (n_tested, n_correct, details)."""
    details, tested, correct = [], 0, 0
    sizes = defaultdict(int)
    for v in ref_clade.values():
        sizes[v] += 1
    names = {n.name for n in leaves(root)}
    for r, lab in ref_clade.items():
        if r not in names or sizes[lab] < 2:
            continue                       # a clade with a single reference cannot be recovered without itself
        res = classify(root, ref_clade, [r], **kw)[r]
        tested += 1
        ok = res['call'] == 'assigned' and res['clade'] == lab
        correct += ok
        details.append((r, lab, res['call'], res['clade'], ok))
    return tested, correct, details
