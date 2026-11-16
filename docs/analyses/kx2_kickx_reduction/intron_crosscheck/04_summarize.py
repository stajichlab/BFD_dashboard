#!/usr/bin/env python3
"""Summarize the miniprot intron calls: per-genome metrics, group tables, effect sizes, identity strata, splice sites, annotation
comparison and a bootstrap over reference proteins. Thresholds are those fixed in README.md before the runs.

Usage (HPCC, /usr/bin/python3.12 with numpy): 04_summarize.py <BENCH>
Reads  <BENCH>/genomes_resolved.tsv, <BENCH>/tables/{A,B,C}_<locustag>.tsv, <BENCH>/annotation_introns.tsv
Writes <BENCH>/results/*.tsv and prints a text summary.
"""
import csv
import os
import statistics
import sys
from collections import defaultdict

import numpy as np

BENCH = sys.argv[1]
OUT = os.path.join(BENCH, 'results')
os.makedirs(OUT, exist_ok=True)
SETTINGS = ['A', 'B', 'C']
THR = [0.40, 0.50, 0.60]
QCOV = 0.80
MIN_N = 50
MIN_BIN = 15
BINS = [(0.4, 0.5), (0.5, 0.6), (0.6, 0.7), (0.7, 1.01)]
CORE = {'Coemansia', 'Kickxella', 'Linderina', 'Dipsacomyces', 'Martensiomyces'}
MYCO = {'Mycoemilia', 'Spiromyces'}


def group_of(g):
    if g['set'] == 'anchor_yeast':
        return 'Saccharomycetes anchors'
    if g['set'] == 'anchor_pezizomycotina':
        return 'Pezizomycotina anchors'
    if g['set'] == 'coemansia_strain':
        return 'Coemansia strains'
    if g['genus'] in CORE:
        return 'Kickxellales core'
    if g['genus'] in MYCO:
        return 'Mycoemilia+Spiromyces'
    if g['ord'] == 'Harpellales':
        return 'Harpellales'
    if g['ord'] in ('Dimargaritales', 'Ramicandelaberales'):
        return 'Dimargaritales+Ramicandelaber'
    if g['ord'] == 'Zoopagales':
        return 'Zoopagales'
    if g['ord'] in ('Entomophthorales', 'Basidiobolales'):
        return 'Entomophthorales+Basidiobolus'
    if g['ord'] == 'Mucorales':
        return 'Mucorales'
    if g['phylum'] in ('Mucoromycota', 'Mortierellomycota', 'Glomeromycota'):
        return 'Other Mucoromycota group'
    return 'Root (chytrids etc.)'


GROUP_ORDER = ['Kickxellales core', 'Mycoemilia+Spiromyces', 'Harpellales', 'Dimargaritales+Ramicandelaber', 'Zoopagales',
               'Entomophthorales+Basidiobolus', 'Mucorales', 'Other Mucoromycota group', 'Root (chytrids etc.)',
               'Saccharomycetes anchors', 'Pezizomycotina anchors', 'Coemansia strains']


def load(setting, loc):
    path = f'{BENCH}/tables/{setting}_{loc}.tsv'
    rows = []
    if not os.path.exists(path):          # e.g. setting C failed for this genome (miniprot abort, see failures.tsv)
        return rows
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            ms = [m for m in r['motifs'].split(',') if m and m != 'NA']
            rows.append({'query': r['query'], 'identity': float(r['identity']), 'qcov': float(r['qcov']) if r['qcov'] else 0.0,
                         'n_introns': int(r['n_introns']), 'cds_nt': int(r['cds_nt']), 'motifs': r['motifs'],
                         # POST HOC (added after seeing that the intron count depends strongly on the aligner penalty):
                         # introns with a canonical GT-AG or GC-AG splice site only
                         'canon': sum(1 for m in ms if m in ('GT-AG', 'GC-AG'))})
    return rows


def qualifying(rows, lo, hi=1.01):
    return [r for r in rows if r['qcov'] >= QCOV and lo <= r['identity'] < hi]


def metrics(q):
    n = len(q)
    if n == 0:
        return None
    ints = [r['n_introns'] for r in q]
    motifs = [m for r in q for m in r['motifs'].split(',') if m and m != 'NA']
    nm = len(motifs)
    return {'n': n, 'median_introns': statistics.median(ints), 'mean_introns': sum(ints) / n,
            'frac_ge1': sum(1 for x in ints if x >= 1) / n,
            'per_kb': 1000 * sum(ints) / max(1, sum(r['cds_nt'] for r in q)),
            'canon_per_kb': 1000 * sum(r['canon'] for r in q) / max(1, sum(r['cds_nt'] for r in q)),
            'frac_canon_ge1': sum(1 for r in q if r['canon'] >= 1) / n,
            'gt_ag': sum(1 for m in motifs if m == 'GT-AG') / nm if nm else float('nan'),
            'gc_ag': sum(1 for m in motifs if m == 'GC-AG') / nm if nm else float('nan'),
            'other_sites': sum(1 for m in motifs if m not in ('GT-AG', 'GC-AG')) / nm if nm else float('nan'),
            'median_identity': statistics.median(r['identity'] for r in q)}


def cliff(x, y):
    x = [v for v in x if v == v]
    y = [v for v in y if v == v]
    if not x or not y:
        return float('nan')
    gt = sum(1 for a in x for b in y if a > b)
    lt = sum(1 for a in x for b in y if a < b)
    return (gt - lt) / (len(x) * len(y))


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(x, y):
    pairs = [(a, b) for a, b in zip(x, y) if a == a and b == b]
    if len(pairs) < 5:
        return float('nan'), len(pairs)
    rx, ry = ranks([p[0] for p in pairs]), ranks([p[1] for p in pairs])
    return float(np.corrcoef(rx, ry)[0, 1]), len(pairs)


genomes = list(csv.DictReader(open(f'{BENCH}/genomes_resolved.tsv'), delimiter='\t'))
for g in genomes:
    g['group'] = group_of(g)
data = {s: {g['locustag']: load(s, g['locustag']) for g in genomes} for s in SETTINGS}
annot = {r['locustag']: float(r['introns_per_transcript']) for r in csv.DictReader(open(f'{BENCH}/annotation_introns.tsv'), delimiter='\t')}
backbone = [g for g in genomes if g['set'] == 'backbone']

# ---- per-genome metrics
per_genome = []
M = {}   # (setting, thr, loc) -> metrics or None
for s in SETTINGS:
    for thr in THR:
        for g in genomes:
            m = metrics(qualifying(data[s][g['locustag']], thr))
            M[(s, thr, g['locustag'])] = m
            per_genome.append((s, thr, g, m))
with open(f'{OUT}/per_genome_metrics.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['setting', 'ident_thr', 'locustag', 'genus', 'group', 'set', 'n_aligned', 'n_qualifying', 'median_introns',
                'mean_introns', 'frac_ge1', 'per_kb', 'gt_ag', 'gc_ag', 'other_sites', 'median_identity', 'annotation_introns_per_tx'])
    for s, thr, g, m in per_genome:
        row = [s, thr, g['locustag'], g['genus'], g['group'], g['set'], len(data[s][g['locustag']])]
        if m:
            row += [m['n'], m['median_introns'], f"{m['mean_introns']:.3f}", f"{m['frac_ge1']:.3f}", f"{m['per_kb']:.3f}",
                    f"{m['gt_ag']:.3f}", f"{m['gc_ag']:.3f}", f"{m['other_sites']:.3f}", f"{m['median_identity']:.3f}"]
        else:
            row += [0] + [''] * 8
        row.append(annot.get(g['locustag'], ''))
        w.writerow(row)

METRICS = ['median_introns', 'frac_ge1', 'per_kb', 'mean_introns', 'canon_per_kb', 'frac_canon_ge1']   # last two are POST HOC


def group_values(s, thr, metric, group, genomes_subset=backbone):
    vals, excluded = [], []
    for g in genomes_subset:
        if g['group'] != group:
            continue
        m = M[(s, thr, g['locustag'])]
        if m is None or m['n'] < MIN_N:
            excluded.append(g['genus'])
        else:
            vals.append(m[metric])
    return vals, excluded


def pooled(s, thr, metric, groups, subset=backbone):
    v = []
    for gp in groups:
        v += group_values(s, thr, metric, gp, subset)[0]
    return v


lines = []
P = lines.append
P('== Group medians of genus values (backbone genera; qualifying = qcov>=0.8, identity>=thr; genomes with <50 qualifying proteins excluded) ==')
group_rows = []
for s in SETTINGS:
    for thr in THR:
        P(f'\n-- setting {s}, identity >= {thr}')
        P(f"{'group':34s} {'n_gen':>5s} {'med_introns':>11s} {'frac>=1':>8s} {'per_kb':>7s} {'mean':>6s} {'canon/kb*':>9s} {'canon>=1*':>9s}   (annotation introns/transcript)   *post hoc")
        for gp in GROUP_ORDER:
            if gp == 'Coemansia strains':
                continue
            subset = [g for g in genomes if g['set'] != 'coemansia_strain']
            vals = {mt: group_values(s, thr, mt, gp, subset)[0] for mt in METRICS}
            ex = group_values(s, thr, 'per_kb', gp, subset)[1]
            n = len(vals['per_kb'])
            if n == 0:
                P(f'{gp:34s} {0:5d}   (no genome with >= {MIN_N} qualifying proteins; excluded: {ex})')
                continue
            an = [annot[g['locustag']] for g in subset if g['group'] == gp and g['locustag'] in annot]
            P(f"{gp:34s} {n:5d} {statistics.median(vals['median_introns']):11.2f} {statistics.median(vals['frac_ge1']):8.2f} "
              f"{statistics.median(vals['per_kb']):7.2f} {statistics.median(vals['mean_introns']):6.2f} "
              f"{statistics.median(vals['canon_per_kb']):9.2f} {statistics.median(vals['frac_canon_ge1']):9.2f}   ({statistics.median(an):.2f})"
              + (f'  excluded: {ex}' if ex else ''))
            group_rows.append([s, thr, gp, n] + [statistics.median(vals[mt]) for mt in METRICS] + [statistics.median(an)])
with open(f'{OUT}/group_table.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['setting', 'ident_thr', 'group', 'n_genera'] + [f'median_{m}' for m in METRICS] + ['annotation_introns_per_tx_median'])
    w.writerows(group_rows)

# ---- effect sizes vs Mucorales (one value per genus)
P('\n== Cliff delta vs Mucorales (genus values). Negative = fewer introns than Mucorales ==')
cmp_groups = {'Kickxellales core (5)': ['Kickxellales core'], 'Kickxellales all (7)': ['Kickxellales core', 'Mycoemilia+Spiromyces'],
              'Harpellales (4)': ['Harpellales'], 'Dimargaritales+Ramicandelaber (4)': ['Dimargaritales+Ramicandelaber'],
              'Zoopagales (8)': ['Zoopagales'], 'Entomophthorales+Basidiobolus (4)': ['Entomophthorales+Basidiobolus']}
eff_rows = []
success = defaultdict(list)
success_ph = defaultdict(list)
subset_all = [g for g in genomes if g['set'] != 'coemansia_strain']
for s in SETTINGS:
    for thr in THR:
        P(f'\n-- setting {s}, identity >= {thr}')
        P(f"{'group':36s} {'med_introns':>11s} {'frac>=1':>8s} {'per_kb':>7s} {'canon/kb*':>9s} {'canon>=1*':>9s}   *post hoc")
        ref = {mt: pooled(s, thr, mt, ['Mucorales'], subset_all) for mt in METRICS}
        for name, gps in cmp_groups.items():
            d = {mt: cliff(pooled(s, thr, mt, gps, subset_all), ref[mt]) for mt in METRICS}
            P(f"{name:36s} {d['median_introns']:11.2f} {d['frac_ge1']:8.2f} {d['per_kb']:7.2f} {d['canon_per_kb']:9.2f} {d['frac_canon_ge1']:9.2f}")
            eff_rows.append([s, thr, name] + [d[mt] for mt in METRICS])
            if name == 'Kickxellales core (5)':
                success[(s, thr)] = [d['median_introns'], d['frac_ge1'], d['per_kb']]
                success_ph[(s, thr)] = [d['canon_per_kb'], d['frac_canon_ge1']]
with open(f'{OUT}/effect_sizes.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['setting', 'ident_thr', 'group'] + [f'cliff_{m}' for m in METRICS])
    w.writerows(eff_rows)
ok = all(all(v <= -0.5 for v in vals) for vals in success.values())
P(f"\nPre-registered success rule (core Kickxellales Cliff delta <= -0.5 vs Mucorales on median introns, share >=1 and per kb, in every "
  f"setting and identity threshold): {'MET' if ok else 'NOT MET'}")
P('  worst case (max) across the 9 setting x threshold cells: ' + ', '.join(
    f"{mt}={max(v[i] for v in success.values()):.2f}" for i, mt in enumerate(['median_introns', 'frac_ge1', 'per_kb'])))
ok_ph = all(all(v <= -0.5 for v in vals) for vals in success_ph.values())
P(f"POST HOC (canonical GT-AG/GC-AG introns only; added after seeing the penalty dependence, so not part of the pre-registered rule): "
  f"core Kickxellales Cliff delta <= -0.5 in every cell: {'MET' if ok_ph else 'NOT MET'}; worst case "
  + ', '.join(f"{mt}={max(v[i] for v in success_ph.values()):.2f}" for i, mt in enumerate(['canon_per_kb', 'frac_canon_ge1'])))

# ---- identity strata (setting A): mean introns per protein in identity bins, genus values
P('\n== Setting A, mean introns per qualifying protein by identity bin (genus values; bin needs >= 15 proteins) ==')
P(f"{'group':34s} " + ' '.join(f'{lo:.1f}-{min(hi, 1):.1f}'.rjust(11) for lo, hi in BINS))
bin_rows = []
bin_vals = {}
for gp in ['Kickxellales core', 'Harpellales', 'Dimargaritales+Ramicandelaber', 'Zoopagales', 'Mucorales', 'Other Mucoromycota group',
           'Saccharomycetes anchors', 'Pezizomycotina anchors']:
    cells = []
    for lo, hi in BINS:
        vals = []
        for g in subset_all:
            if g['group'] != gp:
                continue
            q = qualifying(data['A'][g['locustag']], lo, hi)
            q = [r for r in q if r['qcov'] >= QCOV]
            if len(q) >= MIN_BIN:
                vals.append(sum(r['n_introns'] for r in q) / len(q))
        bin_vals[(gp, lo)] = vals
        cells.append(f'{statistics.median(vals):.2f} (n={len(vals)})'.rjust(11) if vals else '-'.rjust(11))
        bin_rows.append([gp, lo, hi, len(vals), statistics.median(vals) if vals else ''])
    P(f'{gp:34s} ' + ' '.join(cells))
P('Cliff delta, Kickxellales core vs Mucorales, within each bin: ' + ', '.join(
    f"{lo:.1f}-{min(hi, 1):.1f}: {cliff(bin_vals[('Kickxellales core', lo)], bin_vals[('Mucorales', lo)]):.2f}" for lo, hi in BINS))
with open(f'{OUT}/identity_bins.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['group', 'ident_lo', 'ident_hi', 'n_genera', 'median_mean_introns'])
    w.writerows(bin_rows)

# ---- splice sites (setting A and C, ident>=0.5)
P('\n== Splice-site dinucleotides of called introns (group medians of genome shares; ident >= 0.5) ==')
for s in ['A', 'B', 'C']:
    P(f'-- setting {s}')
    for gp in ['Kickxellales core', 'Harpellales', 'Mucorales', 'Saccharomycetes anchors', 'Pezizomycotina anchors']:
        sel = [M[(s, 0.50, g['locustag'])] for g in subset_all if g['group'] == gp and M[(s, 0.50, g['locustag'])]
               and M[(s, 0.50, g['locustag'])]['gt_ag'] == M[(s, 0.50, g['locustag'])]['gt_ag']]
        if sel:
            P(f"{gp:28s} GT-AG {statistics.median(x['gt_ag'] for x in sel):.2f}  GC-AG {statistics.median(x['gc_ag'] for x in sel):.2f}  "
              f"other {statistics.median(x['other_sites'] for x in sel):.2f}   (n genera {len(sel)})")

# ---- Coemansia within-genus spread
P('\n== Within Coemansia (12 contiguous strains + the backbone representative), setting A, ident >= 0.5 ==')
coe = [g for g in genomes if g['genus'] == 'Coemansia']
for mt in ['median_introns', 'frac_ge1', 'per_kb', 'n']:
    v = [M[('A', 0.50, g['locustag'])][mt] for g in coe if M[('A', 0.50, g['locustag'])]]
    P(f'{mt:15s} n={len(v)} min {min(v):.2f} median {statistics.median(v):.2f} max {max(v):.2f}')

# ---- annotation comparison
P('\n== Agreement with funannotate introns per transcript (Spearman across genomes) ==')
cmp_rows = []
for s in SETTINGS:
    for mt in ['per_kb', 'mean_introns', 'frac_ge1']:
        for label, subset in (('83 genomes', genomes), ('63 backbone', backbone)):
            x = [M[(s, 0.50, g['locustag'])][mt] if (M[(s, 0.50, g['locustag'])] and M[(s, 0.50, g['locustag'])]['n'] >= MIN_N) else float('nan')
                 for g in subset]
            y = [annot.get(g['locustag'], float('nan')) for g in subset]
            rho, n = spearman(x, y)
            cmp_rows.append([s, mt, label, n, rho])
            if s in ('A', 'C') and mt == 'per_kb':
                P(f'setting {s}, {mt}, {label}: rho = {rho:.2f} (n = {n})')
with open(f'{OUT}/annotation_comparison.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['setting', 'metric', 'genomes', 'n', 'spearman_rho_vs_annotation_introns_per_tx'])
    w.writerows(cmp_rows)

# ---- calibration
P('\n== Calibration anchors (setting A, ident >= 0.5, per_kb; genome values) ==')
for g in genomes:
    if g['set'].startswith('anchor'):
        m = M[('A', 0.50, g['locustag'])]
        P(f"{g['genus']:18s} {g['set']:22s} n={m['n'] if m else 0:4d} per_kb={m['per_kb'] if m else float('nan'):.2f} annotation={annot.get(g['locustag'], float('nan')):.2f}")

# ---- bootstrap over reference proteins
P('\n== Bootstrap over reference proteins (1,000 resamples; setting A and C, ident >= 0.5; genus-level Cliff delta vs Mucorales) ==')
rng = np.random.default_rng(20261111)
boot_rows = []
for s in ['A', 'C']:
    refs = sorted({r['query'] for g in genomes for r in data[s][g['locustag']]})
    pidx = {p: i for i, p in enumerate(refs)}
    nG, nP = len(genomes), len(refs)
    I = np.zeros((nG, nP)); C = np.zeros((nG, nP)); V = np.zeros((nG, nP))
    for gi, g in enumerate(genomes):
        for r in qualifying(data[s][g['locustag']], 0.50):
            j = pidx[r['query']]
            I[gi, j], C[gi, j], V[gi, j] = r['n_introns'], r['cds_nt'], 1
    gi_by_group = defaultdict(list)
    for gi, g in enumerate(genomes):
        if g['set'] != 'coemansia_strain':
            gi_by_group[g['group']].append(gi)
    comps = {'Kickxellales core': ['Kickxellales core'], 'Kickxellales all': ['Kickxellales core', 'Mycoemilia+Spiromyces'],
             'Harpellales': ['Harpellales']}
    res = defaultdict(list)
    for b in range(1000):
        w = np.bincount(rng.integers(0, nP, nP), minlength=nP).astype(float)
        nq = (V * w).sum(1)
        perkb = np.where((C * w).sum(1) > 0, 1000 * (I * w).sum(1) / np.maximum((C * w).sum(1), 1), np.nan)
        mean_i = np.where(nq > 0, (I * w).sum(1) / np.maximum(nq, 1), np.nan)
        muc = [perkb[i] for i in gi_by_group['Mucorales'] if nq[i] >= MIN_N]
        muc_m = [mean_i[i] for i in gi_by_group['Mucorales'] if nq[i] >= MIN_N]
        for name, gps in comps.items():
            idx = [i for gp in gps for i in gi_by_group[gp] if nq[i] >= MIN_N]
            res[(name, 'per_kb')].append(cliff([perkb[i] for i in idx], muc))
            res[(name, 'mean_introns')].append(cliff([mean_i[i] for i in idx], muc_m))
    for (name, mt), v in sorted(res.items()):
        v = np.array([x for x in v if x == x])
        lo, med, hi = np.percentile(v, [2.5, 50, 97.5])
        share = float((v <= -0.5).mean())
        P(f'setting {s} {name:20s} {mt:13s} delta median {med:.2f} (95% {lo:.2f} to {hi:.2f}); share of resamples <= -0.5: {share:.2f}')
        boot_rows.append([s, name, mt, med, lo, hi, share])
with open(f'{OUT}/bootstrap.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['setting', 'group', 'metric', 'delta_median', 'delta_2.5', 'delta_97.5', 'share_le_minus0.5'])
    w.writerows(boot_rows)

open(f'{OUT}/summary.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
