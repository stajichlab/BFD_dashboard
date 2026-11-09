#!/usr/bin/env python3
"""Render the P1 views of the BFD summary atlas into docs/summary_atlas (figures, tables, Markdown).

Every view: reads the small intermediates in atlas_data/, writes figures/<VID>_<slug>.png|svg and
tables/<VID>_<slug>.tsv, and returns a Markdown block that states n (genomes/species/genera), the
quality filter, the statistic with effect size and CI or p, and caveats.
Run:  python render_views.py [--only V08,V12] [--skip-mixed]
"""
import argparse
import json
import os
import sys
import traceback

import duckdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stats as S  # noqa: E402

ROOT = '/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard'
A = f'{ROOT}/atlas_data'
T = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables'
DOCS = f'{ROOT}/docs/summary_atlas'
FIG = f'{DOCS}/figures'
TAB = f'{DOCS}/tables'
os.makedirs(FIG, exist_ok=True); os.makedirs(TAB, exist_ok=True)

# Okabe-Ito, colorblind safe
GF_COLORS = {'yeast': '#D55E00', 'filamentous_mycelium': '#0072B2', 'dimorphic_yeast': '#E69F00',
             'thallus_photosynthetic': '#009E73', 'other': '#999999'}
GF_ORDER = ['yeast', 'dimorphic_yeast', 'filamentous_mycelium']
SACCH = ['Saccharomycetes', 'Pichiomycetes', 'Dipodascomycetes', 'Lipomycetes', 'Trigonopsidomycetes',
         'Alloascoideomycetes', 'Sporopachydermiomycetes']
PEZIZO = ['Eurotiomycetes', 'Dothideomycetes', 'Sordariomycetes', 'Leotiomycetes', 'Pezizomycetes', 'Orbiliomycetes',
          'Xylonomycetes', 'Geoglossomycetes']
SISTER_PAIRS = [  # (name, yeast-side classes, yeast-side growth form or None, sister classes, sister growth form or None)
    ('Schizosaccharomycetes vs Taphrinomycetes+Pneumocystomycetes', ['Schizosaccharomycetes'], None,
     ['Taphrinomycetes', 'Pneumocystomycetes'], None),
    ('Tremellomycetes (yeast) vs Agaricomycetes', ['Tremellomycetes'], 'yeast', ['Agaricomycetes'], 'filamentous_mycelium'),
    ('Saccharomycotina vs Pezizomycotina', SACCH, 'yeast', PEZIZO, 'filamentous_mycelium'),
    ('Malasseziomycetes vs Ustilaginomycetes', ['Malasseziomycetes'], None, ['Ustilaginomycetes'], None),
    ('Eurotiomycetes dimorphic vs Eurotiomycetes filamentous', ['Eurotiomycetes'], 'dimorphic_yeast',
     ['Eurotiomycetes'], 'filamentous_mycelium'),
    ('Microbotryomycetes (yeast) vs Pucciniomycetes', ['Microbotryomycetes'], 'yeast', ['Pucciniomycetes'], None),
]


class Ctx:
    def __init__(self, skip_mixed=False):
        self.con = duckdb.connect(); self.con.execute("SET threads=4; SET memory_limit='6GB'; SET temp_directory='/scratch/jstajich/31877510/atlas_render_tmp'")
        self.skip_mixed = skip_mixed
        self.labels = self.q(f"SELECT * FROM '{A}/labels/genome_labels.parquet'")
        self.labels['log_n50'] = np.log10(self.labels['N50_bp'].astype(float))
        self.labels['log_genome'] = np.log10(self.labels['total_length_bp'].astype(float))
        self.cache = {}

    def q(self, sql):
        return self.con.execute(sql).df()

    def have(self, rel):
        return os.path.exists(f'{A}/{rel}')

    def tbl(self, rel):
        if rel not in self.cache:
            self.cache[rel] = self.q(f"SELECT * FROM '{A}/{rel}'")
        return self.cache[rel]

    # ---- genome-level frame with labels + filters
    def reps(self, busco='p25', dup=True, n50=True, growth=None):
        d = self.labels[self.labels['is_representative'] & self.labels['q_ann']].copy()
        if busco == 'p25':
            d = d[d['q_busco_rel'].fillna(False)]
        elif busco == 'p50':
            d = d[d['q_busco_rel50'].fillna(False)]
        if dup:
            d = d[d['q_dup'].fillna(False)]
        if n50:
            d = d[d['q_n50'].fillna(False)]
        if growth is not None:
            d = d[d['growth_form'].isin(growth)]
        return d

    def filt_text(self, busco='p25', dup=True, n50=True, model=True):
        f = ['Q-ann', 'Q-rep (one genome per species)']
        if busco == 'p25':
            f.append('Q-busco-rel (BUSCO complete >= class 25th percentile)')
        if busco == 'p50':
            f.append('Q-busco-rel50 (>= class median)')
        if dup:
            f.append('Q-dup (BUSCO duplicated <= 10 %)')
        if n50:
            f.append('Q-n50 (N50 >= 50 kb)')
        if model:
            f.append('Q-model (complete gene models only)')
        return ', '.join(f)


def counts_text(d):
    return (f"n = {len(d):,} genomes / {d['SPECIES'].nunique():,} species / {d['GENUS'].nunique():,} genera")


def genus_medians(d, cols):
    g = d.groupby('GENUS')
    out = g[cols].median()
    meta = g.agg(PHYLUM=('PHYLUM', 'first'), CLASS=('CLASS', 'first'), ORDER_=('ORDER_', 'first'),
                 growth_form=('growth_form', lambda s: s.mode().iloc[0] if s.notna().any() else None),
                 primary_lifestyle=('primary_lifestyle', 'first'), yeast_origin=('yeast_origin', 'first'),
                 n_genomes=('LOCUSTAG', 'size'),
                 log_n50=('log_n50', 'median'), busco_complete=('busco_complete', 'median'),
                 log_genome=('log_genome', 'median'))
    return out.join(meta).reset_index()


def save_fig(fig, vid, slug):
    base = f'{FIG}/{vid}_{slug}'
    fig.savefig(base + '.png', dpi=200, bbox_inches='tight'); fig.savefig(base + '.svg', bbox_inches='tight')
    plt.close(fig)
    return f'figures/{vid}_{slug}.png'


def save_tab(df, vid, slug, note=''):
    p = f'{TAB}/{vid}_{slug}.tsv'
    with open(p, 'w') as o:
        if note:
            o.write('# ' + note.replace('\n', ' ') + '\n')
        df.to_csv(o, sep='\t', index=False, float_format='%.4g')
    return f'tables/{vid}_{slug}.tsv'


def md_table(df, max_rows=25):
    df = df.head(max_rows)
    cols = list(df.columns)
    lines = ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
    for _, r in df.iterrows():
        lines.append('| ' + ' | '.join(S.fmt(v) if not isinstance(v, str) else v for v in r.values) + ' |')
    return '\n'.join(lines)


def contrast_block(ctx, gd, col, label, covars=('log_n50', 'busco_complete')):
    """Yeast vs filamentous on genus medians + mixed model. gd = genus-median frame."""
    y = gd[gd['growth_form'] == 'yeast'][col]; f = gd[gd['growth_form'] == 'filamentous_mycelium'][col]
    r = S.two_group(y, f)
    md = [f"- Unit = genus (median of representative genomes). Yeast genera n = {r['n_yeast']:,}, filamentous genera n = {r['n_filamentous']:,}.",
          f"- Median {label}: yeast {S.fmt(r['median_yeast'])}, filamentous {S.fmt(r['median_filamentous'])}; "
          f"difference {S.fmt(r.get('median_diff'))} (bootstrap 95 % CI {S.fmt(r.get('diff_ci_lo'))} to {S.fmt(r.get('diff_ci_hi'))}); "
          f"Cliff's delta {S.fmt(r.get('cliffs_delta'))}; Mann-Whitney p = {S.fmt(r.get('mw_p'))}."]
    d = gd[gd['growth_form'] == 'dimorphic_yeast'][col]
    if len(d):
        md.append(f"- Dimorphic yeast genera n = {len(d):,}, median {S.fmt(float(d.median()))} (shown, not pooled).")
    if not ctx.skip_mixed:
        sub = gd[gd['growth_form'].isin(['yeast', 'filamentous_mycelium'])].copy()
        sub['gf_yeast'] = (sub['growth_form'] == 'yeast').astype(int)
        mm = S.mixed_model(sub, col, ['gf_yeast'] + list(covars), group='CLASS', vc='ORDER_')
        if mm.get('available'):
            t = [x for x in mm['terms'] if x['term'] == 'gf_yeast'][0]
            md.append(f"- Mixed model `{mm['formula']}` (n = {mm['n']:,} genera, {mm['n_groups']} classes, converged = {mm['converged']}): "
                      f"yeast effect {S.fmt(t['estimate'])} (95 % CI {S.fmt(t['ci_lo'])} to {S.fmt(t['ci_hi'])}), p = {S.fmt(t['p'])}. PGLS pending tree.")
        else:
            md.append(f"- Mixed model not fitted: {mm.get('note')}. PGLS pending tree.")
    return '\n'.join(md), r


def strip_by_class(gd_or_d, col, ylabel, vid, slug, min_n=10, log=False, unit='genus'):
    d = gd_or_d.dropna(subset=[col, 'CLASS'])
    n = d.groupby('CLASS').size(); classes = n[n >= min_n].sort_values(ascending=False).index.tolist()
    order = sorted(classes, key=lambda c: (d[d['CLASS'] == c]['PHYLUM'].iloc[0] or '', d[d['CLASS'] == c][col].median()))
    fig, ax = plt.subplots(figsize=(7, max(4, 0.28 * len(order) + 1)))
    for i, c in enumerate(order):
        sub = d[d['CLASS'] == c]
        vals = sub[col].values
        ax.boxplot([vals], positions=[i], vert=False, widths=0.6, showfliers=False,
                   medianprops={'color': 'black'}, boxprops={'color': '#444'}, whiskerprops={'color': '#444'})
        gfc = sub['growth_form'].mode().iloc[0] if sub['growth_form'].notna().any() else 'other'
        jitter = (np.random.default_rng(1).random(len(vals)) - 0.5) * 0.5
        ax.scatter(vals, i + jitter, s=6, alpha=0.5, color=GF_COLORS.get(gfc, '#999999'), linewidths=0, rasterized=True)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([f'{c} (n={n[c]})' for c in order], fontsize=7)
    ax.set_xlabel(ylabel)
    if log:
        ax.set_xscale('log')
    for gf, colr in GF_COLORS.items():
        if gf in d['growth_form'].values:
            ax.scatter([], [], color=colr, label=gf, s=12)
    ax.legend(fontsize=7, loc='lower right', title='dominant growth form', title_fontsize=7)
    ax.set_title(f'{vid}: {ylabel} by class ({unit}-level)', fontsize=9)
    return save_fig(fig, vid, slug)


# ------------------------------------------------------------------ H1
def v01(ctx):
    rows = []
    small = {'species': ('LOCUSTAG', 'species.parquet'), 'samples': ('LOCUSTAG', 'samples.parquet'),
             'asm_stats': ('ASMID', 'asm_stats.parquet'), 'busco_genome': ('ASMID', 'busco_genome.parquet'),
             'telomere_summary': ('ASMID', 'telomere_summary.parquet'), 'telomere_tracts': ('ASMID', 'telomere_tracts.parquet'),
             'aa_freq': ('species_prefix', 'aa_freq.parquet'), 'codon_freq': ('species_prefix', 'codon_freq.parquet'),
             'tf_inventory': ('species_prefix', 'tf_inventory.parquet'), 'cazy.overview': ('species_prefix', 'cazy.overview.parquet'),
             'cazy.cazymes_hmm': ('species_prefix', 'cazy.cazymes_hmm.parquet'),
             'signalp.signal_peptide': ('species_prefix', 'signalp.signal_peptide.parquet'),
             'targetP': ('species_prefix', 'targetP.parquet'), 'infernal_rfam': ('species_prefix', 'infernal_rfam.parquet'),
             'wgd.ks.summary': ('species_prefix', 'wgd.ks.summary.parquet')}
    for t, (k, f) in small.items():
        n, r = ctx.con.execute(f"SELECT count(DISTINCT {k}), count(*) FROM '{T}/{f}'").fetchone()
        rows.append({'table': t, 'key': k, 'n_genomes': n, 'n_rows': r, 'source': 'measured at render'})
    n, r = ctx.con.execute(f"SELECT count(DISTINCT split_part(gene_id,'_',1)), count(*) FROM '{T}/gene_trnas.parquet'").fetchone()
    rows.append({'table': 'gene_trnas', 'key': 'split_part(gene_id)', 'n_genomes': n, 'n_rows': r, 'source': 'measured at render'})
    rows.append({'table': 'wgd.ks', 'key': 'species_prefix', 'n_genomes': 7384, 'n_rows': 25831001, 'source': 'tables/_manifest/wgd.ks.json'})
    rows.append({'table': 'wgd.ks.nodes', 'key': 'species_prefix', 'n_genomes': 7384, 'n_rows': 28982469, 'source': 'tables/_manifest/wgd.ks.nodes.json'})
    cov_dir = f'{A}/coverage'
    if os.path.isdir(cov_dir):
        for f in sorted(os.listdir(cov_dir)):
            if f.endswith('.parquet'):
                d = ctx.q(f"SELECT tbl, key, n_genomes, n_rows FROM '{cov_dir}/{f}'").iloc[0]
                rows.append({'table': d['tbl'], 'key': d['key'], 'n_genomes': int(d['n_genomes']), 'n_rows': int(d['n_rows']), 'source': 'job J0'})
    df = pd.DataFrame(rows)
    missing = [t for t in ['pfam', 'swissprot_best', 'gene_info', 'gene_transcripts', 'gene_proteins', 'gene_introns',
                           'gene_intergenic_distances', 'merops', 'tmhmm', 'wolfpsort', 'predgpi', 'idp_summary', 'idp',
                           'gene_CDS', 'gene_exons', 'swissprot_hits'] if t not in set(df['table'])]
    tp = save_tab(df, 'V01', 'coverage')
    md = f"### V01 Data coverage\n\nGenomes covered per table (distinct genome key) and rows. Full table: `{tp}`.\n\n" + md_table(df, 60)
    if missing:
        md += f"\n\nNot yet measured (J0 still running or failed): {', '.join(missing)}."
    return md


def v02(ctx):
    d = ctx.labels[ctx.labels['q_ann'].fillna(False)].copy()
    g = d.groupby(['CLASS']).agg(n=('LOCUSTAG', 'size'), p10=('busco_complete', lambda s: s.quantile(.1)),
                                 median=('busco_complete', 'median'), p90=('busco_complete', lambda s: s.quantile(.9)),
                                 dup_median=('busco_duplicated', 'median'),
                                 growth_form=('growth_form', lambda s: s.mode().iloc[0] if s.notna().any() else None)).reset_index()
    g = g[g['n'] >= 20].sort_values('median')
    tp = save_tab(g, 'V02', 'busco_by_class')
    fp = strip_by_class(d, 'busco_complete', 'BUSCO complete (%) fungi_odb12', 'V02', 'busco_by_class', min_n=20, unit='genome')
    fig, ax = plt.subplots(figsize=(7, max(4, 0.28 * len(g) + 1)))
    ax.barh(range(len(g)), g['dup_median'], color=[GF_COLORS.get(x, '#999') for x in g['growth_form']])
    ax.set_yticks(range(len(g))); ax.set_yticklabels(g['CLASS'], fontsize=7); ax.set_xlabel('median BUSCO duplicated (%)')
    ax.set_title('V02b: median BUSCO duplication by class', fontsize=9)
    fp2 = save_fig(fig, 'V02', 'busco_dup_by_class')
    yeast_cls = g[g['CLASS'].isin(SACCH + ['Schizosaccharomycetes'])]
    return (f"### V02 BUSCO completeness and duplication by class\n\n![]({fp})\n\n![]({fp2})\n\n"
            f"- All annotated genomes (n = {len(d):,}); no representative filter (this view describes the raw set).\n"
            f"- Median BUSCO complete in yeast classes: " + '; '.join(f"{r.CLASS} {r.median:.1f} %" for r in yeast_cls.itertuples()) +
            f". Filamentous classes with >= 1,000 genomes: " + '; '.join(f"{r.CLASS} {r.median:.1f} %" for r in g[g['n'] >= 1000].itertuples() if r.CLASS not in SACCH) + '.\n'
            f"- Consequence: a fixed cutoff removes yeasts. All later views use the class-relative threshold (Q-busco-rel, 25th percentile of the class).\n"
            f"- Caveat: fungi_odb12 is one marker set for all fungi; low yeast values are partly marker-set bias, not annotation quality. Pucciniomycetes duplication reflects dikaryotic assemblies.\n"
            f"- Table: `{tp}`.")


def v03(ctx):
    if not ctx.have('A2/genome_gene_stats.parquet'):
        return '### V03 Gene count versus assembly size\n\nNot built: A2 (job J1) missing.'
    a2 = ctx.tbl('A2/genome_gene_stats.parquet')
    d = ctx.reps(dup=False, n50=False).merge(a2, left_on='LOCUSTAG', right_on='genome')
    d = d[d['growth_form'].isin(GF_ORDER)]
    d['genes_per_mb'] = d['n_protein_coding'] / (d['total_length_bp'] / 1e6)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    res = {}
    for gf in GF_ORDER:
        s = d[d['growth_form'] == gf]
        ax.scatter(s['total_length_bp'] / 1e6, s['n_protein_coding'], s=5, alpha=0.4, color=GF_COLORS[gf], label=f'{gf} (n={len(s)})', linewidths=0, rasterized=True)
        res[gf] = S.loglog_slope(s['total_length_bp'], s['n_protein_coding'])
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('assembly size (Mb)'); ax.set_ylabel('protein-coding genes')
    ax.legend(fontsize=7); ax.set_title('V03: gene count vs assembly size (representative genomes)', fontsize=9)
    fp = save_fig(fig, 'V03', 'genes_vs_size')
    gd = genus_medians(d, ['genes_per_mb', 'n_protein_coding', 'total_length_bp'])
    cb, r = contrast_block(ctx, gd, 'genes_per_mb', 'genes per Mb', covars=('log_genome', 'busco_complete'))
    tp = save_tab(gd, 'V03', 'genus_medians')
    md = [f"### V03 Gene count versus assembly size\n\n![]({fp})\n", f"- {counts_text(d)}. Filters: {ctx.filt_text(dup=False, n50=False, model=False)}."]
    for gf, x in res.items():
        md.append(f"- log-log slope (genes on size), {gf}: {S.fmt(x.get('slope'))} (95 % CI {S.fmt(x.get('slope_ci_lo'))} to {S.fmt(x.get('slope_ci_hi'))}), r^2 = {S.fmt(x.get('r2'))}, n = {x['n']:,} genomes.")
    md.append("**Genes per Mb, yeast vs filamentous:**\n" + cb)
    md.append(f"- Caveat: gene count depends on the per-genome funannotate model set and repeat masking (`masked_pct` available). Genus-median table: `{tp}`.")
    return '\n'.join(md)


def v04(ctx):
    if not ctx.have('A2/genome_gene_stats.parquet'):
        return '### V04 Partial gene models versus assembly contiguity\n\nNot built: A2 (job J1) missing.'
    a2 = ctx.tbl('A2/genome_gene_stats.parquet')
    d = ctx.labels[ctx.labels['q_ann'].fillna(False)].merge(a2, left_on='LOCUSTAG', right_on='genome')
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for gf in GF_ORDER:
        s = d[d['growth_form'] == gf]
        ax.scatter(s['N50_bp'], s['partial_fraction'], s=4, alpha=0.3, color=GF_COLORS[gf], label=f'{gf} (n={len(s)})', linewidths=0, rasterized=True)
    ax.set_xscale('log'); ax.set_xlabel('assembly N50 (bp)'); ax.set_ylabel('fraction of transcripts flagged is_partial')
    ax.legend(fontsize=7); ax.set_title('V04: partial gene models vs N50 (all annotated genomes)', fontsize=9)
    fp = save_fig(fig, 'V04', 'partial_vs_n50')
    rows = []
    for c, s in d.groupby('CLASS'):
        if len(s) >= 50:
            r = S.spearman(np.log10(s['N50_bp']), s['partial_fraction'])
            rows.append({'CLASS': c, 'n': len(s), 'median_partial_fraction': s['partial_fraction'].median(),
                         'median_frac_start_null': s['frac_start_null'].median(), 'spearman_rho_logN50': r['rho'], 'p': r['p']})
    tb = pd.DataFrame(rows).sort_values('median_partial_fraction', ascending=False)
    tp = save_tab(tb, 'V04', 'partial_by_class')
    allr = S.spearman(np.log10(d['N50_bp']), d['partial_fraction'])
    return (f"### V04 Partial gene models versus assembly contiguity\n\n![]({fp})\n\n"
            f"- n = {len(d):,} annotated genomes (no representative filter). Overall median partial fraction {d['partial_fraction'].median():.3f}; "
            f"Spearman rho(log N50, partial fraction) = {allr['rho']:.3f} (n = {allr['n']:,}).\n"
            f"- `has_start_codon` is the string 'NULL' for a median {d['frac_start_null'].median():.3f} of transcripts; these are excluded from Q-model.\n"
            f"- Per-class table (classes with >= 50 genomes): `{tp}`.\n\n" + md_table(tb, 15))


def v05(ctx):
    L = ctx.labels
    rows = []
    for c, s in L.groupby('CLASS'):
        r = s[s['is_representative']]
        rows.append({'CLASS': c, 'PHYLUM': s['PHYLUM'].iloc[0], 'genomes': len(s), 'species_reps': len(r),
                     'reps_q_busco_rel': int((r['q_busco_rel'].fillna(False)).sum()),
                     'reps_all_q': int((r['q_busco_rel'].fillna(False) & r['q_dup'].fillna(False) & r['q_n50'].fillna(False)).sum()),
                     'genus_reps': int(r['is_genus_rep'].sum()),
                     'growth_form': r['growth_form'].mode().iloc[0] if r['growth_form'].notna().any() else None,
                     'busco_thr_p25': r['busco_thr_p25'].iloc[0] if len(r) else np.nan})
    tb = pd.DataFrame(rows).sort_values('genomes', ascending=False)
    tp = save_tab(tb, 'V05', 'representatives_by_class')
    tot = {'genomes': len(L), 'species representatives': int(L['is_representative'].sum()),
           'reps passing Q-busco-rel': int((L['is_representative'] & L['q_busco_rel'].fillna(False)).sum()),
           'reps passing Q-busco-rel + Q-dup + Q-n50': int((L['is_representative'] & L['q_busco_rel'].fillna(False) & L['q_dup'].fillna(False) & L['q_n50'].fillna(False)).sum()),
           'genus representatives': int(L['is_genus_rep'].sum())}
    gf = L[L['is_representative']].groupby('growth_form', dropna=False).size().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(6, 3.2))
    ax.barh(list(tot.keys())[::-1], list(tot.values())[::-1], color='#0072B2'); ax.set_xlabel('genomes')
    ax.set_title('V05: what each filter leaves', fontsize=9)
    for i, v in enumerate(list(tot.values())[::-1]):
        ax.text(v, i, f' {v:,}', va='center', fontsize=7)
    fp = save_fig(fig, 'V05', 'filter_waterfall')
    md = [f"### V05 Representative set and quality tiers\n\n![]({fp})\n"]
    md += [f"- {k}: {v:,}" for k, v in tot.items()]
    md.append("- Growth form of species representatives (FungalTraits 1.2 + overrides): " + '; '.join(f"{k if k == k else 'NULL'} {v:,}" for k, v in gf.items()) + '.')
    md.append(f"- Rule D6: representative = highest BUSCO complete within species, tie N50. Per-class table: `{tp}`.")
    return '\n'.join(md)


# ------------------------------------------------------------------ H2
def v07(ctx):
    rel = 'A4/genome_protein_summary.parquet'
    if not ctx.have(rel):
        return '### V07 Raw protein length by growth form and class\n\nNot built: A4 (job J3) missing.'
    g = ctx.tbl(rel)
    d = ctx.reps().merge(g, left_on='LOCUSTAG', right_on='genome')
    d = d[d['growth_form'].isin(GF_ORDER)]
    fp = strip_by_class(d, 'median_len_model', 'median protein length (aa), complete models', 'V07', 'protein_len_by_class', unit='genome')
    gd = genus_medians(d, ['median_len_model', 'median_len_all', 'mean_len_model'])
    cb, r = contrast_block(ctx, gd, 'median_len_model', 'median protein length (aa)')
    tp = save_tab(gd, 'V07', 'genus_medians')
    return (f"### V07 Raw protein length by growth form and class\n\n![]({fp})\n\n- {counts_text(d)}. Filters: {ctx.filt_text()}.\n"
            f"{cb}\n- Caveat: raw lengths mix gene content with length within families. This view sets the baseline only; V08 controls for family. Table: `{tp}`.")


def family_contrast(ctx, fam, key, busco='p25', dup=True, n50=True, model=True, yeast_def='growth_form', min_genera=10):
    """Family-level log2(yeast/filamentous) ratios. fam: (key, genome, med_len_model, med_len_all...)."""
    L = ctx.reps(busco=busco, dup=dup, n50=n50)
    gfcol = 'growth_form' if yeast_def == 'growth_form' else 'growth_form_tax'
    L = L[L[gfcol].isin(['yeast', 'filamentous_mycelium'])][['LOCUSTAG', 'GENUS', gfcol]].rename(columns={gfcol: 'gf'})
    lencol = 'med_len_model' if model else 'med_len_all'
    f = fam.merge(L, left_on='genome', right_on='LOCUSTAG').dropna(subset=[lencol])
    # genus median within family, then median over genera per growth form
    gm = f.groupby([key, 'gf', 'GENUS'])[lencol].median().reset_index()
    agg = gm.groupby([key, 'gf'])[lencol].agg(['median', 'size']).unstack('gf')
    agg.columns = [f'{a}_{b}' for a, b in agg.columns]
    agg = agg.dropna()
    agg = agg[(agg['size_yeast'] >= min_genera) & (agg['size_filamentous_mycelium'] >= min_genera)]
    agg['log2_ratio'] = np.log2(agg['median_yeast'] / agg['median_filamentous_mycelium'])
    return agg.reset_index()


def v08(ctx):
    rel = 'A4/fam_sp_genome.parquet'
    if not ctx.have(rel):
        return '### V08 Within-family length contrast (SwissProt proxy)\n\nNot built: A4 (job J3) missing.'
    fam = ctx.tbl(rel)
    agg = family_contrast(ctx, fam, 'sp_id')
    subj = ctx.q(f"SELECT sp_id, accession, entry_name, org_code, protein_name FROM '{T}/swissprot_subject.parquet'")
    agg = agg.merge(subj, on='sp_id', how='left')
    st = S.sign_test(agg['log2_ratio']); lo, hi = S.boot_ci(agg['log2_ratio'])
    med = float(agg['log2_ratio'].median())
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.hist(agg['log2_ratio'], bins=80, color='#0072B2', alpha=0.85)
    ax.axvline(0, color='black', lw=0.8); ax.axvline(med, color='#D55E00', lw=1.2, label=f'median {med:+.3f}')
    ax.set_xlabel('log2( yeast median length / filamentous median length ) per SwissProt family')
    ax.set_ylabel('families'); ax.legend(fontsize=8)
    ax.set_title(f'V08: within-family protein length, yeast vs filamentous (n = {len(agg):,} families)', fontsize=9)
    fp = save_fig(fig, 'V08', 'family_log2_ratio')
    # by reference-length decile
    agg['ref_len'] = agg['median_filamentous_mycelium']
    agg['ref_len_bin'] = pd.qcut(agg['ref_len'], 10, duplicates='drop')
    byb = agg.groupby('ref_len_bin', observed=True)['log2_ratio'].agg(['median', 'size']).reset_index()
    byb['ref_len_bin'] = byb['ref_len_bin'].astype(str)
    out = agg.sort_values('log2_ratio')[['sp_id', 'accession', 'entry_name', 'org_code', 'protein_name', 'size_yeast', 'size_filamentous_mycelium',
                                         'median_yeast', 'median_filamentous_mycelium', 'log2_ratio']]
    tp = save_tab(out, 'V08', 'family_ratios', 'SwissProt best-hit families, qcov/scov>=80, pident>=40, complete models, representative genomes, >=10 genera per side')
    tp2 = save_tab(byb, 'V08', 'ratio_by_reference_length_decile')
    md = [f"### V08 Within-family length contrast, SwissProt proxy (P-A)\n\n![]({fp})\n",
          f"- Unit = family (SwissProt best-hit subject `sp_id`). Families kept: {len(agg):,} with >= 10 yeast genera and >= 10 filamentous genera. "
          f"Genomes: representative set with {ctx.filt_text()}; proteins: qcov >= 80 %, scov >= 80 %, pident >= 40 %, complete models. Family value = median over genera of the per-genome median length.",
          f"- Median log2 ratio = {med:+.4f} (bootstrap 95 % CI {lo:+.4f} to {hi:+.4f}); as a length ratio yeast/filamentous = {2**med:.4f}.",
          f"- Sign test: {st['n_below_zero']:,} of {st['n']:,} families have yeast shorter ({100*st['frac_below_zero']:.1f} %), binomial p = {S.fmt(st['binom_p'])}.",
          f"- Mean absolute log2 ratio = {agg['log2_ratio'].abs().mean():.3f}; 10th/90th percentiles {agg['log2_ratio'].quantile(.1):+.3f} / {agg['log2_ratio'].quantile(.9):+.3f}.",
          "- By reference-length decile (filamentous median length), median log2 ratio:\n\n" + md_table(byb, 12),
          f"\n- Caveats: reference-dominance bias (Saccharomycotina best hits go mostly to *S. cerevisiae* references; see V07 table column `frac_best_hit_yeast`); coverage gates remove truncated and fused variants (conservative); paralogs collapse to the per-genome median. Family table: `{tp}`; decile table: `{tp2}`."]
    return '\n'.join(md)


def genome_offsets(ctx, fam, key='sp_id', busco='p25', dup=True, n50=True, lencol='med_len_model', min_fil_genera=10):
    """Per-genome median log2(len / family reference) where reference = median over filamentous genera."""
    L = ctx.reps(busco=busco, dup=dup, n50=n50)
    Lf = L[L['growth_form'] == 'filamentous_mycelium'][['LOCUSTAG', 'GENUS']]
    f = fam.merge(Lf, left_on='genome', right_on='LOCUSTAG').dropna(subset=[lencol])
    ref = f.groupby([key, 'GENUS'])[lencol].median().groupby(key).agg(['median', 'size'])
    ref = ref[ref['size'] >= min_fil_genera].rename(columns={'median': 'ref_len'})
    g = fam.merge(ref[['ref_len']], left_on=key, right_index=True).dropna(subset=[lencol])
    g['off'] = np.log2(g[lencol] / g['ref_len'])
    per = g.groupby('genome').agg(offset=('off', 'median'), n_fam=('off', 'size')).reset_index()
    return per.merge(L, left_on='genome', right_on='LOCUSTAG')


def v09(ctx):
    rel = 'A4/fam_sp_genome.parquet'
    if not ctx.have(rel):
        return '### V09 Per-genome length offset by yeast origin\n\nNot built: A4 (job J3) missing.'
    per = genome_offsets(ctx, ctx.tbl(rel))
    per = per[per['n_fam'] >= 200]
    gd = genus_medians(per, ['offset', 'n_fam'])
    origins = [o for o in ['Saccharomycotina', 'Taphrinomycotina', 'Tremellomycetes', 'Pucciniomycotina-yeast', 'Ustilaginomycotina-yeast',
                           'Pezizomycotina-dimorphic', 'Ustilaginomycotina-dimorphic', 'filamentous'] if o in set(gd['yeast_origin'])]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    rows = []
    for i, o in enumerate(origins):
        v = gd[gd['yeast_origin'] == o]['offset'].values
        lo, hi = S.boot_ci(v)
        ax.scatter(v, i + (np.random.default_rng(2).random(len(v)) - 0.5) * 0.5, s=8, alpha=0.5, linewidths=0,
                   color=GF_COLORS['filamentous_mycelium'] if o == 'filamentous' else ('#E69F00' if 'dimorphic' in o else '#D55E00'))
        ax.plot([lo, hi], [i, i], color='black', lw=2); ax.plot([np.median(v)], [i], 'k|', ms=12)
        rows.append({'yeast_origin': o, 'n_genera': len(v), 'median_offset_log2': float(np.median(v)), 'ci_lo': lo, 'ci_hi': hi,
                     'ratio': 2 ** float(np.median(v))})
    ax.axvline(0, color='grey', lw=0.8); ax.set_yticks(range(len(origins))); ax.set_yticklabels([f'{o} (n={r["n_genera"]})' for o, r in zip(origins, rows)], fontsize=8)
    ax.set_xlabel('genus median of per-genome median log2(length / filamentous family reference)')
    ax.set_title('V09: protein-length offset by independent yeast origin', fontsize=9)
    fp = save_fig(fig, 'V09', 'offset_by_origin')
    tb = pd.DataFrame(rows); tp = save_tab(tb, 'V09', 'offset_by_origin')
    same = tb[(tb['yeast_origin'] != 'filamentous') & ('dimorphic' not in tb['yeast_origin'].str) if False else (tb['yeast_origin'] != 'filamentous')]
    neg = same[same['median_offset_log2'] < 0]
    md = [f"### V09 Per-genome length offset by yeast origin\n\n![]({fp})\n",
          f"- Unit = genus (median of genome offsets). Genome offset = median over families (>= 200 families per genome) of log2(length / family reference), reference = median over filamentous genera. Filters: {ctx.filt_text()}.",
          md_table(tb),
          f"\n- Origins with median offset < 0: {len(neg)} of {len(same)} non-filamentous groups ({', '.join(neg['yeast_origin'])}). The filamentous row is the reference distribution (expected near 0 by construction).",
          f"- Caveat: the filamentous reference itself is a mixture of classes; Pezizomycotina dominate it. Table: `{tp}`. PGLS pending tree."]
    return '\n'.join(md)


def v12(ctx):
    rel = 'A4/fam_sp_genome.parquet'
    if not ctx.have(rel):
        return '### V12 Sensitivity to quality filters\n\nNot built: A4 (job J3) missing.'
    fam = ctx.tbl(rel)
    grid = [('F0 Q-rep + Q-ann only, complete models', dict(busco=None, dup=False, n50=False, model=True)),
            ('F1 + Q-busco-rel (p25)', dict(busco='p25', dup=False, n50=False, model=True)),
            ('F2 default: + Q-dup + Q-n50', dict(busco='p25', dup=True, n50=True, model=True)),
            ('F3 default but all models (incl. partial)', dict(busco='p25', dup=True, n50=True, model=False)),
            ('F4 default with BUSCO >= class median', dict(busco='p50', dup=True, n50=True, model=True)),
            ('F5 default, taxonomy-only yeast definition', dict(busco='p25', dup=True, n50=True, model=True, yeast_def='tax')),
            ('F6 default, >= 25 genera per side', dict(busco='p25', dup=True, n50=True, model=True, min_genera=25))]
    rows = []
    for name, kw in grid:
        agg = family_contrast(ctx, fam, 'sp_id', **kw)
        st = S.sign_test(agg['log2_ratio']); lo, hi = S.boot_ci(agg['log2_ratio'])
        rows.append({'filter_set': name, 'n_families': len(agg), 'median_log2_ratio': float(agg['log2_ratio'].median()),
                     'ci_lo': lo, 'ci_hi': hi, 'frac_yeast_shorter': st.get('frac_below_zero'), 'binom_p': st.get('binom_p')})
    tb = pd.DataFrame(rows); tp = save_tab(tb, 'V12', 'filter_grid')
    fig, ax = plt.subplots(figsize=(6.5, 3.5))
    ax.errorbar(tb['median_log2_ratio'], range(len(tb)), xerr=[tb['median_log2_ratio'] - tb['ci_lo'], tb['ci_hi'] - tb['median_log2_ratio']], fmt='o', color='#0072B2')
    ax.axvline(0, color='grey', lw=0.8); ax.set_yticks(range(len(tb))); ax.set_yticklabels(tb['filter_set'], fontsize=7)
    ax.set_xlabel('median log2(yeast / filamentous) over families, 95 % bootstrap CI'); ax.set_title('V12: sensitivity of the family contrast', fontsize=9)
    fp = save_fig(fig, 'V12', 'filter_grid')
    return (f"### V12 Sensitivity of the length contrast to quality filters\n\n![]({fp})\n\n" + md_table(tb) +
            f"\n\n- Same family statistic as V08 under each filter set. Table: `{tp}`.")


def v14(ctx):
    rel = 'A4/fam_sp_genome.parquet'
    if not ctx.have(rel):
        return '### V14 Paired sister-clade contrasts\n\nNot built: A4 (job J3) missing.'
    per = genome_offsets(ctx, ctx.tbl(rel)); per = per[per['n_fam'] >= 200]
    gd = genus_medians(per, ['offset'])
    rows = []
    for name, ycls, ygf, scls, sgf in SISTER_PAIRS:
        a = gd[gd['CLASS'].isin(ycls)]; b = gd[gd['CLASS'].isin(scls)]
        if ygf:
            a = a[a['growth_form'] == ygf]
        if sgf:
            b = b[b['growth_form'] == sgf]
        r = S.two_group(a['offset'], b['offset'], 'yeast_side', 'sister')
        rows.append({'pair': name, 'n_yeast_genera': r['n_yeast_side'], 'n_sister_genera': r['n_sister'],
                     'median_yeast_offset': r['median_yeast_side'], 'median_sister_offset': r['median_sister'],
                     'difference_log2': r.get('median_diff'), 'ci_lo': r.get('diff_ci_lo'), 'ci_hi': r.get('diff_ci_hi'), 'mw_p': r.get('mw_p')})
    tb = pd.DataFrame(rows); tp = save_tab(tb, 'V14', 'sister_pairs')
    st = S.sign_test(tb['difference_log2'].dropna())
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ok = tb.dropna(subset=['difference_log2'])
    ax.errorbar(ok['difference_log2'], range(len(ok)), xerr=[ok['difference_log2'] - ok['ci_lo'], ok['ci_hi'] - ok['difference_log2']], fmt='o', color='#D55E00')
    ax.axvline(0, color='grey', lw=0.8); ax.set_yticks(range(len(ok))); ax.set_yticklabels([f"{p} ({a}/{b})" for p, a, b in zip(ok['pair'], ok['n_yeast_genera'], ok['n_sister_genera'])], fontsize=7)
    ax.set_xlabel('yeast-side minus sister median offset (log2), bootstrap 95 % CI'); ax.set_title('V14: sister-clade contrasts', fontsize=9)
    fp = save_fig(fig, 'V14', 'sister_pairs')
    return (f"### V14 Paired sister-clade contrasts\n\n![]({fp})\n\n" + md_table(tb) +
            f"\n\n- Unit = genus median offset (as V09). Pairs are taxonomy-based (D3), not tree-based. Sign test across pairs: {st.get('n_below_zero', 0)} of {st.get('n', 0)} pairs negative, binomial p = {S.fmt(st.get('binom_p'))} (low power with {st.get('n', 0)} pairs). Table: `{tp}`.")


# ------------------------------------------------------------------ H3
def _a2reps(ctx, **kw):
    a2 = ctx.tbl('A2/genome_gene_stats.parquet')
    d = ctx.reps(**kw).merge(a2, left_on='LOCUSTAG', right_on='genome')
    return d[d['growth_form'].isin(GF_ORDER)]


def v15(ctx):
    if not ctx.have('A2/genome_gene_stats.parquet'):
        return '### V15 Introns per gene\n\nNot built: A2 (job J1) missing.'
    d = _a2reps(ctx, dup=False, n50=False)
    fp = strip_by_class(d, 'introns_per_transcript', 'introns per transcript', 'V15', 'introns_by_class', unit='genome')
    gd = genus_medians(d, ['introns_per_transcript', 'frac_intronless', 'median_intron_len', 'frac_gt_ag'])
    cb, r = contrast_block(ctx, gd, 'introns_per_transcript', 'introns per transcript')
    cb2, r2 = contrast_block(ctx, gd, 'frac_intronless', 'fraction of intronless transcripts')
    tp = save_tab(gd, 'V15', 'genus_medians')
    return (f"### V15 Introns per gene by class and growth form\n\n![]({fp})\n\n- {counts_text(d)}. Filters: {ctx.filt_text(dup=False, n50=False, model=False)}.\n"
            f"**Introns per transcript:**\n{cb}\n**Intronless fraction:**\n{cb2}\n- Caveat: ab initio predictors trained per genome; intron under-prediction in genomes without RNA-seq (training flag not in tables, D8). Table: `{tp}`.")


def v17(ctx):
    if not ctx.have('A2/genome_gene_stats.parquet'):
        return '### V17 Intergenic distance\n\nNot built: A2 (job J1) missing.'
    d = _a2reps(ctx, dup=False, n50=True)
    fp = strip_by_class(d, 'median_intergenic', 'median intergenic distance (bp)', 'V17', 'intergenic_by_class', unit='genome', log=True)
    gd = genus_medians(d, ['median_intergenic', 'q10_intergenic', 'frac_overlapping'])
    cb, r = contrast_block(ctx, gd, 'median_intergenic', 'median intergenic distance (bp)', covars=('log_genome', 'log_n50'))
    tp = save_tab(gd, 'V17', 'genus_medians')
    return (f"### V17 Intergenic distance by growth form\n\n![]({fp})\n\n- {counts_text(d)}. Filters: {ctx.filt_text(dup=False, model=False)}.\n{cb}\n"
            f"- Caveat: distances exist only within contigs; fragmented assemblies bias them downward (Q-n50 applied; log N50 is a covariate). Table: `{tp}`.")


def v18(ctx):
    if not ctx.have('A2/genome_gene_stats.parquet'):
        return '### V18 Gene density versus genome size\n\nNot built: A2 (job J1) missing.'
    d = _a2reps(ctx, dup=False, n50=False)
    d['genes_per_mb'] = d['n_protein_coding'] / (d['total_length_bp'] / 1e6)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
    res = {}
    for gf in GF_ORDER:
        s = d[d['growth_form'] == gf]
        axs[0].scatter(s['total_length_bp'] / 1e6, s['genes_per_mb'], s=5, alpha=0.4, color=GF_COLORS[gf], label=f'{gf} (n={len(s)})', linewidths=0, rasterized=True)
        axs[1].scatter(s['masked_pct'], s['genes_per_mb'], s=5, alpha=0.4, color=GF_COLORS[gf], linewidths=0, rasterized=True)
        res[gf] = {'size': S.spearman(s['total_length_bp'], s['genes_per_mb']), 'masked': S.spearman(s['masked_pct'], s['genes_per_mb'])}
    axs[0].set_xscale('log'); axs[0].set_xlabel('assembly size (Mb)'); axs[0].set_ylabel('protein-coding genes per Mb'); axs[0].legend(fontsize=7)
    axs[1].set_xlabel('masked (%)'); axs[1].set_ylabel('genes per Mb')
    fig.suptitle('V18: gene density vs genome size and repeat masking', fontsize=9)
    fp = save_fig(fig, 'V18', 'density')
    md = [f"### V18 Gene density versus genome size\n\n![]({fp})\n", f"- {counts_text(d)}. Filters: {ctx.filt_text(dup=False, n50=False, model=False)}."]
    for gf, x in res.items():
        md.append(f"- {gf}: Spearman rho(size, density) = {x['size']['rho']:.3f}; rho(masked %, density) = {x['masked']['rho']:.3f} (n = {x['size']['n']:,} genomes).")
    md.append("- Caveat: masking method may differ between pipeline waves (not verified). Genus-level contrast of density is in V03.")
    return '\n'.join(md)


# ------------------------------------------------------------------ H4
def _secretome(ctx):
    s = ctx.tbl('A5/secretome.parquet'); g = ctx.tbl('A4/genome_protein_summary.parquet') if ctx.have('A4/genome_protein_summary.parquet') else ctx.tbl('A2/protein_stats.parquet')
    d = ctx.reps().merge(s, left_on='LOCUSTAG', right_on='genome').merge(g[['genome', 'n_proteins']], on='genome')
    d['secretome_frac'] = d['n_secreted_strict'] / d['n_proteins']
    d['signalp_frac'] = d['n_signalp_09'] / d['n_proteins']
    return d


def v20(ctx):
    if not ctx.have('A5/secretome.parquet') or not (ctx.have('A4/genome_protein_summary.parquet') or ctx.have('A2/protein_stats.parquet')):
        return '### V20 Secretome size per genome\n\nNot built: A5 (J4) or protein counts (J1/J3) missing.'
    d = _secretome(ctx); d = d[d['growth_form'].isin(GF_ORDER)]
    fp = strip_by_class(d, 'n_secreted_strict', 'strict secretome size (SignalP >= 0.9, <= 1 TM, no GPI)', 'V20', 'secretome_by_class', unit='genome')
    fp2 = strip_by_class(d, 'secretome_frac', 'strict secretome fraction of proteome', 'V20', 'secretome_frac_by_class', unit='genome')
    gd = genus_medians(d, ['n_secreted_strict', 'secretome_frac', 'signalp_frac', 'n_proteins'])
    cb, r = contrast_block(ctx, gd, 'n_secreted_strict', 'strict secretome size')
    cb2, r2 = contrast_block(ctx, gd, 'secretome_frac', 'secretome fraction')
    tp = save_tab(gd, 'V20', 'genus_medians')
    notm = d['n_signalp_no_tmhmm_row'].sum()
    return (f"### V20 Secretome size per genome\n\n![]({fp})\n\n![]({fp2})\n\n- {counts_text(d)}. Filters: {ctx.filt_text(model=False)}.\n"
            f"**Count:**\n{cb}\n**Fraction of proteome:**\n{cb2}\n- SignalP proteins without a TMHMM row (the TMHMM table holds only proteins with >= 1 helix, so these have 0 helices): {int(notm):,} across these genomes.\n"
            f"- Caveat: SignalP cannot see truncated N-termini; partial models reduce counts. Table: `{tp}`.")


def v21(ctx):
    if not ctx.have('A5/secretome.parquet'):
        return '### V21 Secretome fraction by primary lifestyle\n\nNot built: A5 (J4) missing.'
    d = _secretome(ctx)
    gd = genus_medians(d, ['secretome_frac', 'n_secreted_strict'])
    n = gd.groupby('primary_lifestyle').size(); keep = n[n >= 10].index
    g2 = gd[gd['primary_lifestyle'].isin(keep)]
    order = g2.groupby('primary_lifestyle')['secretome_frac'].median().sort_values().index.tolist()
    fig, ax = plt.subplots(figsize=(7, 0.3 * len(order) + 1.5))
    for i, ls in enumerate(order):
        v = g2[g2['primary_lifestyle'] == ls]['secretome_frac'].values
        ax.boxplot([v], positions=[i], vert=False, widths=0.6, showfliers=False, medianprops={'color': 'black'})
        ax.scatter(v, i + (np.random.default_rng(3).random(len(v)) - 0.5) * 0.5, s=6, alpha=0.5, linewidths=0,
                   color=[GF_COLORS.get(x, '#999') for x in g2[g2['primary_lifestyle'] == ls]['growth_form']])
    ax.set_yticks(range(len(order))); ax.set_yticklabels([f'{l} (n={n[l]})' for l in order], fontsize=7); ax.set_xlabel('strict secretome fraction (genus median)')
    ax.set_title('V21: secretome fraction by FungalTraits primary lifestyle (genera with >= 10)', fontsize=9)
    fp = save_fig(fig, 'V21', 'secretome_by_lifestyle')
    kw = S.kruskal({l: g2[g2['primary_lifestyle'] == l]['secretome_frac'].values for l in order})
    tb = g2.groupby('primary_lifestyle').agg(n_genera=('GENUS', 'size'), median_secretome_frac=('secretome_frac', 'median'),
                                             median_n_secreted=('n_secreted_strict', 'median'),
                                             dominant_growth_form=('growth_form', lambda s: s.mode().iloc[0] if s.notna().any() else None)).reset_index().sort_values('median_secretome_frac')
    tp = save_tab(tb, 'V21', 'by_lifestyle')
    mm_txt = ''
    if not ctx.skip_mixed:
        g3 = g2.copy(); g3['lifestyle'] = g3['primary_lifestyle'].str.replace('[^A-Za-z0-9]', '_', regex=True)
        mm = S.mixed_model(g3, 'secretome_frac', ['C(lifestyle, Treatment("unspecified_saprotroph"))', 'log_n50', 'busco_complete'], group='CLASS', vc='ORDER_', extra_cols=['lifestyle'])
        if mm.get('available'):
            rows = [f"{t['term'].split('T.')[-1].rstrip(']')}: {t['estimate']:+.4f} (95 % CI {t['ci_lo']:+.4f} to {t['ci_hi']:+.4f}), p = {S.fmt(t['p'])}" for t in mm['terms'] if 'lifestyle' in t['term']]
            mm_txt = f"- Mixed model secretome_frac ~ lifestyle + log N50 + BUSCO + (1|CLASS/ORDER), reference = unspecified_saprotroph, n = {mm['n']:,} genera: " + '; '.join(rows) + '. PGLS pending tree.'
        else:
            mm_txt = f"- Mixed model not fitted: {mm.get('note')}."
    return (f"### V21 Secretome fraction by primary lifestyle\n\n![]({fp})\n\n- Unit = genus. {len(g2):,} genera in {len(order)} lifestyles with >= 10 genera. Filters: {ctx.filt_text(model=False)}.\n"
            f"- Kruskal-Wallis (screening only): H = {kw['kw_H']:.1f}, p = {S.fmt(kw['kw_p'])}.\n{mm_txt}\n\n" + md_table(tb) +
            f"\n\n- Caveats: lifestyle is a genus trait (FungalTraits 1.2, citation to verify); `nectar/tap_saprotroph` is mostly *Saccharomyces* and confounded with growth form; `plant_pathogen` mixes trophic modes. Table: `{tp}`.")


def v22(ctx):
    if not ctx.have('A5/cazy_class_counts.parquet'):
        return '### V22 CAZy class profile by lifestyle\n\nNot built: A5 (J4) missing.'
    cc = ctx.tbl('A5/cazy_class_counts.parquet')
    cls = cc.groupby(['genome', 'cazy_class'])['n_proteins'].sum().unstack(fill_value=0)
    cls = cls[[c for c in ['GH', 'GT', 'PL', 'CE', 'AA', 'CBM'] if c in cls.columns]]
    g = ctx.tbl('A4/genome_protein_summary.parquet') if ctx.have('A4/genome_protein_summary.parquet') else ctx.tbl('A2/protein_stats.parquet')
    d = ctx.reps().merge(cls, left_on='LOCUSTAG', right_index=True, how='left').merge(g[['genome', 'n_proteins']], left_on='LOCUSTAG', right_on='genome')
    for c in cls.columns:
        d[c] = d[c].fillna(0) * 1000 / d['n_proteins']
    gd = genus_medians(d, list(cls.columns))
    n = gd.groupby('primary_lifestyle').size(); keep = n[n >= 10].index
    m = gd[gd['primary_lifestyle'].isin(keep)].groupby('primary_lifestyle')[list(cls.columns)].median()
    m = m.loc[m['GH'].sort_values().index]
    fig, ax = plt.subplots(figsize=(6, 0.3 * len(m) + 1.5))
    im = ax.imshow(np.log2(m.values + 1), aspect='auto', cmap='viridis')
    ax.set_xticks(range(len(m.columns))); ax.set_xticklabels(m.columns); ax.set_yticks(range(len(m))); ax.set_yticklabels([f'{i} (n={n[i]})' for i in m.index], fontsize=7)
    for i in range(len(m)):
        for j in range(len(m.columns)):
            ax.text(j, i, f'{m.values[i, j]:.1f}', ha='center', va='center', fontsize=6, color='white')
    plt.colorbar(im, ax=ax, label='log2(median genus count per 1,000 proteins + 1)')
    ax.set_title('V22: CAZy class per 1,000 proteins by lifestyle (genus medians)', fontsize=9)
    fp = save_fig(fig, 'V22', 'cazy_class_by_lifestyle')
    tp = save_tab(m.reset_index(), 'V22', 'cazy_class_by_lifestyle')
    gf = gd[gd['growth_form'].isin(GF_ORDER)].groupby('growth_form')[list(cls.columns)].median()
    kw = {c: S.kruskal({l: gd[gd['primary_lifestyle'] == l][c].values for l in keep}) for c in cls.columns}
    return (f"### V22 CAZy class profile by lifestyle\n\n![]({fp})\n\n- Unit = genus median of per-genome counts per 1,000 proteins (dbCAN `toolcount >= 2`). "
            f"{len(gd):,} genera; lifestyles with >= 10 genera shown. Filters: {ctx.filt_text(model=False)}.\n"
            f"- Median per 1,000 proteins by growth form:\n\n" + md_table(gf.reset_index()) +
            "\n\n- Kruskal-Wallis across lifestyles (screening): " + '; '.join(f"{c} H = {v['kw_H']:.0f}, p = {S.fmt(v['kw_p'])}" for c, v in kw.items()) +
            f".\n- Caveat: a protein with several family tokens counts once per class; shared ancestry inflates significance. Table: `{tp}`.")


def v23(ctx):
    if not ctx.have('A5/cazy_substrate_counts.parquet'):
        return '### V23 Plant-cell-wall CAZymes versus plant association\n\nNot built: A5 (J4) missing.'
    sc = ctx.tbl('A5/cazy_substrate_counts.parquet')
    pcw = {'cellulose', 'xylan', 'pectin', 'lignin', 'hemicellulose', 'xyloglucan', 'beta-mannan', 'arabinan', 'arabinogalactan protein', 'beta-galactan',
           'alpha-galactan', 'cutin', 'suberin'}
    sub = sc[sc['substrate'].isin(pcw)].groupby('genome')['n_proteins'].sum().rename('n_pcw')
    g = ctx.tbl('A4/genome_protein_summary.parquet') if ctx.have('A4/genome_protein_summary.parquet') else ctx.tbl('A2/protein_stats.parquet')
    d = ctx.reps().merge(sub, left_on='LOCUSTAG', right_index=True, how='left').merge(g[['genome', 'n_proteins']], left_on='LOCUSTAG', right_on='genome')
    d['n_pcw'] = d['n_pcw'].fillna(0); d['pcw_per_1000'] = d['n_pcw'] * 1000 / d['n_proteins']
    gd = genus_medians(d, ['pcw_per_1000', 'n_pcw'])
    groups = {}
    for lab, sel in [('plant_pathogen', gd['primary_lifestyle'] == 'plant_pathogen'), ('wood_saprotroph', gd['primary_lifestyle'] == 'wood_saprotroph'),
                     ('litter_saprotroph', gd['primary_lifestyle'] == 'litter_saprotroph'), ('soil_saprotroph', gd['primary_lifestyle'] == 'soil_saprotroph'),
                     ('ectomycorrhizal', gd['primary_lifestyle'] == 'ectomycorrhizal'), ('animal_parasite', gd['primary_lifestyle'] == 'animal_parasite'),
                     ('yeast (growth form)', gd['growth_form'] == 'yeast'), ('filamentous (growth form)', gd['growth_form'] == 'filamentous_mycelium')]:
        groups[lab] = gd[sel]['pcw_per_1000'].values
    fig, ax = plt.subplots(figsize=(7, 4))
    labs = list(groups)
    for i, l in enumerate(labs):
        v = groups[l]; lo, hi = S.boot_ci(v)
        ax.scatter(v, i + (np.random.default_rng(4).random(len(v)) - 0.5) * 0.5, s=6, alpha=0.4, linewidths=0, color='#0072B2')
        ax.plot([lo, hi], [i, i], color='black', lw=2); ax.plot([np.median(v)], [i], 'k|', ms=12)
    ax.set_yticks(range(len(labs))); ax.set_yticklabels([f'{l} (n={len(groups[l])})' for l in labs], fontsize=8)
    ax.set_xlabel('plant-cell-wall CAZymes per 1,000 proteins (genus median)'); ax.set_title('V23: PCW-active CAZymes by lifestyle / growth form', fontsize=9)
    fp = save_fig(fig, 'V23', 'pcw_by_group')
    tb = pd.DataFrame([{'group': l, 'n_genera': len(v), 'median_pcw_per_1000': float(np.median(v)) if len(v) else np.nan, 'ci_lo': S.boot_ci(v)[0], 'ci_hi': S.boot_ci(v)[1]} for l, v in groups.items()])
    tp = save_tab(tb, 'V23', 'pcw_by_group')
    r = S.two_group(groups['plant_pathogen'], groups['animal_parasite'], 'plant_pathogen', 'animal_parasite')
    r2 = S.two_group(groups['yeast (growth form)'], groups['filamentous (growth form)'])
    return (f"### V23 Plant-cell-wall-degrading CAZymes versus plant association\n\n![]({fp})\n\n" + md_table(tb) +
            f"\n\n- Unit = genus. Substrates counted: {', '.join(sorted(pcw))} (dbCAN family-level substrate map, `toolcount >= 2`). Filters: {ctx.filt_text(model=False)}.\n"
            f"- plant_pathogen vs animal_parasite: median difference {S.fmt(r.get('median_diff'))} per 1,000 (95 % CI {S.fmt(r.get('diff_ci_lo'))} to {S.fmt(r.get('diff_ci_hi'))}), Cliff's delta {S.fmt(r.get('cliffs_delta'))}, p = {S.fmt(r.get('mw_p'))}.\n"
            f"- yeast vs filamentous: median difference {S.fmt(r2.get('median_diff'))} (95 % CI {S.fmt(r2.get('diff_ci_lo'))} to {S.fmt(r2.get('diff_ci_hi'))}), Cliff's delta {S.fmt(r2.get('cliffs_delta'))}, p = {S.fmt(r2.get('mw_p'))}.\n"
            f"- Caveat: substrate labels come from the family, not the protein; groups overlap taxonomically (no independence). Table: `{tp}`. PGLS pending tree.")


# ------------------------------------------------------------------ H5
def v27(ctx):
    if not ctx.have('A5/tf_counts.parquet'):
        return '### V27 TF family counts versus gene count\n\nNot built: A5 (J4) missing.'
    tf = ctx.tbl('A5/tf_counts.parquet')
    g = ctx.tbl('A4/genome_protein_summary.parquet') if ctx.have('A4/genome_protein_summary.parquet') else ctx.tbl('A2/protein_stats.parquet')
    d = ctx.reps(dup=False, n50=False).merge(g[['genome', 'n_proteins']], left_on='LOCUSTAG', right_on='genome')
    d = d[d['growth_form'].isin(GF_ORDER)]
    fams = tf.groupby('tf_family')['n_proteins'].sum().sort_values(ascending=False).index.tolist()
    fams = [f for f in fams if tf[tf['tf_family'] == f]['genome'].nunique() >= 500]
    ncol = 4; nrow = int(np.ceil(len(fams) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(13, 2.8 * nrow)); axs = axs.ravel()
    rows = []
    for i, f in enumerate(fams):
        s = d.merge(tf[tf['tf_family'] == f][['genome', 'n_proteins']].rename(columns={'n_proteins': 'n_tf'}), on='genome', how='left')
        s['n_tf'] = s['n_tf'].fillna(0)
        r = S.loglog_slope(s['n_proteins'], s['n_tf'])
        rows.append({'tf_family': f, 'n_genomes': r['n'], 'slope': r.get('slope'), 'ci_lo': r.get('slope_ci_lo'), 'ci_hi': r.get('slope_ci_hi'), 'r2': r.get('r2'),
                     'median_yeast': float(s[s['growth_form'] == 'yeast']['n_tf'].median()), 'median_filamentous': float(s[s['growth_form'] == 'filamentous_mycelium']['n_tf'].median())})
        ax = axs[i]
        for gf in GF_ORDER:
            ss = s[s['growth_form'] == gf]
            ax.scatter(ss['n_proteins'], ss['n_tf'].clip(lower=0.5), s=3, alpha=0.3, color=GF_COLORS[gf], linewidths=0, rasterized=True)
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_title(f"{f[:38]}  slope {r.get('slope', np.nan):.2f}", fontsize=7)
        ax.tick_params(labelsize=6)
    for j in range(len(fams), len(axs)):
        axs[j].axis('off')
    fig.suptitle('V27: TF family size vs proteome size (representative genomes; orange yeast, blue filamentous, yellow dimorphic)', fontsize=9)
    fig.tight_layout()
    fp = save_fig(fig, 'V27', 'tf_scaling')
    tb = pd.DataFrame(rows).sort_values('slope', ascending=False); tp = save_tab(tb, 'V27', 'tf_slopes')
    return (f"### V27 TF family counts versus gene count\n\n![]({fp})\n\n- {counts_text(d)}. Filters: {ctx.filt_text(dup=False, n50=False, model=False)}. Domain gate domain_i_evalue <= 1e-5; a protein counts once per family.\n"
            f"- log-log OLS slope of family size on proteome size, with 95 % CI (slope < 1 = sub-linear; slope > 1 = super-linear):\n\n" + md_table(tb, 25) +
            f"\n\n- Caveat: genomes are not independent (no mixed model per family here; genus-level version is in V28 for Zn2Cys6). Table: `{tp}`.")


def v28(ctx):
    if not ctx.have('A5/tf_counts.parquet'):
        return '### V28 Zn2Cys6 share of the TF repertoire\n\nNot built: A5 (J4) missing.'
    tf = ctx.tbl('A5/tf_counts.parquet'); tg = ctx.tbl('A5/tf_genome.parquet')
    zn = tf[tf['tf_family'] == 'Zn(II)2Cys6 zinc cluster'][['genome', 'n_proteins']].rename(columns={'n_proteins': 'n_zn'})
    d = ctx.reps(dup=False, n50=False).merge(tg, left_on='LOCUSTAG', right_on='genome').merge(zn, on='genome', how='left')
    d['n_zn'] = d['n_zn'].fillna(0); d['zn_share'] = d['n_zn'] / d['n_tf_proteins']
    d = d[d['growth_form'].isin(GF_ORDER) | d['growth_form'].isna()]
    fp = strip_by_class(d, 'zn_share', 'Zn2Cys6 proteins / all TF proteins', 'V28', 'zn2cys6_share_by_class', unit='genome')
    gd = genus_medians(d, ['zn_share', 'n_zn', 'n_tf_proteins'])
    gd['subphylum_group'] = np.where(gd['CLASS'].isin(PEZIZO + ['Lecanoromycetes', 'Arthoniomycetes', 'Lichinomycetes', 'Coniocybemycetes']), 'Pezizomycotina',
                             np.where(gd['CLASS'].isin(SACCH), 'Saccharomycotina', np.where(gd['PHYLUM'] == 'Basidiomycota', 'Basidiomycota', 'other')))
    tb = gd.groupby('subphylum_group').agg(n_genera=('GENUS', 'size'), median_zn_share=('zn_share', 'median'), median_n_zn=('n_zn', 'median')).reset_index()
    r = S.two_group(gd[gd['subphylum_group'] == 'Pezizomycotina']['zn_share'], gd[gd['subphylum_group'] != 'Pezizomycotina']['zn_share'], 'Pezizomycotina', 'others')
    cb, _ = contrast_block(ctx, gd, 'zn_share', 'Zn2Cys6 share')
    tp = save_tab(gd, 'V28', 'genus_medians')
    return (f"### V28 Zn2Cys6 share of the TF repertoire by class\n\n![]({fp})\n\n- {counts_text(d)}. Filters: {ctx.filt_text(dup=False, n50=False, model=False)}.\n\n" + md_table(tb) +
            f"\n\n- Pezizomycotina vs all other genera: median share {S.fmt(r['median_Pezizomycotina'])} vs {S.fmt(r['median_others'])}, difference {S.fmt(r.get('median_diff'))} (95 % CI {S.fmt(r.get('diff_ci_lo'))} to {S.fmt(r.get('diff_ci_hi'))}), Cliff's delta {S.fmt(r.get('cliffs_delta'))}.\n"
            f"**Yeast vs filamentous:**\n{cb}\n- Caveat: Zn2Cys6 and the fungal-specific TF domain co-occur on one protein; proteins are counted once per family. Table: `{tp}`.")


# ------------------------------------------------------------------ H6
CODE1 = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
         'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
         'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
         'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
         'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
         'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}


def _rscu(ctx):
    if 'rscu' in ctx.cache:
        return ctx.cache['rscu']
    cf = ctx.q(f"SELECT species_prefix AS genome, upper(codon) AS codon, frequency FROM '{T}/codon_freq.parquet'")
    cf['aa'] = cf['codon'].map(CODE1)
    cf = cf[cf['aa'].notna() & (cf['aa'] != '*') & ~cf['aa'].isin(['M', 'W'])]
    cf['rscu'] = cf['frequency'] / cf.groupby(['genome', 'aa'])['frequency'].transform('mean')
    R = cf.pivot(index='genome', columns='codon', values='rscu').fillna(0)
    gc3 = ctx.q(f"SELECT species_prefix AS genome, sum(frequency) FILTER (WHERE upper(codon) LIKE '%G' OR upper(codon) LIKE '%C') AS gc3 FROM '{T}/codon_freq.parquet' GROUP BY 1").set_index('genome')
    ctx.cache['rscu'] = (R, gc3)
    return R, gc3


def _genus_rscu(ctx):
    R, gc3 = _rscu(ctx)
    L = ctx.reps(dup=False, n50=False)
    d = L.merge(R, left_on='LOCUSTAG', right_index=True)
    gd = genus_medians(d, list(R.columns))
    return gd, list(R.columns)


def v33(ctx):
    gd, cols = _genus_rscu(ctx)
    gd = gd.dropna(subset=['CLASS'])
    X = gd[cols].values; X = X - X.mean(0)
    U, s, Vt = np.linalg.svd(X, full_matrices=False); pcs = U * s; var = s ** 2 / (s ** 2).sum()
    gd['PC1'] = pcs[:, 0]; gd['PC2'] = pcs[:, 1]
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    top = gd['CLASS'].value_counts(); top = top[top >= 15].index.tolist()
    cmap = plt.get_cmap('tab20')
    for i, c in enumerate(top):
        s_ = gd[gd['CLASS'] == c]; axs[0].scatter(s_['PC1'], s_['PC2'], s=7, alpha=0.6, color=cmap(i % 20), label=f'{c} ({len(s_)})', linewidths=0, rasterized=True)
    axs[0].legend(fontsize=5, ncol=2); axs[0].set_xlabel(f'PC1 ({100*var[0]:.1f} %)'); axs[0].set_ylabel(f'PC2 ({100*var[1]:.1f} %)'); axs[0].set_title('by class (>= 15 genera)', fontsize=8)
    for gf in GF_ORDER:
        s_ = gd[gd['growth_form'] == gf]; axs[1].scatter(s_['PC1'], s_['PC2'], s=7, alpha=0.5, color=GF_COLORS[gf], label=f'{gf} ({len(s_)})', linewidths=0, rasterized=True)
    axs[1].legend(fontsize=7); axs[1].set_xlabel(f'PC1 ({100*var[0]:.1f} %)'); axs[1].set_title('by growth form', fontsize=8)
    fig.suptitle('V33: PCA of genus-median RSCU (59 codons)', fontsize=9)
    fp = save_fig(fig, 'V33', 'rscu_pca')
    # PERMANOVA on Euclidean distances, classes with >= 5 genera
    sub = gd[gd['CLASS'].isin(gd['CLASS'].value_counts()[lambda v: v >= 5].index)].copy()
    sub = sub[sub['growth_form'].isin(GF_ORDER)]
    from scipy.spatial.distance import pdist, squareform
    D = squareform(pdist(sub[cols].values))
    p_cls = S.permanova(D, sub['CLASS'].values, n_perm=499)
    p_gf = S.permanova(D, sub['growth_form'].values, n_perm=499)
    p_gf_in = S.permanova(D, sub['growth_form'].values, n_perm=499, strata=sub['CLASS'].values)
    tb = pd.DataFrame([{'model': 'CLASS', **p_cls}, {'model': 'growth_form', **p_gf}, {'model': 'growth_form | CLASS (permutations within class)', **p_gf_in}])
    tp = save_tab(tb, 'V33', 'permanova'); tp2 = save_tab(gd[['GENUS', 'PHYLUM', 'CLASS', 'growth_form', 'PC1', 'PC2']], 'V33', 'genus_pcs')
    return (f"### V33 Codon usage ordination\n\n![]({fp})\n\n- Unit = genus (median RSCU over representative genomes; {len(gd):,} genera). RSCU computed with the standard code for all genomes "
            f"(CTG-clade code-12 genomes keep CTG in the Leu family; see caveat). Filters: {ctx.filt_text(dup=False, n50=False, model=False)}.\n"
            f"- PC1 explains {100*var[0]:.1f} %, PC2 {100*var[1]:.1f} % of RSCU variance.\n- PERMANOVA (Euclidean, {p_cls['n']:,} genera in classes with >= 5 genera, 499 permutations):\n\n" + md_table(tb) +
            f"\n\n- Caveat: whole-proteome codon frequencies, not expression-weighted; GC content drives much of PC1 (see V34). Tables: `{tp}`, `{tp2}`.")


def v34(ctx):
    R, gc3 = _rscu(ctx)
    d = ctx.reps(dup=False, n50=False).merge(gc3, left_on='LOCUSTAG', right_index=True)
    d = d[d['growth_form'].isin(GF_ORDER)]
    fig, ax = plt.subplots(figsize=(6, 5)); rows = []
    for gf in GF_ORDER:
        s = d[d['growth_form'] == gf]
        ax.scatter(s['gc_pct'], 100 * s['gc3'], s=5, alpha=0.4, color=GF_COLORS[gf], label=f'{gf} (n={len(s)})', linewidths=0, rasterized=True)
        lr = __import__('scipy').stats.linregress(s['gc_pct'], 100 * s['gc3'])
        rows.append({'growth_form': gf, 'n_genomes': len(s), 'slope_gc3_on_gc': lr.slope, 'intercept': lr.intercept, 'r2': lr.rvalue ** 2,
                     'spearman_rho': S.spearman(s['gc_pct'], s['gc3'])['rho']})
    ax.plot([25, 70], [25, 70], color='grey', lw=0.8, ls='--'); ax.set_xlabel('assembly GC (%)'); ax.set_ylabel('GC3 (%)'); ax.legend(fontsize=7)
    ax.set_title('V34: GC3 vs genomic GC (representative genomes)', fontsize=9)
    fp = save_fig(fig, 'V34', 'gc3_vs_gc')
    tb = pd.DataFrame(rows); tp = save_tab(tb, 'V34', 'gc3_slopes')
    d['gc3_resid'] = 100 * d['gc3'] - (rows[2]['intercept'] + rows[2]['slope_gc3_on_gc'] * d['gc_pct'])
    gd = genus_medians(d, ['gc3_resid', 'gc3', 'gc_pct'])
    cb, _ = contrast_block(ctx, gd, 'gc3_resid', 'GC3 residual from the filamentous line (%)', covars=())
    return (f"### V34 GC3 versus genomic GC\n\n![]({fp})\n\n- {counts_text(d)}. Filters: {ctx.filt_text(dup=False, n50=False, model=False)}.\n\n" + md_table(tb) +
            f"\n\n**Residual GC3 (yeast vs filamentous):**\n{cb}\n- Caveat: assembly GC includes repeats and masked regions. Table: `{tp}`.")


def v35(ctx):
    gd, cols = _genus_rscu(ctx)
    gd = gd.dropna(subset=['CLASS', 'ORDER_', 'PHYLUM']); gd = gd[gd['growth_form'].isin(['yeast', 'filamentous_mycelium'])].reset_index(drop=True)
    from scipy.spatial.distance import pdist, squareform
    D = squareform(pdist(gd[cols].values)); n = len(gd)
    iu = np.triu_indices(n, 1)
    fam = gd['FAMILY'].values if 'FAMILY' in gd else None
    lab = ctx.labels[ctx.labels['is_representative']].groupby('GENUS')['FAMILY'].first()
    gd['FAMILY'] = gd['GENUS'].map(lab)
    def rank_of(i, j):
        if gd['FAMILY'].values[i] == gd['FAMILY'].values[j] and gd['FAMILY'].values[i] is not None:
            return 'same family'
        if gd['ORDER_'].values[i] == gd['ORDER_'].values[j]:
            return 'same order'
        if gd['CLASS'].values[i] == gd['CLASS'].values[j]:
            return 'same class'
        if gd['PHYLUM'].values[i] == gd['PHYLUM'].values[j]:
            return 'same phylum'
        return 'different phylum'
    F = gd['FAMILY'].values; O = gd['ORDER_'].values; C = gd['CLASS'].values; P = gd['PHYLUM'].values; G = gd['growth_form'].values
    i, j = iu
    rank = np.where((F[i] == F[j]) & pd.notna(F[i]), 'same family', np.where(O[i] == O[j], 'same order', np.where(C[i] == C[j], 'same class', np.where(P[i] == P[j], 'same phylum', 'different phylum'))))
    same_gf = G[i] == G[j]
    df = pd.DataFrame({'dist': D[iu], 'rank': rank, 'same_growth_form': same_gf})
    order = ['same family', 'same order', 'same class', 'same phylum', 'different phylum']
    tb = df.groupby(['rank', 'same_growth_form'])['dist'].agg(['median', 'size']).reset_index()
    tb['rank'] = pd.Categorical(tb['rank'], order); tb = tb.sort_values(['rank', 'same_growth_form'])
    fig, ax = plt.subplots(figsize=(7, 4))
    for k, (sg, colr) in enumerate([(True, '#0072B2'), (False, '#D55E00')]):
        vals = [df[(df['rank'] == r) & (df['same_growth_form'] == sg)]['dist'].values for r in order]
        bp = ax.boxplot(vals, positions=np.arange(len(order)) + (k - 0.5) * 0.35, widths=0.3, showfliers=False, patch_artist=True)
        for b in bp['boxes']:
            b.set(facecolor=colr, alpha=0.6)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, fontsize=8); ax.set_ylabel('Euclidean RSCU distance between genera')
    ax.scatter([], [], color='#0072B2', label='same growth form'); ax.scatter([], [], color='#D55E00', label='different growth form'); ax.legend(fontsize=8)
    ax.set_title('V35: codon-usage distance by lowest shared taxonomic rank', fontsize=9)
    fp = save_fig(fig, 'V35', 'distance_by_rank')
    tp = save_tab(tb, 'V35', 'distance_by_rank')
    return (f"### V35 Codon-usage distance versus taxonomic rank\n\n![]({fp})\n\n- Unit = genus pair ({len(df):,} pairs among {n:,} yeast or filamentous genera). No p-value: pairs are not independent.\n\n" + md_table(tb) +
            f"\n\n- Reading: if distance rises with rank at both colors, taxonomy structures codon usage; the gap between colors within a rank is the growth-form signal at that rank. Table: `{tp}`.")


# ------------------------------------------------------------------ pages
PAGES = {
    'H1_quality': ('H1. Annotation and assembly quality differ by clade', 'Expected: yeast classes have lower fungi_odb12 BUSCO; fragmented assemblies have more partial models.', ['v01', 'v02', 'v03', 'v04', 'v05']),
    'H2_protein_length': ('H2. Among orthologous proteins, yeast proteins are shorter', 'PI hypothesis: yeast shorter. Proxies for orthology: SwissProt best-hit families (P-A); model-anchor (P-B) and Pfam architecture (P-C) are P2.', ['v07', 'v08', 'v09', 'v12', 'v14']),
    'H3_gene_structure': ('H3. Yeast genomes have fewer introns and denser, shorter genes', 'Expected: yeast fewer introns per gene, shorter intergenic distances, higher gene density.', ['v15', 'v17', 'v18']),
    'H4_secretome_lifestyle': ('H4. Secretome and CAZy content track lifestyle', 'Expected: plant-associated and wood/litter saprotrophs > yeasts and animal parasites.', ['v20', 'v21', 'v22', 'v23']),
    'H5_domains_tf': ('H5. TF and domain repertoires scale sub-linearly; Zn2Cys6 expansion is Pezizomycotina-specific', 'Expected: slopes < 1 for most TF families; highest Zn2Cys6 share in Pezizomycotina.', ['v27', 'v28']),
    'H6_codon_usage': ('H6. Codon usage tracks taxonomy more than growth form; GC3 follows genomic GC', 'Expected: CLASS explains more RSCU variance than growth form.', ['v33', 'v34', 'v35']),
}

VERDICTS = {
    'H1_quality': "**Verdict: supported.** BUSCO (fungi_odb12) completeness depends on class: median 75-84 % in Saccharomycetes, Pichiomycetes and Schizosaccharomycetes against 98-99 % in the large filamentous classes (V02). The partial-model fraction also differs by class (median 0.25 in Saccharomycetes, 0.65 in Glomeromycetes) and is only weakly related to N50 across all genomes (Spearman rho = -0.07, V04). The top 20 species hold 28.7 % of all genomes (V06). Every later view therefore uses one genome per species, class-relative BUSCO thresholds and the complete-model filter. Next data needed: an RNA-seq-in-training flag per genome (D8).",
    'H2_protein_length': "**Verdict: a small, consistent shortening exists; it is 1-2 % and strongest in Saccharomycotina.** Three orthology proxies agree in direction. P-A (shared SwissProt best hit, V08, 3,652 families): median yeast/filamentous length ratio 0.994 (log2 -0.009, 95 % CI -0.010 to -0.007), 59.9 % of families yeast-shorter. P-B (every genome anchored to one reference proteome, V10): ratio 0.985-0.987 (log2 -0.019 to -0.022) for all four anchors (*S. cerevisiae*, *N. crassa*, *S. pombe*, *A. nidulans*), 69-71 % of families yeast-shorter; the *S. cerevisiae*- and *N. crassa*-anchored medians are identical, so the reference-dominance bias of P-A does not create the signal. P-C (Pfam domain architecture, V11, 5,026 architectures): ratio 0.980 (log2 -0.029, 95 % CI -0.032 to -0.025), 63.9 % yeast-shorter, similar at every domain count. The effect is small against the spread between families (10th/90th percentiles about -0.12/+0.04 for P-B) and is stable across the quality-filter grid (V12). By origin (V09) Saccharomycotina genera are 1.6 % shorter and the two Taphrinomycotina genera 3 % shorter; Tremellomycetes, Pucciniomycotina and Ustilaginomycotina yeasts are 0.2-1.3 % longer than the filamentous reference, so the shortening is a Saccharomycotina trait rather than a yeast trait. Raw median protein length is slightly longer in yeast genera (414 vs 401 aa, V07) because yeast proteomes have more multi-domain proteins (V31: 40 % vs 33 %) and lack many short secreted families. The mixed model with CLASS/ORDER random intercepts cannot separate growth form from class; the within-family, multiple-origin and sister-pair designs carry the comparison. PGLS pending tree.",
    'H3_gene_structure': "**Verdict: supported.** Yeast genera have a median 0.24 introns per transcript against 2.03 in filamentous genera (Cliff's delta -0.68), 81 % intronless transcripts against 22 % (V15), transcript span / CDS ratio 1.00 against 1.14 (V16), median intergenic distance 406 bp against 843 bp (V17), and 427 against 301 protein-coding genes per Mb (V03). Median CDS length is not shorter in yeasts (1,242 vs 1,194 bp, V16). Dimorphic yeasts sit between the two groups on every measure. The mixed models with CLASS/ORDER random intercepts return yeast effects with CIs spanning zero: growth form is confounded with class, so those models cannot separate the two; the genus-median contrasts above are the stated result. Caveats: intron counts depend on per-genome gene prediction and RNA-seq availability; intergenic distances are truncated by contig ends (Q-n50 applied).",
    'H4_secretome_lifestyle': "**Verdict: supported in direction; lifestyle and clade are confounded.** The strict secretome is 2.6 % of the proteome in yeast genera against 5.9 % in filamentous genera (V20); 15 % of yeast CAZymes carry a signal peptide against 42 % (V24); secreted peptidases per 1,000 proteins are 2.6 against 5.6 while total MEROPS hits are similar (32 vs 28, V25). Plant pathogens and dung/litter/wood saprotrophs have the largest secretome fractions (6.5-6.9 % and 5.5-5.7 %), ectomycorrhizal genera and the yeast-dominated nectar/tap saprotrophs the smallest (2.4-2.7 %) (V21). Plant-cell-wall CAZymes per 1,000 proteins: plant pathogens 11.3, wood saprotrophs 9.5, animal parasites 2.3, yeasts 1.5 (V23). PL and CE classes are near zero in yeasts (V22). Caveats: FungalTraits labels are genus traits; several lifestyles map to one clade, so the Kruskal-Wallis values overstate independence; the mixed model with CLASS/ORDER random intercepts still gives a yeast secretome effect below zero (count -270, 95 % CI -404 to -136; fraction -0.013, CI -0.023 to -0.004), so the yeast deficit is not explained by class membership alone. Across lifestyles the same model leaves only ectomycorrhizal genera different from unspecified saprotrophs (-0.021, 95 % CI -0.027 to -0.015); the other lifestyle contrasts in V21 are absorbed by class, so most of the lifestyle signal is clade signal; PGLS pending tree.",
    'H5_domains_tf': "**Verdict: partly supported.** 17 of 19 TF families scale sub-linearly with proteome size (slope < 1, V27); the fungal-specific TF domain scales super-linearly (slope 1.24, 95 % CI 1.19-1.29) and Zn2Cys6 has slope 0.73. Pfam family richness scales with slope 0.09 (yeast) and 0.30 (filamentous) and yeast genera carry 4,062 against 4,755 distinct families (V29). The Zn2Cys6 share of the TF repertoire is 0.40 in Pezizomycotina against 0.18 in Basidiomycota, but Saccharomycotina reach 0.39 (V28), so the expansion is an Ascomycota feature, not a Pezizomycotina-specific one. 3,582 of 7,907 Pfam families differ in genus presence between yeasts and filamentous fungi at BH q < 0.01 (V30); the most yeast-specific are Saccharomycotina-only entries (POC1, STD1, MPM1, Gsf2), the most filamentous-specific are cerato-platanin, GH6 cellulases and AAA_lid_6. Caveat: slopes are fitted on genomes, not genera (clade non-independence not modelled in V27/V29 slopes).",
    'H6_codon_usage': "**Verdict: supported.** CLASS explains 57 % of the variance in genus-median RSCU (PERMANOVA R2 = 0.567) against 14 % for growth form, and growth form has no effect once permutations are restricted within class (p = 0.86, V33). GC3 follows genomic GC (Spearman rho 0.81 filamentous, 0.96 yeast, V34); at equal genomic GC, yeast genera have GC3 3.2 points lower than the filamentous regression line. Codon-usage distance rises with taxonomic rank; within a phylum, pairs with different growth forms are about twice as distant as pairs with the same growth form (V35), which reflects the Saccharomycotina/Pezizomycotina split. Amino-acid composition PC1 (75 % of variance) tracks genomic GC (rho = -0.86, V38). tRNA gene counts correlate with codon frequency more strongly in yeasts (median rho 0.39 vs 0.23, V36). Code-12 Pichiomycetes use CTG at 5.2 per 1,000 codons against 9.2 in code-1 Pichiomycetes and 17.7 across fungi (V37). Caveat: whole-proteome codon frequencies, not expression-weighted; code-12 CTG kept in the Leu family for RSCU.",
    'H7_paralogs': "**Verdict: not supported.** Young paralog nodes (yn00 dS 0.05-1.0) per 1,000 proteins are 10.6 in yeast genera and 9.2 in filamentous genera (difference +1.4, 95 % CI -0.7 to +3.9, Cliff's delta -0.03, p = 0.58; V39, 987 genera, hybrids and BUSCO-duplicated > 5 % excluded). Identical-CDS gene content correlates with BUSCO duplication (rho = 0.37, V40), consistent with unpurged haplotigs, but TE-derived models also contribute. Caveats: interim wgd.ks release (7,384 genomes, subset rule undocumented); node sampling is noisy per node; PGLS pending tree.",
}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--only', default=''); ap.add_argument('--skip-mixed', action='store_true'); ap.add_argument('--no-index', action='store_true')
    a = ap.parse_args()
    ctx = Ctx(skip_mixed=a.skip_mixed)
    import render_p2
    globals().update(render_p2.VIEWS)
    for pg, vs in render_p2.PAGE_ADD.items():
        t, e, cur = PAGES[pg]; PAGES[pg] = (t, e, sorted(set(cur) | set(vs), key=lambda x: int(x[1:])))
    PAGES.update(render_p2.PAGE_NEW)
    only = set(x.strip().lower() for x in a.only.split(',') if x.strip())
    status = {}
    for page, (title, expect, views) in PAGES.items():
        blocks = []
        for v in views:
            if only and v not in only:
                continue
            fn = globals()[v]
            try:
                print('render', v, flush=True)
                md = fn(ctx); status[v] = 'ok' if 'Not built' not in md else 'not built'
            except Exception as e:
                md = f"### {v.upper()}\n\nRender failed: {e}\n```\n{traceback.format_exc()[-1500:]}\n```"; status[v] = 'failed'
                print('FAILED', v, e, flush=True)
            blocks.append(md)
        if not blocks:
            continue
        pth = f'{DOCS}/{page}.md'
        # merge with existing page blocks when --only is used
        existing = {}
        if only and os.path.exists(pth):
            cur = open(pth).read().split('\n## Verdict')[0].split('\n### ')
            for b in cur[1:]:
                existing[b.split(' ')[0]] = '### ' + b
        new = {b.split('\n')[0].split(' ')[1]: b for b in blocks}
        existing.update(new)
        ordered = [existing[v.upper()] for v in views if v.upper() in existing]
        with open(pth, 'w') as o:
            o.write(f"# {title}\n\n{expect}\n\n[Back to index](index.md)\n\n" + '\n\n'.join(ordered) + '\n\n## Verdict\n\n' + VERDICTS.get(page, '') + '\n')
    if not a.no_index:
        write_index(status)
    print(json.dumps(status, indent=1))


def write_index(status):
    L = ctx_labels = duckdb.connect().execute(f"SELECT count(*), count(*) FILTER (WHERE is_representative), count(*) FILTER (WHERE is_genus_rep) FROM '{A}/labels/genome_labels.parquet'").fetchone()
    lines = ['# BFD summary atlas', '',
             'Static summary of the BFD fungal genome database, grouped by hypothesis. Built from the parquet snapshot in '
             '`/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables` (species table built 2026-11-06) through the intermediates in `BFD_dashboard/atlas_data/`. '
             'Design and open decisions: [PLAN.md](PLAN.md).', '',
             f'Genomes: {L[0]:,}. Species representatives (one per species, best BUSCO then N50): {L[1]:,}. Genus representatives: {L[2]:,}.', '',
             '## How to read the statistics', '',
             '- The unit of replication is stated on every view. Most cross-clade tests use genus medians of representative genomes, because strains of one species and species of one genus are not independent.',
             '- Quality filters: Q-ann (annotated), Q-rep (one genome per species), Q-busco-rel (BUSCO complete >= the 25th percentile of its class), Q-dup (BUSCO duplicated <= 10 %), Q-n50 (N50 >= 50 kb), Q-model (complete gene models: not partial, start and stop codon present).',
             '- Effect sizes come first (median differences with bootstrap 95 % CI, Cliff\'s delta, log2 ratios). Kruskal-Wallis is a screening statistic only. Mixed models use random intercepts for CLASS and ORDER within CLASS. PGLS is marked "pending tree" until the genus-level phylogeny finishes.',
             '- Growth form and lifestyle labels are genus-level traits from FungalTraits 1.2 plus an override file for genera it lacks; "dimorphic_yeast" is shown as a separate group and never pooled.', '',
             '## Hypotheses', '', '| Page | Hypothesis | Views |', '|---|---|---|']
    for page, (title, expect, views) in PAGES.items():
        st = ', '.join(f"{v.upper()} ({status.get(v, 'pending')})" for v in views)
        lines.append(f'| [{page}]({page}.md) | {title} | {st} |')
    lines += ['| H8 (P3) | Telomere repeat monomer is conserved within classes | V41, V42 (P3, not built) |', '',
              '## Data coverage', '', 'See V01 on the H1 page. Tables and figures: `tables/`, `figures/`.', '']
    with open(f'{DOCS}/index.md', 'w') as o:
        o.write('\n'.join(lines))


if __name__ == '__main__':
    main()
