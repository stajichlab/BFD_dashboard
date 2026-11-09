#!/usr/bin/env python3
"""Render the size-transition hypothesis pages (HS_G1, HS_G3, HS_P2, HS_S1, HS_C1, HS_C3) under docs/summary_atlas/
from atlas_data/A11 (job J9) plus the small A2/A4/A5/A7 intermediates and labels/genome_labels_ext.parquet (hs_labels.py).

Run under SLURM with atlas_env/venv (statsmodels). Usage: render_hs.py [--only g1,g3,...]
"""
import argparse
import os
import sys
import traceback

import duckdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stats as S
import render_views as RV

A, T, DOCS = RV.A, RV.T, RV.DOCS
save_fig, save_tab, md_table = RV.save_fig, RV.save_tab, RV.md_table
TMP = os.environ.get('SCRATCH', '/tmp') + '/atlas_render_tmp'
os.makedirs(TMP, exist_ok=True)
con = duckdb.connect(); con.execute(f"SET threads=4; SET memory_limit='20GB'; SET temp_directory='{TMP}'")
q = lambda s: con.execute(s).df()

GROUPS = ['zoosporic', 'Mucoromycota_sl', 'Zoopagomycota_sl', 'filamentous_other', 'lichen', 'Agaricomycetes', 'dimorphic', 'yeast', 'Microsporidia']
GCOL = {'zoosporic': '#009E73', 'Mucoromycota_sl': '#56B4E9', 'Zoopagomycota_sl': '#0072B2', 'filamentous_other': '#999999', 'lichen': '#CC79A7',
        'Agaricomycetes': '#E69F00', 'dimorphic': '#F0E442', 'yeast': '#D55E00', 'Microsporidia': '#000000', 'other': '#dddddd'}
PEZ = ['Eurotiomycetes', 'Dothideomycetes', 'Sordariomycetes', 'Leotiomycetes', 'Pezizomycetes', 'Orbiliomycetes', 'Xylonomycetes', 'Geoglossomycetes']
# (pair name, derived mask fn, sister mask fn) on a genus-median frame
PAIRS = [
    ('Saccharomycotina vs Pezizomycotina', lambda g: (g.yeast_origin == 'Saccharomycotina') & (g.growth_form == 'yeast'), lambda g: g.CLASS.isin(PEZ) & (g.growth_form == 'filamentous_mycelium')),
    ('Taphrinomycotina vs Pezizomycotina', lambda g: (g.yeast_origin == 'Taphrinomycotina') & (g.growth_form == 'yeast'), lambda g: g.CLASS.isin(PEZ) & (g.growth_form == 'filamentous_mycelium')),
    ('Tremellomycetes vs Agaricomycetes', lambda g: (g.yeast_origin == 'Tremellomycetes') & (g.growth_form == 'yeast'), lambda g: (g.CLASS == 'Agaricomycetes') & (g.growth_form == 'filamentous_mycelium')),
    ('Pucciniomycotina yeasts vs Pucciniomycetes', lambda g: (g.yeast_origin == 'Pucciniomycotina-yeast') & (g.growth_form == 'yeast'), lambda g: g.CLASS == 'Pucciniomycetes'),
    ('Malasseziomycetes vs Ustilaginomycetes', lambda g: (g.CLASS == 'Malasseziomycetes') & (g.growth_form == 'yeast'), lambda g: g.CLASS == 'Ustilaginomycetes'),
    ('Eurotiomycetes dimorphic vs filamentous', lambda g: (g.CLASS == 'Eurotiomycetes') & (g.growth_form == 'dimorphic_yeast'), lambda g: (g.CLASS == 'Eurotiomycetes') & (g.growth_form == 'filamentous_mycelium')),
]
YEAST_PAIRS = 5  # first five pairs are yeast origins; the sixth is the dimorphic within-class contrast
FB_ORDER = ['agaricoid', 'polyporoid', 'corticioid', 'gasteroid', 'clavarioid', 'gasteroid-hypogeous']
HY_ORDER = ['gills', 'poroid', 'hydnoid', 'smooth', 'closed']
PENDING = 'PGLS, phylogenetic signal (lambda, K) and ancestral states: pending tree (atlas_data/phylo).'


# ------------------------------------------------------------------ frames
def genome_frame():
    d = q(f"""SELECT l.*, e.fruitbody_type, e.hymenium_type, e.ecm_exploration_type, e.primary_photobiont, e.hs_group, e.septation, e.obligate_mode, e.obligate_lineage,
                     e.micro_size_ok, e.hs_keep, e.busco_ge90,
                     g.n_proteins, g.n_protein_coding, g.introns_per_transcript, g.frac_intronless, g.median_intron_len, g.median_intergenic, g.median_protein_len, g.frac_short_lt100,
                     p.frac_with_pfam, p.frac_multidomain, p.frac_with_swissprot, p.n_strict_hits,
                     i.median_idp_fraction AS idp_a5, i.frac_idp_ge50 AS idp50_a5
              FROM '{A}/labels/genome_labels.parquet' l JOIN '{A}/labels/genome_labels_ext.parquet' e USING (LOCUSTAG)
              LEFT JOIN '{A}/A2/genome_gene_stats.parquet' g ON g.genome = l.LOCUSTAG
              LEFT JOIN '{A}/A4/genome_protein_summary.parquet' p ON p.genome = l.LOCUSTAG
              LEFT JOIN '{A}/A5/idp_genome.parquet' i ON i.genome = l.LOCUSTAG
              WHERE l.is_representative""")
    d['log_n50'] = np.log10(d['N50_bp'].astype(float)); d['log_genome'] = np.log10(d['total_length_bp'].astype(float))
    d['genome_mb'] = d['total_length_bp'] / 1e6; d['genes_per_mb'] = d['n_protein_coding'] / d['genome_mb']
    d['log_nprot'] = np.log10(d['n_proteins'].astype(float))
    return d


META = dict(PHYLUM=('PHYLUM', 'first'), CLASS=('CLASS', 'first'), ORDER_=('ORDER_', 'first'), FAMILY=('FAMILY', 'first'),
            growth_form=('growth_form', lambda s: s.mode().iloc[0] if s.notna().any() else None), yeast_origin=('yeast_origin', 'first'),
            hs_group=('hs_group', lambda s: s.mode().iloc[0]), septation=('septation', 'first'), obligate_mode=('obligate_mode', 'first'), obligate_lineage=('obligate_lineage', 'first'),
            fruitbody_type=('fruitbody_type', 'first'), hymenium_type=('hymenium_type', 'first'), primary_lifestyle=('primary_lifestyle', 'first'),
            n_genomes=('LOCUSTAG', 'size'), busco_complete=('busco_complete', 'median'), log_n50=('log_n50', 'median'), log_genome=('log_genome', 'median'),
            gc_pct=('gc_pct', 'median'), masked_pct=('masked_pct', 'median'))


def genus_frame(d, cols):
    g = d.groupby('GENUS')
    cols = [c for c in dict.fromkeys(cols) if c not in META]
    out = g[cols].median().join(g.agg(**META)).reset_index()
    for c in cols + ['busco_complete', 'log_n50', 'log_genome', 'gc_pct', 'masked_pct']:
        out[c] = pd.to_numeric(out[c], errors='coerce').astype(float)  # nullable Int/Float -> float with NaN (no pd.NA)
    out['gf_yeast'] = (out['growth_form'] == 'yeast').astype(int)
    return out


def counts(d):
    return f"n = {len(d):,} genomes / {d['SPECIES'].nunique():,} species / {d['GENUS'].nunique():,} genera"


def group_counts(gd, col='hs_group'):
    return gd.groupby(col).size().to_dict()


# ------------------------------------------------------------------ generic stats/plots
def pair_table(gd, col, label):
    rows = []
    for name, fd, fs in PAIRS:
        a = gd[fd(gd)][col].dropna(); b = gd[fs(gd)][col].dropna()
        r = S.two_group(a, b, 'derived', 'sister')
        rows.append({'pair': name, 'n_derived': r['n_derived'], 'n_sister': r['n_sister'], 'median_derived': r['median_derived'], 'median_sister': r['median_sister'],
                     'diff': r.get('median_diff', np.nan), 'ci_lo': r.get('diff_ci_lo', np.nan), 'ci_hi': r.get('diff_ci_hi', np.nan), 'cliffs_delta': r.get('cliffs_delta', np.nan), 'mw_p': r.get('mw_p', np.nan)})
    tb = pd.DataFrame(rows)
    yp = tb.head(YEAST_PAIRS).dropna(subset=['diff']); st = S.sign_test(yp['diff'].values) if len(yp) else {'n': 0}
    txt = (f"- Sister pairs ({label}; unit = genus; difference = derived median minus sister median, bootstrap 95 % CI). Sign across the {len(yp)} yeast-origin pairs with data: "
           f"{st.get('n_below_zero', 0)} of {st.get('n', 0)} below zero (binomial p = {S.fmt(st.get('binom_p', np.nan))}). Pairs with < 3 genera on one side have no CI.")
    return tb, txt


def forest(tb, vid, slug, title, xlabel, est='diff', lo='ci_lo', hi='ci_hi', lab='pair', vline=0):
    fig, ax = plt.subplots(figsize=(7, 0.45 * len(tb) + 1.2))
    y = np.arange(len(tb))[::-1]
    ax.errorbar(tb[est], y, xerr=[tb[est] - tb[lo], tb[hi] - tb[est]], fmt='o', color='#0072B2', ecolor='#444', capsize=3)
    ax.axvline(vline, color='black', lw=0.8); ax.set_yticks(y); ax.set_yticklabels([f"{p} ({a}/{b})" for p, a, b in zip(tb[lab], tb.get('n_derived', ''), tb.get('n_sister', ''))], fontsize=8)
    ax.set_xlabel(xlabel); ax.set_title(f'{vid}: {title}', fontsize=9)
    return save_fig(fig, vid, slug)


def strip_groups(gd, col, xlabel, vid, slug, groups=None, log=False, col2=None, title=None):
    groups = [g for g in (groups or GROUPS) if (gd['hs_group'] == g).sum() >= 2]
    fig, ax = plt.subplots(figsize=(7, 0.42 * len(groups) + 1.2))
    for i, g in enumerate(groups):
        v = gd[gd['hs_group'] == g][col].dropna().values
        ax.boxplot([v], positions=[i], vert=False, widths=0.6, showfliers=False, medianprops={'color': 'black'}, boxprops={'color': '#444'}, whiskerprops={'color': '#444'})
        jit = (np.random.default_rng(1).random(len(v)) - 0.5) * 0.5
        ax.scatter(v, i + jit, s=7, alpha=0.5, color=GCOL.get(g, '#999'), linewidths=0, rasterized=True)
    ax.set_yticks(range(len(groups))); ax.set_yticklabels([f"{g} (n={(gd['hs_group'] == g).sum()})" for g in groups], fontsize=8); ax.set_xlabel(xlabel)
    if log:
        ax.set_xscale('log')
    ax.set_title(f'{vid}: {title or xlabel} by organism group (genus medians)', fontsize=9)
    return save_fig(fig, vid, slug)


def allometry(gd, ycol, label, vid, slug, groups=None):
    """log-log OLS of genus-median count on proteome size; returns residual column name and summary dict."""
    m = (gd[ycol] > 0) & gd['n_proteins'].notna()
    r = S.loglog_slope(gd.loc[m, 'n_proteins'], gd.loc[m, ycol])
    gd[f'res_{ycol}'] = np.log10(gd[ycol].clip(lower=0.5)) - (r['intercept'] + r['slope'] * np.log10(gd['n_proteins']))
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    for g in (groups or GROUPS):
        s = gd[(gd['hs_group'] == g) & m]
        if len(s):
            ax.scatter(s['n_proteins'], s[ycol], s=8, alpha=0.6, color=GCOL[g], label=f'{g} ({len(s)})', linewidths=0, rasterized=True)
    xx = np.logspace(np.log10(gd.loc[m, 'n_proteins'].min()), np.log10(gd.loc[m, 'n_proteins'].max()), 50)
    ax.plot(xx, 10 ** (r['intercept'] + r['slope'] * np.log10(xx)), color='black', lw=1)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('proteins per genome (genus median)'); ax.set_ylabel(label)
    ax.set_title(f"{vid}: {label}; slope {r['slope']:.2f} ({r['slope_ci_lo']:.2f}-{r['slope_ci_hi']:.2f}), n = {r['n']}", fontsize=8); ax.legend(fontsize=6)
    return save_fig(fig, vid, slug), r


def residual_by(gd, rescol, bycol, order, vid, slug, title, min_n=3):
    cats = [c for c in order if (gd[bycol] == c).sum() >= min_n]
    fig, ax = plt.subplots(figsize=(6, 0.45 * len(cats) + 1.2))
    rows = []
    for i, c in enumerate(cats):
        v = gd[gd[bycol] == c][rescol].dropna().values; lo, hi = S.boot_ci(v)
        rows.append({bycol: c, 'n_genera': len(v), 'median_residual': float(np.median(v)), 'ci_lo': lo, 'ci_hi': hi})
        ax.boxplot([v], positions=[i], vert=False, widths=0.6, showfliers=False, medianprops={'color': 'black'})
        ax.scatter(v, i + (np.random.default_rng(2).random(len(v)) - 0.5) * 0.5, s=7, alpha=0.5, color='#0072B2', linewidths=0, rasterized=True)
    ax.axvline(0, color='black', lw=0.8); ax.set_yticks(range(len(cats))); ax.set_yticklabels([f'{c} (n={n})' for c, n in zip(cats, [r['n_genera'] for r in rows])], fontsize=8)
    ax.set_xlabel('residual log10 count after proteome-size allometry'); ax.set_title(f'{vid}: {title}', fontsize=9)
    fp = save_fig(fig, vid, slug)
    tb = pd.DataFrame(rows)
    kw = S.kruskal({r[bycol]: gd[gd[bycol] == r[bycol]][rescol].dropna().values for r in rows}) if len(rows) >= 2 else {}
    return fp, tb, kw


def mixed_txt(gd, response, fixed, label):
    mm = S.mixed_model(gd, response, fixed, group='CLASS', vc='ORDER_')
    if not mm.get('available'):
        return f"- Mixed model not fitted: {mm.get('note')}. {PENDING}"
    t = [x for x in mm['terms'] if x['term'] == fixed[0]][0]
    return (f"- Mixed model `{mm['formula']}` (n = {mm['n']:,} genera, {mm['n_groups']} classes, converged = {mm['converged']}): {label} effect {S.fmt(t['estimate'])} "
            f"(95 % CI {S.fmt(t['ci_lo'])} to {S.fmt(t['ci_hi'])}), p = {S.fmt(t['p'])}. {PENDING}")


def outgroup_text(gd):
    c = group_counts(gd)
    return (f"Outgroups available in this frame (genera): zoosporic {c.get('zoosporic', 0)}, Mucoromycota s.l. {c.get('Mucoromycota_sl', 0)}, Zoopagomycota s.l. {c.get('Zoopagomycota_sl', 0)}; "
            f"Dikarya: filamentous non-Agaricomycetes {c.get('filamentous_other', 0)}, Agaricomycetes {c.get('Agaricomycetes', 0)}, lichen {c.get('lichen', 0)}, dimorphic {c.get('dimorphic', 0)}, yeast {c.get('yeast', 0)}. "
            f"Microsporidia are excluded (PI decision 2026-11-08) unless the page says otherwise.")


def page(name, title, intro, blocks, verdict):
    with open(f'{DOCS}/{name}.md', 'w') as o:
        o.write(f"# {title}\n\n{intro}\n\n[Back to index](index.md) · [Hypothesis text](HYPOTHESES_size_transitions.md)\n\n" + '\n\n'.join(blocks) + f"\n\n## Verdict\n\n{verdict}\n")
    print('wrote', name, flush=True)


# ------------------------------------------------------------------ G1
G1_SISTER = {  # lineage -> mask on genome frame for the free-living sister set
    'Microsporidia': lambda d: d.hs_group == 'zoosporic', 'Rozellomycota': lambda d: d.hs_group == 'zoosporic', 'Aphelidiomycota': lambda d: d.hs_group == 'zoosporic',
    'Olpidiomycota': lambda d: d.hs_group == 'zoosporic', 'Synchytrium': lambda d: (d.CLASS == 'Chytridiomycetes'),
    'Pneumocystis': lambda d: d.CLASS.isin(['Schizosaccharomycetes', 'Taphrinomycetes']),
    'Pucciniales': lambda d: d.CLASS.isin(['Microbotryomycetes', 'Cystobasidiomycetes', 'Agaricostilbomycetes', 'Mixiomycetes']) | ((d.CLASS == 'Pucciniomycetes') & (d.ORDER_ != 'Pucciniales')),
    'Erysiphaceae': lambda d: (d.ORDER_ == 'Helotiales') & (d.FAMILY != 'Erysiphaceae'),
    'Glomeromycota': lambda d: d.PHYLUM.isin(['Mortierellomycota', 'Mucoromycota']),
    'Entomophthorales': lambda d: d.PHYLUM.isin(['Basidiobolomycota', 'Kickxellomycota', 'Zoopagomycota']),
}
G1_TRAITS = [('log_genome', 'log10 genome size'), ('genes_per_mb', 'genes per Mb'), ('masked_pct', 'masked %'), ('introns_per_transcript', 'introns per transcript'),
             ('median_intergenic', 'median intergenic bp'), ('n_pfam_families', 'Pfam families'), ('n_proteins', 'proteins')]


def g1(d):
    # G1 filter: representatives, annotated, class-relative BUSCO OR (obligate lineage with BUSCO >= 90), Microsporidia < 10 Mb
    b = lambda c: d[c].fillna(False).astype(bool)
    keep = b('q_ann') & (b('q_busco_rel') | ((d['obligate_mode'] != 'none') & b('busco_ge90')) | (d['PHYLUM'] == 'Microsporidia')) & b('micro_size_ok')
    dd = d[keep].copy()
    fam = q(f"SELECT genome, count(*) AS n_pfam_families FROM '{A}/A7/pfam_genome_matrix.parquet' WHERE genome IN (SELECT LOCUSTAG FROM '{A}/labels/genome_labels.parquet' WHERE is_representative) GROUP BY 1")
    dd = dd.merge(fam, left_on='LOCUSTAG', right_on='genome', how='left')
    cols = ['log_genome', 'genome_mb', 'genes_per_mb', 'masked_pct', 'introns_per_transcript', 'median_intergenic', 'n_pfam_families', 'n_proteins', 'busco_duplicated']
    gd = genus_frame(dd, cols)
    gl = gd[gd['obligate_mode'] != 'none'].copy()
    # per-lineage contrast to its sister set (genus medians)
    rows = []
    for lin, grp in gl.groupby('obligate_lineage'):
        sis = gd[G1_SISTER[lin](gd) & (gd['obligate_mode'] == 'none')]
        r = {'lineage': lin, 'mode': grp['obligate_mode'].iloc[0], 'n_genera': len(grp), 'n_sister_genera': len(sis),
             'genome_mb_lineage': float(grp['genome_mb'].median()), 'genome_mb_sister': float(sis['genome_mb'].median()) if len(sis) else np.nan,
             'busco_duplicated_lineage': float(grp['busco_duplicated'].median())}
        for c, _ in G1_TRAITS:
            r[f'{c}_lineage'] = float(grp[c].median()); r[f'{c}_sister'] = float(sis[c].median()) if len(sis) else np.nan
            if c == 'log_genome':
                r['size_ratio'] = 10 ** (r[f'{c}_lineage'] - r[f'{c}_sister']) if len(sis) else np.nan
        rows.append(r)
    tb = pd.DataFrame(rows).sort_values(['mode', 'lineage'])
    tp = save_tab(tb, 'HS_G1', 'lineage_vs_sister', 'genus medians; sister sets defined in render_hs.G1_SISTER')
    # Figure 1: genome size vs gene density, all genera grey, obligate lineages coloured
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    bg = gd[gd['obligate_mode'] == 'none']; ax.scatter(bg['genome_mb'], bg['genes_per_mb'], s=5, color='#bbbbbb', alpha=0.5, linewidths=0, rasterized=True, label=f'other genera ({len(bg)})')
    mk = {'intracellular': ('o', '#0072B2'), 'extracellular_obligate': ('s', '#D55E00')}
    for mode, (m, c) in mk.items():
        s = gl[gl['obligate_mode'] == mode]; ax.scatter(s['genome_mb'], s['genes_per_mb'], s=28, marker=m, color=c, edgecolor='black', linewidths=0.4, label=f'{mode} ({len(s)} genera)')
    for _, r in tb.iterrows():
        s = gl[gl['obligate_lineage'] == r['lineage']]; ax.annotate(r['lineage'], (s['genome_mb'].median(), s['genes_per_mb'].median()), fontsize=6, xytext=(3, 3), textcoords='offset points')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('genome size (Mb, genus median)'); ax.set_ylabel('protein-coding genes per Mb'); ax.legend(fontsize=7)
    ax.set_title('HS_G1: genome size and gene density of obligate parasite lineages', fontsize=9)
    f1 = save_fig(fig, 'HS_G1', 'size_vs_density')
    # Figure 2: per-lineage log2 ratio to sister for each trait
    fig, axs = plt.subplots(1, len(G1_TRAITS), figsize=(2.3 * len(G1_TRAITS), 0.4 * len(tb) + 1.5), sharey=True)
    y = np.arange(len(tb))[::-1]
    for ax, (c, lab) in zip(axs, G1_TRAITS):
        if c == 'log_genome':
            v = (tb[f'{c}_lineage'] - tb[f'{c}_sister']) / np.log10(2)
        elif c == 'masked_pct':
            v = tb[f'{c}_lineage'] - tb[f'{c}_sister']
        else:
            v = np.log2(tb[f'{c}_lineage'].clip(lower=0.01) / tb[f'{c}_sister'].clip(lower=0.01))
        cols_ = ['#0072B2' if m == 'intracellular' else '#D55E00' for m in tb['mode']]
        ax.barh(y, v, color=cols_); ax.axvline(0, color='black', lw=0.8); ax.set_title(lab, fontsize=8)
        ax.set_xlabel('difference (points)' if c == 'masked_pct' else 'log2(lineage / sister)', fontsize=7)
    axs[0].set_yticks(y); axs[0].set_yticklabels([f"{l} ({n}/{m})" for l, n, m in zip(tb['lineage'], tb['n_genera'], tb['n_sister_genera'])], fontsize=7)
    fig.suptitle('HS_G1: obligate lineage vs free-living sister set (genus medians); blue = intracellular, orange = extracellular obligate', fontsize=9)
    f2 = save_fig(fig, 'HS_G1', 'lineage_ratios')
    # Figure 3: intron and intergenic strips by lineage + sister
    fig, axs = plt.subplots(1, 2, figsize=(10, 0.4 * len(tb) + 1.5))
    for ax, c, lab in zip(axs, ['introns_per_transcript', 'median_intergenic'], ['introns per transcript', 'median intergenic distance (bp)']):
        for i, (_, r) in enumerate(tb.iterrows()):
            s = gl[gl['obligate_lineage'] == r['lineage']][c].dropna().values; sis = gd[G1_SISTER[r['lineage']](gd) & (gd['obligate_mode'] == 'none')][c].dropna().values
            ax.scatter(sis, np.full(len(sis), len(tb) - 1 - i) + 0.15, s=6, color='#bbbbbb', linewidths=0, rasterized=True)
            ax.scatter(s, np.full(len(s), len(tb) - 1 - i) - 0.15, s=14, color='#0072B2' if r['mode'] == 'intracellular' else '#D55E00', linewidths=0)
        ax.set_yticks(y); ax.set_yticklabels(tb['lineage'], fontsize=7); ax.set_xlabel(lab)
        if c == 'median_intergenic':
            ax.set_xscale('log')
    fig.suptitle('HS_G1: lineage genera (colour) against sister genera (grey)', fontsize=9)
    f3 = save_fig(fig, 'HS_G1', 'intron_intergenic')
    # Figure 4: Pfam richness vs proteome size allometry with lineages marked
    gd2 = gd.copy(); f4, allo = allometry(gd2, 'n_pfam_families', 'distinct Pfam families per genome', 'HS_G1', 'pfam_richness_allometry')
    res = gd2.set_index('GENUS')['res_n_pfam_families']
    tb['pfam_richness_residual'] = [float(res.reindex(gl[gl['obligate_lineage'] == l]['GENUS']).median()) for l in tb['lineage']]
    # Test: intracellular vs extracellular lineages, unit = lineage
    def modes_test(t):
        a = t[t['mode'] == 'intracellular']; b = t[t['mode'] == 'extracellular_obligate']
        out = {}
        for c in ['log_genome', 'genes_per_mb', 'masked_pct', 'introns_per_transcript', 'n_pfam_families']:
            r = S.two_group(a[f'{c}_lineage'], b[f'{c}_lineage'], 'intracellular', 'extracellular')
            out[c] = r
        return out
    t_all = modes_test(tb); t_nomicro = modes_test(tb[tb['lineage'] != 'Microsporidia'])
    rows2 = []
    for c, lab in G1_TRAITS[:6]:
        if c in t_all:
            ra, rb = t_all[c], t_nomicro[c]
            rows2.append({'trait': lab, 'median_intracellular': ra['median_intracellular'], 'median_extracellular': ra['median_extracellular'], 'mw_p_all': ra.get('mw_p', np.nan),
                          'median_intracellular_noMicro': rb['median_intracellular'], 'mw_p_noMicro': rb.get('mw_p', np.nan), 'n_lineages': f"{ra['n_intracellular']}+{ra['n_extracellular']}"})
    tb2 = pd.DataFrame(rows2); tp2 = save_tab(tb2, 'HS_G1', 'mode_test', 'unit = lineage (genus-median of lineage); with and without Microsporidia')
    # gene-content loss in both directions?
    both_loss = tb[(tb['n_pfam_families_lineage'] < tb['n_pfam_families_sister'])]
    intro = ("Hypothesis G1: obligate parasitism moves genome size in two directions set by the host-cell boundary. Intracellular parasites shrink (small genomes, high gene density, few introns); "
             "haustorial or hemocoel-living extracellular obligate biotrophs expand through repeats; both lose gene content (Pfam family richness). Expected: a bimodal genome-size distribution among obligate lineages.")
    micro_n = int((gl['obligate_lineage'] == 'Microsporidia').sum())
    b1 = (f"### HS_G1a Genome size and gene density\n\n![]({f1})\n\n- {counts(dd)} after the G1 filter: representative, annotated, class-relative BUSCO or (obligate lineage with BUSCO >= 90 %), Microsporidia only < 10 Mb ({micro_n} genera). "
          f"Obligate lineages (hand list, `hs_labels.py`): intracellular = Microsporidia, Rozellomycota, Aphelidiomycota, Olpidiomycota, *Pneumocystis*, *Synchytrium*; extracellular obligate = Pucciniales, Erysiphaceae, Glomeromycota, Entomophthorales. "
          f"Unit = genus (median of representative genomes). Lineages present after the filter: {len(tb)}.\n\n" + md_table(tb[['lineage', 'mode', 'n_genera', 'n_sister_genera', 'genome_mb_lineage', 'genome_mb_sister', 'size_ratio', 'busco_duplicated_lineage', 'genes_per_mb_lineage', 'genes_per_mb_sister', 'masked_pct_lineage', 'masked_pct_sister', 'introns_per_transcript_lineage', 'introns_per_transcript_sister', 'n_pfam_families_lineage', 'n_pfam_families_sister', 'pfam_richness_residual']]) +
          f"\n\n- `size_ratio` = lineage genome size / sister genome size (genus medians). Table: `{tp}`.")
    b2 = (f"### HS_G1b Lineage versus sister, trait by trait\n\n![]({f2})\n\n- Bars are log2 ratios of genus medians (lineage / its free-living sister set); masked % is a difference in points. Sister sets: Microsporidia, Rozellomycota, Aphelidiomycota, Olpidiomycota -> zoosporic genera; "
          f"*Synchytrium* -> other Chytridiomycetes; *Pneumocystis* -> *Schizosaccharomyces* + *Taphrina*; Pucciniales -> other Pucciniomycotina; Erysiphaceae -> other Helotiales; Glomeromycota -> Mortierellomycota + Mucoromycota; Entomophthorales -> Basidiobolomycota + Kickxellomycota + Zoopagomycota.\n"
          f"- Intracellular vs extracellular lineages (unit = lineage; Mann-Whitney has very low power at these n and is a screen only):\n\n" + md_table(tb2) + f"\n\n- Table: `{tp2}`.")
    b3 = (f"### HS_G1c Introns and intergenic distance\n\n![]({f3})\n\n- Coloured points = genera of the obligate lineage; grey = genera of its sister set. Q-n50 is not applied here (it would remove most rust and AM assemblies); intergenic distances in fragmented assemblies are biased downward.")
    b4 = (f"### HS_G1d Pfam family richness after proteome-size allometry\n\n![]({f4})\n\n- Log-log slope of distinct Pfam families on proteins per genome, all genera in the frame: {allo['slope']:.3f} (95 % CI {allo['slope_ci_lo']:.3f} to {allo['slope_ci_hi']:.3f}, n = {allo['n']}). "
          f"Per-lineage median residual is in the table above (`pfam_richness_residual`, log10 units; negative = fewer families than expected for the proteome size). "
          f"Lineages with fewer families than their sister (raw genus medians): {len(both_loss)} of {len(tb)} ({', '.join(both_loss['lineage'])}).")
    evo = ("### Evolutionary framing\n\n- Polarity: the ancestral fungal genome in this frame is 25-40 Mb (genus medians: zoosporic " + f"{gd[gd.hs_group == 'zoosporic']['genome_mb'].median():.1f} Mb, Mucoromycota s.l. {gd[gd.hs_group == 'Mucoromycota_sl']['genome_mb'].median():.1f} Mb, filamentous Pezizomycotina-type {gd[gd.hs_group == 'filamentous_other']['genome_mb'].median():.1f} Mb). " +
           "Reduction (intracellular lineages) and expansion (extracellular obligate lineages) are both derived. " + outgroup_text(gd).replace('Microsporidia are excluded (PI decision 2026-11-08) unless the page says otherwise.', f'Microsporidia are kept here as one replicate ({micro_n} genera, < 10 Mb), and the mode test is reported with and without them.') +
           "\n- Placement: each lineage is one origin (single origin inherited by its genera); the replication is across lineages (" + f"{int((tb['mode'] == 'intracellular').sum())} intracellular, {int((tb['mode'] == 'extracellular_obligate').sum())} extracellular). " +
           "Rozellomycota and Aphelidiomycota are already small, so part of the Microsporidia reduction may predate Microsporidia; a 1-2 genus outgroup gives weak polarization.\n"
           "- Mode: reduction = gene loss plus intergenic deletion; expansion = repeat accumulation (masked %) with gene count roughly kept; family richness loss in both. wgd.ks is not used here (rust and AM duplication is an assembly property).\n"
           f"- What the tree adds: `bayou`/`OUwie` shift placement on log genome size, `fitPagel` between obligate mode and size class, ancestral size and richness along each stem to order gene loss and size change. {PENDING}\n"
           "- Caveats: Erysiphaceae pass only with the relaxed BUSCO filter; rust assemblies may be dikaryotic (`busco_duplicated` is in the genus table); Microsporidia fungi_odb12 BUSCO is not a valid completeness measure; *Synchytrium* and Entomophthorales classification as intracellular / extracellular is a hand assignment (to check).")
    # verdict
    intra = tb[tb['mode'] == 'intracellular']; extra = tb[tb['mode'] == 'extracellular_obligate']
    bimodal = (intra['size_ratio'] < 1).sum(), (extra['size_ratio'] > 1).sum()
    verdict = (f"**Direction: {'supported' if bimodal[0] == len(intra) and bimodal[1] == len(extra) else 'partly supported'}.** {bimodal[0]} of {len(intra)} intracellular lineages are smaller than their sister set and {bimodal[1]} of {len(extra)} extracellular obligate lineages are larger "
               f"(genome size ratios: intracellular {', '.join(f'{l} {r:.2f}' for l, r in zip(intra['lineage'], intra['size_ratio']))}; extracellular {', '.join(f'{l} {r:.2f}' for l, r in zip(extra['lineage'], extra['size_ratio']))}). "
               f"Gene density: intracellular median {t_all['genes_per_mb']['median_intracellular']:.0f} vs extracellular {t_all['genes_per_mb']['median_extracellular']:.0f} genes per Mb (unit = lineage, Mann-Whitney p = {S.fmt(t_all['genes_per_mb'].get('mw_p', np.nan))}; without Microsporidia p = {S.fmt(t_nomicro['genes_per_mb'].get('mw_p', np.nan))}). "
               f"Pfam family richness is below the sister in {len(both_loss)} of {len(tb)} lineages. Low lineage counts (n = {len(intra)} + {len(extra)}) make every p-value a screen; the effect sizes are the result. {PENDING}")
    page('HS_G1_genome_size_parasitism', 'HS_G1. Two directions of genome size under obligate parasitism', intro, [b1, b2, b3, b4, evo], verdict)


# ------------------------------------------------------------------ G3
def g3(d):
    dd = d[d['hs_keep'] & (d['hs_group'] != 'Microsporidia')].copy()
    gd = genus_frame(dd, ['introns_per_transcript', 'frac_intronless', 'median_intron_len', 'genome_mb', 'n_proteins'])
    tb, txt = pair_table(gd, 'introns_per_transcript', 'introns per transcript'); tp = save_tab(tb, 'HS_G3', 'pairs_introns')
    f1 = forest(tb, 'HS_G3', 'forest_introns', 'introns per transcript, derived minus sister (genus medians)', 'difference in introns per transcript')
    f2 = strip_groups(gd, 'introns_per_transcript', 'introns per transcript', 'HS_G3', 'introns_by_group')
    f3 = strip_groups(gd, 'median_intron_len', 'median intron length (bp)', 'HS_G3', 'intron_len_by_group')
    tb2, _ = pair_table(gd, 'frac_intronless', 'fraction of intronless transcripts'); tp2 = save_tab(tb2, 'HS_G3', 'pairs_intronless')
    tb3, _ = pair_table(gd, 'median_intron_len', 'median intron length'); tp3 = save_tab(tb3, 'HS_G3', 'pairs_intron_len')
    grp = gd.groupby('hs_group').agg(n_genera=('GENUS', 'size'), introns_per_transcript=('introns_per_transcript', 'median'), frac_intronless=('frac_intronless', 'median'), median_intron_len=('median_intron_len', 'median')).reindex([g for g in GROUPS if g != 'Microsporidia']).dropna().reset_index()
    tpg = save_tab(grp, 'HS_G3', 'group_medians')
    # Basidiomycete yeast origins only (the falsification set)
    bas = tb.iloc[2:5]
    n_bas_drop = int((bas['diff'] < 0).sum()); n_bas_ci = int(((bas['ci_hi'] < 0)).sum())
    mm = mixed_txt(gd[gd['growth_form'].isin(['yeast', 'filamentous_mycelium'])], 'introns_per_transcript', ['gf_yeast', 'log_n50', 'busco_complete'], 'yeast')
    intro = ("Hypothesis G3 (falsification test): unicellular growth does not by itself drive intron loss. Prediction that supports the popular reading of H3: every yeast origin has fewer introns per transcript than its sister. "
             "Prediction that falsifies it: the basidiomycete yeast origins (Tremellomycetes, Pucciniomycotina yeasts, Malasseziomycetes) and the Pezizomycotina dimorphic genera keep sister-like intron densities, so the loss is a Saccharomycotina (and Taphrinomycotina) event.")
    b1 = (f"### HS_G3a Intron density by yeast origin, derived minus sister\n\n![]({f1})\n\n- {counts(dd)}. Filters: Q-ann, Q-rep, Q-busco-rel, Microsporidia excluded (hs_keep). Unit = genus.\n\n" + md_table(tb) + f"\n\n{txt}\n{mm}\n- Table: `{tp}`.")
    b2 = (f"### HS_G3b Intron density and intron length by organism group\n\n![]({f2})\n\n![]({f3})\n\n" + md_table(grp) + f"\n\n- Table: `{tpg}`. Extends V15 (which is by class) with the zoosporic and coenocytic groups.")
    b3 = (f"### HS_G3c Intronless fraction and intron length by pair\n\n" + md_table(tb2) + "\n\n" + md_table(tb3) + f"\n\n- Tables: `{tp2}`, `{tp3}`.")
    evo = ("### Evolutionary framing\n\n- Polarity: the ancestral fungal gene is intron-rich. " + f"Genus-median introns per transcript: zoosporic {grp.set_index('hs_group').loc['zoosporic', 'introns_per_transcript']:.2f}, Mucoromycota s.l. {grp.set_index('hs_group').loc['Mucoromycota_sl', 'introns_per_transcript']:.2f}, Agaricomycetes {grp.set_index('hs_group').loc['Agaricomycetes', 'introns_per_transcript']:.2f}, filamentous non-Agaricomycetes {grp.set_index('hs_group').loc['filamentous_other', 'introns_per_transcript']:.2f}, yeast {grp.set_index('hs_group').loc['yeast', 'introns_per_transcript']:.2f}. " + outgroup_text(gd) +
           f"\n- Placement: the pairs test whether the drop sits on the Saccharomycotina stem only (one event inherited by {int(tb.loc[0, 'n_derived'])} genera) or repeats on the basidiomycete yeast stems. Result above: {n_bas_drop} of 3 basidiomycete yeast origins have a lower median than the sister, {n_bas_ci} with a CI that excludes zero.\n"
           "- Mode: intron loss; a stem burst vs crown decline needs `bayou` on the genus tree. No WGD link tested.\n"
           f"- What the tree adds: `phylosig` (lambda expected near 1), PGLS on growth form, `fitPagel` for growth form vs intron class, ancestral density at the Ascomycota, Pezizomycotina and Saccharomycotina nodes. {PENDING}\n"
           "- Caveats: ab initio intron prediction without RNA-seq (D8); Taphrinomycotina is 3 genera; the Saccharomycotina branch is long; Pezizomycotina (1.8 per transcript) are already reduced relative to the early-diverging groups, so the filamentous sister is not the ancestral state.")
    sac = tb.iloc[0]
    verdict = (f"**Falsification test: {'the popular reading fails' if n_bas_ci == 0 else ('partly fails' if n_bas_ci < 3 else 'the popular reading holds')}.** Saccharomycotina genera have {sac['median_derived']:.2f} introns per transcript against {sac['median_sister']:.2f} in Pezizomycotina (difference {sac['diff']:.2f}, 95 % CI {sac['ci_lo']:.2f} to {sac['ci_hi']:.2f}). "
               f"Basidiomycete yeast origins: " + '; '.join(f"{r['pair']} {r['median_derived']:.2f} vs {r['median_sister']:.2f} (diff {S.fmt(r['diff'], 2)}, CI {S.fmt(r['ci_lo'], 2)} to {S.fmt(r['ci_hi'], 2)}; n = {r['n_derived']} / {r['n_sister']})" for _, r in bas.iterrows()) +
               f". Eurotiomycetes dimorphic vs filamentous: {tb.iloc[5]['median_derived']:.2f} vs {tb.iloc[5]['median_sister']:.2f} (diff {tb.iloc[5]['diff']:.2f}). "
               f"{n_bas_ci} of 3 basidiomycete origins show a loss with a CI excluding zero, so intron loss is not a general property of the yeast transition in this data. {PENDING}")
    page('HS_G3_intron_loss_by_origin', 'HS_G3. Intron loss by independent yeast origin (falsification test)', intro, [b1, b2, b3, evo], verdict)


# ------------------------------------------------------------------ P2
def p2(d):
    dd = d[d['hs_keep'] & (d['hs_group'] != 'Microsporidia')].copy()
    if os.path.exists(f'{A}/A11/idp_genome_model.parquet'):
        im = q(f"SELECT * FROM '{A}/A11/idp_genome_model.parquet'"); dd = dd.merge(im, left_on='LOCUSTAG', right_on='genome', how='left')
        idpcol, idpnote = 'med_idp_model', 'median per-protein disorder fraction over complete gene models (A11)'
    else:
        dd['med_idp_model'] = dd['idp_a5']; dd['med_idp_nonsecreted'] = np.nan; idpcol, idpnote = 'med_idp_model', 'median per-protein disorder fraction over all proteins (A5)'
    gd = genus_frame(dd, [idpcol, 'med_idp_nonsecreted', 'idp_a5', 'n_proteins', 'frac_with_pfam'])
    gd = gd.dropna(subset=[idpcol])
    # GC control: OLS of idp on gc_pct over all genera, residual
    m = gd['gc_pct'].notna(); r = sps.linregress(gd.loc[m, 'gc_pct'], gd.loc[m, idpcol])
    gd['idp_res_gc'] = gd[idpcol] - (r.intercept + r.slope * gd['gc_pct'])
    f1 = strip_groups(gd, idpcol, 'median IDP fraction (complete models)', 'HS_P2', 'idp_by_group')
    tb, txt = pair_table(gd, idpcol, 'median IDP fraction'); tp = save_tab(tb, 'HS_P2', 'pairs_idp')
    tbr, txtr = pair_table(gd, 'idp_res_gc', 'IDP residual after GC'); tpr = save_tab(tbr, 'HS_P2', 'pairs_idp_gc_residual')
    f2 = forest(tb, 'HS_P2', 'forest_idp', 'IDP fraction, derived minus sister', 'difference in median IDP fraction')
    f2b = forest(tbr, 'HS_P2', 'forest_idp_gcres', 'IDP residual after GC regression, derived minus sister', 'difference in GC-adjusted IDP fraction')
    fig, ax = plt.subplots(figsize=(6, 4.3))
    for g in GROUPS:
        s = gd[gd['hs_group'] == g]
        if len(s):
            ax.scatter(s['gc_pct'], s[idpcol], s=8, alpha=0.6, color=GCOL[g], label=f'{g} ({len(s)})', linewidths=0, rasterized=True)
    xx = np.linspace(gd['gc_pct'].min(), gd['gc_pct'].max(), 20); ax.plot(xx, r.intercept + r.slope * xx, color='black', lw=1)
    ax.set_xlabel('genomic GC (%)'); ax.set_ylabel('median IDP fraction'); ax.legend(fontsize=6); ax.set_title(f'HS_P2: disorder vs GC; OLS slope {r.slope:.4f} per GC point, r = {r.rvalue:.2f}, n = {int(m.sum())}', fontsize=8)
    f3 = save_fig(fig, 'HS_P2', 'idp_vs_gc')
    mm1 = mixed_txt(gd[gd['growth_form'].isin(['yeast', 'filamentous_mycelium'])], idpcol, ['gf_yeast', 'gc_pct', 'busco_complete', 'frac_with_pfam'], 'yeast')
    # within-family contrast
    b4 = '### HS_P2d Within-family disorder contrast\n\nNot built: A11/idp_family_genome.parquet missing.'
    famtxt = ''
    if os.path.exists(f'{A}/A11/idp_family_genome.parquet'):
        lab = dd[['LOCUSTAG', 'GENUS', 'growth_form', 'yeast_origin', 'CLASS']].rename(columns={'LOCUSTAG': 'genome'})
        fam = q(f"SELECT genome, sp_id, med_idp, med_len, n_prot FROM '{A}/A11/idp_family_genome.parquet' WHERE genome IN (SELECT LOCUSTAG FROM '{A}/labels/genome_labels_ext.parquet' WHERE hs_keep)")
        fam = fam.merge(lab, on='genome')
        rows = []
        for name, sel_y, sel_f in [('all yeast origins vs filamentous', lambda f: f.growth_form == 'yeast', lambda f: f.growth_form == 'filamentous_mycelium'),
                                   ('Saccharomycotina vs Pezizomycotina', lambda f: (f.yeast_origin == 'Saccharomycotina') & (f.growth_form == 'yeast'), lambda f: f.CLASS.isin(PEZ) & (f.growth_form == 'filamentous_mycelium')),
                                   ('Tremellomycetes vs Agaricomycetes', lambda f: (f.yeast_origin == 'Tremellomycetes') & (f.growth_form == 'yeast'), lambda f: (f.CLASS == 'Agaricomycetes') & (f.growth_form == 'filamentous_mycelium')),
                                   ('Pucciniomycotina yeasts vs Pucciniomycetes', lambda f: (f.yeast_origin == 'Pucciniomycotina-yeast') & (f.growth_form == 'yeast'), lambda f: f.CLASS == 'Pucciniomycetes'),
                                   ('Eurotiomycetes dimorphic vs filamentous', lambda f: (f.CLASS == 'Eurotiomycetes') & (f.growth_form == 'dimorphic_yeast'), lambda f: (f.CLASS == 'Eurotiomycetes') & (f.growth_form == 'filamentous_mycelium'))]:
            fy = fam[sel_y(fam)].groupby(['sp_id', 'GENUS'])['med_idp'].median().groupby('sp_id').agg(['median', 'size'])
            ff = fam[sel_f(fam)].groupby(['sp_id', 'GENUS'])['med_idp'].median().groupby('sp_id').agg(['median', 'size'])
            minn = 10 if name.startswith('all') or name.startswith('Sacch') else 3
            j = fy.join(ff, lsuffix='_y', rsuffix='_f', how='inner'); j = j[(j['size_y'] >= minn) & (j['size_f'] >= minn)]
            diff = (j['median_y'] - j['median_f']).values
            st = S.sign_test(diff); lo, hi = S.boot_ci(diff)
            rows.append({'contrast': name, 'min_genera_per_side': minn, 'n_families': len(j), 'median_diff_idp': float(np.median(diff)) if len(diff) else np.nan, 'ci_lo': lo, 'ci_hi': hi,
                         'frac_derived_lower': st.get('frac_below_zero', np.nan), 'binom_p': st.get('binom_p', np.nan)})
            if name.startswith('all'):
                fig, ax = plt.subplots(figsize=(6, 3.8)); ax.hist(diff, bins=80, range=(-0.3, 0.3), color='#0072B2'); ax.axvline(0, color='black', lw=0.8); ax.axvline(np.median(diff), color='#D55E00', label=f'median {np.median(diff):+.4f}')
                ax.set_xlabel('family median IDP fraction: yeast minus filamentous (genus-level medians inside each family)'); ax.set_ylabel('SwissProt families'); ax.legend(fontsize=8); ax.set_title(f'HS_P2: within-family disorder, n = {len(j):,} families', fontsize=9)
                f4 = save_fig(fig, 'HS_P2', 'family_idp_hist')
        tbf = pd.DataFrame(rows); tpf = save_tab(tbf, 'HS_P2', 'family_contrast')
        b4 = (f"### HS_P2d Within-family disorder contrast\n\n![]({f4})\n\n- Family = SwissProt best-hit family (P-A proxy as V08: qcov, scov >= 80 %, pident >= 40 %, complete models). Per family: genus-median of the per-genome median IDP fraction; difference derived minus sister; families need the stated minimum genera per side. Unit = family.\n\n"
              + md_table(tbf) + f"\n\n- This test removes gene content: it compares disorder in the same protein family. Table: `{tpf}`.")
        famtxt = f" Within families (P-A proxy, {int(tbf.iloc[0]['n_families']):,} families), the yeast copy has a median IDP difference of {tbf.iloc[0]['median_diff_idp']:+.4f} (95 % CI {tbf.iloc[0]['ci_lo']:+.4f} to {tbf.iloc[0]['ci_hi']:+.4f}; {100 * tbf.iloc[0]['frac_derived_lower']:.0f} % of families yeast-lower)."
    grp = gd.groupby('hs_group').agg(n_genera=('GENUS', 'size'), med_idp=(idpcol, 'median'), med_idp_nonsecreted=('med_idp_nonsecreted', 'median'), gc_pct=('gc_pct', 'median'), idp_res_gc=('idp_res_gc', 'median')).reindex([g for g in GROUPS if g != 'Microsporidia']).dropna(subset=['n_genera']).reset_index()
    tpg = save_tab(grp, 'HS_P2', 'group_medians')
    intro = ("Hypothesis P2: intrinsic disorder is lower in every independent yeast origin and in reduced parasites than in filamentous sisters, after GC control. "
             "Disorder predictors respond to amino-acid composition, which tracks GC (H6), so the GC-adjusted residual and the within-family contrast are the stated tests.")
    b1 = (f"### HS_P2a Disorder by organism group\n\n![]({f1})\n\n- {counts(dd)}. Filters: hs_keep (Q-ann, Q-rep, Q-busco-rel, Microsporidia excluded). Trait: {idpnote}; `med_idp_nonsecreted` excludes the strict secreted set. Unit = genus.\n\n" + md_table(grp) + f"\n\n- Table: `{tpg}`.")
    b2 = (f"### HS_P2b Disorder by yeast origin, derived minus sister\n\n![]({f2})\n\n" + md_table(tb) + f"\n\n{txt}\n{mm1}\n- Table: `{tp}`.")
    b3 = (f"### HS_P2c GC control\n\n![]({f3})\n\n![]({f2b})\n\n- OLS of the genus-median IDP fraction on genomic GC over all genera: slope {r.slope:.4f} per GC point (r = {r.rvalue:.2f}, n = {int(m.sum())}). Residuals from this line, derived minus sister:\n\n" + md_table(tbr) + f"\n\n{txtr}\n- Table: `{tpr}`.")
    evo = ("### Evolutionary framing\n\n- Polarity: " + f"genus-median IDP fraction in the early-diverging groups (zoosporic {grp.set_index('hs_group').loc['zoosporic', 'med_idp']:.3f}, Mucoromycota s.l. {grp.set_index('hs_group').loc['Mucoromycota_sl', 'med_idp']:.3f}, Zoopagomycota s.l. {grp.set_index('hs_group').loc['Zoopagomycota_sl', 'med_idp']:.3f}) sets the ancestral value; yeasts ({grp.set_index('hs_group').loc['yeast', 'med_idp']:.3f}) are derived decreases where the pairs show it, dimorphic Pezizomycotina ({grp.set_index('hs_group').loc['dimorphic', 'med_idp']:.3f}) a derived increase. " + outgroup_text(gd) +
           f"\n- Placement: {YEAST_PAIRS} yeast-origin stems plus the Onygenales dimorphic branches; the sign across origins above is the no-tree placement test.\n"
           "- Mode: loss of disordered sequence inside proteins (H2b: 188 of 301 Saccharomycotina shrinker families lose >= 10 disordered residues), tested here within families (HS_P2d), not gene loss.\n"
           f"- What the tree adds: `phylosig` on IDP and on GC, PGLS with GC as covariate, `OUwie` two-regime, ancestral IDP and GC on the same branches to see whether disorder changes where GC does not. {PENDING}\n"
           "- Caveats: predictor composition bias (GC residual applied; the predictor is AIUPred-type per `BFD_function/results/function/aiupred`, to check); Saccharomycotina annotation completeness; truncated models lose termini (complete models only).")
    sac = tb.iloc[0]; sacr = tbr.iloc[0]
    n_ci = int((tb.head(YEAST_PAIRS)['ci_hi'] < 0).sum()); n_cir = int((tbr.head(YEAST_PAIRS)['ci_hi'] < 0).sum())
    n_up = int((tb.head(YEAST_PAIRS)['ci_lo'] > 0).sum())
    verdict = (f"**Not a general yeast trait: {n_ci} of {YEAST_PAIRS} yeast origins are lower with a CI excluding zero ({n_cir} after GC adjustment), {n_up} is higher; the drop is a Saccharomycotina (and Taphrinomycotina, 2 genera) result.** "
               f"Saccharomycotina genera: median IDP fraction {sac['median_derived']:.3f} vs {sac['median_sister']:.3f} in Pezizomycotina (difference {sac['diff']:+.3f}, 95 % CI {sac['ci_lo']:+.3f} to {sac['ci_hi']:+.3f}; GC-adjusted {sacr['diff']:+.3f}, CI {sacr['ci_lo']:+.3f} to {sacr['ci_hi']:+.3f}). "
               f"Basidiomycete yeast origins: " + '; '.join(f"{r['pair']} {S.fmt(r['diff'])} (CI {S.fmt(r['ci_lo'])} to {S.fmt(r['ci_hi'])}; n = {r['n_derived']} / {r['n_sister']})" for _, r in tb.iloc[2:5].iterrows()) + f". Eurotiomycetes dimorphic vs filamentous {S.fmt(tb.iloc[5]['diff'])} (CI {S.fmt(tb.iloc[5]['ci_lo'])} to {S.fmt(tb.iloc[5]['ci_hi'])}).{famtxt} "
               f"The within-family contrasts (HS_P2d) show the same split: Saccharomycotina lower, Tremellomycetes and Pucciniomycotina yeasts higher. The reduced-parasite part of P2 is not tested (Microsporidia excluded by PI decision). {PENDING}")
    page('HS_P2_disorder', 'HS_P2. Intrinsic disorder across yeast origins, with GC control', intro, [b1, b2, b3, b4, evo], verdict)


# ------------------------------------------------------------------ S1
S1_SETS = [('kinase', 'protein kinases (Pkinase + PK_Tyr_Ser-Thr)'), ('ras_arf', 'small GTPases (Ras + Arf)'), ('rho_reg', 'RhoGEF + RhoGAP'), ('gpcr', 'GPCRs (7tm_1, STE2, STE3, opsin)'),
           ('two_comp', 'two-component (HisKA + Response_reg)'), ('g_alpha', 'G-alpha'), ('calcium', 'EF-hand + C2'), ('sh3', 'SH3'), ('ph', 'PH')]


def trait_counts(d):
    tc = q(f"SELECT genome, set_name, n_proteins FROM '{A}/A11/pfam_trait_counts.parquet'").pivot(index='genome', columns='set_name', values='n_proteins').fillna(0)
    return d.merge(tc, left_on='LOCUSTAG', right_index=True, how='left').fillna({c: 0 for c in tc.columns})


def s1(d):
    dd = trait_counts(d[d['hs_keep'] & (d['hs_group'] != 'Microsporidia')].copy())
    sets = [s for s, _ in S1_SETS if s in dd.columns]
    for s in sets:
        dd[f'{s}_per1000'] = 1000 * dd[s] / dd['n_proteins']
    gd = genus_frame(dd, sets + [f'{s}_per1000' for s in sets] + ['n_proteins'])
    # allometry per set
    figs, rows = [], []
    for s, lab in S1_SETS:
        if s not in sets:
            continue
        fp, r = allometry(gd, s, lab, 'HS_S1', f'allometry_{s}'); figs.append(fp)
        rows.append({'set': s, 'label': lab, 'n_genera': r['n'], 'slope': r['slope'], 'slope_ci_lo': r['slope_ci_lo'], 'slope_ci_hi': r['slope_ci_hi'], 'r2': r['r2']})
    tbs = pd.DataFrame(rows); tps = save_tab(tbs, 'HS_S1', 'allometry_slopes')
    # absolute counts and per-1000 by group
    grp = gd.groupby('hs_group').agg(n_genera=('GENUS', 'size'), n_proteins=('n_proteins', 'median'), **{s: (s, 'median') for s in sets}, **{f'{s}_per1000': (f'{s}_per1000', 'median') for s in sets}).reindex([g for g in GROUPS if g != 'Microsporidia']).dropna(subset=['n_genera']).reset_index()
    tpg = save_tab(grp, 'HS_S1', 'group_medians')
    # residuals by fruit-body type and hymenium (Agaricomycetes), by Ascomycota fruit body within class, by yeast pairs
    aga = gd[gd['CLASS'] == 'Agaricomycetes'].copy()
    blocks_res, tabs_res = [], []
    for s in ['kinase', 'ras_arf', 'gpcr', 'two_comp', 'g_alpha']:
        if s not in sets:
            continue
        fpa, ta, kwa = residual_by(aga, f'res_{s}', 'fruitbody_type', FB_ORDER, 'HS_S1', f'res_{s}_fruitbody', f'{s} residual by Agaricomycetes fruit-body type')
        fph, th, kwh = residual_by(aga, f'res_{s}', 'hymenium_type', HY_ORDER, 'HS_S1', f'res_{s}_hymenium', f'{s} residual by Agaricomycetes hymenium type')
        ag_vs_co = S.two_group(aga[aga['fruitbody_type'] == 'agaricoid'][f'res_{s}'], aga[aga['fruitbody_type'] == 'corticioid'][f'res_{s}'], 'agaricoid', 'corticioid')
        gi_vs_sm = S.two_group(aga[aga['hymenium_type'] == 'gills'][f'res_{s}'], aga[aga['hymenium_type'] == 'smooth'][f'res_{s}'], 'gills', 'smooth')
        mm = S.mixed_model(aga.assign(fb_complex=aga['fruitbody_type'].isin(['agaricoid', 'polyporoid', 'gasteroid', 'clavarioid', 'gasteroid-hypogeous']).astype(int)), f'res_{s}', ['fb_complex', 'busco_complete'], group='ORDER_', vc='FAMILY')
        mmt = (f"mixed model residual ~ fb_complex (agaricoid/polyporoid/gasteroid/clavarioid = 1, corticioid = 0) + busco (1|ORDER/FAMILY), n = {mm['n']}: effect {S.fmt([t for t in mm['terms'] if t['term'] == 'fb_complex'][0]['estimate'])} "
               f"(95 % CI {S.fmt([t for t in mm['terms'] if t['term'] == 'fb_complex'][0]['ci_lo'])} to {S.fmt([t for t in mm['terms'] if t['term'] == 'fb_complex'][0]['ci_hi'])})" if mm.get('available') else f"mixed model not fitted: {mm.get('note')}")
        ta['set'] = s; th['set'] = s; tabs_res.append(ta.rename(columns={'fruitbody_type': 'category'}).assign(axis='fruitbody')); tabs_res.append(th.rename(columns={'hymenium_type': 'category'}).assign(axis='hymenium'))
        if s == 'kinase':
            kin_tests = {'ag_vs_co': ag_vs_co, 'gi_vs_sm': gi_vs_sm, 'mm': mm}
        blocks_res.append(f"**{s}.** ![]({fpa}) ![]({fph})\n\n" + md_table(ta) + "\n\n" + md_table(th) +
                          f"\n\n- Kruskal-Wallis across fruit-body types (screen): H = {S.fmt(kwa.get('kw_H'))}, p = {S.fmt(kwa.get('kw_p'))}; across hymenium types: p = {S.fmt(kwh.get('kw_p'))}. "
                          f"Agaricoid vs corticioid residual: {S.fmt(ag_vs_co.get('median_diff'))} (95 % CI {S.fmt(ag_vs_co.get('diff_ci_lo'))} to {S.fmt(ag_vs_co.get('diff_ci_hi'))}, n = {ag_vs_co['n_agaricoid']} / {ag_vs_co['n_corticioid']}); "
                          f"gills vs smooth: {S.fmt(gi_vs_sm.get('median_diff'))} (CI {S.fmt(gi_vs_sm.get('diff_ci_lo'))} to {S.fmt(gi_vs_sm.get('diff_ci_hi'))}). {mmt}.")
    tpr = save_tab(pd.concat(tabs_res), 'HS_S1', 'residuals_agaricomycetes')
    # Ascomycota fruit body within class
    asc = gd[gd['PHYLUM'] == 'Ascomycota'].copy()
    rows = []
    for cls, a_, b_ in [('Sordariomycetes', 'perithecium_(hymenium_hidden,_narrow_opening)', 'none'), ('Leotiomycetes', 'apothecium_(hymenium_on_surface)', 'none'), ('Dothideomycetes', 'perithecium_(hymenium_hidden,_narrow_opening)', 'none')]:
        for s in ['kinase', 'ras_arf', 'gpcr', 'two_comp']:
            if s not in sets:
                continue
            x = asc[(asc['CLASS'] == cls) & (asc['fruitbody_type'] == a_)][f'res_{s}']; y = asc[(asc['CLASS'] == cls) & (asc['fruitbody_type'] == b_)][f'res_{s}']
            r = S.two_group(x, y, 'fruitbody', 'none'); rows.append({'class': cls, 'fruitbody': a_.split('_')[0], 'set': s, 'n_fruitbody': r['n_fruitbody'], 'n_none': r['n_none'], 'median_res_fruitbody': r['median_fruitbody'], 'median_res_none': r['median_none'], 'diff': r.get('median_diff', np.nan), 'ci_lo': r.get('diff_ci_lo', np.nan), 'ci_hi': r.get('diff_ci_hi', np.nan)})
    tba = pd.DataFrame(rows); tpa = save_tab(tba, 'HS_S1', 'residuals_ascomycota_fruitbody')
    # yeast pairs on the kinase residual and on absolute kinase count
    tbk, txtk = pair_table(gd, 'res_kinase', 'kinase residual'); tpk = save_tab(tbk, 'HS_S1', 'pairs_kinase_residual')
    tbc, txtc = pair_table(gd, 'kinase', 'absolute kinase count'); tpc = save_tab(tbc, 'HS_S1', 'pairs_kinase_count')
    fk = forest(tbk, 'HS_S1', 'forest_kinase_residual', 'kinase residual after allometry, derived minus sister', 'difference in residual log10 kinase count')
    intro = ("Hypothesis S1 (falsification test): complex multicellularity does not require more signaling genes. The popular idea predicts positive residuals of kinase, small-GTPase, GPCR and two-component counts "
             "(after proteome-size allometry) in fruit-body-forming Agaricomycetes and Pezizomycotina relative to genera without a fruit body in the same class. Residuals near zero or negative falsify it. "
             "Counts are proteins with >= 1 domain of the set (A3 Pfam architectures, domain e-value <= 1e-5), one count per protein.")
    b1 = (f"### HS_S1a Allometry of signaling-gene counts on proteome size\n\n" + ' '.join(f'![]({f})' for f in figs) + f"\n\n- {counts(dd)}. Filters: hs_keep. Unit = genus. Log-log OLS of genus-median count on genus-median proteins per genome; slope < 1 = sub-linear.\n\n" + md_table(tbs) +
          f"\n\n- Absolute counts and copies per 1,000 proteins by organism group (genus medians):\n\n" + md_table(grp) + f"\n\n- Tables: `{tps}`, `{tpg}`.")
    b2 = (f"### HS_S1b Residuals by Agaricomycetes fruit-body and hymenium type\n\n- Residual = log10(count) minus the allometric prediction from HS_S1a. Fruit-body and hymenium types are FungalTraits genus labels. Unit = genus.\n\n" + '\n\n'.join(blocks_res) + f"\n\n- Table: `{tpr}`.")
    b3 = (f"### HS_S1c Residuals by Ascomycota fruit-body type within class\n\n" + md_table(tba) + f"\n\n- Table: `{tpa}`. 'none' is the FungalTraits label; it mixes asexual-only genera and unknowns.")
    b4 = (f"### HS_S1d Kinase residual and absolute count by yeast origin\n\n![]({fk})\n\n" + md_table(tbk) + f"\n\n{txtk}\n\n" + md_table(tbc) + f"\n\n{txtc}\n- Tables: `{tpk}`, `{tpc}`.")
    g = grp.set_index('hs_group')
    evo = ("### Evolutionary framing\n\n- Polarity: " + f"kinase density is highest in the early-diverging groups (per 1,000 proteins: zoosporic {g.loc['zoosporic', 'kinase_per1000']:.1f}, Mucoromycota s.l. {g.loc['Mucoromycota_sl', 'kinase_per1000']:.1f}, Zoopagomycota s.l. {g.loc['Zoopagomycota_sl', 'kinase_per1000']:.1f}) and lowest in filamentous non-Agaricomycetes ({g.loc['filamentous_other', 'kinase_per1000']:.1f}); absolute kinase counts (genus medians): zoosporic {g.loc['zoosporic', 'kinase']:.0f}, Mucoromycota s.l. {g.loc['Mucoromycota_sl', 'kinase']:.0f}, filamentous non-Agaricomycetes {g.loc['filamentous_other', 'kinase']:.0f}, Agaricomycetes {g.loc['Agaricomycetes', 'kinase']:.0f}, yeast {g.loc['yeast', 'kinase']:.0f}. "
           "The polarized question is therefore whether Dikarya lineages diluted or lost kinases, not whether complex forms gained them. " + outgroup_text(gd) +
           "\n- Placement: no-tree version = residual contrasts by fruit-body type (several Agaricomycetes orders, HS_S1b) and by class (HS_S1c); a single Dikarya-wide shift vs separate shifts at fruit-body origins needs ancestral counts at the Dikarya, Agaricomycetes and fruit-body nodes.\n"
           "- Mode: expected dilution (proteome growth by non-signaling families) rather than kinome change; Mucoromycota kinome expansion may be duplication-driven (wgd.ks young nodes in kinase families: not tested here).\n"
           f"- What the tree adds: PGLS of log count on log proteome size with fruit-body type; `phylosig`; `bayou` on the residual; CAFE-style on kinase family size after dating. {PENDING}\n"
           "- Caveats: subfamily turnover is invisible to counts; the deep order of Mucoromycota and Zoopagomycota is uncertain; fragmented models split kinases (complete-model filter not applied to Pfam counts); FungalTraits 'none' fruit body includes unknowns.")
    # verdict text
    kin_tb = pd.concat(tabs_res); kf = kin_tb[(kin_tb['set'] == 'kinase') & (kin_tb['axis'] == 'fruitbody')].set_index('category')
    kh = kin_tb[(kin_tb['set'] == 'kinase') & (kin_tb['axis'] == 'hymenium')].set_index('category')
    kslope = tbs.set_index('set').loc['kinase']
    ac = kin_tests['ag_vs_co']; gs_ = kin_tests['gi_vs_sm']; mmk = kin_tests['mm']
    mmk_t = [t for t in mmk['terms'] if t['term'] == 'fb_complex'][0] if mmk.get('available') else None
    supported = (ac.get('diff_ci_lo', np.nan) > 0) or (mmk_t is not None and mmk_t['ci_lo'] > 0)
    verdict = (f"**Falsification test: the popular idea is {'supported' if supported else 'not supported'}: inside Agaricomycetes, genera with complex fruit bodies do not carry more kinases than corticioid genera after proteome-size allometry.** "
               f"Kinase count scales with slope {kslope['slope']:.2f} (95 % CI {kslope['slope_ci_lo']:.2f} to {kslope['slope_ci_hi']:.2f}) on proteome size (sub-linear). "
               f"Agaricomycetes as a class sit above the all-fungi allometric line (residuals by fruit-body type, median and 95 % CI: " + '; '.join(f"{c} {kf.loc[c, 'median_residual']:+.3f} ({kf.loc[c, 'ci_lo']:+.3f} to {kf.loc[c, 'ci_hi']:+.3f}, n = {int(kf.loc[c, 'n_genera'])})" for c in kf.index) +
               f"), so the test is the within-class contrast: agaricoid minus corticioid residual {S.fmt(ac.get('median_diff'))} (95 % CI {S.fmt(ac.get('diff_ci_lo'))} to {S.fmt(ac.get('diff_ci_hi'))}; n = {ac['n_agaricoid']} / {ac['n_corticioid']}); gills minus smooth hymenium {S.fmt(gs_.get('median_diff'))} (CI {S.fmt(gs_.get('diff_ci_lo'))} to {S.fmt(gs_.get('diff_ci_hi'))}); "
               + (f"mixed model fb_complex effect {S.fmt(mmk_t['estimate'])} (CI {S.fmt(mmk_t['ci_lo'])} to {S.fmt(mmk_t['ci_hi'])}, ORDER/FAMILY random intercepts). " if mmk_t else "mixed model not fitted. ") +
               f"Absolute kinase counts are {g.loc['yeast', 'kinase']:.0f} in yeast genera and {g.loc['filamentous_other', 'kinase']:.0f} in filamentous non-Agaricomycetes genera; the Saccharomycotina residual is above its sister ({tbk.iloc[0]['diff']:+.3f}, CI {S.fmt(tbk.iloc[0]['ci_lo'])} to {S.fmt(tbk.iloc[0]['ci_hi'])}), so the unicellular proteome keeps more kinases than its size predicts. "
               f"Other signaling classes: see HS_S1b and HS_S1c. {PENDING}")
    page('HS_S1_signaling_scaling', 'HS_S1. Signaling-gene counts scale with proteome size, not with fruit-body complexity (falsification test)', intro, [b1, b2, b3, b4, evo], verdict)


# ------------------------------------------------------------------ C1
def c1(d):
    dd = trait_counts(d[d['hs_keep'] & (d['hs_group'] != 'Microsporidia')].copy())
    nl = q(f"SELECT * FROM '{A}/A11/nlr_arch_counts.parquet'"); dd = dd.merge(nl, left_on='LOCUSTAG', right_on='genome', how='left').fillna({c: 0 for c in nl.columns if c != 'genome'})
    cols = ['het', 'nacht', 'nbarc', 'n_nlr_arch', 'n_het_nlr', 'n_het_only', 'n_sesa', 'n_hets']
    for c in cols:
        dd[f'{c}_per1000'] = 1000 * dd[c] / dd['n_proteins']; dd[f'{c}_present'] = (dd[c] > 0).astype(float)
    gd = genus_frame(dd, cols + [f'{c}_per1000' for c in cols] + [f'{c}_present' for c in cols] + ['n_proteins', 'frac_with_pfam'])
    gd['nlr_sum_per1000'] = gd['het_per1000'] + gd['nacht_per1000'] + gd['nbarc_per1000']
    figs = [strip_groups(gd, f'{c}_per1000', lab, 'HS_C1', f'{c}_by_group') for c, lab in [('het', 'HET-domain proteins per 1,000'), ('nacht', 'NACHT proteins per 1,000'), ('nbarc', 'NB-ARC proteins per 1,000'), ('n_nlr_arch', 'NLR architectures (NACHT/NB-ARC + Ank/TPR/WD40) per 1,000')]]
    grp = gd.groupby('hs_group').agg(n_genera=('GENUS', 'size'), **{f'{c}_per1000': (f'{c}_per1000', 'median') for c in cols}, **{f'{c}_present': (f'{c}_present', 'mean') for c in ['het', 'nacht', 'nbarc', 'n_nlr_arch']}).reindex([g for g in GROUPS if g != 'Microsporidia']).dropna(subset=['n_genera']).reset_index()
    tpg = save_tab(grp, 'HS_C1', 'group_medians')
    pairs, ptxt, fps = {}, {}, {}
    for c in ['het', 'nacht', 'n_nlr_arch']:
        gd[f'log_{c}'] = np.log10(gd[f'{c}_per1000'] + 0.01)
        pairs[c], ptxt[c] = pair_table(gd, f'log_{c}', f'log10({c} per 1,000 + 0.01)'); save_tab(pairs[c], 'HS_C1', f'pairs_{c}')
        fps[c] = forest(pairs[c], 'HS_C1', f'forest_{c}', f'log10 {c} density, derived minus sister', 'difference in log10(per 1,000 + 0.01)')
    aga = gd[gd['CLASS'] == 'Agaricomycetes']
    fa, ta, kwa = residual_by(aga.assign(res_=aga['log_nacht']), 'res_', 'fruitbody_type', FB_ORDER, 'HS_C1', 'nacht_fruitbody', 'log10 NACHT density by Agaricomycetes fruit-body type')
    ta = ta.rename(columns={'median_residual': 'median_log10_density'}); tpa = save_tab(ta, 'HS_C1', 'agaricomycetes_fruitbody_nacht')
    mm = mixed_txt(gd[gd['growth_form'].isin(['yeast', 'filamentous_mycelium'])], 'log_het', ['gf_yeast', 'busco_complete', 'frac_with_pfam'], 'yeast')
    mm2 = mixed_txt(gd[gd['growth_form'].isin(['yeast', 'filamentous_mycelium'])], 'log_nacht', ['gf_yeast', 'busco_complete', 'frac_with_pfam'], 'yeast')
    comp = gd.groupby('hs_group').agg(het_only=('n_het_only_per1000', 'median'), het_nlr=('n_het_nlr_per1000', 'median'), nlr_arch=('n_nlr_arch_per1000', 'median'), sesa=('n_sesa_per1000', 'median'), hets=('n_hets_per1000', 'median')).reindex([g for g in GROUPS if g != 'Microsporidia']).dropna().reset_index()
    tpc = save_tab(comp, 'HS_C1', 'composition')
    intro = ("Hypothesis C1: NLR-type allorecognition repertoires (HET, NACHT, NB-ARC) track septate heterokaryotic hyphae and collapse in every yeast origin, in coenocytic lineages and in dikaryon-forming Agaricomycetes. "
             "Expected: filamentous Pezizomycotina > lichens > Agaricomycetes > Mucoromycota s.l. ~ yeasts ~ 0; within Pezizomycotina, dimorphic Eurotiomycetes below filamentous Eurotiomycetes.")
    b1 = (f"### HS_C1a Density by organism group\n\n" + ' '.join(f'![]({f})' for f in figs) + f"\n\n- {counts(dd)}. Filters: hs_keep. Counts from A3 Pfam architectures (one count per protein). `present` = fraction of genera (genus-level mean of per-genome presence). Unit = genus.\n\n" + md_table(grp) +
          f"\n\n- Composition per 1,000 proteins (genus medians): HET without an NTPase, HET fused to NACHT/NB-ARC, NLR architectures (NTPase + repeat), SesA, HET-s:\n\n" + md_table(comp) + f"\n\n- Tables: `{tpg}`, `{tpc}`.")
    b2 = (f"### HS_C1b Yeast origins and the dimorphic within-class contrast\n\n![]({fps['het']})\n\n![]({fps['nacht']})\n\n![]({fps['n_nlr_arch']})\n\n**HET.**\n\n" + md_table(pairs['het']) + f"\n\n{ptxt['het']}\n{mm}\n\n**NACHT.**\n\n" + md_table(pairs['nacht']) + f"\n\n{ptxt['nacht']}\n{mm2}\n\n**NLR architectures.**\n\n" + md_table(pairs['n_nlr_arch']) + f"\n\n{ptxt['n_nlr_arch']}")
    b3 = (f"### HS_C1c Agaricomycetes by fruit-body type\n\n![]({fa})\n\n" + md_table(ta) + f"\n\n- Kruskal-Wallis (screen): p = {S.fmt(kwa.get('kw_p'))}. Table: `{tpa}`.")
    g = grp.set_index('hs_group')
    evo = ("### Evolutionary framing\n\n- Polarity: " + f"HET density per 1,000 (genus medians): zoosporic {g.loc['zoosporic', 'het_per1000']:.2f} (present in {100 * g.loc['zoosporic', 'het_present']:.0f} % of genera), Mucoromycota s.l. {g.loc['Mucoromycota_sl', 'het_per1000']:.2f} ({100 * g.loc['Mucoromycota_sl', 'het_present']:.0f} %), Zoopagomycota s.l. {g.loc['Zoopagomycota_sl', 'het_per1000']:.2f}, filamentous non-Agaricomycetes {g.loc['filamentous_other', 'het_per1000']:.2f}, lichen {g.loc['lichen', 'het_per1000']:.2f}, Agaricomycetes {g.loc['Agaricomycetes', 'het_per1000']:.2f}, yeast {g.loc['yeast', 'het_per1000']:.2f}. "
           f"NACHT: zoosporic {g.loc['zoosporic', 'nacht_per1000']:.2f} ({100 * g.loc['zoosporic', 'nacht_present']:.0f} %), Mucoromycota s.l. {g.loc['Mucoromycota_sl', 'nacht_per1000']:.2f} ({100 * g.loc['Mucoromycota_sl', 'nacht_present']:.0f} %), filamentous non-Agaricomycetes {g.loc['filamentous_other', 'nacht_per1000']:.2f}, Agaricomycetes {g.loc['Agaricomycetes', 'nacht_per1000']:.2f}, yeast {g.loc['yeast', 'nacht_per1000']:.2f}. "
           "So the HET-domain expansion is a derived gain on the Pezizomycotina stem (not a loss elsewhere), NACHT is present at low copy outside Dikarya and expanded in both Dikarya subphyla, and the yeast values are derived losses within Dikarya. " + outgroup_text(gd) +
           "\n- Placement: HET gain = one event inherited by Pezizomycotina (no replication); NACHT expansion at the Dikarya node or twice (Ascomycota and Basidiomycota) needs ancestral counts; losses on the 5 yeast stems and the Onygenales dimorphic branches are the replicated part (HS_C1b).\n"
           "- Mode: family expansion by duplication and domain shuffling along Pezizomycotina branches (gradual; genus counts vary widely), loss as a burst on yeast stems; wgd.ks young nodes inside HET/NACHT families not tested here.\n"
           f"- What the tree adds: `phylosig` (high lambda expected), PGLS on growth form, `bayou` shift branches on log density, CAFE-style after dating, ancestral HET count along the Eurotiomycetes path to dimorphic Onygenales to order NLR contraction and the yeast phase. {PENDING}\n"
           "- Caveats: HET and NACHT HMMs are clade-sensitive (`frac_with_pfam` is a covariate in the mixed model); strain variation in *Neurospora* and *Podospora* (genus medians); Pezizomycotina stem is one branch.")
    ph = pairs['het']; pn = pairs['nacht']
    dim_fail = ph.iloc[5]['ci_lo'] > 0
    verdict = (f"**Supported across organism groups and in {int((ph.head(5)['ci_hi'] < 0).sum())} of 5 yeast origins (2 origins have zero HET on both sides); the within-class dimorphic prediction {'fails: dimorphic Eurotiomycetes carry more HET proteins than filamentous Eurotiomycetes' if dim_fail else 'holds'}.** HET per 1,000 proteins (genus medians): filamentous non-Agaricomycetes {g.loc['filamentous_other', 'het_per1000']:.2f}, lichen {g.loc['lichen', 'het_per1000']:.2f}, Agaricomycetes {g.loc['Agaricomycetes', 'het_per1000']:.2f}, Mucoromycota s.l. {g.loc['Mucoromycota_sl', 'het_per1000']:.2f}, yeast {g.loc['yeast', 'het_per1000']:.2f}. "
               f"Yeast origins (log10 density difference to sister, 95 % CI): " + '; '.join(f"{r['pair']} {S.fmt(r['diff'], 2)} ({S.fmt(r['ci_lo'], 2)} to {S.fmt(r['ci_hi'], 2)})" for _, r in ph.head(5).iterrows()) +
               f". Eurotiomycetes dimorphic vs filamentous: HET {S.fmt(ph.iloc[5]['diff'], 2)} ({S.fmt(ph.iloc[5]['ci_lo'], 2)} to {S.fmt(ph.iloc[5]['ci_hi'], 2)}), NACHT {S.fmt(pn.iloc[5]['diff'], 2)} ({S.fmt(pn.iloc[5]['ci_lo'], 2)} to {S.fmt(pn.iloc[5]['ci_hi'], 2)}). "
               f"NACHT is not lower in Agaricomycetes ({g.loc['Agaricomycetes', 'nacht_per1000']:.2f} vs {g.loc['filamentous_other', 'nacht_per1000']:.2f} per 1,000), so the Agaricomycetes collapse holds for HET, not for NACHT-type NLRs. {PENDING}")
    page('HS_C1_nlr_allorecognition', 'HS_C1. NLR allorecognition repertoires and the septate heterokaryotic habit', intro, [b1, b2, b3, evo], verdict)


# ------------------------------------------------------------------ C3
C3_GENES = [('HEX1', 'woronin'), ('WSC', 'woronin'), ('PRM1', 'fusion'), ('FIG1', 'fusion'), ('HAM7', 'fusion'), ('FUS3_MAPK', 'fusion_mapk'), ('STE50', 'fusion_mapk'), ('STE11', 'fusion_mapk'), ('STE7', 'fusion_mapk'), ('STE12', 'fusion_mapk'), ('SLT2_MAPK', 'cwi_control'), ('HOG1_MAPK', 'hog_control')]


def c3(d):
    dd = d[d['hs_keep'] & (d['hs_group'] != 'Microsporidia')].copy()
    genes = [g for g, _ in C3_GENES]
    # detector 1 (primary): A9 anchors = best hit inside one reference proteome (YEAST, NEUCR, SCHPO, EMENI); a gene is present when any reference anchor hits it
    pr = q(f"SELECT genome, gene, sum(n_relaxed) AS n_relaxed, sum(n_strict) AS n_strict, max(max_pident_x10) AS max_pident_x10 FROM '{A}/A11/c3_anchor_presence.parquet' GROUP BY 1,2")
    # detector 2 (coverage note only): best SwissProt hit over all organisms equals one of the gene-set entries
    pr2 = q(f"SELECT genome, gene, n_relaxed FROM '{A}/A11/c3_gene_presence.parquet'")
    rel2 = (pr2.pivot(index='genome', columns='gene', values='n_relaxed').reindex(columns=genes).fillna(0) > 0).astype(float)
    bh = dd[['LOCUSTAG', 'GENUS', 'hs_group']].merge(rel2, left_on='LOCUSTAG', right_index=True, how='left').fillna(0)
    bh_hog = bh.groupby('GENUS').agg(hs_group=('hs_group', 'first'), hog=('HOG1_MAPK', 'mean'), hex1=('HEX1', 'mean')).groupby('hs_group')[['hog', 'hex1']].mean()
    rel = pr.pivot(index='genome', columns='gene', values='n_relaxed').reindex(columns=genes).fillna(0); strict = pr.pivot(index='genome', columns='gene', values='n_strict').reindex(columns=genes).fillna(0)
    dd = dd.merge((rel > 0).astype(float).add_prefix('rel_'), left_on='LOCUSTAG', right_index=True, how='left').merge((strict > 0).astype(float).add_prefix('str_'), left_on='LOCUSTAG', right_index=True, how='left').fillna({f'rel_{g}': 0 for g in genes} | {f'str_{g}': 0 for g in genes})
    for mod in ['woronin', 'fusion', 'fusion_mapk']:
        gs = [g for g, m in C3_GENES if m == mod]; dd[f'score_{mod}'] = dd[[f'rel_{g}' for g in gs]].mean(axis=1)
    gd = genus_frame(dd, [f'rel_{g}' for g in genes] + [f'str_{g}' for g in genes] + ['score_woronin', 'score_fusion', 'score_fusion_mapk', 'frac_with_swissprot', 'n_strict_hits'])
    # genus-level presence = fraction of representative genomes in the genus with >= 1 hit, then present if >= 0.5
    for g in genes:
        gd[f'pres_{g}'] = (gd[f'rel_{g}'] >= 0.5).astype(float); gd[f'pres_str_{g}'] = (gd[f'str_{g}'] >= 0.5).astype(float)
    grp_order = [g for g in GROUPS if g != 'Microsporidia']
    mat = gd.groupby('hs_group')[[f'pres_{g}' for g in genes]].mean().reindex(grp_order).dropna(); mat.columns = genes
    mats = gd.groupby('hs_group')[[f'pres_str_{g}' for g in genes]].mean().reindex(grp_order).dropna(); mats.columns = genes
    ng = gd.groupby('hs_group').size().reindex(grp_order)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax, m_, t_ in zip(axs, [mat, mats], ['relaxed gates (pident >= 25 %, qcov and scov >= 50 %)', 'strict gates (pident >= 40 %, qcov and scov >= 80 %)']):
        im = ax.imshow(m_.values, cmap='Blues', vmin=0, vmax=1, aspect='auto'); ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=60, fontsize=7, ha='right'); ax.set_yticks(range(len(m_))); ax.set_yticklabels([f'{g} (n={int(ng[g])})' for g in m_.index], fontsize=7)
        for i in range(m_.shape[0]):
            for j in range(m_.shape[1]):
                ax.text(j, i, f'{m_.values[i, j]:.2f}', ha='center', va='center', fontsize=6, color='white' if m_.values[i, j] > 0.6 else 'black')
        ax.set_title(t_, fontsize=8)
    fig.colorbar(im, ax=axs, fraction=0.02, label='fraction of genera with the gene'); fig.suptitle('HS_C3: Woronin body, fusion and MAPK gene presence by organism group (best SwissProt hit in the gene set)', fontsize=9)
    f1 = save_fig(fig, 'HS_C3', 'presence_matrix')
    tpm = save_tab(mat.reset_index(), 'HS_C3', 'presence_relaxed'); tpms = save_tab(mats.reset_index(), 'HS_C3', 'presence_strict')
    # module score strips
    f2 = [strip_groups(gd, f'score_{m}', f'{m} module score (fraction of genes present)', 'HS_C3', f'score_{m}') for m in ['woronin', 'fusion', 'fusion_mapk']]
    # pairs: presence fraction derived vs sister per gene
    rows = []
    for name, fd, fs in PAIRS:
        a = gd[fd(gd)]; b = gd[fs(gd)]
        for g in genes:
            rows.append({'pair': name, 'gene': g, 'n_derived': len(a), 'n_sister': len(b), 'frac_derived': float(a[f'pres_{g}'].mean()) if len(a) else np.nan, 'frac_sister': float(b[f'pres_{g}'].mean()) if len(b) else np.nan})
    tbp = pd.DataFrame(rows); tbp['diff'] = tbp['frac_derived'] - tbp['frac_sister']; tpp = save_tab(tbp, 'HS_C3', 'pairs_presence')
    piv = tbp.pivot(index='gene', columns='pair', values='diff').reindex(genes)
    fig, ax = plt.subplots(figsize=(9, 4.5)); im = ax.imshow(piv.values, cmap='RdBu', vmin=-1, vmax=1, aspect='auto'); ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=40, ha='right', fontsize=7); ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=7)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, f'{piv.values[i, j]:+.2f}', ha='center', va='center', fontsize=6)
    fig.colorbar(im, ax=ax, fraction=0.03, label='presence fraction, derived minus sister'); ax.set_title('HS_C3: gene presence by sister pair (genus level)', fontsize=9)
    f3 = save_fig(fig, 'HS_C3', 'pairs_heatmap')
    # anchor-rate covariate: swissprot coverage by group
    cov = gd.groupby('hs_group').agg(n_genera=('GENUS', 'size'), frac_with_swissprot=('frac_with_swissprot', 'median'), n_strict_hits=('n_strict_hits', 'median'), score_woronin=('score_woronin', 'median'), score_fusion=('score_fusion', 'median'), score_fusion_mapk=('score_fusion_mapk', 'median')).reindex(grp_order).dropna().reset_index()
    tpc = save_tab(cov, 'HS_C3', 'coverage_and_scores')
    intro = ("Hypothesis C3: Woronin-body and cell-fusion gene sets are present only in lineages with septate, fusing hyphae; Agaricomycetes replace the Pezizomycotina Woronin set rather than lack fusion machinery. "
             "Design change from the hypothesis text: SwissProt holds no *N. crassa* septal-pore (SPA) or HAM-5/HAM-6/HAM-8 proteins and only three of the planned Woronin/fusion proteins for NEUCR (HEX-1, WSC, PRM-1), so the gene set uses SwissProt entries from several fungi "
             "(HEX1 from *N. crassa* and *A. nidulans*, WSC, PRM1 from four reference species, FIG1, HAM-7, the FUS3/KSS1/SPK1 MAPK with STE50/STE11/STE7/STE12) and two control MAPKs (SLT2/SPM1 cell-wall integrity, HOG1). "
             "Presence = a protein whose best hit inside one of the four reference proteomes (A9 anchors: *S. cerevisiae*, *N. crassa*, *S. pombe*, *A. nidulans*) is the gene's entry, with relaxed gates (pident >= 25 %, reference coverage >= 50 %; no query-coverage gate, because HEX-1 orthologs carry N-terminal extensions: mean query coverage 58 % on the 176-aa HEX1_NEUCR) and strict gates (pident >= 40 %, reference coverage >= 80 %). "
             "The best-SwissProt-hit detector (A11/c3_gene_presence) was tried first and rejected: it finds HOG1 in only " + f"{100 * bh_hog.loc['filamentous_other', 'hog']:.0f} % of filamentous non-Agaricomycetes genera, because the best hit of an ortholog is usually an entry from another species. "
             "The septal-pore module is not testable with the tables on disk.")
    b1 = (f"### HS_C3a Presence matrix\n\n![]({f1})\n\n- {counts(dd)}. Filters: hs_keep. Genus presence = >= 50 % of the genus's representative genomes have >= 1 anchored hit. Cell = fraction of genera.\n\n" + md_table(mat.reset_index()) + f"\n\n- SwissProt coverage and module scores by group (a low `n_strict_hits` means absence may be detection failure; `frac_with_swissprot` = fraction of proteins with any SwissProt hit):\n\n" + md_table(cov) + f"\n\n- Tables: `{tpm}`, `{tpms}`, `{tpc}`.")
    b2 = (f"### HS_C3b Module scores\n\n" + ' '.join(f'![]({f})' for f in f2) + "\n\n- Module score = fraction of the module's genes with >= 1 relaxed hit in the genome, then genus median.")
    b3 = (f"### HS_C3c Sister pairs\n\n![]({f3})\n\n" + md_table(tbp[tbp['gene'].isin(['HEX1', 'WSC', 'PRM1', 'FIG1', 'HAM7', 'FUS3_MAPK', 'STE12'])].sort_values(['gene', 'pair']), max_rows=50) + f"\n\n- Table: `{tpp}`.")
    evo = ("### Evolutionary framing\n\n- Polarity: " + f"HEX1 presence (fraction of genera): filamentous non-Agaricomycetes {mat.loc['filamentous_other', 'HEX1']:.2f}, lichen {mat.loc['lichen', 'HEX1']:.2f}, dimorphic {mat.loc['dimorphic', 'HEX1']:.2f}, yeast {mat.loc['yeast', 'HEX1']:.2f}, Agaricomycetes {mat.loc['Agaricomycetes', 'HEX1']:.2f}, Mucoromycota s.l. {mat.loc['Mucoromycota_sl', 'HEX1']:.2f}, zoosporic {mat.loc['zoosporic', 'HEX1']:.2f}; "
           f"PRM1: filamentous non-Agaricomycetes {mat.loc['filamentous_other', 'PRM1']:.2f}, Agaricomycetes {mat.loc['Agaricomycetes', 'PRM1']:.2f}, yeast {mat.loc['yeast', 'PRM1']:.2f}, Mucoromycota s.l. {mat.loc['Mucoromycota_sl', 'PRM1']:.2f}, zoosporic {mat.loc['zoosporic', 'PRM1']:.2f}. "
           "The Woronin set is a Pezizomycotina gain if absent from all outgroups; the fusion/MAPK set is ancestral if present in the early-diverging groups. " + outgroup_text(gd) +
           "\n- Placement: Woronin gain on the Pezizomycotina stem (one event); fusion-gene losses on yeast stems are the replicated part (HS_C3c).\n"
           "- Mode: gene birth (Woronin) and gene loss (fusion genes in yeasts); no expansion expected.\n"
           f"- What the tree adds: `make.simmap` per gene, `fitPagel` growth form vs module presence, phylogenetic logistic regression with `n_strict_hits` as covariate; ancestral presence along the Saccharomycotina path to order fusion-gene loss against growth-form change. {PENDING}\n"
           "- Caveats: anchor presence fades with distance from the reference organisms (coverage table above); absence in Basidiomycota and early-diverging lineages can be divergence; SwissProt lacks the septal-pore proteins; the gene set is small. "
           "STE50, STE7 and STE12 have anchors only in *S. cerevisiae* and *S. pombe*, and their relaxed-gate presence outside yeasts is near zero even in Pezizomycotina, so those three rows measure anchor reach, not gene loss; HOG1 and SLT2 (anchors in all four references) are the detection controls.")
    hx = mat['HEX1']; pm = mat['PRM1']; f3m = mat['FUS3_MAPK']
    ham = mat['HAM7']
    verdict = (f"**Woronin module: {'supported' if hx['filamentous_other'] >= 0.8 and max(hx.get('Mucoromycota_sl', 0), hx.get('zoosporic', 0), hx.get('Agaricomycetes', 0)) <= 0.2 else 'partly supported'} (Pezizomycotina gain, absent from Agaricomycetes, yeasts and all outgroups). "
               f"Fusion module: PRM1 is shared by Dikarya hyphal forms and yeasts and is rare outside Dikarya; HAM-7 is Pezizomycotina-restricted ({100 * ham['filamentous_other']:.0f} % filamentous non-Agaricomycetes, {100 * ham['Agaricomycetes']:.0f} % Agaricomycetes, {100 * ham['yeast']:.0f} % yeast); the FUS3-type MAPK module is universal where detection allows.** "
               f"HEX1 present in {100 * hx['filamentous_other']:.0f} % of filamentous non-Agaricomycetes genera, {100 * hx['lichen']:.0f} % of lichen, {100 * hx['dimorphic']:.0f} % of dimorphic, {100 * hx['yeast']:.0f} % of yeast, {100 * hx['Agaricomycetes']:.0f} % of Agaricomycetes, {100 * hx['Mucoromycota_sl']:.0f} % of Mucoromycota s.l. and {100 * hx['zoosporic']:.0f} % of zoosporic genera. "
               f"PRM1: {100 * pm['filamentous_other']:.0f} % / {100 * pm['Agaricomycetes']:.0f} % / {100 * pm['yeast']:.0f} % / {100 * pm['Mucoromycota_sl']:.0f} % / {100 * pm['zoosporic']:.0f} % (filamentous / Agaricomycetes / yeast / Mucoromycota s.l. / zoosporic). FUS3-type MAPK: {100 * f3m['filamentous_other']:.0f} % / {100 * f3m['Agaricomycetes']:.0f} % / {100 * f3m['yeast']:.0f} % / {100 * f3m['Mucoromycota_sl']:.0f} % / {100 * f3m['zoosporic']:.0f} %. "
               f"Control MAPKs (HOG1): filamentous {100 * mat.loc['filamentous_other', 'HOG1_MAPK']:.0f} %, zoosporic {100 * mat.loc['zoosporic', 'HOG1_MAPK']:.0f} %, which bounds the detection rate in distant clades. The septal-pore module is not testable with the tables on disk. {PENDING}")
    page('HS_C3_woronin_fusion', 'HS_C3. Woronin-body and cell-fusion gene sets across organization shifts', intro, [b1, b2, b3, evo], verdict)


VIEWS = {'g1': g1, 'g3': g3, 'p2': p2, 's1': s1, 'c1': c1, 'c3': c3}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--only', default=''); a = ap.parse_args()
    only = [x.strip().lower() for x in a.only.split(',') if x.strip()] or list(VIEWS)
    d = genome_frame()
    status = {}
    for v in only:
        try:
            print('render', v, flush=True); VIEWS[v](d.copy()); status[v] = 'ok'
        except Exception as e:
            status[v] = f'failed: {e}'; print('FAILED', v, e, flush=True); traceback.print_exc()
    print(status)


if __name__ == '__main__':
    main()
