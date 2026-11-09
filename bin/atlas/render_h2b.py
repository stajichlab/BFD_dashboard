#!/usr/bin/env python3
"""Render docs/summary_atlas/H2b_compaction_extremes.md from atlas_data/A10 (built by h2b_build.py)."""
import os
import sys

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
A10 = f'{A}/A10'
save_fig, save_tab, md_table = RV.save_fig, RV.save_tab, RV.md_table
GROUP_COL = {'Saccharomycotina': '#D55E00', 'Pezizomycotina-filamentous': '#0072B2', 'Tremellomycetes': '#CC79A7', 'filamentous-other': '#56B4E9',
             'Pucciniomycotina-yeast': '#009E73', 'Ustilaginomycotina-yeast': '#E69F00'}
TRAITS = [('log_genome', 'log10 genome size (bp)'), ('genes_per_mb', 'protein-coding genes per Mb'), ('median_intergenic', 'median intergenic distance (bp)'),
          ('introns_per_transcript', 'introns per transcript')]
SHRINK_RULE = 'n_sacch >= 10, n_fil >= 10, >= 80 % of Saccharomycotina genera below the filamentous 10th percentile, loss >= 20 aa or >= 5 %'
GAIN_RULE = 'n_sacch >= 10, n_fil >= 10, >= 80 % of Saccharomycotina genera above the filamentous 90th percentile, gain >= 20 aa or >= 5 %'

con = duckdb.connect(); con.execute("SET threads=4; SET memory_limit='20GB'")
q = lambda s: con.execute(s).df()


def rho_ci(x, y, n=1000):
    x = np.asarray(x, float); y = np.asarray(y, float); m = ~np.isnan(x) & ~np.isnan(y); x, y = x[m], y[m]
    if len(x) < 5:
        return {'n': int(len(x)), 'rho': np.nan, 'ci_lo': np.nan, 'ci_hi': np.nan, 'p': np.nan}
    r = sps.spearmanr(x, y); idx = S.RNG.integers(0, len(x), size=(n, len(x)))
    bs = [sps.spearmanr(x[i], y[i]).statistic for i in idx]
    return {'n': int(len(x)), 'rho': float(r.statistic), 'ci_lo': float(np.quantile(bs, .025)), 'ci_hi': float(np.quantile(bs, .975)), 'p': float(r.pvalue)}


def ols(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); m = ~np.isnan(x) & ~np.isnan(y)
    if m.sum() < 5:
        return {'slope': np.nan, 'slope_lo': np.nan, 'slope_hi': np.nan}
    r = sps.linregress(x[m], y[m]); t = sps.t.ppf(.975, m.sum() - 2)
    return {'slope': float(r.slope), 'slope_lo': float(r.slope - t * r.stderr), 'slope_hi': float(r.slope + t * r.stderr), 'r2': float(r.rvalue ** 2)}


def genus_frame():
    g = q(f"SELECT * FROM '{A10}/genome_offset_traits.parquet'")
    lab = q(f"SELECT LOCUSTAG AS genome, FAMILY FROM '{A}/labels/genome_labels.parquet'")
    g = g.merge(lab, on='genome')
    g['log_intergenic'] = np.log10(g['median_intergenic'].clip(lower=1))
    agg = g.groupby('GENUS').agg(offset=('offset', 'median'), n_fam=('n_fam', 'median'), log_genome=('log_genome', 'median'), genes_per_mb=('genes_per_mb', 'median'),
                                 median_intergenic=('median_intergenic', 'median'), log_intergenic=('log_intergenic', 'median'), introns_per_transcript=('introns_per_transcript', 'median'),
                                 busco_complete=('busco_complete', 'median'), log_n50=('log_n50', 'median'), partial_fraction=('partial_fraction', 'median'),
                                 rnaseq_frac=('rnaseq_trained', 'mean'), n_genomes=('genome', 'size'),
                                 clade_group=('clade_group', lambda s: s.mode().iloc[0]), growth_form=('growth_form', 'first'), CLASS=('CLASS', 'first'),
                                 ORDER_=('ORDER_', 'first'), FAMILY=('FAMILY', 'first'), PHYLUM=('PHYLUM', 'first')).reset_index()
    agg['gf_yeast'] = (agg['growth_form'] == 'yeast').astype(int)
    agg['is_sacch'] = (agg['clade_group'] == 'Saccharomycotina').astype(int)
    return g, agg


# ------------------------------------------------------------------ Q1
def v50(g, agg):
    groups = [c for c in ['Saccharomycotina', 'Pezizomycotina-filamentous', 'Tremellomycetes', 'filamentous-other', 'Pucciniomycotina-yeast', 'Ustilaginomycotina-yeast']
              if (agg['clade_group'] == c).sum() >= 10]
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5)); axs = axs.ravel(); rows = []
    for ax, (tr, lab) in zip(axs, TRAITS):
        xcol = 'log_intergenic' if tr == 'median_intergenic' else tr
        for grp in groups:
            s = agg[agg['clade_group'] == grp]
            ax.scatter(s[xcol], s['offset'], s=9, alpha=0.6, color=GROUP_COL.get(grp, '#999'), label=f'{grp} (n={len(s)})', linewidths=0, rasterized=True)
            r = rho_ci(s[xcol], s['offset']); o = ols(s[xcol], s['offset'])
            rows.append({'trait': tr, 'clade_group': grp, **r, **o})
        r = rho_ci(agg[xcol], agg['offset']); o = ols(agg[xcol], agg['offset']); rows.append({'trait': tr, 'clade_group': 'all genera', **r, **o})
        ax.set_xlabel(('log10 ' if tr == 'median_intergenic' else '') + lab); ax.set_ylabel('genus median length offset (log2)'); ax.axhline(0, color='grey', lw=0.6)
    axs[0].legend(fontsize=7)
    fig.suptitle('V50: per-genome ortholog length offset (P-A) vs genome compaction traits, genus medians', fontsize=9); fig.tight_layout()
    fp = save_fig(fig, 'V50', 'offset_vs_compaction')
    tb = pd.DataFrame(rows); tp = save_tab(tb, 'V50', 'spearman_by_clade')
    tb2 = tb[['trait', 'clade_group', 'n', 'rho', 'ci_lo', 'ci_hi', 'slope', 'slope_lo', 'slope_hi']]
    return (f"### V50 Length offset versus genome compaction, within and between clades\n\n![]({fp})\n\n"
            f"- Unit = genus (median over representative genomes; {len(agg):,} genera, {len(g):,} genomes). Offset = median over SwissProt families (>= 200 per genome) of log2(genome length / filamentous family median), as V09. "
            f"Filters: {RV.Ctx.filt_text(None)}. Compaction traits from A2 and asm_stats. Spearman rho with bootstrap 95 % CI; OLS slope of offset on the trait (intergenic distance on log10 scale).\n\n" + md_table(tb2, 40) +
            f"\n\n- Table: `{tp}`. PGLS pending tree.")


def v51(agg):
    rows = []
    for tr, lab in TRAITS:
        xcol = 'log_intergenic' if tr == 'median_intergenic' else tr
        d = agg.dropna(subset=[xcol, 'offset']).copy(); d['x'] = (d[xcol] - d[xcol].mean()) / d[xcol].std()
        mm = S.mixed_model(d, 'offset', ['x', 'gf_yeast'], group='CLASS', vc='ORDER_')
        if mm.get('available'):
            tx = [t for t in mm['terms'] if t['term'] == 'x'][0]; tg = [t for t in mm['terms'] if t['term'] == 'gf_yeast'][0]
            rows.append({'model': 'all genera: offset ~ trait + yeast + (1|CLASS/ORDER)', 'trait': tr, 'n_genera': mm['n'], 'trait_effect_per_SD': tx['estimate'], 'ci_lo': tx['ci_lo'], 'ci_hi': tx['ci_hi'], 'p': tx['p'],
                         'yeast_effect': tg['estimate'], 'yeast_ci_lo': tg['ci_lo'], 'yeast_ci_hi': tg['ci_hi'], 'converged': mm['converged']})
        else:
            rows.append({'model': 'all genera', 'trait': tr, 'note': mm.get('note')})
        for grp in ['Saccharomycotina', 'Pezizomycotina-filamentous']:
            s = d[d['clade_group'] == grp].copy()
            mm = S.mixed_model(s, 'offset', ['x'], group='ORDER_', vc='FAMILY')
            if mm.get('available'):
                tx = [t for t in mm['terms'] if t['term'] == 'x'][0]
                rows.append({'model': f'{grp}: offset ~ trait + (1|ORDER/FAMILY)', 'trait': tr, 'n_genera': mm['n'], 'trait_effect_per_SD': tx['estimate'], 'ci_lo': tx['ci_lo'], 'ci_hi': tx['ci_hi'], 'p': tx['p'], 'converged': mm['converged']})
            else:
                rows.append({'model': f'{grp}: within-clade', 'trait': tr, 'note': mm.get('note')})
    tb = pd.DataFrame(rows); tp = save_tab(tb, 'V51', 'mixed_models')
    cols = [c for c in ['model', 'trait', 'n_genera', 'trait_effect_per_SD', 'ci_lo', 'ci_hi', 'p', 'yeast_effect', 'yeast_ci_lo', 'yeast_ci_hi', 'converged', 'note'] if c in tb]
    return (f"### V51 Partial association of the offset with compaction after growth form and taxonomy\n\n- Mixed models on genus medians (statsmodels MixedLM, REML). Traits standardized (effect = change in log2 offset per 1 SD of the trait). "
            f"Across all genera: random intercepts CLASS and ORDER within CLASS, growth form (yeast = 1) as fixed covariate. Within a clade: random intercepts ORDER and FAMILY within ORDER.\n\n" + md_table(tb[cols], 20) +
            f"\n\n- Reading: a trait effect with a CI excluding zero inside Saccharomycotina means compaction predicts shorter orthologs within the clade; the yeast effect in the all-genera model is the growth-form difference that the trait does not absorb. Table: `{tp}`. PGLS pending tree.")


def v52(g, agg):
    s = agg[agg['clade_group'] == 'Saccharomycotina'].copy(); rows = []
    for col, lab in [('busco_complete', 'BUSCO complete (%)'), ('log_n50', 'log10 N50'), ('partial_fraction', 'partial-model fraction'), ('log_genome', 'log10 genome size'), ('genes_per_mb', 'genes per Mb')]:
        rows.append({'covariate': lab, **rho_ci(s[col], s['offset'])})
    tb = pd.DataFrame(rows)
    # RNA-seq-trained: genome level within Saccharomycotina
    gs = g[g['clade_group'] == 'Saccharomycotina']; rt = gs[gs['rnaseq_trained'] == True]['offset']; rn = gs[gs['rnaseq_trained'] == False]['offset']
    r2 = S.two_group(rt, rn, 'rnaseq', 'no_rnaseq')
    # joint OLS: offset ~ log_genome + busco + log_n50 + partial (standardized)
    X = s[['log_genome', 'busco_complete', 'log_n50', 'partial_fraction']].dropna(); y = s.loc[X.index, 'offset']
    Xs = (X - X.mean()) / X.std(); Xm = np.column_stack([np.ones(len(Xs)), Xs.values])
    beta, res, rk, sv = np.linalg.lstsq(Xm, y.values, rcond=None); resid = y.values - Xm @ beta; dof = len(y) - Xm.shape[1]
    sigma2 = (resid ** 2).sum() / dof; cov = sigma2 * np.linalg.inv(Xm.T @ Xm); se = np.sqrt(np.diag(cov)); t = sps.t.ppf(.975, dof)
    joint = pd.DataFrame({'term': ['intercept'] + list(X.columns), 'coef_per_SD': beta, 'ci_lo': beta - t * se, 'ci_hi': beta + t * se, 'p': 2 * sps.t.sf(np.abs(beta / se), dof)})
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.6))
    for ax, (col, lab) in zip(axs, [('busco_complete', 'BUSCO complete (%)'), ('log_n50', 'log10 N50'), ('partial_fraction', 'partial-model fraction')]):
        ax.scatter(s[col], s['offset'], s=10, alpha=0.6, color='#D55E00', linewidths=0, rasterized=True); ax.set_xlabel(lab); ax.set_ylabel('offset (log2)'); ax.axhline(0, color='grey', lw=0.6)
    fig.suptitle(f'V52: Saccharomycotina genus offsets vs annotation quality (n = {len(s)} genera)', fontsize=9); fig.tight_layout()
    fp = save_fig(fig, 'V52', 'sacch_quality'); tp = save_tab(tb, 'V52', 'quality_rho'); tp2 = save_tab(joint, 'V52', 'joint_ols')
    return (f"### V52 Does annotation quality explain the Saccharomycotina offset?\n\n![]({fp})\n\n- Within Saccharomycotina: {len(s)} genera ({len(gs):,} representative genomes). Spearman rho of the genus offset with each covariate:\n\n" + md_table(tb) +
            f"\n\n- RNA-seq-trained genomes (a `training/funannotate_train.pasa.gff3` or `kallisto.tsv` file exists): n = {r2['n_rnaseq']:,} genomes, median offset {S.fmt(r2['median_rnaseq'])}; not trained: n = {r2['n_no_rnaseq']:,}, median {S.fmt(r2['median_no_rnaseq'])}; difference {S.fmt(r2.get('median_diff'))} (95 % CI {S.fmt(r2.get('diff_ci_lo'))} to {S.fmt(r2.get('diff_ci_hi'))}), p = {S.fmt(r2.get('mw_p'))} (genome level, not independent).\n"
            f"- Joint OLS on genus medians, standardized covariates (offset ~ log genome size + BUSCO + log N50 + partial fraction):\n\n" + md_table(joint) +
            f"\n\n- Caveat: the offset uses complete models only, so partial-model fraction acts through which genes remain, not through truncated lengths. Tables: `{tp}`, `{tp2}`.")


# ------------------------------------------------------------------ Q2
def annot():
    s = q(f"SELECT sp_id::VARCHAR AS fam_key, entry_name, org_code, protein_name, go_terms FROM '{T}/swissprot_subject.parquet'")
    s['go_short'] = s['go_terms'].fillna('').str.split(';').str[:3].str.join(';').str.replace(r'GO:\d+ \(', '', regex=True).str.replace(')', '')
    return s


def consistency_tables():
    out = {}
    for name in ['P-A', 'P-B-YEAST', 'P-B-NEUCR', 'P-C']:
        d = q(f"SELECT * FROM '{A10}/fam_consistency_{name}.parquet'")
        d['shrinker'] = (d['n_sacch'] >= 10) & (d['sacch_frac_below_p10'] >= 0.8) & ((d['loss_aa'] >= 20) | (d['loss_frac'] >= 0.05))
        d['shrinker_p25_90'] = (d['n_sacch'] >= 10) & (d['sacch_frac_below_p25'] >= 0.9) & ((d['loss_aa'] >= 20) | (d['loss_frac'] >= 0.05))
        d['gainer'] = (d['n_sacch'] >= 10) & (d['sacch_frac_above_p90'] >= 0.8) & ((-d['loss_aa'] >= 20) | (-d['loss_frac'] >= 0.05))
        out[name] = d
    return out


def v54(cons):
    rows = []
    for name, d in cons.items():
        e = d[d['n_sacch'] >= 10]
        rows.append({'proxy': name, 'families (n_sacch >= 10)': len(e), 'median loss aa (fil - sacch)': float(e['loss_aa'].median()),
                     '>= 80 % sacch below fil p10': int((e['sacch_frac_below_p10'] >= 0.8).sum()), 'consistent shrinkers (rule)': int(e['shrinker'].sum()),
                     'shrinkers (p25 rule, >= 90 %)': int(e['shrinker_p25_90'].sum()), '>= 80 % sacch above fil p90': int((e['sacch_frac_above_p90'] >= 0.8).sum()),
                     'consistent gainers (rule)': int(e['gainer'].sum()),
                     'median sacch_frac_below_p25': float(e['sacch_frac_below_p25'].median())})
    tb = pd.DataFrame(rows); tp = save_tab(tb, 'V54', 'consistency_summary')
    d = cons['P-A']; e = d[d['n_sacch'] >= 10]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    axs[0].hist(e['sacch_frac_below_p10'], bins=20, color='#0072B2'); axs[0].set_xlabel('fraction of Saccharomycotina genera below the filamentous 10th percentile'); axs[0].set_ylabel('P-A families')
    axs[0].axvline(0.8, color='#D55E00'); axs[0].set_title(f'expected 0.10 if no shift; median {e["sacch_frac_below_p10"].median():.2f}', fontsize=8)
    axs[1].scatter(e['loss_aa'], e['sacch_frac_below_p10'], s=5, alpha=0.4, color='#0072B2', linewidths=0, rasterized=True); axs[1].axhline(0.8, color='#D55E00', lw=0.8); axs[1].axvline(20, color='#D55E00', lw=0.8)
    axs[1].set_xlabel('loss in aa (filamentous median - Saccharomycotina median)'); axs[1].set_ylabel('fraction of Sacch. genera below fil p10'); axs[1].set_xlim(-200, 300)
    fig.suptitle('V54: consistency of shortening across Saccharomycotina genera (P-A families, n_sacch >= 10)', fontsize=9); fig.tight_layout()
    fp = save_fig(fig, 'V54', 'consistency')
    return (f"### V54 How consistent is the shortening across families?\n\n![]({fp})\n\n- Per family: genus medians of the per-genome median length (complete models, strict gates, representative genomes with default filters). "
            f"Filamentous percentiles are taken over filamentous genera. Shrinker rule: {SHRINK_RULE}. Gainer rule: {GAIN_RULE}.\n\n" + md_table(tb) +
            f"\n\n- Reading: with no shift, 10 % of Saccharomycotina genera would fall below the filamentous 10th percentile. Table: `{tp}`.")


def loc_class(r):
    comps = {'N-terminal': r['d_nterm'], 'C-terminal': r['d_cterm'], 'internal': r['d_internal']}
    comps = {k: v for k, v in comps.items() if pd.notna(v)}
    if not comps:
        return 'NA'
    tot = sum(max(v, 0) for v in comps.values())
    if tot <= 0:
        return 'none'
    k, v = max(comps.items(), key=lambda kv: kv[1])
    return k if v / tot >= 0.6 else 'mixed'


def v55(cons, an):
    loc = q(f"SELECT * FROM '{A10}/top_family_location.parquet'"); dom = q(f"SELECT * FROM '{A10}/top_family_domains.parquet'")
    lost = dom[(dom['fil_frac'] >= 0.5) & (dom['sacch_frac'] <= 0.2)].groupby('fam_key')['dom'].agg(lambda s: ';'.join(sorted(s))).rename('domains_lost')
    gained = dom[(dom['sacch_frac'] >= 0.5) & (dom['fil_frac'] <= 0.2)].groupby('fam_key')['dom'].agg(lambda s: ';'.join(sorted(s))).rename('domains_gained')
    loc['location'] = loc.apply(loc_class, axis=1); loc['d_idp_res'] = loc['sacch_idp_res'] - loc['fil_idp_res']; loc['d_ndom'] = loc['sacch_ndom'] - loc['fil_ndom']
    loc = loc.merge(lost, on='fam_key', how='left').merge(gained, on='fam_key', how='left')
    d = cons['P-A']; sh = d[d['shrinker']].merge(an, on='fam_key', how='left').merge(loc, on='fam_key', how='left')
    sh = sh.sort_values(['sacch_frac_below_p10', 'loss_frac'], ascending=False)
    cols = ['fam_key', 'entry_name', 'org_code', 'protein_name', 'n_sacch', 'n_fil', 'fil_median', 'sacch_median', 'loss_aa', 'loss_frac', 'sacch_frac_below_p10', 'sacch_frac_below_p25', 'yeast_frac_below_p10',
            'n_origins_shorter', 'n_origins_tested', 'origin_deltas', 'location', 'd_nterm', 'd_cterm', 'd_internal', 'd_ndom', 'domains_lost', 'd_idp_res', 'sacch_idp_frac', 'fil_idp_frac', 'go_short']
    tp = save_tab(sh[cols], 'V55', 'consistent_shrinkers_PA', f'rule: {SHRINK_RULE}; location from swissprot_best alignment coordinates (Sacch. genus median minus filamentous genus median of the unaligned reference residues at each end and the internal gap difference)')
    for name in ['P-B-YEAST', 'P-B-NEUCR']:
        e = cons[name]; e = e[e['shrinker']].merge(an, on='fam_key', how='left').merge(loc, on='fam_key', how='left').sort_values(['sacch_frac_below_p10', 'loss_frac'], ascending=False)
        save_tab(e[[c for c in cols if c in e]], 'V55', f'consistent_shrinkers_{name}', f'rule: {SHRINK_RULE}; location columns only where the family is also a P-A family')
    pc = cons['P-C']; pc = pc[pc['shrinker']].sort_values(['sacch_frac_below_p10', 'loss_frac'], ascending=False)
    save_tab(pc[['fam_key', 'n_sacch', 'n_fil', 'fil_median', 'sacch_median', 'loss_aa', 'loss_frac', 'sacch_frac_below_p10', 'n_origins_shorter', 'n_origins_tested', 'origin_deltas']], 'V55', 'consistent_shrinkers_PC', f'rule: {SHRINK_RULE}; fam_key = Pfam architecture')
    # overlap across proxies
    setA = set(cons['P-A'][cons['P-A']['shrinker']]['fam_key']); setB = set(cons['P-B-YEAST'][cons['P-B-YEAST']['shrinker']]['fam_key'])
    # location summary
    ls = sh['location'].value_counts(); n_dl = int(sh['domains_lost'].notna().sum()); n_idp = int((sh['d_idp_res'] <= -10).sum())
    top = sh.head(10)[['entry_name', 'protein_name', 'n_sacch', 'fil_median', 'sacch_median', 'loss_aa', 'loss_frac', 'sacch_frac_below_p10', 'n_origins_shorter', 'n_origins_tested', 'location', 'd_nterm', 'd_cterm', 'd_internal', 'domains_lost', 'd_idp_res']]
    top = top.copy(); top['protein_name'] = top['protein_name'].str.slice(0, 55)
    # figure: stacked bars of loss components for top 20
    t20 = sh.head(20)
    fig, ax = plt.subplots(figsize=(8, 0.33 * len(t20) + 1.5))
    left = np.zeros(len(t20))
    for comp, colr, lab in [('d_nterm', '#0072B2', 'N-terminal'), ('d_internal', '#009E73', 'internal'), ('d_cterm', '#D55E00', 'C-terminal')]:
        v = t20[comp].fillna(0).clip(lower=0).values; ax.barh(range(len(t20)), v, left=left, color=colr, label=lab); left += v
    ax.set_yticks(range(len(t20))); ax.set_yticklabels([f"{e} ({l:.0f} aa)" for e, l in zip(t20['entry_name'], t20['loss_aa'])], fontsize=7); ax.invert_yaxis()
    ax.set_xlabel('reference residues missing in Saccharomycotina but present in filamentous copies (genus-median difference, aa)'); ax.legend(fontsize=7)
    ax.set_title('V55: where the length is lost, top 20 consistent shrinkers (P-A)', fontsize=9)
    fp = save_fig(fig, 'V55', 'loss_location')
    return (f"### V55 Consistent extreme shrinkers\n\n![]({fp})\n\n- Rule: {SHRINK_RULE}. P-A: {len(sh)} families; P-B (S. cerevisiae anchor): {int(cons['P-B-YEAST']['shrinker'].sum())}; P-B (N. crassa anchor): {int(cons['P-B-NEUCR']['shrinker'].sum())}; P-C: {len(pc)} architectures. "
            f"{len(setA & setB)} families are consistent shrinkers in both P-A and P-B (S. cerevisiae anchor).\n"
            f"- Location (P-A shrinkers, from alignment coordinates against the SwissProt reference): " + '; '.join(f"{k} {v}" for k, v in ls.items()) + f". {n_dl} of {len(sh)} families lose a Pfam domain (present in >= 50 % of filamentous genera, <= 20 % of Saccharomycotina genera); {n_idp} lose >= 10 disordered residues (idp_summary genus medians).\n"
            f"- Top 10 by consistency then relative loss:\n\n" + md_table(top, 10) +
            f"\n\n- Full ranked tables: `{tp}`, `tables/V55_consistent_shrinkers_P-B-YEAST.tsv`, `tables/V55_consistent_shrinkers_P-B-NEUCR.tsv`, `tables/V55_consistent_shrinkers_PC.tsv`. `n_origins_shorter` counts independent yeast origins (>= 3 genera) whose median is below the filamentous median.\n"
            f"- Caveats: alignment coordinates come from the best SwissProt hit, so 'internal' includes alignment gaps, not only deletions; a domain counted as lost may be present below the e-value gate; paralogs can be mixed in a family; protein names and GO terms are those of the reference, not of the fungal protein.")


def v56(cons, an):
    """Gainers and domain gain."""
    d = cons['P-A'].copy(); d = d[(d['n_sacch'] >= 10) & d['sacch_ndom'].notna() & d['fil_ndom'].notna()]
    d['class'] = np.where(d['loss_frac'] >= 0.05, 'yeast shorter (>= 5 %)', np.where(d['loss_frac'] <= -0.05, 'yeast longer (>= 5 %)', np.where(d['loss_frac'].abs() < 0.02, 'unchanged (< 2 %)', 'intermediate')))
    d['dom_gain'] = d['ndom_gain_sacch'] >= 1; d['dom_loss'] = d['ndom_gain_sacch'] <= -1; d['dom_gain_half'] = d['ndom_gain_sacch'] >= 0.5
    tb = d.groupby('class').agg(n_families=('fam_key', 'size'), share_domain_gain_ge1=('dom_gain', 'mean'), share_domain_gain_ge05=('dom_gain_half', 'mean'), share_domain_loss_ge1=('dom_loss', 'mean'),
                                median_ndom_fil=('fil_ndom', 'median'), median_ndom_sacch=('sacch_ndom', 'median')).reset_index()
    order = ['yeast longer (>= 5 %)', 'intermediate', 'unchanged (< 2 %)', 'yeast shorter (>= 5 %)']; tb['class'] = pd.Categorical(tb['class'], order); tb = tb.sort_values('class')
    tp = save_tab(tb, 'V56', 'domain_gain_by_class')
    g = cons['P-A']; ga = g[g['gainer']].merge(an, on='fam_key', how='left')
    dom = q(f"SELECT * FROM '{A10}/top_family_domains.parquet'"); loc = q(f"SELECT fam_key, d_nterm, d_cterm, d_internal, sacch_ndom AS sacch_ndom_top, fil_ndom AS fil_ndom_top, sacch_idp_res, fil_idp_res FROM '{A10}/top_family_location.parquet'")
    gained = dom[(dom['sacch_frac'] >= 0.5) & (dom['fil_frac'] <= 0.2)].groupby('fam_key')['dom'].agg(lambda s: ';'.join(sorted(s))).rename('domains_gained')
    ga = ga.merge(gained, on='fam_key', how='left').merge(loc, on='fam_key', how='left'); ga['gain_aa'] = -ga['loss_aa']; ga['gain_frac'] = -ga['loss_frac']
    ga['location'] = ga.apply(lambda r: loc_class({'d_nterm': -r['d_nterm'] if pd.notna(r['d_nterm']) else np.nan, 'd_cterm': -r['d_cterm'] if pd.notna(r['d_cterm']) else np.nan, 'd_internal': -r['d_internal'] if pd.notna(r['d_internal']) else np.nan}), axis=1)
    ga['d_idp_res'] = ga['sacch_idp_res'] - ga['fil_idp_res']
    ga = ga.sort_values(['sacch_frac_above_p90', 'gain_frac'], ascending=False)
    cols = ['fam_key', 'entry_name', 'org_code', 'protein_name', 'n_sacch', 'n_fil', 'fil_median', 'sacch_median', 'gain_aa', 'gain_frac', 'sacch_frac_above_p90', 'sacch_frac_above_p75', 'n_origins_longer', 'n_origins_tested', 'origin_deltas',
            'ndom_gain_sacch', 'domains_gained', 'location', 'd_idp_res', 'go_short']
    tp2 = save_tab(ga[cols], 'V56', 'consistent_gainers_PA', f'rule: {GAIN_RULE}')
    for name in ['P-B-YEAST', 'P-B-NEUCR']:
        e = cons[name]; e = e[e['gainer']].merge(an, on='fam_key', how='left').sort_values(['sacch_frac_above_p90', 'loss_frac'])
        save_tab(e[[c for c in ['fam_key', 'entry_name', 'org_code', 'protein_name', 'n_sacch', 'n_fil', 'fil_median', 'sacch_median', 'loss_aa', 'loss_frac', 'sacch_frac_above_p90', 'n_origins_longer', 'n_origins_tested', 'origin_deltas'] if c in e]], 'V56', f'consistent_gainers_{name}', f'rule: {GAIN_RULE}')
    top = ga.head(10)[['entry_name', 'protein_name', 'n_sacch', 'fil_median', 'sacch_median', 'gain_aa', 'gain_frac', 'sacch_frac_above_p90', 'n_origins_longer', 'n_origins_tested', 'ndom_gain_sacch', 'domains_gained', 'location', 'd_idp_res']].copy()
    top['protein_name'] = top['protein_name'].str.slice(0, 55)
    share_g = float(ga['ndom_gain_sacch'].ge(1).mean()) if len(ga) else np.nan; n_dg = int(ga['domains_gained'].notna().sum())
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.bar(range(len(tb)), tb['share_domain_gain_ge1'], color='#0072B2', label='median domain count +1 or more in Saccharomycotina'); ax.bar(range(len(tb)), -tb['share_domain_loss_ge1'], color='#D55E00', label='-1 or less')
    ax.set_xticks(range(len(tb))); ax.set_xticklabels([f"{c}\n(n={n})" for c, n in zip(tb['class'], tb['n_families'])], fontsize=7); ax.axhline(0, color='black', lw=0.6); ax.set_ylabel('share of families'); ax.legend(fontsize=7)
    ax.set_title('V56: Pfam domain-count change by family length class (P-A)', fontsize=9)
    fp = save_fig(fig, 'V56', 'domain_gain_by_class')
    return (f"### V56 Families where Saccharomycotina proteins are longer: is it extra domains?\n\n![]({fp})\n\n- P-A families with n_sacch >= 10 and domain counts on both sides ({len(d):,}); class by relative length change of the Saccharomycotina median vs the filamentous median; domain change = Saccharomycotina genus-median Pfam domain count minus filamentous.\n\n" + md_table(tb) +
            f"\n\n- Consistent gainers (rule: {GAIN_RULE}): {len(ga)} P-A families; {int(cons['P-B-YEAST']['gainer'].sum())} P-B (S. cerevisiae anchor); {int(cons['P-B-NEUCR']['gainer'].sum())} P-B (N. crassa anchor); {int(cons['P-C']['gainer'].sum())} P-C architectures. "
            f"Among the P-A gainers, {100*share_g:.0f} % have a median domain count higher by >= 1, and {n_dg} carry a Pfam domain present in >= 50 % of Saccharomycotina genera but <= 20 % of filamentous genera.\n- Top 10 consistent gainers:\n\n" + md_table(top, 10) +
            f"\n\n- Tables: `{tp}`, `{tp2}`, `tables/V56_consistent_gainers_P-B-YEAST.tsv`, `tables/V56_consistent_gainers_P-B-NEUCR.tsv`. Location for gainers = where the extra residues sit relative to the reference (sign reversed). Caveats as V55.")


def main():
    g, agg = genus_frame(); an = annot(); cons = consistency_tables()
    blocks = []
    for fn, args in [(v50, (g, agg)), (v51, (agg,)), (v52, (g, agg)), (v54, (cons,)), (v55, (cons, an)), (v56, (cons, an))]:
        print('render', fn.__name__, flush=True)
        try:
            blocks.append(fn(*args))
        except Exception as e:
            import traceback; blocks.append(f"### {fn.__name__.upper()}\n\nRender failed: {e}\n```\n{traceback.format_exc()[-1500:]}\n```"); print('FAILED', fn.__name__, e, flush=True)
    head = ("# H2b. Genome compaction, consistent shrinkers and gainers\n\nFollow-up to [H2](H2_protein_length.md) and the [H2 report](H2_report.md). Questions: (1) is the Saccharomycotina ortholog shortening explained by genome compaction, within clades and not only between them, "
            "and does annotation quality explain it; (2) which families shrink consistently across Saccharomycotina genera, where in the protein, and what is lost; (3) which families are consistently longer, and is that extra domains.\n\n"
            "[Back to index](index.md)\n\n## Q1. Genome compaction\n\n")
    txt = head + '\n\n'.join(blocks[:3]) + "\n\n## Q2. Consistent extreme shrinkers\n\n" + '\n\n'.join(blocks[3:5]) + "\n\n## Q3. Consistent gainers and domain gain\n\n" + blocks[5] + "\n\n## Verdict\n\n" + VERDICT + '\n'
    with open(f'{DOCS}/H2b_compaction_extremes.md', 'w') as o:
        o.write(txt)
    print('wrote page', flush=True)


VERDICT = os.environ.get('H2B_VERDICT', '_Verdict written after the numbers are reviewed._')

if __name__ == '__main__':
    main()
