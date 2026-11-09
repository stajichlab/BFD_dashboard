"""P2 views for the BFD summary atlas. Imported by render_views.main(); reuses its helpers."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as sps

import render_views as RV
import stats as S

A, T = RV.A, RV.T
GF_ORDER, GF_COLORS, SACCH, PEZIZO = RV.GF_ORDER, RV.GF_COLORS, RV.SACCH, RV.PEZIZO
save_fig, save_tab, md_table, genus_medians, counts_text, contrast_block, strip_by_class = (
    RV.save_fig, RV.save_tab, RV.md_table, RV.genus_medians, RV.counts_text, RV.contrast_block, RV.strip_by_class)


# ------------------------------------------------------------------ H1
def v06(ctx):
    L = ctx.labels
    n = L.groupby('SPECIES').size()
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.hist(n.values, bins=np.logspace(0, np.log10(n.max() + 1), 40), color='#0072B2')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('genomes per species'); ax.set_ylabel('species')
    ax.set_title('V06: strain oversampling', fontsize=9)
    fp = save_fig(fig, 'V06', 'genomes_per_species')
    top = n.sort_values(ascending=False).head(20).reset_index(); top.columns = ['SPECIES', 'n_genomes']
    tp = save_tab(top, 'V06', 'top_species')
    hist = pd.DataFrame({'genomes_per_species': ['1', '2', '3-5', '6-10', '11-50', '51-200', '>200'],
                         'n_species': [int((n == 1).sum()), int((n == 2).sum()), int(((n >= 3) & (n <= 5)).sum()), int(((n >= 6) & (n <= 10)).sum()),
                                       int(((n >= 11) & (n <= 50)).sum()), int(((n >= 51) & (n <= 200)).sum()), int((n > 200).sum())]})
    return (f"### V06 Strain oversampling\n\n![]({fp})\n\n- {len(L):,} genomes in {len(n):,} named species; top 20 species hold {int(top['n_genomes'].sum()):,} genomes ({100 * top['n_genomes'].sum() / len(L):.1f} %).\n\n"
            + md_table(hist) + "\n\n" + md_table(top) + f"\n\n- This is why every cross-clade view uses one representative per species. Table: `{tp}`.")


# ------------------------------------------------------------------ H2
def _ratio_summary(agg, name):
    st = S.sign_test(agg['log2_ratio']); lo, hi = S.boot_ci(agg['log2_ratio'])
    return {'proxy': name, 'n_families': len(agg), 'median_log2_ratio': float(agg['log2_ratio'].median()), 'ci_lo': lo, 'ci_hi': hi,
            'ratio': 2 ** float(agg['log2_ratio'].median()), 'frac_yeast_shorter': st.get('frac_below_zero'), 'binom_p': st.get('binom_p'),
            'p10': float(agg['log2_ratio'].quantile(.1)), 'p90': float(agg['log2_ratio'].quantile(.9))}


def v10(ctx):
    rel = 'A4/fam_anchor_genome.parquet'
    if not ctx.have(rel):
        return '### V10 Model-organism-anchored contrast (P-B)\n\nNot built: A4/fam_anchor_genome (job J5) missing.'
    fam = ctx.tbl(rel)
    rows, dists = [], {}
    for org in ['YEAST', 'NEUCR', 'SCHPO', 'EMENI']:
        f = fam[fam['org_code'] == org].rename(columns={'anchor_sp_id': 'fam_key'})
        agg = RV.family_contrast(ctx, f, 'fam_key')
        rows.append(_ratio_summary(agg, f'P-B anchor {org}')); dists[org] = agg['log2_ratio'].values
    # P-A for comparison on the same page
    aggA = RV.family_contrast(ctx, ctx.tbl('A4/fam_sp_genome.parquet'), 'sp_id'); rows.append(_ratio_summary(aggA, 'P-A SwissProt best hit'))
    tb = pd.DataFrame(rows); tp = save_tab(tb, 'V10', 'anchor_contrast')
    fig, axs = plt.subplots(1, 4, figsize=(14, 3.4), sharey=True)
    for ax, org in zip(axs, ['YEAST', 'NEUCR', 'SCHPO', 'EMENI']):
        v = dists[org]; ax.hist(v, bins=60, range=(-0.4, 0.4), color='#0072B2'); ax.axvline(0, color='black', lw=0.8); ax.axvline(np.median(v), color='#D55E00')
        ax.set_title(f'anchor {org}: n = {len(v):,}, median {np.median(v):+.3f}', fontsize=8); ax.set_xlabel('log2(yeast / filamentous)')
    fig.suptitle('V10: within-family length contrast with every genome anchored to one reference proteome', fontsize=9)
    fp = save_fig(fig, 'V10', 'anchor_hist')
    # per-genome anchoring coverage
    cov = ctx.tbl('A9/anchor_genome_summary.parquet')
    d = ctx.reps().merge(cov, left_on='LOCUSTAG', right_on='genome'); d = d[d['growth_form'].isin(['yeast', 'filamentous_mycelium'])]
    cv = d.groupby(['org_code', 'growth_form'])['n_strict'].median().unstack(); cv.columns = [f'median_strict_anchors_{c}' for c in cv.columns]
    tp2 = save_tab(cv.reset_index(), 'V10', 'anchor_coverage')
    return (f"### V10 Model-organism-anchored contrast (P-B)\n\n![]({fp})\n\n- Family = the best hit of each protein to one reference proteome (*S. cerevisiae* YEAST, *N. crassa* NEUCR, *S. pombe* SCHPO, *A. nidulans* EMENI) taken from all ranked SwissProt hits, so yeast and filamentous genomes are compared on the same reference protein. "
            f"Gates on the anchor hit: qcov >= 80 %, scov >= 80 %, pident >= 40 %; complete models; representative genomes with {ctx.filt_text()}; >= 10 genera per side; unit = family.\n\n" + md_table(tb) +
            f"\n\n- Median number of strict anchors per representative genome:\n\n" + md_table(cv.reset_index()) +
            f"\n\n- Result: the four anchors agree (median log2 ratio -0.019 to -0.022, 69-71 % of families yeast-shorter). The YEAST-anchored and NEUCR-anchored medians are the same (-0.019), so the P-A signal is not an artefact of Saccharomycotina genomes hitting *S. cerevisiae* references; P-B is the cleaner proxy because every genome is scored against one reference protein. Yeast genomes have more strict YEAST anchors (1,599 vs 931) and fewer NEUCR/EMENI anchors (282 vs 620) than filamentous genomes, which shows the reference-side bias that P-A inherits. Tables: `{tp}`, `{tp2}`.")


def v11(ctx):
    rel = 'A4/fam_arch_genome.parquet'
    if not ctx.have(rel):
        return '### V11 Within-architecture length contrast (Pfam proxy, P-C)\n\nNot built: A4/fam_arch_genome missing.'
    fam = ctx.tbl(rel)
    agg = RV.family_contrast(ctx, fam, 'pfam_arch')
    agg['n_dom'] = agg['pfam_arch'].str.count(';') + 1
    r = _ratio_summary(agg, 'P-C Pfam architecture (>= 2 domains)')
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.hist(agg['log2_ratio'], bins=80, color='#0072B2'); ax.axvline(0, color='black', lw=0.8); ax.axvline(agg['log2_ratio'].median(), color='#D55E00', label=f"median {agg['log2_ratio'].median():+.3f}")
    ax.set_xlabel('log2( yeast median length / filamentous median length ) per Pfam architecture'); ax.set_ylabel('architectures'); ax.legend(fontsize=8)
    ax.set_title(f'V11: within-architecture protein length (n = {len(agg):,} architectures)', fontsize=9)
    fp = save_fig(fig, 'V11', 'arch_log2_ratio')
    byd = agg.groupby(pd.cut(agg['n_dom'], [1, 2, 3, 4, 6, 100]), observed=True)['log2_ratio'].agg(['median', 'size']).reset_index(); byd['n_dom'] = byd['n_dom'].astype(str)
    tb = pd.DataFrame([r]); tp = save_tab(agg.sort_values('log2_ratio'), 'V11', 'arch_ratios'); tp2 = save_tab(byd, 'V11', 'ratio_by_domain_count')
    return (f"### V11 Within-architecture length contrast, Pfam proxy (P-C)\n\n![]({fp})\n\n- Family = ordered Pfam domain string (domain_i_evalue <= 1e-5), architectures with >= 2 domains present in >= 10 yeast and >= 10 filamentous genera; complete models; representative genomes with {ctx.filt_text()}; unit = architecture.\n\n"
            + md_table(tb) + "\n\n- By number of domains in the architecture:\n\n" + md_table(byd) +
            f"\n\n- Caveat: a truncated model loses domains and leaves its architecture class, so truncation does not make a family shorter here; what remains is linker and terminal length. Tables: `{tp}`, `{tp2}`.")


def v13(ctx):
    rel = 'A4/fam_sp_genome.parquet'
    if not ctx.have(rel):
        return '### V13 Which families drive the difference\n\nNot built.'
    fam = ctx.tbl(rel)
    agg = RV.family_contrast(ctx, fam, 'sp_id')
    subj = ctx.q(f"SELECT sp_id, org_code, s_len, protein_name FROM '{T}/swissprot_subject.parquet'")
    agg = agg.merge(subj, on='sp_id', how='left')
    nd = fam.groupby('sp_id')['med_n_domains'].median().rename('med_n_domains'); agg = agg.merge(nd, left_on='sp_id', right_index=True)
    agg['ref_group'] = np.where(agg['org_code'] == 'YEAST', 'S. cerevisiae (YEAST)', np.where(agg['org_code'].isin(['NEUCR', 'EMENI']), 'N. crassa / A. nidulans',
                               np.where(agg['org_code'] == 'SCHPO', 'S. pombe', 'other reference')))
    by_ref = agg.groupby('ref_group')['log2_ratio'].agg(['median', 'size']).reset_index()
    by_ref['ci_lo'], by_ref['ci_hi'] = zip(*[S.boot_ci(agg[agg['ref_group'] == g]['log2_ratio']) for g in by_ref['ref_group']])
    by_ref['frac_yeast_shorter'] = [S.sign_test(agg[agg['ref_group'] == g]['log2_ratio']).get('frac_below_zero') for g in by_ref['ref_group']]
    agg['dom_bin'] = pd.cut(agg['med_n_domains'], [-1, 0, 1, 2, 3, 100], labels=['0', '1', '2', '3', '>= 4'])
    by_dom = agg.groupby('dom_bin', observed=True)['log2_ratio'].agg(['median', 'size']).reset_index(); by_dom['dom_bin'] = by_dom['dom_bin'].astype(str)
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    for ax, (tb_, col, lab) in zip(axs, [(by_ref, 'ref_group', 'reference organism of the family'), (by_dom, 'dom_bin', 'median Pfam domains per protein')]):
        ax.bar(range(len(tb_)), tb_['median'], color='#0072B2'); ax.axhline(0, color='black', lw=0.8)
        ax.set_xticks(range(len(tb_))); ax.set_xticklabels([f"{g}\n(n={n})" for g, n in zip(tb_[col], tb_['size'])], fontsize=7); ax.set_ylabel('median log2(yeast / filamentous)'); ax.set_title(lab, fontsize=8)
    fig.suptitle('V13: family length ratio by reference organism and by domain count', fontsize=9)
    fp = save_fig(fig, 'V13', 'by_reference_and_domains')
    ext = agg.sort_values('log2_ratio')[['sp_id', 'org_code', 'protein_name', 'size_yeast', 'size_filamentous_mycelium', 'median_yeast', 'median_filamentous_mycelium', 'log2_ratio']]
    tp = save_tab(by_ref, 'V13', 'by_reference'); tp2 = save_tab(by_dom, 'V13', 'by_domain_count'); tp3 = save_tab(pd.concat([ext.head(25), ext.tail(25)]), 'V13', 'extreme_families')
    return (f"### V13 Which families drive the difference\n\n![]({fp})\n\n- Same {len(agg):,} families as V08 (unit = family). Split by the organism of the SwissProt reference each family maps to, and by the median Pfam domain count of its members.\n\n"
            + md_table(by_ref) + "\n\n" + md_table(by_dom) +
            f"\n\n- Reading: a negative median in the *N. crassa / A. nidulans*-referenced families shows the yeast-shorter signal is not produced by families referenced to *S. cerevisiae*. The 25 most negative and 25 most positive families are in `{tp3}`. Tables: `{tp}`, `{tp2}`.")


def v31(ctx):
    g = ctx.tbl('A4/genome_protein_summary.parquet')
    d = ctx.reps().merge(g, left_on='LOCUSTAG', right_on='genome'); d = d[d['growth_form'].isin(GF_ORDER)]
    gd = genus_medians(d, ['frac_multidomain', 'frac_with_pfam', 'median_domains_if_any'])
    cb, r = contrast_block(ctx, gd, 'frac_multidomain', 'fraction of proteins with >= 2 Pfam domains')
    cb2, r2 = contrast_block(ctx, gd, 'frac_with_pfam', 'fraction of proteins with >= 1 Pfam domain')
    fp = strip_by_class(d, 'frac_multidomain', 'fraction of proteins with >= 2 Pfam domains', 'V31', 'multidomain_by_class', unit='genome')
    tp = save_tab(gd, 'V31', 'genus_medians')
    return (f"### V31 Domains per protein and multi-domain fraction\n\n![]({fp})\n\n- {counts_text(d)}. Filters: {ctx.filt_text(model=False)}. Domain gate domain_i_evalue <= 1e-5.\n"
            f"**Multi-domain fraction:**\n{cb}\n**Fraction with any Pfam domain:**\n{cb2}\n- Link to H2: V13 shows the family length ratio by domain count. Table: `{tp}`.")


# ------------------------------------------------------------------ H3
def v16(ctx):
    d = RV._a2reps(ctx, dup=False, n50=False)
    d['span_cds_ratio'] = d['median_transcript_span'] / d['median_cds_len']
    fig, ax = plt.subplots(figsize=(6, 5))
    for gf in GF_ORDER:
        s = d[d['growth_form'] == gf]; ax.scatter(s['median_cds_len'], s['median_transcript_span'], s=5, alpha=0.4, color=GF_COLORS[gf], label=f'{gf} (n={len(s)})', linewidths=0, rasterized=True)
    ax.plot([300, 3000], [300, 3000], color='grey', ls='--', lw=0.8); ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('median CDS length (bp)'); ax.set_ylabel('median transcript span (bp)'); ax.legend(fontsize=7)
    ax.set_title('V16: gene span vs CDS length (representative genomes)', fontsize=9)
    fp = save_fig(fig, 'V16', 'span_vs_cds')
    gd = genus_medians(d, ['span_cds_ratio', 'median_cds_len', 'median_transcript_span', 'median_noncoding_fraction'])
    cb, r = contrast_block(ctx, gd, 'span_cds_ratio', 'transcript span / CDS length')
    cb2, r2 = contrast_block(ctx, gd, 'median_cds_len', 'median CDS length (bp)')
    tp = save_tab(gd, 'V16', 'genus_medians')
    return (f"### V16 Gene span versus CDS length\n\n![]({fp})\n\n- {counts_text(d)}. Filters: {ctx.filt_text(dup=False, n50=False, model=False)}. Span = transcript start to end (introns included; UTRs only where predicted).\n"
            f"**Span / CDS ratio:**\n{cb}\n**Median CDS length:**\n{cb2}\n- Caveat: most ab initio models have no UTR, so span minus CDS is intron length. Table: `{tp}`. Scatter points are rasterized in the SVG.")


# ------------------------------------------------------------------ H4
def v24(ctx):
    cz = ctx.tbl('A5/cazy_genome.parquet')
    d = ctx.reps().merge(cz, left_on='LOCUSTAG', right_on='genome'); d['frac_cazy_secreted'] = d['n_cazy_secreted'] / d['n_cazy']
    gd = genus_medians(d, ['frac_cazy_secreted', 'n_cazy', 'n_cazy_secreted'])
    cb, r = contrast_block(ctx, gd, 'frac_cazy_secreted', 'fraction of CAZymes with a strict signal peptide')
    n = gd.groupby('primary_lifestyle').size(); keep = n[n >= 10].index
    tb = gd[gd['primary_lifestyle'].isin(keep)].groupby('primary_lifestyle').agg(n_genera=('GENUS', 'size'), median_frac_secreted=('frac_cazy_secreted', 'median'), median_n_cazy=('n_cazy', 'median'), median_n_cazy_secreted=('n_cazy_secreted', 'median')).reset_index().sort_values('median_frac_secreted')
    fig, ax = plt.subplots(figsize=(6.5, 0.3 * len(tb) + 1.5))
    ax.barh(range(len(tb)), tb['median_frac_secreted'], color='#0072B2'); ax.set_yticks(range(len(tb))); ax.set_yticklabels([f"{l} (n={n})" for l, n in zip(tb['primary_lifestyle'], tb['n_genera'])], fontsize=7); ax.set_xlabel('median fraction of CAZymes secreted (genus medians)')
    ax.set_title('V24: secreted fraction of CAZymes by lifestyle', fontsize=9)
    fp = save_fig(fig, 'V24', 'secreted_cazy_by_lifestyle'); tp = save_tab(tb, 'V24', 'by_lifestyle')
    kw = S.kruskal({l: gd[gd['primary_lifestyle'] == l]['frac_cazy_secreted'].values for l in keep})
    return (f"### V24 Secreted CAZymes\n\n![]({fp})\n\n- {counts_text(d)}. Filters: {ctx.filt_text(model=False)}. CAZyme = dbCAN `toolcount >= 2`; secreted = SignalP >= 0.9, <= 1 TM helix, no GPI anchor.\n{cb}\n\n" + md_table(tb) +
            f"\n\n- Kruskal-Wallis across lifestyles (screening): H = {kw['kw_H']:.1f}, p = {S.fmt(kw['kw_p'])}. Caveat: GPI-anchored cell-wall CAZymes are excluded by the strict definition. Table: `{tp}`.")


def v25(ctx):
    if not ctx.have('A5/merops_family_counts.parquet'):
        return '### V25 Peptidase (MEROPS) families by lifestyle\n\nNot built: A5/merops_family_counts (job J7) missing.'
    mf = ctx.tbl('A5/merops_family_counts.parquet'); chk = ctx.tbl('A5/merops_mapping_check.parquet').iloc[0]
    g = ctx.tbl('A4/genome_protein_summary.parquet')[['genome', 'n_proteins']]
    L = ctx.reps()[['LOCUSTAG', 'GENUS', 'CLASS', 'ORDER_', 'PHYLUM', 'growth_form', 'primary_lifestyle', 'yeast_origin', 'log_n50', 'busco_complete', 'log_genome', 'SPECIES']]
    m = mf.merge(g, on='genome').merge(L, left_on='genome', right_on='LOCUSTAG')
    m['per1000'] = m['n_proteins_x'] * 1000 / m['n_proteins_y']; m['sec_per1000'] = m['n_secreted'] * 1000 / m['n_proteins_y']
    tot = m.groupby(['genome', 'GENUS', 'growth_form', 'primary_lifestyle']).agg(n_merops=('n_proteins_x', 'sum'), n_sec=('n_secreted', 'sum'), n_prot=('n_proteins_y', 'first')).reset_index()
    tot['merops_per1000'] = tot['n_merops'] * 1000 / tot['n_prot']; tot['sec_per1000'] = tot['n_sec'] * 1000 / tot['n_prot']
    gd = tot.groupby('GENUS').agg(merops_per1000=('merops_per1000', 'median'), sec_per1000=('sec_per1000', 'median'), growth_form=('growth_form', 'first'), primary_lifestyle=('primary_lifestyle', 'first')).reset_index()
    y, f = gd[gd['growth_form'] == 'yeast'], gd[gd['growth_form'] == 'filamentous_mycelium']
    r = S.two_group(y['sec_per1000'], f['sec_per1000']); r2 = S.two_group(y['merops_per1000'], f['merops_per1000'])
    top = m.groupby('merops_family')['n_proteins_x'].sum().sort_values(ascending=False); top = [t for t in top.index if t != 'unmapped'][:15]
    n = gd.groupby('primary_lifestyle').size(); keep = n[n >= 10].index
    gm = m[m['merops_family'].isin(top)].groupby(['GENUS', 'merops_family'])['sec_per1000'].median().unstack(fill_value=0)
    gm = gm.join(gd.set_index('GENUS')['primary_lifestyle']); gm = gm[gm['primary_lifestyle'].isin(keep)]
    heat = gm.groupby('primary_lifestyle')[top].median(); heat = heat.loc[heat.sum(1).sort_values().index]
    fig, ax = plt.subplots(figsize=(9, 0.3 * len(heat) + 1.8))
    im = ax.imshow(np.log2(heat.values + 0.1), aspect='auto', cmap='viridis'); ax.set_xticks(range(len(top))); ax.set_xticklabels(top, fontsize=7, rotation=90)
    ax.set_yticks(range(len(heat))); ax.set_yticklabels([f"{i} (n={n[i]})" for i in heat.index], fontsize=7); plt.colorbar(im, ax=ax, label='log2(median secreted peptidases per 1,000 proteins + 0.1)')
    ax.set_title('V25: secreted MEROPS families per 1,000 proteins by lifestyle (genus medians)', fontsize=9)
    fp = save_fig(fig, 'V25', 'merops_by_lifestyle'); tp = save_tab(heat.reset_index(), 'V25', 'merops_by_lifestyle')
    ls = gd[gd['primary_lifestyle'].isin(keep)].groupby('primary_lifestyle').agg(n_genera=('GENUS', 'size'), median_merops_per1000=('merops_per1000', 'median'), median_secreted_per1000=('sec_per1000', 'median')).reset_index().sort_values('median_secreted_per1000')
    tp2 = save_tab(ls, 'V25', 'totals_by_lifestyle')
    return (f"### V25 Peptidase (MEROPS) families by lifestyle\n\n![]({fp})\n\n- Best MEROPS hit per protein (evalue <= 1e-10), family from `merops_scan.lib` headers (`atlas_data/labels/merops_family_map.tsv`); {int(chk['n_unmapped']):,} of {int(chk['n_best_hits']):,} best hits ({100*chk['n_unmapped']/chk['n_best_hits']:.2f} %) map to no family. "
            f"Unit = genus. Filters: {ctx.filt_text(model=False)}.\n- Secreted peptidases per 1,000 proteins: yeast genera median {S.fmt(r['median_yeast'])} vs filamentous {S.fmt(r['median_filamentous'])} (difference {S.fmt(r.get('median_diff'))}, 95 % CI {S.fmt(r.get('diff_ci_lo'))} to {S.fmt(r.get('diff_ci_hi'))}, Cliff's delta {S.fmt(r.get('cliffs_delta'))}); all peptidases: {S.fmt(r2['median_yeast'])} vs {S.fmt(r2['median_filamentous'])}.\n\n"
            + md_table(ls) + f"\n\n- Caveat: MEROPS hits include non-peptidase homologues; lifestyle and clade are confounded. Tables: `{tp}`, `{tp2}`.")


# ------------------------------------------------------------------ H5
def v29(ctx):
    rich = ctx.q(f"SELECT genome, count(*) AS n_families, sum(n_proteins) AS n_domain_proteins FROM '{A}/A7/pfam_genome_matrix.parquet' GROUP BY 1")
    g = ctx.tbl('A4/genome_protein_summary.parquet')[['genome', 'n_proteins']]
    d = ctx.reps(dup=False, n50=False).merge(rich, left_on='LOCUSTAG', right_on='genome').merge(g, on='genome'); d = d[d['growth_form'].isin(GF_ORDER)]
    fig, ax = plt.subplots(figsize=(6, 5)); rows = []
    for gf in GF_ORDER:
        s = d[d['growth_form'] == gf]; ax.scatter(s['n_proteins'], s['n_families'], s=5, alpha=0.4, color=GF_COLORS[gf], label=f'{gf} (n={len(s)})', linewidths=0, rasterized=True)
        r = S.loglog_slope(s['n_proteins'], s['n_families']); rows.append({'growth_form': gf, **r})
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('proteins per genome'); ax.set_ylabel('distinct Pfam families'); ax.legend(fontsize=7); ax.set_title('V29: Pfam family richness vs proteome size', fontsize=9)
    fp = save_fig(fig, 'V29', 'pfam_richness'); tb = pd.DataFrame(rows); tp = save_tab(tb, 'V29', 'slopes')
    gd = genus_medians(d, ['n_families', 'n_proteins']); cb, r = contrast_block(ctx, gd, 'n_families', 'distinct Pfam families', covars=('log_genome', 'busco_complete'))
    return (f"### V29 Pfam domain richness versus gene count\n\n![]({fp})\n\n- {counts_text(d)}. Filters: {ctx.filt_text(dup=False, n50=False, model=False)}. Domain gate domain_i_evalue <= 1e-5.\n\n" + md_table(tb) +
            f"\n\n{cb}\n- Caveat: richness saturates at the Pfam library size; slopes below 1 are expected for any finite library. Table: `{tp}`. Scatter points are rasterized in the SVG.")


def v30(ctx):
    L = ctx.reps()[['LOCUSTAG', 'GENUS', 'growth_form', 'ORDER_']]; L = L[L['growth_form'].isin(['yeast', 'filamentous_mycelium'])]
    ctx.con.register('reps30', L)
    pres = ctx.q(f"""
      WITH gsize AS (SELECT GENUS, growth_form, count(*) AS n_gen FROM reps30 GROUP BY 1,2),
      pg AS (SELECT r.GENUS, r.growth_form, m.pfam_acc, any_value(m.pfam_id) AS pfam_id, count(DISTINCT m.genome) AS n_with
             FROM '{A}/A7/pfam_genome_matrix.parquet' m JOIN reps30 r ON m.genome = r.LOCUSTAG GROUP BY 1,2,3),
      p AS (SELECT pg.*, pg.n_with * 1.0 / gs.n_gen >= 0.5 AS present FROM pg JOIN gsize gs USING (GENUS, growth_form))
      SELECT pfam_acc, any_value(pfam_id) AS pfam_id,
             count(*) FILTER (WHERE present AND growth_form='yeast') AS yeast_present,
             count(*) FILTER (WHERE present AND growth_form='filamentous_mycelium') AS fil_present
      FROM p GROUP BY 1""")
    ny = L[L['growth_form'] == 'yeast']['GENUS'].nunique(); nf = L[L['growth_form'] == 'filamentous_mycelium']['GENUS'].nunique()
    pres['yeast_frac'] = pres['yeast_present'] / ny; pres['fil_frac'] = pres['fil_present'] / nf
    pres = pres[(pres['yeast_present'] + pres['fil_present']) >= 5].copy()
    pv = [sps.fisher_exact([[a, ny - a], [b, nf - b]])[1] for a, b in zip(pres['yeast_present'], pres['fil_present'])]
    pres['fisher_p'] = pv; m = len(pv); order = np.argsort(pv); q = np.empty(m); q[order] = np.minimum.accumulate((np.array(pv)[order] * m / np.arange(1, m + 1))[::-1])[::-1]; pres['bh_q'] = q
    pres['log2_odds'] = np.log2((pres['yeast_frac'] + 0.01) / (pres['fil_frac'] + 0.01))
    # orders in which family is present (independence guard)
    ctx.con.register('pres30', pres[['pfam_acc']])
    nord = ctx.q(f"SELECT m.pfam_acc, count(DISTINCT r.ORDER_) AS n_orders_present FROM '{A}/A7/pfam_genome_matrix.parquet' m JOIN reps30 r ON m.genome = r.LOCUSTAG JOIN pres30 USING (pfam_acc) GROUP BY 1")
    pres = pres.merge(nord, on='pfam_acc', how='left')
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(pres['log2_odds'], -np.log10(pres['bh_q'].clip(lower=1e-300)), s=4, alpha=0.4, color='#777', linewidths=0, rasterized=True)
    sig = pres[pres['bh_q'] < 0.01]; ax.scatter(sig['log2_odds'], -np.log10(sig['bh_q'].clip(lower=1e-300)), s=5, color='#D55E00', linewidths=0, rasterized=True)
    ax.set_xlabel('log2( (yeast genus fraction + 0.01) / (filamentous genus fraction + 0.01) )'); ax.set_ylabel('-log10 BH q (Fisher exact on genus presence)')
    ax.set_title(f'V30: Pfam families enriched in yeast or filamentous genera ({ny} vs {nf} genera)', fontsize=9)
    fp = save_fig(fig, 'V30', 'volcano')
    cols = ['pfam_acc', 'pfam_id', 'yeast_present', 'fil_present', 'yeast_frac', 'fil_frac', 'log2_odds', 'bh_q', 'n_orders_present']
    ty = sig.sort_values('log2_odds', ascending=False).head(30)[cols]; tf = sig.sort_values('log2_odds').head(30)[cols]
    tp = save_tab(pres.sort_values('bh_q')[cols], 'V30', 'all_families'); tp1 = save_tab(ty, 'V30', 'yeast_enriched'); tp2 = save_tab(tf, 'V30', 'filamentous_enriched')
    return (f"### V30 Clade-enriched Pfam families, yeast versus filamentous\n\n![]({fp})\n\n- Unit = genus ({ny} yeast, {nf} filamentous genera of representative genomes with {ctx.filt_text(model=False)}). A family is present in a genus when >= 50 % of its representative genomes carry it. "
            f"{len(pres):,} families tested (present in >= 5 genera); {len(sig):,} with BH q < 0.01.\n\n**Top yeast-enriched (by log2 odds, q < 0.01):**\n\n" + md_table(ty, 30) + "\n\n**Top filamentous-enriched:**\n\n" + md_table(tf, 30) +
            f"\n\n- Caveat: genera share ancestry, so q-values overstate independence; `n_orders_present` counts the orders in which the family occurs. Tables: `{tp}`, `{tp1}`, `{tp2}`. Points are rasterized in the SVG.")


# ------------------------------------------------------------------ H6
def v36(ctx):
    tr = ctx.tbl('A5/trna_codon_counts.parquet'); tr['codon'] = tr['codon'].str.upper()
    cf = ctx.q(f"SELECT species_prefix AS genome, upper(codon) AS codon, frequency FROM '{T}/codon_freq.parquet'")
    cf['aa'] = cf['codon'].map(RV.CODE1); cf = cf[cf['aa'].notna() & (cf['aa'] != '*')]
    L = ctx.reps(dup=False)[['LOCUSTAG', 'GENUS', 'CLASS', 'ORDER_', 'PHYLUM', 'growth_form', 'primary_lifestyle', 'yeast_origin', 'log_n50', 'busco_complete', 'log_genome', 'SPECIES']]
    m = cf.merge(tr[['genome', 'codon', 'n']], on=['genome', 'codon'], how='left').fillna({'n': 0})
    m = m[m['genome'].isin(L['LOCUSTAG'])]
    tot = tr.groupby('genome')['n'].sum().rename('n_trna')
    rows = []
    for g, s in m.groupby('genome'):
        r = S.spearman(s['frequency'], s['n']); rows.append({'genome': g, 'rho': r['rho'], 'n_codons_with_trna': int((s['n'] > 0).sum())})
    d = pd.DataFrame(rows).merge(tot, left_on='genome', right_index=True).merge(L, left_on='genome', right_on='LOCUSTAG')
    fp = strip_by_class(d, 'rho', 'Spearman rho(codon frequency, tRNA gene count) per genome', 'V36', 'trna_codon_rho_by_class', unit='genome')
    gd = genus_medians(d, ['rho', 'n_trna', 'n_codons_with_trna']); cb, r = contrast_block(ctx, gd, 'rho', 'Spearman rho', covars=('log_n50',))
    tp = save_tab(gd, 'V36', 'genus_medians')
    anti = tr.groupby('codon')['n'].sum().sort_values(ascending=False).head(5)
    return (f"### V36 tRNA gene copy number versus codon usage\n\n![]({fp})\n\n- {counts_text(d)}. Filters: {ctx.filt_text(dup=False, model=False)}. Per genome: Spearman correlation across 61 sense codons between codon frequency (codon_freq) and the number of tRNA genes whose `codon` field matches.\n"
            f"- Median tRNA genes per representative genome: {d['n_trna'].median():.0f}; median codons with >= 1 tRNA gene: {d['n_codons_with_trna'].median():.0f} of 61. Most frequent `codon` values in gene_trnas: {', '.join(f'{c} ({n:,})' for c, n in anti.items())}.\n{cb}\n"
            f"- Caveat: the `codon` column of gene_trnas was taken as the decoded codon (the sampled rows are consistent with codons, not anticodons; not verified against tRNAscan output); wobble decoding is ignored, so rho is bounded below 1. Table: `{tp}`.")


def v37(ctx):
    cf = ctx.q(f"SELECT species_prefix AS genome, upper(codon) AS codon, frequency FROM '{T}/codon_freq.parquet' WHERE upper(codon) IN ('CTG','CTT','CTC','CTA','TTA','TTG','TCT','TCC','TCA','TCG','AGT','AGC')")
    piv = cf.pivot(index='genome', columns='codon', values='frequency')
    piv['ctg_frac_leu'] = piv['CTG'] / piv[['CTG', 'CTT', 'CTC', 'CTA', 'TTA', 'TTG']].sum(1); piv['ctg_freq'] = piv['CTG']
    d = ctx.reps(dup=False, n50=False).merge(piv[['ctg_freq', 'ctg_frac_leu']], left_on='LOCUSTAG', right_index=True)
    p = d[d['CLASS'] == 'Pichiomycetes'].copy(); p['code'] = p['transl_table'].astype(str)
    gd = genus_medians(p.assign(code=p['code']), ['ctg_freq', 'ctg_frac_leu']).merge(p.groupby('GENUS')['code'].agg(lambda s: s.mode().iloc[0]).rename('code'), left_on='GENUS', right_index=True)
    r = S.two_group(gd[gd['code'] == '12']['ctg_freq'], gd[gd['code'] == '1']['ctg_freq'], 'code12', 'code1')
    fig, ax = plt.subplots(figsize=(7, 4))
    orders = p.groupby('ORDER_').size().sort_values(ascending=False).index.tolist()
    for i, o in enumerate(orders):
        for code, colr, off in [('1', '#0072B2', -0.15), ('12', '#D55E00', 0.15)]:
            v = p[(p['ORDER_'] == o) & (p['code'] == code)]['ctg_freq'].values
            if len(v):
                ax.scatter(1000 * v, i + off + (np.random.default_rng(5).random(len(v)) - 0.5) * 0.2, s=6, alpha=0.5, color=colr, linewidths=0, rasterized=True)
    ax.set_yticks(range(len(orders))); ax.set_yticklabels(orders, fontsize=8); ax.set_xlabel('CTG per 1,000 codons'); ax.scatter([], [], color='#0072B2', label='code 1'); ax.scatter([], [], color='#D55E00', label='code 12 (CTG = Ser)'); ax.legend(fontsize=8)
    ax.set_title('V37: CTG frequency in Pichiomycetes by genetic code (representative genomes)', fontsize=9)
    fp = save_fig(fig, 'V37', 'ctg_by_code'); tp = save_tab(gd, 'V37', 'genus_medians')
    allg = genus_medians(d, ['ctg_freq'])
    return (f"### V37 CTG-clade codon usage\n\n![]({fp})\n\n- Pichiomycetes representative genomes: {len(p):,} ({int((p['code']=='12').sum()):,} with TRANSL_TABLE 12, {int((p['code']=='1').sum()):,} with code 1); unit = genus ({r['n_code12']} code-12 genera, {r['n_code1']} code-1 genera). Filters: {ctx.filt_text(dup=False, n50=False, model=False)}.\n"
            f"- CTG per 1,000 codons, genus medians: code 12 {1000*r['median_code12']:.2f}, code 1 {1000*r['median_code1']:.2f}; difference {1000*r.get('median_diff', np.nan):.2f} (95 % CI {1000*r.get('diff_ci_lo', np.nan):.2f} to {1000*r.get('diff_ci_hi', np.nan):.2f}), Cliff's delta {S.fmt(r.get('cliffs_delta'))}, p = {S.fmt(r.get('mw_p'))}. All-fungi genus median: {1000*allg['ctg_freq'].median():.2f} per 1,000.\n"
            f"- Caveat: the genetic code comes from the samples sheet, not from the sequence; the single code-26 genome (*Pachysolen*) is not shown. Table: `{tp}`.")


def v38(ctx):
    aa = ctx.q(f"SELECT species_prefix AS genome, amino_acid, frequency FROM '{T}/aa_freq.parquet' WHERE amino_acid NOT IN ('*','X')")
    P = aa.pivot(index='genome', columns='amino_acid', values='frequency').fillna(0); P = P.div(P.sum(1), axis=0)
    d = ctx.reps(dup=False, n50=False).merge(P, left_on='LOCUSTAG', right_index=True)
    cols = list(P.columns); gd = genus_medians(d, cols + ['gc_pct']).dropna(subset=['CLASS'])
    X = gd[cols].values; X = X - X.mean(0); U, s, Vt = np.linalg.svd(X, full_matrices=False); pcs = U * s; var = s ** 2 / (s ** 2).sum()
    gd['PC1'], gd['PC2'] = pcs[:, 0], pcs[:, 1]
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    sc = axs[0].scatter(gd['PC1'], gd['PC2'], c=gd['gc_pct'], s=7, cmap='viridis', linewidths=0, rasterized=True); plt.colorbar(sc, ax=axs[0], label='genomic GC (%)')
    axs[0].set_xlabel(f'PC1 ({100*var[0]:.1f} %)'); axs[0].set_ylabel(f'PC2 ({100*var[1]:.1f} %)'); axs[0].set_title('genus-median amino-acid composition, colored by GC', fontsize=8)
    for gf in GF_ORDER:
        s_ = gd[gd['growth_form'] == gf]; axs[1].scatter(s_['PC1'], s_['PC2'], s=7, alpha=0.5, color=GF_COLORS[gf], label=f'{gf} ({len(s_)})', linewidths=0, rasterized=True)
    axs[1].legend(fontsize=7); axs[1].set_xlabel(f'PC1 ({100*var[0]:.1f} %)'); axs[1].set_title('by growth form', fontsize=8)
    fig.suptitle('V38: PCA of amino-acid composition (20 aa, genus medians)', fontsize=9)
    fp = save_fig(fig, 'V38', 'aa_pca')
    load = pd.DataFrame({'amino_acid': cols, 'PC1_loading': Vt[0], 'PC2_loading': Vt[1]}).sort_values('PC1_loading')
    r1 = S.spearman(gd['PC1'], gd['gc_pct']); gd['pc1_score'] = gd['PC1']; gd['genomic_gc'] = gd['gc_pct']
    cb, r = contrast_block(ctx, gd.drop(columns=cols), 'pc1_score', 'PC1 score', covars=('genomic_gc',))
    tp = save_tab(gd[['GENUS', 'PHYLUM', 'CLASS', 'growth_form', 'gc_pct', 'PC1', 'PC2']], 'V38', 'genus_pcs'); tp2 = save_tab(load, 'V38', 'loadings')
    return (f"### V38 Amino-acid composition ordination\n\n![]({fp})\n\n- Unit = genus ({len(gd):,} genera). Filters: {ctx.filt_text(dup=False, n50=False, model=False)}. `*` and `X` symbols excluded; frequencies renormalized to the 20 amino acids.\n"
            f"- PC1 {100*var[0]:.1f} %, PC2 {100*var[1]:.1f} %. Spearman rho(PC1, genomic GC) = {r1['rho']:.3f}. Strongest PC1 loadings: " + ', '.join(f"{a} {v:+.2f}" for a, v in zip(load['amino_acid'].head(3), load['PC1_loading'].head(3))) + ' ... ' + ', '.join(f"{a} {v:+.2f}" for a, v in zip(load['amino_acid'].tail(3), load['PC1_loading'].tail(3))) +
            f".\n**PC1 score, yeast vs filamentous (GC as covariate in the mixed model):**\n{cb}\n- Tables: `{tp}`, `{tp2}`. Points are rasterized in the SVG.")


# ------------------------------------------------------------------ H7
def v39(ctx):
    w = ctx.tbl('A5/wgd_genome.parquet'); ws = ctx.q(f"SELECT species_prefix AS genome, hybrid, n_identical_groups, n_genes_in_identical_groups, largest_identical_group FROM '{T}/wgd.ks.summary.parquet'")
    g = ctx.tbl('A4/genome_protein_summary.parquet')[['genome', 'n_proteins']]
    L = ctx.labels[ctx.labels['is_representative'] & ctx.labels['q_ann'] & ctx.labels['q_busco_rel'].fillna(False) & (ctx.labels['busco_duplicated'] <= 5) & ctx.labels['q_n50'].fillna(False)]
    d = L.merge(w, left_on='LOCUSTAG', right_on='genome').merge(ws, on='genome').merge(g, on='genome'); d = d[(d['hybrid'] == 0) & d['growth_form'].isin(GF_ORDER)]
    d['young_per1000'] = d['n_young_nodes'] * 1000 / d['n_proteins']; d['nodes_per1000'] = d['n_tree_nodes'] * 1000 / d['n_proteins']
    fp = strip_by_class(d, 'young_per1000', 'paralog nodes with 0.05 <= dS <= 1.0 per 1,000 proteins', 'V39', 'young_paralogs_by_class', unit='genome')
    gd = genus_medians(d, ['young_per1000', 'nodes_per1000']); cb, r = contrast_block(ctx, gd, 'young_per1000', 'young paralog nodes per 1,000 proteins', covars=('log_genome', 'busco_complete'))
    tp = save_tab(gd, 'V39', 'genus_medians')
    return (f"### V39 Young paralog nodes per 1,000 genes by growth form\n\n![]({fp})\n\n- {counts_text(d)} of the 7,384 genomes in the interim wgd.ks release. Filters: Q-ann, Q-rep, Q-busco-rel, BUSCO duplicated <= 5 %, Q-n50, `hybrid = 0`. Young = tree node with mean yn00 dS in 0.05-1.0 (D9).\n{cb}\n"
            f"- Caveat: per fable_context.md, dS < 0.05 is heterozygosity/haplotigs and dS > 1.5 is saturated; the wgd.ks subset selection is not documented in the manifest; TE-derived models inflate node counts in repeat-rich genomes. Table: `{tp}`.")


def v40(ctx):
    ws = ctx.q(f"SELECT species_prefix AS genome, n_identical_groups, n_genes_in_identical_groups, largest_identical_group, identical_pair_share, hybrid FROM '{T}/wgd.ks.summary.parquet'")
    g = ctx.tbl('A4/genome_protein_summary.parquet')[['genome', 'n_proteins']]
    d = ctx.labels[ctx.labels['q_ann'].fillna(False)].merge(ws, left_on='LOCUSTAG', right_on='genome').merge(g, on='genome')
    d['ident_per1000'] = d['n_genes_in_identical_groups'] * 1000 / d['n_proteins']
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for ph, colr in [('Ascomycota', '#0072B2'), ('Basidiomycota', '#D55E00'), ('Mucoromycota', '#009E73')]:
        s = d[d['PHYLUM'] == ph]; ax.scatter(s['busco_duplicated'], s['ident_per1000'] + 0.1, s=5, alpha=0.4, color=colr, label=f'{ph} (n={len(s)})', linewidths=0, rasterized=True)
    s = d[~d['PHYLUM'].isin(['Ascomycota', 'Basidiomycota', 'Mucoromycota'])]; ax.scatter(s['busco_duplicated'], s['ident_per1000'] + 0.1, s=5, alpha=0.4, color='#999', label=f'other (n={len(s)})', linewidths=0, rasterized=True)
    ax.set_yscale('log'); ax.set_xlabel('BUSCO duplicated (%)'); ax.set_ylabel('genes in byte-identical CDS groups per 1,000 proteins (+0.1)'); ax.legend(fontsize=7)
    ax.set_title('V40: identical-copy genes vs BUSCO duplication (all annotated genomes in wgd.ks)', fontsize=9)
    fp = save_fig(fig, 'V40', 'identical_vs_dup')
    r = S.spearman(d['busco_duplicated'], d['ident_per1000'])
    out = d[(d['ident_per1000'] >= 50) | (d['largest_identical_group'] >= 20)].sort_values('ident_per1000', ascending=False)[['LOCUSTAG', 'SPECIES', 'PHYLUM', 'CLASS', 'busco_duplicated', 'ident_per1000', 'largest_identical_group', 'n_identical_groups', 'hybrid', 'N50_bp']]
    tp = save_tab(out, 'V40', 'outliers')
    return (f"### V40 Identical-copy groups versus BUSCO duplication\n\n![]({fp})\n\n- n = {len(d):,} annotated genomes with a wgd.ks summary (no representative filter: this is an artefact screen). Spearman rho(BUSCO duplicated, identical genes per 1,000) = {r['rho']:.3f}.\n"
            f"- {len(out):,} genomes have >= 50 identical-group genes per 1,000 proteins or a largest identical group >= 20; {int((d['busco_duplicated'] > 10).sum()):,} genomes have BUSCO duplication > 10 %. Outlier table: `{tp}`.\n\n" + md_table(out.head(15)) +
            "\n\n- Caveat: identical groups also arise from TE-derived gene models and annotation-transfer duplicates (fable_context.md Q1 ranks 2 and 8), not only unpurged haplotigs.")


VIEWS = {k: v for k, v in globals().items() if k.startswith('v') and k[1:].isdigit()}
PAGE_ADD = {'H1_quality': ['v06'], 'H2_protein_length': ['v10', 'v11', 'v13', 'v31'], 'H3_gene_structure': ['v16'], 'H4_secretome_lifestyle': ['v24', 'v25'],
            'H5_domains_tf': ['v29', 'v30'], 'H6_codon_usage': ['v36', 'v37', 'v38']}
PAGE_NEW = {'H7_paralogs': ('H7. Recent paralog retention is lower in yeasts; large identical-copy groups mark assembly artefacts',
                            'Expected: fewer young paralog nodes per 1,000 genes in yeasts. 7,384 genomes only (interim wgd.ks release).', ['v39', 'v40'])}
