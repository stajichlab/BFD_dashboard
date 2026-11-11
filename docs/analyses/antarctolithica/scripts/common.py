"""Shared helpers: comparison set, block transforms, colours."""
import os, sys, numpy as np, pandas as pd
# Run from AO_BASE (the working dir holding data/); repo root found relative to this file.
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[4]))
FIG = os.environ.get('AO_BASE', '.')  # where figures are written
U = pd.read_parquet('data/ao_set.parquet').set_index('LOCUSTAG')
HQ = ['FF1FC3B3', 'F46C92EF', 'F7B4B958', 'F7E47E53']
AO = list(U.index[U.role == 'AO'])
MAGNAME = {i: U.STRAIN[i] for i in AO}
BIG = ['Sordariomycetes', 'Dothideomycetes', 'Eurotiomycetes', 'Lecanoromycetes', 'Leotiomycetes']
def plot_group(r):
    if r.role == 'AO': return 'Antarctolithica'
    if r.role == 'outgroup': return 'Basidiomycota (outgroup)'
    if r.role == 'unplaced': return 'unplaced Ascomycota'
    if r.subphylum in ('Saccharomycotina', 'Taphrinomycotina'): return r.subphylum
    return r.CLASS if r.CLASS in BIG else 'other Pezizomycotina'
U['pg'] = U.apply(plot_group, axis=1)
COL = {'Sordariomycetes': '#2a78d6', 'Dothideomycetes': '#eb6834', 'Eurotiomycetes': '#1baf7a', 'Lecanoromycetes': '#eda100',
       'Leotiomycetes': '#e87ba4', 'Saccharomycotina': '#008300', 'Taphrinomycotina': '#6250d6',
       'other Pezizomycotina': '#8a8a84', 'unplaced Ascomycota': '#55554f', 'Basidiomycota (outgroup)': '#c9c8c0', 'Antarctolithica': '#111111'}
def feat(n): return pd.read_parquet(f'data/feat/{n}.parquet').reindex(U.index).astype(float)
def clr(X, pseudo):
    X = X + pseudo; X = X.div(X.sum(1), axis=0); L = np.log(X); return L.sub(L.mean(1), axis=0)
def zscore(X, ref):
    return (X - X.loc[ref].mean()) / X.loc[ref].std(ddof=0).replace(0, 1)
LOGS = ['median_gene_span','median_cds_len','median_protein_len','mean_protein_len','q10_protein_len','q90_protein_len','median_intron_len','q90_intron_len',
        'median_intergenic','q10_intergenic','q90_intergenic','mean_intergenic','genes_per_mb','trna_per_mb']
EXTENSIVE = ['log10_genome_mb','log10_n_genes','n_trna','masked_pct','partial_fraction']
def blocks():
    """Return dict name -> (matrix, metric, description). Matrices are already transformed."""
    nprot = 10 ** feat('struct').log10_n_genes
    s = feat('struct').copy()
    for c in LOGS: s[c] = np.log10(s[c].clip(lower=1))
    cell = feat('cell').copy()
    pf = feat('pfam'); pf = pf.loc[:, (pf > 0).mean() >= 0.01]
    B = {}
    B['struct'] = (s.drop(columns=EXTENSIVE), 'z-euclidean', 'gene structure/density (intensive scalars)')
    B['struct_all'] = (s.drop(columns=['partial_fraction']), 'z-euclidean', 'gene structure incl. genome size and gene count')
    B['codon'] = (clr(feat('codon'), 1e-5), 'euclidean', 'codon usage, CLR (Aitchison)')
    B['aa'] = (clr(feat('aa'), 1e-6), 'euclidean', 'amino-acid composition, CLR (Aitchison)')
    B['trna'] = (clr(feat('trna'), 0.5), 'euclidean', 'tRNA gene pool by anticodon class, CLR (Aitchison)')
    B['pfam_pa'] = ((pf > 0).astype(float), 'jaccard', 'Pfam presence/absence (Jaccard)')
    B['pfam_ab'] = (np.sqrt(pf.div(nprot, axis=0)), 'braycurtis', 'Pfam proteins per proteome, sqrt (Bray-Curtis)')
    for n in ['cazy', 'merops']:
        x = feat(n); x = x.loc[:, (x > 0).mean() >= 0.01]
        B[n] = (np.sqrt(x.div(nprot, axis=0) * 1000), 'braycurtis', f'{n} families per 1k proteins, sqrt (Bray-Curtis)')
    B['tf'] = (clr(feat('tf'), 0.5), 'euclidean', 'TF family composition, CLR (Aitchison)')
    B['cell'] = (cell, 'z-euclidean', 'secretome/TM/GPI/disorder/localization/CAZy-MEROPS-TF fractions/signaling')
    return B
