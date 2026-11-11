"""SwissProt best-hit taxonomy as placement signal + per-contig contamination screen of the MAGs."""
import duckdb, gzip, glob, re, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from common import *
d = duckdb.connect()
om = pd.read_csv('data/orgcode_map0.csv')
euk = set('Arabidopsis Homo Mus Bos Dictyostelium Danio Rattus Gallus Drosophila Xenopus Pongo Oryza Caenorhabditis Oryctolagus Aedes Solanum Macaca Medicago Equus Sus Ovis Canis Felis Pan Zea Nicotiana Glycine Triticum Hordeum Chlamydomonas Plasmodium Toxoplasma Trypanosoma Leishmania Tetrahymena Paramecium Anopheles Bombyx Apis Mesocricetus Cricetulus Cavia Physcomitrium Pisum Brassica Spinacia Vitis Populus Schistosoma Strongylocentrotus Ciona Takifugu Oncorhynchus Salmo Oryzias Gorilla Callithrix Capra Entamoeba Giardia Trichomonas Phytophthora Acanthamoeba Naegleria Cyanidioschyzon Ostreococcus Volvox Selaginella Sorghum Manihot Ricinus Daucus Petunia Cucumis Gossypium Lotus Phaseolus Vigna Beta Helianthus Catharanthus Papaver Coffea Camellia Citrus Malus Prunus Fragaria Rosa Mentha Ocimum Artemisia Taxus Pinus Picea Ginkgo Marchantia'.split())
om['kingdom'] = np.where(om.fungal, 'Fungi', np.where(om.genus.isin(euk), 'other Eukaryota', 'prokaryote/virus/unclassified'))
om.to_csv('results/orgcode_kingdom_map.csv', index=False)
L = pd.read_parquet('data/feat/sprot_org_long.parquet').merge(om[['org_code', 'kingdom']], on='org_code', how='left')
tot = L.groupby('g').n.sum()
def frac(mask): return (L[mask].groupby('g').n.sum() / tot).reindex(U.index).fillna(0)
F = pd.DataFrame(index=U.index)
F['frac_SCHPO'] = frac(L.org_code == 'SCHPO'); F['frac_YEAST'] = frac(L.org_code == 'YEAST')
F['frac_NEUCR'] = frac(L.org_code == 'NEUCR'); F['frac_EMENI_ASP'] = frac(L.org_code.str.startswith(('EMEN', 'ASP')))
F['frac_Fungi'] = frac(L.kingdom == 'Fungi'); F['frac_otherEuk'] = frac(L.kingdom == 'other Eukaryota'); F['frac_prok_like'] = frac(L.kingdom == 'prokaryote/virus/unclassified')
for oc in ['SCHPO', 'YEAST', 'NEUCR', 'EMENI', 'YARLI', 'CANAL']:
    x = L[L.org_code == oc].set_index('g'); F['pid_' + oc] = x.pid.reindex(U.index); F['pidcov_' + oc] = x.pid_cov.reindex(U.index)
F['pid_NEUCR_minus_SCHPO'] = F.pid_NEUCR - F.pid_SCHPO; F['pid_NEUCR_minus_YEAST'] = F.pid_NEUCR - F.pid_YEAST; F['pid_SCHPO_minus_YEAST'] = F.pid_SCHPO - F.pid_YEAST
F.to_parquet('data/feat/sprot_scalars.parquet')
# composition block: top organisms, CLR
W = L.pivot_table(index='g', columns='org_code', values='n', aggfunc='sum', fill_value=0).reindex(U.index).fillna(0)
W = W.loc[:, W.sum().sort_values(ascending=False).index[:80]]; W.to_parquet('data/feat/sprot_org.parquet')
tg = U.index[U.is_genus_unit & (U.role == 'train')]
g = U.loc[tg].copy(); g['grp'] = np.where(g.subphylum == 'Pezizomycotina', g.CLASS, g.CLASS)
tab = F.loc[tg].join(g[['subphylum', 'CLASS']]).groupby(['subphylum', 'CLASS']).median(numeric_only=True)
tab.insert(0, 'n_genera', g.groupby(['subphylum', 'CLASS']).size())
ao = F.loc[HQ].median().to_frame().T; ao.index = pd.MultiIndex.from_tuples([('AO', 'Antarctolithica (4 HQ MAGs)')]); ao.insert(0, 'n_genera', 4)
tab = pd.concat([ao, tab]); tab.to_csv('results/sprot_class_table.csv')
pd.set_option('display.width', 300); pd.set_option('display.max_columns', 40)
print(tab.round(3).to_string())
print('AO percentile of frac_prok_like among Ascomycota genus units:', (F.loc[tg].frac_prok_like < F.loc[HQ].frac_prok_like.median()).mean(), F.loc[AO].frac_prok_like.round(4).to_dict())
# ---------- per-contig screen
rows = []
sb = pd.read_parquet('data/ao_swissprot_best.parquet').merge(om[['org_code', 'kingdom']], on='org_code', how='left')
for dpath in sorted(glob.glob('mag/*/predict_results')):
    fa = glob.glob(dpath + '/*.scaffolds.fa.gz')[0]; gff = glob.glob(dpath + '/*.gff3.gz')[0]
    seqs = {}; name = None
    for line in gzip.open(fa, 'rt'):
        if line[0] == '>': name = line[1:].split()[0]; seqs[name] = []
        else: seqs[name].append(line.strip().upper())
    seqs = {k: ''.join(v) for k, v in seqs.items()}
    g2c = {}
    for line in gzip.open(gff, 'rt'):
        if line[0] == '#': continue
        p = line.split('\t')
        if len(p) > 8 and p[2] == 'gene':
            m = re.search(r'ID=([^;\n]+)', p[8]); g2c[m.group(1)] = p[0]
    tag = next(iter(g2c)).split('_')[0]
    s = sb[sb.species_prefix == tag].copy(); s['contig'] = s.protein_id.str.replace(r'-T\d+$', '', regex=True).map(g2c)
    genes = pd.Series(g2c).value_counts()
    for c, sq in seqs.items():
        acgt = sum(sq.count(x) for x in 'ACGT'); x = s[s.contig == c]
        rows.append(dict(mag=tag, strain=MAGNAME[tag], contig=c, length=len(sq), gc=100 * (sq.count('G') + sq.count('C')) / max(acgt, 1), n_genes=int(genes.get(c, 0)), n_hits=len(x),
                         n_fungi=int((x.kingdom == 'Fungi').sum()), n_prok_like=int((x.kingdom == 'prokaryote/virus/unclassified').sum()), n_otherEuk=int((x.kingdom == 'other Eukaryota').sum()), mean_pid=x.pident_x10.mean() / 10 if len(x) else np.nan))
C = pd.DataFrame(rows)
wm = C.groupby('mag').apply(lambda x: np.average(x.gc, weights=x.length)); C['gc_dev'] = C.gc - C.mag.map(wm)
C['genes_per_kb'] = C.n_genes / C.length * 1000
C['flag'] = ((C.gc_dev.abs() > 6) & (C.length >= 3000)) | ((C.n_hits >= 4) & (C.n_prok_like / C.n_hits.clip(lower=1) >= 0.5)) | ((C.length >= 10000) & (C.n_genes == 0))
C.to_csv('results/mag_contigs.csv', index=False)
summ = C.groupby('strain').apply(lambda x: pd.Series(dict(contigs=len(x), total_mb=x.length.sum() / 1e6, gc_wmean=np.average(x.gc, weights=x.length), gc_sd_contigs=np.sqrt(np.average((x.gc - np.average(x.gc, weights=x.length)) ** 2, weights=x.length)),
            flagged_contigs=int(x.flag.sum()), flagged_kb=x[x.flag].length.sum() / 1e3, flagged_pct_bp=100 * x[x.flag].length.sum() / x.length.sum(), genes_on_flagged=int(x[x.flag].n_genes.sum()),
            hits_fungi_pct=100 * x.n_fungi.sum() / x.n_hits.sum(), hits_prok_like_pct=100 * x.n_prok_like.sum() / x.n_hits.sum())))
summ.to_csv('results/mag_contam_summary.csv'); print(summ.round(2).to_string())
print(C[C.flag].sort_values('length', ascending=False).head(25).round(2).to_string())
