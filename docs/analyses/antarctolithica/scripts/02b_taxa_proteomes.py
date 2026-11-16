"""Write results/taxa_proteomes_hpcc.tsv: LOCUSTAG, label and the absolute HPCC protein FASTA path of every phylogeny taxon
(the 159 taxa of the published FastTree design in the committed results/taxa.tsv, which include the 6 annotated MAGs),
for downstream IQ-TREE / PhyloBayes runs.

Usage: python 02b_taxa_proteomes.py <taxa.tsv> [out.tsv]   (run from AO_BASE; needs phylo/locustag_to_pep.tsv from 00_prep_data.py)
Columns: LOCUSTAG, label (tip name used in the trees), CLASS, SPECIES, STRAIN, is_AO, pep_path (path in Fungi_BFD_runs/input/pep,
usually a symlink), pep_realpath (resolved file), exists."""
import os, sys, pandas as pd


def proteome_table(taxa, l2p):
    T = taxa[['LOCUSTAG', 'label', 'CLASS', 'SPECIES', 'STRAIN']].merge(l2p[['LOCUSTAG', 'pep_path']], on='LOCUSTAG', how='left')
    T['is_AO'] = T.SPECIES.str.contains('Antarctolith', na=False)
    T['pep_realpath'] = [os.path.realpath(p) if isinstance(p, str) else None for p in T.pep_path]
    T['exists'] = [isinstance(p, str) and os.path.exists(p) for p in T.pep_realpath]
    return T[['LOCUSTAG', 'label', 'CLASS', 'SPECIES', 'STRAIN', 'is_AO', 'pep_path', 'pep_realpath', 'exists']]


if __name__ == '__main__':
    taxa = pd.read_csv(sys.argv[1], sep='\t')
    out = sys.argv[2] if len(sys.argv) > 2 else 'results/taxa_proteomes_hpcc.tsv'
    l2p = pd.read_csv('phylo/locustag_to_pep.tsv', sep='\t')
    T = proteome_table(taxa, l2p)
    T.to_csv(out, sep='\t', index=False)
    print(f'{len(T)} taxa ({int(T.is_AO.sum())} Antarctolithica MAGs); missing proteome: {int((~T.exists).sum())} -> {out}')
    if not T.exists.all():
        print(T[~T.exists].to_string()); sys.exit(1)
