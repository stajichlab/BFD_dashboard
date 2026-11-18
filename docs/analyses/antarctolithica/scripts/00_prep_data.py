"""Stage the inputs of analysis C into the working directory AO_BASE (step 0, run before 01_build_set.py).

Everything here is a symlink or a small copy; nothing under BFD_DATA is written.

Environment variables
  AO_BASE   working directory that will hold data/, phylo/, mag/, report/, results/ (default: current directory)
  BFD_DATA  root of the shared BFD tree (default /bigdata/stajichlab/shared/projects/BFD). Derived locations:
              atlas   $BFD_DATA/BFD_dashboard/atlas_data          (A2/A4/A5/A7 intermediates, labels/, phylo/)
              tables  $BFD_DATA/Fungi_BFD_runs/tables              (raw BFD parquet snapshot)
              annot   $BFD_DATA/Fungi_BFD_runs/genome_annotation   (funannotate predict_results of the MAGs)
              odb10   $atlas/phylo/phyling_db/fungi_odb10          (BUSCO HMMs + scores_cutoff)
            Each can be overridden with AO_ATLAS, AO_TABLES, AO_ANNOT, AO_ODB10.

Layout produced (the paths the numbered scripts read):
  data/labels/{genome_labels.parquet,pfam_trait_sets.tsv}
  data/A2/genome_gene_stats.parquet, data/A4/genome_protein_summary.parquet, data/A5/*.parquet, data/A7/pfam_genome_matrix.parquet
  data/raw/{codon_freq,aa_freq,telomere_summary,wgd.ks.summary}.parquet
  phylo/locustag_to_pep.tsv, phylo/db/{all.hmm,scores_cutoff}
  mag/<MAG name>/predict_results -> annotation dir of each annotated Antarctolithica MAG
  report/{01_main.txt,02_methods.txt} (prose fragments for 17_assemble_report.py), results/
"""
import glob, os, shutil, sys
from pathlib import Path

ATLAS_FILES = {
    'A2': ['genome_gene_stats.parquet'],
    'A4': ['genome_protein_summary.parquet'],
    'A5': ['cazy_class_counts.parquet', 'cazy_substrate_counts.parquet', 'merops_family_counts.parquet', 'tf_counts.parquet',
           'secretome.parquet', 'wolfpsort_counts.parquet', 'idp_genome.parquet', 'trna_codon_counts.parquet', 'wgd_genome.parquet'],
    'A7': ['pfam_genome_matrix.parquet'],
    'labels': ['genome_labels.parquet', 'pfam_trait_sets.tsv'],
}
RAW_TABLES = ['codon_freq.parquet', 'aa_freq.parquet', 'telomere_summary.parquet', 'wgd.ks.summary.parquet']
MAG_GLOB = '*Antarctolithica*'


def sources(env=os.environ):
    bfd = Path(env.get('BFD_DATA', '/bigdata/stajichlab/shared/projects/BFD'))
    atlas = Path(env.get('AO_ATLAS', bfd / 'BFD_dashboard' / 'atlas_data'))
    return dict(bfd=bfd, atlas=atlas,
                tables=Path(env.get('AO_TABLES', bfd / 'Fungi_BFD_runs' / 'tables')),
                annot=Path(env.get('AO_ANNOT', bfd / 'Fungi_BFD_runs' / 'genome_annotation')),
                odb10=Path(env.get('AO_ODB10', atlas / 'phylo' / 'phyling_db' / 'fungi_odb10')))


def link(src, dst):
    """Symlink dst -> src (absolute), replacing an existing link; fail loudly if src is missing."""
    src, dst = Path(src).resolve(), Path(dst)
    if not src.exists():
        raise FileNotFoundError(f'missing input: {src}')
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    dst.symlink_to(src)


def build_hmm_db(odb10, out):
    """Concatenate the BUSCO marker HMMs into phylo/db/all.hmm and copy scores_cutoff (as 03_phylo_search.sh expects)."""
    out.mkdir(parents=True, exist_ok=True)
    hmms = sorted(glob.glob(str(odb10 / 'hmms' / '*.hmm')))
    if not hmms:
        raise FileNotFoundError(f'no HMMs in {odb10}/hmms')
    with open(out / 'all.hmm', 'w') as o:
        for h in hmms:
            o.write(open(h).read())
    shutil.copy(odb10 / 'scores_cutoff', out / 'scores_cutoff')
    return len(hmms)


def stage(base, src, script_dir):
    base = Path(base)
    for sub, files in ATLAS_FILES.items():
        for f in files:
            link(src['atlas'] / sub / f, base / 'data' / sub / f)
    for f in RAW_TABLES:
        link(src['tables'] / f, base / 'data' / 'raw' / f)
    link(src['atlas'] / 'phylo' / 'locustag_to_pep.tsv', base / 'phylo' / 'locustag_to_pep.tsv')
    n_hmm = build_hmm_db(src['odb10'], base / 'phylo' / 'db')
    mags = []
    for d in sorted(glob.glob(str(src['annot'] / MAG_GLOB))):
        pr = Path(d) / 'predict_results'
        if pr.is_dir():                     # MAGs skipped by the annotation pipeline (too small) have no predict_results
            link(pr, base / 'mag' / Path(d).name / 'predict_results'); mags.append(Path(d).name)
    (base / 'report').mkdir(parents=True, exist_ok=True); (base / 'results').mkdir(exist_ok=True); (base / 'data' / 'feat').mkdir(exist_ok=True)
    for f in ['01_main.txt', '02_methods.txt']:
        shutil.copy(Path(script_dir) / f, base / 'report' / f)
    return dict(n_hmm=n_hmm, mags=mags)


if __name__ == '__main__':
    base = Path(os.environ.get('AO_BASE', '.')).resolve()
    src = sources()
    info = stage(base, src, Path(__file__).resolve().parent)
    for k, v in src.items():
        print(f'{k:7s} {v}')
    print('AO_BASE', base, '| HMMs:', info['n_hmm'], '| annotated MAGs:', len(info['mags']), info['mags'])
    sys.exit(0)
