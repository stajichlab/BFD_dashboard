"""Fixture tests for the analysis-C data-prep helpers (00_prep_data.py, 02b_taxa_proteomes.py, 08b_mag_aai.py).
Run: pixi run --manifest-path docs/analyses/antarctolithica/pixi.toml pytest docs/analyses/antarctolithica/tests -q"""
import importlib.util, os
from pathlib import Path
import pandas as pd, pytest

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace('.py', '').lstrip('0123456789_b'), SCRIPTS / name)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


prep = load('00_prep_data.py')


def fake_bfd(root):
    """Minimal BFD tree with every file 00_prep_data.py links, two annotated MAGs and one skipped MAG."""
    atlas = root / 'BFD_dashboard' / 'atlas_data'
    for sub, files in prep.ATLAS_FILES.items():
        for f in files:
            (atlas / sub).mkdir(parents=True, exist_ok=True); (atlas / sub / f).write_text(f'{sub}/{f}')
    (atlas / 'phylo').mkdir(parents=True, exist_ok=True)
    (atlas / 'phylo' / 'locustag_to_pep.tsv').write_text('LOCUSTAG\tout\tpep_path\nF1\tA_b\t/x/A_b.proteins.fa.gz\n')
    odb = atlas / 'phylo' / 'phyling_db' / 'fungi_odb10'; (odb / 'hmms').mkdir(parents=True)
    for m in ['m2', 'm1']:
        (odb / 'hmms' / f'{m}.hmm').write_text(f'HMMER3/f\nNAME  {m}\n//\n')
    (odb / 'scores_cutoff').write_text('m1\t10\nm2\t20\n')
    tables = root / 'Fungi_BFD_runs' / 'tables'; tables.mkdir(parents=True)
    for f in prep.RAW_TABLES:
        (tables / f).write_text(f)
    ann = root / 'Fungi_BFD_runs' / 'genome_annotation'
    for m in ['Mars-82.9', 'Mars-17.30']:
        (ann / f'Ascomycota_sp._Antarctolithica_obscura_{m}' / 'predict_results').mkdir(parents=True)
    (ann / 'Ascomycota_sp._Antarctolithica_obscura_Mars-9.1').mkdir()          # skipped_too_small: no predict_results
    (ann / 'Other_fungus' / 'predict_results').mkdir(parents=True)
    return root


def test_sources_env_override(tmp_path):
    s = prep.sources({'BFD_DATA': str(tmp_path)})
    assert s['atlas'] == tmp_path / 'BFD_dashboard' / 'atlas_data'
    assert s['tables'] == tmp_path / 'Fungi_BFD_runs' / 'tables'
    assert s['odb10'] == s['atlas'] / 'phylo' / 'phyling_db' / 'fungi_odb10'
    assert prep.sources({'BFD_DATA': str(tmp_path), 'AO_TABLES': '/t'})['tables'] == Path('/t')


def test_stage_builds_layout(tmp_path):
    bfd = fake_bfd(tmp_path / 'bfd'); base = tmp_path / 'work'
    info = prep.stage(base, prep.sources({'BFD_DATA': str(bfd)}), SCRIPTS)
    assert (base / 'data/A7/pfam_genome_matrix.parquet').read_text() == 'A7/pfam_genome_matrix.parquet'
    assert (base / 'data/labels/pfam_trait_sets.tsv').is_symlink()
    assert (base / 'data/raw/wgd.ks.summary.parquet').read_text() == 'wgd.ks.summary.parquet'
    assert info['n_hmm'] == 2 and (base / 'phylo/db/all.hmm').read_text().count('//') == 2
    assert (base / 'phylo/db/scores_cutoff').read_text().startswith('m1')
    assert info['mags'] == ['Ascomycota_sp._Antarctolithica_obscura_Mars-17.30', 'Ascomycota_sp._Antarctolithica_obscura_Mars-82.9']
    assert (base / 'mag/Ascomycota_sp._Antarctolithica_obscura_Mars-82.9/predict_results').is_dir()
    assert (base / 'report/01_main.txt').exists() and (base / 'results').is_dir()
    # idempotent: a second run replaces links instead of failing
    prep.stage(base, prep.sources({'BFD_DATA': str(bfd)}), SCRIPTS)
    # and nothing is written into the source tree
    assert not any(p.name in ('all.hmm', 'data') for p in bfd.rglob('*'))


def test_stage_missing_input_fails(tmp_path):
    bfd = fake_bfd(tmp_path / 'bfd'); os.remove(bfd / 'Fungi_BFD_runs/tables/codon_freq.parquet')
    with pytest.raises(FileNotFoundError, match='codon_freq'):
        prep.stage(tmp_path / 'work', prep.sources({'BFD_DATA': str(bfd)}), SCRIPTS)


def test_taxa_proteomes(tmp_path):
    tp = load('02b_taxa_proteomes.py')
    f = tmp_path / 'A.proteins.fa.gz'; f.write_text('>p\nM\n'); ln = tmp_path / 'link.fa.gz'; ln.symlink_to(f)
    taxa = pd.DataFrame(dict(LOCUSTAG=['F1', 'F2'], label=['a', 'b'], CLASS=['X', None], SPECIES=['Ascomycota sp. Antarctolithica obscura', 'Y y'], STRAIN=['Mars-82.9', 's']))
    l2p = pd.DataFrame(dict(LOCUSTAG=['F1', 'F2'], out=['a', 'b'], pep_path=[str(ln), str(tmp_path / 'missing.fa.gz')]))
    T = tp.proteome_table(taxa, l2p)
    assert T.is_AO.tolist() == [True, False]
    assert T.pep_realpath[0] == str(f.resolve()) and T.exists.tolist() == [True, False]


def test_aai_summary():
    aai = load('08b_mag_aai.py')
    h = pd.DataFrame([['q1', 's1', 100.0, 100, 100, 100, 0, 200], ['q1', 's2', 50.0, 100, 100, 100, 0, 50],    # best hit kept
                      ['q2', 's2', 98.0, 95, 100, 100, 0, 150],                                               # passes 90 % coverage
                      ['q3', 's3', 99.0, 50, 100, 100, 0, 90]], columns=aai.COLS)                             # fails coverage
    r = aai.summarise(h, 'Mars-17.30')
    assert r['n_pairs'] == 2 and r['median_AAI'] == 99.0 and r['frac_identical'] == 0.5 and r['frac_lt90'] == 0
    assert aai.strain('/x/Ascomycota_sp._Antarctolithica_obscura_Mars-17.30.proteins.fa.gz') == 'Mars-17.30'


def test_swissprot_aggregation():
    import duckdb
    sx = load('01b_swissprot_extract.py')
    df = pd.DataFrame(dict(species_prefix=['G1'] * 3 + ['G2'], org_code=['SCHPO'] * 3 + ['YEAST'],
                           pident_x10=pd.array([400, 600, 500, 700], dtype='uint16'), qcov_x10=pd.array([900, 900, 500, 1000], dtype='uint16'),
                           scov_x10=pd.array([850, 700, 900, 1000], dtype='uint16')))
    con = duckdb.connect(); con.register('df', df)
    r = con.sql(sx.org_long_sql('df')).df().set_index('g')
    assert r.at['G1', 'n'] == 3 and r.at['G1', 'pid'] == 50.0 and r.at['G1', 'pid_cov'] == 40.0   # only the first hit passes 80 %/80 %
    assert r.at['G2', 'pid_cov'] == 70.0
    m = sx.orgcode_map(pd.DataFrame(dict(org_code=['SCHPO', 'HUMAN', 'EMENI'], organism=['Schizosaccharomyces pombe (strain 972)', 'Homo sapiens', 'Emericella nidulans'])),
                       ['Schizosaccharomyces'], ['Emericella'])
    assert m.set_index('org_code').fungal.to_dict() == {'EMENI': True, 'HUMAN': False, 'SCHPO': True}
