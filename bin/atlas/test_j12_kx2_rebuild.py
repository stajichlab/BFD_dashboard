"""Fixture test for the J12 driver (small synthetic data, no HPCC data, runs in a second).

Builds a tiny stand-in for BFD.duckdb (species, asm_stats, busco_genome with the LOCAL column names), runs the real driver on two toy SQL
stages and checks: the shim column names work, statements are split correctly (semicolons inside comments and strings), result queries
are written as TSV, DDL echoes are not, scratch tables are fingerprinted, and `bfd` was not written.
Also checks statically that every bfd.<table> read by the saved H1-H4 SQL exists in the compatibility view list.
"""
import csv
import re
import sys
from pathlib import Path

import duckdb
import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import j12_kx2_rebuild as J  # noqa: E402


@pytest.fixture()
def toy(tmp_path):
    bfd = tmp_path / 'BFD.duckdb'
    c = duckdb.connect(str(bfd))
    c.execute('''CREATE TABLE species(ASMID VARCHAR, LOCUSTAG VARCHAR, PHYLUM VARCHAR, SUBPHYLUM VARCHAR, CLASS VARCHAR, "ORDER" VARCHAR,
                 FAMILY VARCHAR, GENUS VARCHAR, SPECIES VARCHAR, BUSCO_LINEAGE VARCHAR, SUBCLASS VARCHAR, NCBI_TAXONID BIGINT, STRAIN VARCHAR)''')
    c.execute('''CREATE TABLE asm_stats(LOCUSTAG VARCHAR, ASMID VARCHAR, contig_count BIGINT, TOTAL_LENGTH BIGINT, N50_bp BIGINT,
                 GC_PERCENT DOUBLE, annotation_status VARCHAR)''')
    c.execute('CREATE TABLE busco_genome(LOCUSTAG VARCHAR, ASMID VARCHAR, complete_pct DOUBLE, duplicated_pct DOUBLE)')
    c.execute("INSERT INTO species VALUES ('A1','AAA','P','', 'C','Kickxellales','F','Coemansia','Coemansia sp','',NULL,1,'s1'),"
              "('A2','AAB','P','', 'C','Kickxellales','F','Coemansia','Coemansia sp','',NULL,2,'s2'),"
              "('A3','CCC','P','', 'C','Mucorales','F','Mucor','Mucor sp','',NULL,3,'s3')")
    c.execute("INSERT INTO asm_stats VALUES ('AAA','A1',10,1000000,50000,40.0,'annotated'),('AAB','A2',5,1100000,90000,41.0,'annotated'),"
              "('CCC','A3',7,2000000,80000,42.0,'annotated')")
    c.execute("INSERT INTO busco_genome VALUES ('AAA','A1',80,1),('AAB','A2',80,1),('CCC','A3',90,2)")
    c.close()
    root = tmp_path / 'root'
    for s, files in {
        'sql_upstream': {'00_x.sql': "-- a comment; with a semicolon\nCREATE OR REPLACE TABLE bfd_scratch.kw_toy AS SELECT s.LOCUSTAG, s.ORDER_, a.total_length_bp FROM bfd.species s JOIN bfd.asm_stats a USING (ASMID);\nSELECT count(*) AS n FROM bfd_scratch.kw_toy;"},
        'sql_h1': {'01_y.sql': "CREATE OR REPLACE TABLE kx2_h1_t AS SELECT ORDER_, count(*) n, 'a;b' AS s FROM kw_toy GROUP BY 1;\nSELECT * FROM kx2_h1_t ORDER BY 1;"},
    }.items():
        (root / s).mkdir(parents=True)
        for fn, body in files.items():
            (root / s / fn).write_text(body)
    return bfd, root, tmp_path


def test_split_and_rewrite():
    con = duckdb.connect()
    parts = J.split_statements(con, "-- c; d\nSELECT 1; SELECT ';' AS x;\n-- tail;")
    # DuckDB's splitter keeps comments attached to the neighbouring statement; semicolons in comments and strings do not split
    assert len(parts) == 2 and parts[0].endswith('SELECT 1') and "SELECT ';' AS x" in parts[1]
    assert [con.execute(p).fetchall() for p in parts] == [[(1,)], [(';',)]]
    assert J.rewrite('FROM bfd.species s JOIN bfd.asm_stats a') == 'FROM bfd_scratch.j12_species s JOIN bfd_scratch.j12_asm_stats a'
    assert J.rewrite('bfd.species_foo') == 'bfd.species_foo'


def test_driver_end_to_end(toy):
    bfd, root, tmp = toy
    out = tmp / 'out'
    J.main(['--sql-root', str(root), '--out', str(out), '--scratch-db', str(tmp / 's.duckdb'), '--bfd-db', str(bfd),
            '--lib', str(ROOT / 'lib'), '--stages', 'upstream,h1', '--threads', '1', '--memory', '1GB'])
    res = sorted(p.name for p in (out / 'results').glob('*.tsv'))
    assert res == ['h1__01_y__q2.tsv', 'upstream__00_x__q2.tsv']       # only SELECTs, not the DDL echo
    rows = list(csv.reader(open(out / 'results' / 'h1__01_y__q2.tsv'), delimiter='\t'))
    assert rows == [['ORDER_', 'n', 's'], ['Kickxellales', '2', 'a;b'], ['Mucorales', '1', 'a;b']]
    counts = dict(csv.reader(open(out / 'table_counts.tsv'), delimiter='\t'))
    assert counts['kw_toy'] == '3' and counts['kx2_h1_t'] == '2'
    fp = {(r[0], r[1]): r for r in csv.reader(open(out / 'table_fingerprint.tsv'), delimiter='\t')}
    assert fp[('kw_toy', 'total_length_bp')][4] == '4.1e+06'
    # bfd must be untouched (read-only attach)
    c = duckdb.connect(str(bfd), read_only=True)
    assert c.execute('SELECT count(*) FROM duckdb_tables() WHERE table_name LIKE \'kw_%\' OR table_name LIKE \'j12%\'').fetchone()[0] == 0


def test_unknown_stage(toy):
    bfd, root, tmp = toy
    with pytest.raises(SystemExit):
        J.main(['--sql-root', str(root), '--out', str(tmp / 'o'), '--scratch-db', str(tmp / 's.duckdb'), '--bfd-db', str(bfd),
                '--lib', str(ROOT / 'lib'), '--stages', 'h9'])


def test_saved_sql_reads_only_known_bfd_tables():
    """Static check: every bfd.<name> in the saved SQL is a table or view that the local BFD.duckdb builder defines."""
    defined = set()
    for f in [ROOT / 'sql' / 'bfd_compat_views.sql', ROOT / 'sql' / 'bfd_extra_views.sql']:
        defined |= set(re.findall(r'CREATE OR REPLACE (?:MACRO|VIEW|TABLE)\s+([a-zA-Z_0-9]+)', f.read_text()))
    # base views come from workflows/build_bfd_db.sh (one view per parquet file in tables/)
    base = {'species', 'asm_stats', 'busco_genome', 'pfam', 'gene_proteins', 'gene_transcripts', 'gene_exons', 'gene_info', 'gene_introns',
            'gene_intergenic_distances', 'idp_summary', 'swissprot_best', 'swissprot_subject', 'tf_inventory', 'merops', 'cazy_overview',
            'codon_freq', 'aa_freq', 'signalp_signal_peptide', 'gene_trnas'}
    used = set()
    for d in ['sql_upstream', 'sql_h1', 'sql_h2', 'sql_h3', 'sql_h4']:
        for p in (ROOT / 'docs/analyses/kx2_kickx_reduction' / d).glob('*.sql'):
            used |= set(re.findall(r'(?<![A-Za-z_])bfd\.([a-z_0-9]+)', p.read_text()))
    assert used, 'no bfd.* reads found'
    assert used <= (base | defined), sorted(used - base - defined)
