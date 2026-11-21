"""sql/dashboard_views.sql writes the labels into its own file; the result must equal the old views."""
import sys
from pathlib import Path

import duckdb
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_traits_sql as T  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DASH_SQL = ROOT / "sql" / "dashboard_views.sql"

SPECIES = [  # LOCUSTAG, ASMID, PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES, STRAIN, NCBI_TAXONID, SUBPHYLUM
    ("A1", "G1", "Mucoromycota", "Mucoromycetes", "Mucorales", "Mucoraceae", "Mucor", "Mucor a", "s1", "1", None),
    ("A2", "G2", "Mucoromycota", "Mucoromycetes", "Mucorales", "Mucoraceae", "Mucor", "Mucor a", "s2", "1", None),  # same species
    ("A3", "G3", "Mucoromycota", "Mucoromycetes", "Mucorales", "Mucoraceae", "Mucor", "Mucor b", "s3", "2", None),
    ("A4", "G4", "Ascomycota", None, None, None, None, "Ascomycota sp. Placeholder x", "s4", "3", None),             # species table
    ("A5", "G5", "Ascomycota", "Eurotiomycetes", "Eurotiales", "Aspergillaceae", "Aspergillus", "Aspergillus c", "s5", "4", "Given"),
    ("A6", "G6", None, None, None, None, None, None, None, None, None),                                           # unannotated, empty
]
FT = [T.ft("j1", "Mucor", "saprotroph", "filamentous_mycelium"), T.ft("j2", "Aspergillus", "saprotroph", "filamentous_mycelium")]
CLS = [("Mucoromycetes", "Mucoromycotina", "proposed", ""), ("Eurotiomycetes", "Pezizomycotina", "ncbi", "")]
SP = [("Ascomycota sp. Placeholder x", "PlaceholderSub", "proposed", "")]


def make_source(con):
    con.execute("""CREATE TABLE species (LOCUSTAG VARCHAR, ASMID VARCHAR, PHYLUM VARCHAR, CLASS VARCHAR, "ORDER" VARCHAR,
                   FAMILY VARCHAR, GENUS VARCHAR, SPECIES VARCHAR, STRAIN VARCHAR, NCBI_TAXONID VARCHAR, SUBPHYLUM VARCHAR)""")
    con.executemany("INSERT INTO species VALUES (?,?,?,?,?,?,?,?,?,?,?)", SPECIES)
    con.execute("CREATE TABLE samples (LOCUSTAG VARCHAR, TRANSL_TABLE INTEGER)")
    con.executemany("INSERT INTO samples VALUES (?,1)", [(r[0],) for r in SPECIES])
    con.execute("""CREATE TABLE asm_stats (ASMID VARCHAR, annotation_status VARCHAR, TOTAL_LENGTH BIGINT, contig_count INTEGER,
                   N50_bp BIGINT, GC_PERCENT DOUBLE, masked_pct DOUBLE)""")
    con.executemany("INSERT INTO asm_stats VALUES (?,?,?,?,?,?,?)", [
        ("G1", "annotated", 10, 1, 60000, 50, 1), ("G2", "annotated", 10, 1, 70000, 50, 1), ("G3", "annotated", 10, 1, 1, 50, 1),
        ("G4", "annotated", 10, 1, 60000, 50, 1), ("G5", "annotated", 10, 1, 60000, 50, 1), ("G6", "unannotated", 10, 1, 1, 50, 1)])
    con.execute("""CREATE TABLE busco_genome (ASMID VARCHAR, complete_pct DOUBLE, single_pct DOUBLE, duplicated_pct DOUBLE,
                   fragmented_pct DOUBLE, missing_pct DOUBLE)""")
    con.executemany("INSERT INTO busco_genome VALUES (?,?,?,?,?,?)", [
        ("G1", 90, 89, 1, 5, 5), ("G2", 95, 94, 1, 3, 2), ("G3", 80, 79, 1, 10, 10),
        ("G4", 70, 60, 10, 15, 15), ("G5", 85, 84, 1, 7, 8), ("G6", 50, 50, 0, 25, 25)])


def build_labels(tmp_path, source_path, labels_path, mtime=123):
    T.write_refs(tmp_path, T_ft(), (), CLS, [("Ascomycota", "Phy-other", "proposed", "")], SP)
    sql = DASH_SQL.read_text()
    first = sql.splitlines()[0]
    assert first.startswith("-- REF_DIR="), "line 1 of dashboard_views.sql must be '-- REF_DIR=<path>'"
    sql = sql.replace(first[len("-- REF_DIR="):], str(tmp_path)).replace("@SOURCE_DB@", str(source_path)).replace("@SOURCE_MTIME@", str(mtime))
    con = duckdb.connect(str(labels_path))
    con.execute(f"ATTACH '{source_path}' AS bfd (READ_ONLY)")
    con.execute(sql)
    con.close()


def T_ft():
    return FT


def test_labels_file_equals_the_old_views(tmp_path):
    src, labels = tmp_path / "src.duckdb", tmp_path / "labels.duckdb"
    s = duckdb.connect(str(src))
    make_source(s)
    s.close()
    build_labels(tmp_path, src, labels)
    new = duckdb.connect(str(labels), read_only=True)
    got = new.execute("SELECT * FROM v_genome_labels ORDER BY LOCUSTAG").df()
    got_tr = new.execute("SELECT * FROM v_genome_traits ORDER BY LOCUSTAG").df()
    old = duckdb.connect()
    make_source(old)
    T.write_refs(tmp_path, FT, (), CLS, [("Ascomycota", "Phy-other", "proposed", "")], SP)
    T.apply_traits(old, tmp_path)
    old.execute(T._labels_statement())
    want = old.execute("SELECT * FROM v_genome_labels ORDER BY LOCUSTAG").df()
    want_tr = old.execute("SELECT * FROM v_genome_traits ORDER BY LOCUSTAG").df()
    pd.testing.assert_frame_equal(got, want)
    pd.testing.assert_frame_equal(got_tr, want_tr)
    assert len(got) == 6 and got["is_representative"].sum() == 4   # A2 (Mucor a), A3, A4, A5; A1 loses to A2, A6 is unannotated


def test_species_placeholder_and_representatives(tmp_path):
    src, labels = tmp_path / "src.duckdb", tmp_path / "labels.duckdb"
    s = duckdb.connect(str(src))
    make_source(s)
    s.close()
    build_labels(tmp_path, src, labels)
    new = duckdb.connect(str(labels), read_only=True)
    r = new.execute("SELECT SUBPHYLUM, subphylum_source FROM v_genome_labels WHERE LOCUSTAG = 'A4'").fetchone()
    assert r == ("PlaceholderSub", "species_table")
    reps = {x[0] for x in new.execute("SELECT LOCUSTAG FROM v_genome_labels WHERE is_representative").fetchall()}
    assert "A2" in reps and "A1" not in reps          # best BUSCO per species


def test_labels_file_stands_alone_and_records_its_source(tmp_path):
    src, labels = tmp_path / "src.duckdb", tmp_path / "labels.duckdb"
    s = duckdb.connect(str(src))
    make_source(s)
    s.close()
    build_labels(tmp_path, src, labels, mtime=4242)
    src.unlink()                                    # no pipeline file: the labels must not need it
    new = duckdb.connect(str(labels), read_only=True)
    assert new.execute("SELECT count(*) FROM v_genome_labels").fetchone()[0] == 6
    meta = new.execute("SELECT source_db, source_mtime_epoch, n_genomes FROM labels_meta").fetchone()
    assert meta == (str(src), 4242, 6)
    kinds = {r[0] for r in new.execute("SELECT DISTINCT table_type FROM information_schema.tables WHERE table_schema='main'").fetchall()}
    assert kinds == {"BASE TABLE"}
