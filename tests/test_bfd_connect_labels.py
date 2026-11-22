"""lib/bfd_connect.py: optional labels file attached next to the pipeline database."""
import os
import sys
import warnings
from pathlib import Path

import duckdb
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import bfd_connect as bc  # noqa: E402


def _dbs(tmp_path, source_mtime=None):
    bfd = tmp_path / "BFD.duckdb"
    c = duckdb.connect(str(bfd))
    c.execute("CREATE TABLE species (LOCUSTAG VARCHAR)")
    c.execute("INSERT INTO species VALUES ('A1'),('A2')")
    c.close()
    lab = tmp_path / "labels.duckdb"
    c = duckdb.connect(str(lab))
    c.execute("CREATE TABLE v_genome_labels (LOCUSTAG VARCHAR, primary_lifestyle VARCHAR)")
    c.execute("INSERT INTO v_genome_labels VALUES ('A1','x'),('A2','y')")
    c.execute("CREATE TABLE labels_meta (source_db VARCHAR, source_mtime_epoch BIGINT, n_genomes BIGINT)")
    mt = int(bfd.stat().st_mtime) if source_mtime is None else source_mtime
    c.execute("INSERT INTO labels_meta VALUES (?, ?, 2)", [str(bfd), mt])
    c.close()
    return bfd, lab


def test_default_connection_has_no_labels(tmp_path):
    bfd, _ = _dbs(tmp_path)
    con = bc.connect(bfd, scratch=False)
    assert con.execute("SELECT count(*) FROM species").fetchone()[0] == 2
    with pytest.raises(duckdb.Error):
        con.execute("SELECT * FROM v_genome_labels")


def test_labels_resolve_unqualified_next_to_bfd(tmp_path):
    bfd, lab = _dbs(tmp_path)
    con = bc.connect(bfd, scratch=False, labels_db=lab)
    assert con.execute("SELECT count(*) FROM v_genome_labels").fetchone()[0] == 2     # from labels
    assert con.execute("SELECT count(*) FROM species").fetchone()[0] == 2             # from bfd
    assert con.execute("SELECT count(*) FROM labels.v_genome_labels").fetchone()[0] == 2


def test_labels_file_is_read_only(tmp_path):
    bfd, lab = _dbs(tmp_path)
    con = bc.connect(bfd, scratch=False, labels_db=lab)
    with pytest.raises(duckdb.Error):
        con.execute("INSERT INTO labels.v_genome_labels VALUES ('A3','z')")


def test_env_variable_selects_the_labels_file(tmp_path, monkeypatch):
    bfd, lab = _dbs(tmp_path)
    monkeypatch.setenv("BFD_LABELS_DB", str(lab))
    con = bc.connect(bfd, scratch=False)
    assert con.execute("SELECT count(*) FROM v_genome_labels").fetchone()[0] == 2


def test_missing_labels_file_raises_and_names_the_build_script(tmp_path):
    bfd, _ = _dbs(tmp_path)
    with pytest.raises(FileNotFoundError, match="build_labels_db.sh"):
        bc.connect(bfd, scratch=False, labels_db=tmp_path / "nope.duckdb")


def test_stale_labels_warn(tmp_path):
    bfd, lab = _dbs(tmp_path, source_mtime=1)           # labels built from a much older pipeline file
    with pytest.warns(UserWarning, match="older than"):
        bc.connect(bfd, scratch=False, labels_db=lab)


def test_fresh_labels_do_not_warn(tmp_path):
    bfd, lab = _dbs(tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        bc.connect(bfd, scratch=False, labels_db=lab)


def test_default_pipeline_db_attaches_default_labels(tmp_path, monkeypatch):
    bfd, lab = _dbs(tmp_path)
    monkeypatch.delenv("BFD_DB", raising=False)
    monkeypatch.delenv("BFD_LABELS_DB", raising=False)
    monkeypatch.setattr(bc, "DEFAULT_BFD_DB", bfd)
    monkeypatch.setattr(bc, "DEFAULT_LABELS_DB", lab)
    con = bc.connect(scratch=False)
    assert con.execute("SELECT primary_lifestyle FROM v_genome_labels WHERE LOCUSTAG='A1'").fetchone()[0] == "x"


def test_named_bfd_db_does_not_attach_default_labels(tmp_path, monkeypatch):
    bfd, lab = _dbs(tmp_path)
    monkeypatch.delenv("BFD_LABELS_DB", raising=False)
    monkeypatch.setattr(bc, "DEFAULT_LABELS_DB", lab)
    con = bc.connect(bfd, scratch=False)
    with pytest.raises(duckdb.Error):
        con.execute("SELECT * FROM v_genome_labels")


def test_empty_env_turns_labels_off(tmp_path, monkeypatch):
    bfd, lab = _dbs(tmp_path)
    monkeypatch.delenv("BFD_DB", raising=False)
    monkeypatch.setenv("BFD_LABELS_DB", "")
    monkeypatch.setattr(bc, "DEFAULT_BFD_DB", bfd)
    monkeypatch.setattr(bc, "DEFAULT_LABELS_DB", lab)
    con = bc.connect(scratch=False)
    with pytest.raises(duckdb.Error):
        con.execute("SELECT * FROM v_genome_labels")
