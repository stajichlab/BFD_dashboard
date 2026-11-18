"""Tests for the summary suppress list (bin/summary/suppress.py and its use in make_figures / export_clade_json)."""
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_figures as mf  # noqa: E402
from export_clade_json import export_clade_json  # noqa: E402
from export_clade_json import main as export_main  # noqa: E402
from suppress import drop_suppressed, load_suppress  # noqa: E402
from test_export_clade_json import _frames  # noqa: E402
from test_make_figures import synth  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
HEADER = "species_prefix\tspecies\treason\tadded\tuntil\n"


def _write(p, body):
    p.write_text(body)
    return p


def test_load_comments_blank_whitespace(tmp_path):
    p = _write(tmp_path / "s.tsv", "# comment\n" + HEADER + "\n#another\n  G0003 \t Sp a \t why \t2026-11-17\t re-run\n\n")
    df = load_suppress(p)
    assert list(df["species_prefix"]) == ["G0003"]
    assert df.loc[0, "species"] == "Sp a" and df.loc[0, "until"] == "re-run"


def test_load_missing_file_is_empty(tmp_path):
    df = load_suppress(tmp_path / "none.tsv")
    assert len(df) == 0 and "species_prefix" in df.columns


def test_load_empty_id_is_error(tmp_path):
    p = _write(tmp_path / "s.tsv", HEADER + "\tSp\twhy\t2026-11-17\tx\n")
    with pytest.raises(ValueError, match="species_prefix"):
        load_suppress(p)


def test_load_duplicates_ok(tmp_path):
    p = _write(tmp_path / "s.tsv", HEADER + "A\ts\tr\td\tu\nA\ts\tr\td\tu\n")
    assert len(load_suppress(p)) == 2


IDS3 = ["FCF75579", "F5C00A22", "F5CD90C7"]


def test_default_path_works_from_any_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert list(load_suppress()["species_prefix"]) == IDS3


def test_drop_suppressed_counts():
    df = pd.DataFrame({"species_prefix": ["a", "b", "c", "b"], "x": range(4)})
    sup = pd.DataFrame({"species_prefix": ["b", "zz", "b"]})
    out, n = drop_suppressed(df, sup)
    assert list(out["species_prefix"]) == ["a", "c"] and n == 2
    out, n = drop_suppressed(df, sup.iloc[0:0])
    assert len(out) == 4 and n == 0


def test_drop_suppressed_without_key_is_noop():
    df = pd.DataFrame({"pfam_acc": ["a"]})
    out, n = drop_suppressed(df, pd.DataFrame({"species_prefix": ["a"]}))
    assert out is df and n == 0


def _tables(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    return data


def _sup(ids):
    return pd.DataFrame({"species_prefix": ids})


def test_main_suppresses_everywhere(tmp_path, capsys):
    data = _tables(tmp_path)
    mf.main(data, tmp_path / "base", ("phylum",), suppress=_sup([]))
    capsys.readouterr()
    mf.main(data, tmp_path / "sup", ("phylum",), suppress=_sup(["G0003"]))
    assert "suppressed 1 genomes: G0003" in capsys.readouterr().out
    idx0 = pd.read_csv(tmp_path / "base" / "figures_index.tsv", sep="\t")
    idx1 = pd.read_csv(tmp_path / "sup" / "figures_index.tsv", sep="\t")
    assert len(idx0) == len(idx1)
    m = idx0.merge(idx1, on=["name", "level"], suffixes=("_0", "_1"))
    d = (m["n_genomes_0"] - m["n_genomes_1"]).dropna()
    assert len(d) > 0 and d.max() <= 1 and d.min() >= 0 and (d == 1).any()
    # no output table lists the suppressed id
    for t in (tmp_path / "sup").rglob("*.tsv"):
        assert "G0003" not in t.read_text(), t


def test_main_a_figure_n_reduced(tmp_path):
    data = _tables(tmp_path)
    mf.main(data, tmp_path / "base", ("phylum",), suppress=_sup([]))
    mf.main(data, tmp_path / "sup", ("phylum",), suppress=_sup(["G0003"]))
    a0 = pd.read_csv(tmp_path / "base" / "figures_index.tsv", sep="\t")
    a1 = pd.read_csv(tmp_path / "sup" / "figures_index.tsv", sep="\t")
    s0 = a0[a0["name"].str.startswith("A")].sort_values("name")
    s1 = a1[a1["name"].str.startswith("A")].sort_values("name")
    assert len(s0) > 0 and list(s0["name"]) == list(s1["name"])
    d = s0["n_genomes"].values - s1["n_genomes"].values
    assert d.max() == 1 and d.min() >= 0


def test_main_default_reads_repo_file(tmp_path, capsys):
    data = _tables(tmp_path)
    mf.main(data, tmp_path / "o", ("phylum",))
    assert "suppressed 0 genomes" in capsys.readouterr().out    # the 3 suppressed ids are not in the synthetic tables


def test_run_cli_suppress_option(tmp_path, capsys):
    data = _tables(tmp_path)
    f = _write(tmp_path / "s.tsv", HEADER + "G0003\ts\tr\td\tu\n")
    assert mf.run_cli(["x", str(data), str(tmp_path / "o1"), "phylum", "--suppress", str(f)]) == 0
    assert "suppressed 1 genomes: G0003" in capsys.readouterr().out
    assert mf.run_cli(["x", str(data), str(tmp_path / "o2"), "phylum", "--suppress", "none"]) == 0
    assert "suppressed 0 genomes" in capsys.readouterr().out


def test_export_drops_and_meta(tmp_path):
    core, st = _frames(50)
    out = export_clade_json(core, st, tmp_path / "g.json", suppress=_sup(["sp_7"]))
    j = json.loads(out.read_text())
    ids = [r[j["fields"].index("id")] for r in j["rows"]]
    assert "sp_7" not in ids and len(ids) == 49
    assert j["meta"]["n_genomes"] == 49 and j["meta"]["n_suppressed"] == 1 and j["meta"]["suppressed"] == ["sp_7"]


def test_export_no_suppress_keeps(tmp_path):
    core, st = _frames(50)
    j = json.loads(export_clade_json(core, st, tmp_path / "g.json", suppress=_sup([])).read_text())
    assert len(j["rows"]) == 50 and j["meta"]["n_suppressed"] == 0 and j["meta"]["suppressed"] == []


def test_export_cli_option(tmp_path):
    data = _tables(tmp_path)
    pd.read_parquet(data / "genome_all.parquet")
    f = _write(tmp_path / "s.tsv", HEADER + "G0003\ts\tr\td\tu\n")
    out = tmp_path / "g.json"
    export_main(["x", str(data), str(out), "--suppress", str(f)])
    assert json.loads(out.read_text())["meta"]["n_suppressed"] == 1
    export_main(["x", str(data), str(out), "--suppress", "none"])
    assert json.loads(out.read_text())["meta"]["n_suppressed"] == 0


def test_committed_config_and_docs_copy():
    cfg = REPO / "config" / "summary_suppress.tsv"
    doc = REPO / "docs" / "quarto_summary" / "data" / "summary_suppress.tsv"
    assert cfg.read_text() == doc.read_text()
    df = load_suppress(cfg)
    assert list(df["species_prefix"]) == IDS3
    assert df.loc[0, "species"] == "Massospora cicadina"
    assert df.loc[2, "species"] == "Massospora platypediae"


def test_committed_genomes_json_has_no_suppressed():
    j = json.loads((REPO / "docs" / "quarto_summary" / "data" / "genomes.json").read_text())
    ids = {r[j["fields"].index("id")] for r in j["rows"]}
    assert not (set(IDS3) & ids) and len(j["rows"]) == 23677
    assert j["meta"]["n_suppressed"] == 3
