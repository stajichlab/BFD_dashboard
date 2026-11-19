"""Tests for the four TF pages (13..16), the navbar, the family-map copy and the TF section of the data notes."""
import re
import shutil
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
SITE = ROOT / "docs" / "quarto_summary"
sys.path.insert(0, str(SITE))
import _helpers  # noqa: E402

PAGES = {"13_tf_landscape.qmd": "TF landscape", "14_tf_contrasts.qmd": "TF claims",
         "15_tf_orthologs.qmd": "TF orthologs", "16_tf_expansions.qmd": "TF expansions"}

REQUIRED = {
    "13_tf_landscape.qmd": ["Pfam", "distinct proteins", "protein-coding", "STAT", "SBP", "WOPR", "APSES", "QLQ",
                            "Fungal_trans_2", "unmapped"],
    "14_tf_contrasts.qmd": ["Lax et al. 2025, Nature Communications 16:10277", "one value per genus", "best BUSCO",
                            "BUSCO >= 80", "Benjamini-Hochberg", "phylogeny is not otherwise controlled",
                            "only compares the sign", "InterPro", "62 species", "not comparable", "distinct proteins",
                            "18 of 62", "13 of 19", "sensitivity only", "confound"],
    "15_tf_orthologs.qmd": ["Lax et al.", "T21_target_background_tests", "1e-5", "evidence for orthology, not proof", "function is not inferred", "reciprocal best hit",
                            "MOESM6", "JGI", "UniProt", "MMseqs", "has not been built"],
    "16_tf_expansions.qmd": ["Lax et al.", "distinct proteins", "hypotheses", "p_adj < 0.05", "|effect| >= 0.5", "at least 3 genera", "lineage-wide",
                             "per phylum", "annotation completeness", "not independent", "phylogeny is not controlled"],
}


MIN_NOT_BUILT = {"15_tf_orthologs.qmd": 8, "14_tf_contrasts.qmd": 3, "16_tf_expansions.qmd": 5}


def _read(name):
    return (SITE / name).read_text()


def test_navbar_has_tf_pages_after_explore():
    cfg = yaml.safe_load(_read("_quarto.yml"))
    items = cfg["website"]["navbar"]["left"]
    hrefs = [i["href"] for i in items]
    start = hrefs.index("12_explore.qmd")
    assert hrefs[start + 1:start + 5] == list(PAGES)
    for page, text in PAGES.items():
        assert (SITE / page).is_file()
        assert [i["text"] for i in items if i["href"] == page] == [text]
    assert "CNAME" in cfg["project"]["resources"] and "data/genomes.json" in cfg["project"]["resources"]


@pytest.mark.parametrize("page", list(REQUIRED))
def test_required_statements(page):
    text = _read(page)
    for s in REQUIRED[page]:
        assert s in text, f"{page}: {s!r}"


def test_family_map_copy_identical():
    assert (SITE / "data" / "tf_family_map.tsv").read_bytes() == (ROOT / "config" / "tf_family_map.tsv").read_bytes()


def test_species_map_copy_identical():
    assert (SITE / "data" / "paper_species_map.tsv").read_bytes() == (ROOT / "config" / "paper_species_map.tsv").read_bytes()


def test_data_notes_tf_section():
    t = _read("06_data_notes.qmd")
    for s in ("Transcription-factor tables", "sbatch workflows/tf_tables.sh", "tf_family_genome", "tf_summary_genome",
              "tf_paper_family_genome", "data/tf_family_map.tsv", "38.2"):
        assert s in t, s


def _cells(text):
    return re.findall(r"```\{python\}\n(.*?)```", text, flags=re.S)


@pytest.mark.parametrize("page", list(PAGES))
def test_pages_run_without_tf_outputs(page, tmp_path, monkeypatch, capsys):
    """Copy the tracked tree only (no figures or tables): every python cell must run and show 'not built' notes."""
    pytest.importorskip("tabulate")   # DataFrame.to_markdown, as in the existing helpers
    dst = tmp_path / "quarto_summary"
    shutil.copytree(SITE, dst, ignore=shutil.ignore_patterns("figures", "tables", "figures_index.tsv", "__pycache__",
                                                             "_site", ".quarto"))
    monkeypatch.chdir(dst)
    monkeypatch.syspath_prepend(str(dst))
    for m in [k for k in sys.modules if k == "_helpers"]:
        del sys.modules[m]
    ns = {"__name__": "__main__"}
    for cell in _cells((dst / page).read_text()):
        exec(compile(cell, page, "exec"), ns)
    out = capsys.readouterr().out
    sys.modules.pop("_helpers", None)
    sys.path.insert(0, str(SITE))
    import _helpers as _h  # noqa: F401  restore for the other tests
    n = out.count("has not been built")
    assert n >= MIN_NOT_BUILT.get(page, 1), (page, n)


def test_flat_table_skips_leading_comment_lines(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(_helpers, "HERE", tmp_path)
    d = tmp_path / "tables" / "tf"
    d.mkdir(parents=True)
    (d / "T_c.tsv").write_text("# min_cov=0.5 evalue=1e-05\n# second\na\tb\n1\t2\n")
    df = _helpers.read_flat("T_c", "tf")
    assert list(df.columns) == ["a", "b"] and len(df) == 1
    assert _helpers.read_comments("T_c", "tf") == ["min_cov=0.5 evalue=1e-05", "second"]
    _helpers.flat_table("T_c", "tf")
    out = capsys.readouterr().out
    assert "min_cov=0.5" in out and "|   a |" in out
    assert _helpers.read_comments("absent", "tf") == []


def test_gallery_lists_tf_figures(tmp_path):
    idx = tmp_path / "figures_index.tsv"
    head = ("name\tlevel\ttopic\tkind\tscope\tn_genomes\tn_groups\tpooled_share\teffect_size\tscore\tsignal\tsignal_kind\n")
    rows = ["T01_tf_fraction\tphylum\ttf\tfigure\trepresentatives\t7000\t9\t0.1\t\t\t\t",
            "T06_tf_vs_proteome\tall\ttf\tfigure\tall_genomes\t20000\t\t0\t\t\t\t",
            "T10_tf_claims_forest\ttf\ttf\tfigure\trepresentatives\t3000\t\t\t\t\t\t"]
    idx.write_text(head + "\n".join(rows) + "\n")
    h = _helpers.gallery_html(index_path=idx)
    for n in ("T01_tf_fraction", "T06_tf_vs_proteome", "T10_tf_claims_forest"):
        assert n in h
    assert "figures/all/T06_tf_vs_proteome.png" in h and "figures/tf/T10_tf_claims_forest.png" in h
