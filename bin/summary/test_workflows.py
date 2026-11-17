"""Tests for the SLURM scripts and the publish script in workflows/."""
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "workflows"
PUBLISH = WF / "publish_figures_to_site.sh"


def _run(root, out, *args):
    env = dict(os.environ, ROOT=str(root), OUT=str(out))
    return subprocess.run(["bash", str(PUBLISH), *args], env=env, capture_output=True, text=True)


def _make_out(tmp_path, index=True):
    out = tmp_path / "out"
    (out / "figures" / "phylum").mkdir(parents=True)
    (out / "tables").mkdir()
    (out / "figures" / "phylum" / "A01.png").write_bytes(b"x" * 2048)
    (out / "tables" / "t.tsv").write_text("a\tb\n1\t2\n")
    if index:
        (out / "figures_index.tsv").write_text("id\nA01\n")
    root = tmp_path / "repo"
    (root / "docs" / "quarto_summary").mkdir(parents=True)
    return root, out


def test_publish_copies_and_keeps_extras(tmp_path):
    root, out = _make_out(tmp_path)
    site = root / "docs" / "quarto_summary"
    (site / "figures").mkdir()
    (site / "figures" / "extra.png").write_bytes(b"keep")
    r = _run(root, out)
    assert r.returncode == 0, r.stderr
    assert (site / "figures" / "phylum" / "A01.png").stat().st_size == 2048
    assert (site / "tables" / "t.tsv").exists()
    assert (site / "figures_index.tsv").read_text() == "id\nA01\n"
    assert (site / "figures" / "extra.png").exists()
    assert "copied 3 files" in r.stdout and "MB" in r.stdout


def test_publish_refuses_without_index(tmp_path):
    root, out = _make_out(tmp_path, index=False)
    r = _run(root, out)
    assert r.returncode != 0
    assert "figures_index.tsv" in r.stderr
    assert not (root / "docs" / "quarto_summary" / "figures").exists()


def test_publish_dry_run_copies_nothing(tmp_path):
    root, out = _make_out(tmp_path)
    r = _run(root, out, "--dry-run")
    assert r.returncode == 0, r.stderr
    site = root / "docs" / "quarto_summary"
    assert not (site / "figures_index.tsv").exists()
    assert not list((site / "figures").rglob("*.png"))
    assert "dry run" in r.stdout


def test_slurm_scripts_defaults():
    for name in ("summary_figures.sh", "summary_tables.sh", "publish_figures_to_site.sh"):
        p = WF / name
        text = p.read_text()
        assert os.access(p, os.X_OK) and (p.stat().st_mode & 0o755) == 0o755
        assert "BASH_SOURCE" not in text
        assert ".claude/worktrees" not in text
    for name in ("summary_figures.sh", "summary_tables.sh"):
        text = (WF / name).read_text()
        assert "ROOT=${ROOT:-${SLURM_SUBMIT_DIR:?" in text
        assert "#SBATCH -o %x.%j.log" in text
