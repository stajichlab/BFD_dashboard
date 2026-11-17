"""Tests for the release-asset scripts and the workflow wiring (figures are not kept in git)."""
import hashlib
import json
import os
import shutil
import subprocess
import tarfile
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "workflows"
RELEASE = WF / "release_site_assets.sh"
FETCH = WF / "fetch_site_assets.sh"
PAGES = ROOT / ".github" / "workflows" / "pages.yml"
CONCURRENCY = "pages-${{ github.event_name == 'pull_request' && github.ref || 'deploy' }}"
KEYS = {"tag", "asset", "sha256", "bytes", "created", "n_png", "n_files", "source"}


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _make_out(tmp_path, index=True):
    out = tmp_path / "out"
    (out / "figures" / "phylum").mkdir(parents=True)
    (out / "tables").mkdir()
    (out / "figures" / "phylum" / "A01.png").write_bytes(b"png1" * 100)
    (out / "figures" / "phylum" / "A02.png").write_bytes(b"png2" * 100)
    (out / "figures" / "phylum" / "A01.svg").write_text("<svg/>")
    (out / "figures" / "phylum" / "A01.json").write_text("{}")
    (out / "tables" / "t1.tsv").write_text("a\tb\n1\t2\n")
    (out / "run.log").write_text("not in archive")
    if index:
        (out / "figures_index.tsv").write_text("id\nA01\n")
    return out


def _make_root(tmp_path):
    root = tmp_path / "root"
    (root / "docs" / "quarto_summary").mkdir(parents=True)
    return root


def _release(root, out, dest, *args, date="2026-01-02"):
    env = dict(os.environ, ROOT=str(root), OUT=str(out))
    env["PATH"] = os.environ["PATH"]
    return subprocess.run(["bash", str(RELEASE), "--out", str(dest), "--date", date, *args],
                          env=env, capture_output=True, text=True)


def _fetch(root, *args):
    env = dict(os.environ, ROOT=str(root))
    return subprocess.run(["bash", str(FETCH), *args], env=env, capture_output=True, text=True)


@pytest.fixture
def built(tmp_path):
    out = _make_out(tmp_path)
    root = _make_root(tmp_path)
    dest = tmp_path / "dest"
    r = _release(root, out, dest)
    assert r.returncode == 0, r.stderr
    manifest = json.loads((root / "docs" / "quarto_summary" / "figures_release.json").read_text())
    return tmp_path, out, root, dest, manifest, r


def test_scripts_are_executable_and_sparing():
    for p in (RELEASE, FETCH):
        text = p.read_text()
        assert os.access(p, os.X_OK) and (p.stat().st_mode & 0o755) == 0o755
        assert "BASH_SOURCE" not in text
        assert "stajichlab/BFD_dashboard" in text


def test_release_builds_archive_and_manifest(built):
    tmp_path, out, root, dest, m, r = built
    assert KEYS <= set(m)
    assert m["tag"] == "site-figures-2026-01-02"
    assert m["asset"].startswith("site-figures-2026-01-02.tar.")
    arch = dest / m["asset"]
    assert arch.is_file()
    assert m["sha256"] == _sha(arch)
    assert m["bytes"] == arch.stat().st_size
    assert m["n_png"] == 2
    assert m["n_files"] == 6
    assert m["created"] == "2026-01-02"
    assert set(m["source"]) == {"summary_suppress_sha256", "branch_commit"}
    assert "gh release create site-figures-2026-01-02" in r.stdout
    assert "--repo stajichlab/BFD_dashboard" in r.stdout
    assert str(dest) in r.stdout


def test_archive_contents(built):
    tmp_path, out, root, dest, m, r = built
    x = tmp_path / "x"
    x.mkdir()
    if m["asset"].endswith(".zst"):
        subprocess.run(["tar", "-I", "zstd -d", "-xf", str(dest / m["asset"]), "-C", str(x)], check=True)
    else:
        subprocess.run(["tar", "-xzf", str(dest / m["asset"]), "-C", str(x)], check=True)
    assert sorted(p.name for p in x.iterdir()) == ["figures", "figures_index.tsv", "tables"]
    assert not (x / "run.log").exists()


def test_release_refuses_without_index(tmp_path):
    out = _make_out(tmp_path, index=False)
    root = _make_root(tmp_path)
    r = _release(root, out, tmp_path / "dest")
    assert r.returncode == 1
    assert "figures_index.tsv" in r.stderr
    assert not (root / "docs" / "quarto_summary" / "figures_release.json").exists()


def test_release_refuses_out_inside_repo(tmp_path):
    out = _make_out(tmp_path)
    root = _make_root(tmp_path)
    r = _release(root, out, root / "docs" / "dest")
    assert r.returncode != 0
    assert "inside" in r.stderr


def test_release_dry_run_writes_nothing(tmp_path):
    out = _make_out(tmp_path)
    root = _make_root(tmp_path)
    dest = tmp_path / "dest"
    r = _release(root, out, dest, "--dry-run")
    assert r.returncode == 0, r.stderr
    assert "gh release create site-figures-2026-01-02" in r.stdout
    assert not dest.exists()
    assert not (root / "docs" / "quarto_summary" / "figures_release.json").exists()


def test_release_is_deterministic(tmp_path):
    out = _make_out(tmp_path)
    root = _make_root(tmp_path)
    shas = []
    for i in range(2):
        r = _release(root, out, tmp_path / f"d{i}")
        assert r.returncode == 0, r.stderr
        m = json.loads((root / "docs" / "quarto_summary" / "figures_release.json").read_text())
        shas.append(m["sha256"])
        os.utime(out / "tables" / "t1.tsv", (1, 1 + i))  # change mtime only
    assert shas[0] == shas[1]


def test_upload_requires_gh_auth(tmp_path):
    out = _make_out(tmp_path)
    root = _make_root(tmp_path)
    fake = tmp_path / "bin"
    fake.mkdir()
    (fake / "gh").write_text("#!/bin/sh\nexit 1\n")
    (fake / "gh").chmod(0o755)
    env = dict(os.environ, ROOT=str(root), OUT=str(out), PATH=f"{fake}:{os.environ['PATH']}")
    r = subprocess.run(["bash", str(RELEASE), "--out", str(tmp_path / "d"), "--date", "2026-01-02", "--upload"],
                       env=env, capture_output=True, text=True)
    assert r.returncode != 0
    assert "gh auth" in r.stderr


def test_fetch_extracts_and_removes_stale(built):
    tmp_path, out, root, dest, m, r = built
    site = root / "docs" / "quarto_summary"
    (site / "figures" / "old").mkdir(parents=True)
    (site / "figures" / "old" / "stale.png").write_bytes(b"s")
    (site / "tables").mkdir()
    (site / "tables" / "stale.tsv").write_text("s")
    f = _fetch(root, "--from-archive", str(dest / m["asset"]))
    assert f.returncode == 0, f.stderr
    assert (site / "figures" / "phylum" / "A01.png").read_bytes() == b"png1" * 100
    assert (site / "tables" / "t1.tsv").exists()
    assert (site / "figures_index.tsv").read_text() == "id\nA01\n"
    assert not (site / "figures" / "old").exists()
    assert not (site / "tables" / "stale.tsv").exists()
    assert "2 PNG" in f.stdout


def test_fetch_wrong_sha_exit_3(built):
    tmp_path, out, root, dest, m, r = built
    mp = root / "docs" / "quarto_summary" / "figures_release.json"
    m["sha256"] = "0" * 64
    mp.write_text(json.dumps(m))
    f = _fetch(root, "--from-archive", str(dest / m["asset"]))
    assert f.returncode == 3
    assert "sha256" in f.stderr
    assert not (root / "docs" / "quarto_summary" / "figures").exists()
    f2 = _fetch(root, "--from-archive", str(dest / m["asset"]), "--no-verify")
    assert f2.returncode == 0
    assert "WARNING" in f2.stderr


def test_fetch_corrupted_archive_exit_3(built):
    tmp_path, out, root, dest, m, r = built
    bad = tmp_path / m["asset"]
    data = bytearray((dest / m["asset"]).read_bytes())
    data[len(data) // 2] ^= 0xFF
    bad.write_bytes(bytes(data))
    f = _fetch(root, "--from-archive", str(bad))
    assert f.returncode == 3


def test_fetch_missing_manifest_exit_2(tmp_path):
    root = _make_root(tmp_path)
    f = _fetch(root)
    assert f.returncode == 2
    assert "figures_release.json" in f.stderr


# ---- workflow and repo wiring ----

def _wf():
    return yaml.safe_load(PAGES.read_text())


def _build_steps():
    return _wf()["jobs"]["build"]["steps"]


def _idx(name_part):
    for i, s in enumerate(_build_steps()):
        if name_part in s.get("name", ""):
            return i
    raise AssertionError(name_part)


def test_workflow_fetch_step_before_render():
    i, j = _idx("Fetch figures"), _idx("Render site")
    assert i < j
    assert _idx("Install Python packages") < i
    s = _build_steps()[i]
    assert "workflows/fetch_site_assets.sh" in s["run"]
    assert s["env"]["GH_TOKEN"] == "${{ github.token }}"
    assert "zstd" in s["run"] or "zstd" in "".join(x.get("run", "") for x in _build_steps()[:i])


def test_workflow_asserts_png_count_and_index():
    run = _build_steps()[_idx("Assert site is complete")]["run"]
    assert "figures" in run and "*.png" in run and "200" in run
    assert "figures_index.tsv" in run


def test_workflow_paths_and_concurrency_kept():
    wf = _wf()
    on = wf.get("on", wf.get(True))
    for ev in ("push", "pull_request"):
        assert "docs/quarto_summary/figures_release.json" in on[ev]["paths"]
        assert "docs/quarto_summary/**" in on[ev]["paths"]
    assert wf["concurrency"]["group"] == CONCURRENCY
    assert wf["concurrency"]["cancel-in-progress"] is False
    assert wf["permissions"] == {"contents": "read", "pages": "write", "id-token": "write"}


def test_gitignore_patterns():
    lines = (ROOT / ".gitignore").read_text().split()
    for p in ("docs/quarto_summary/figures/", "docs/quarto_summary/tables/",
              "docs/quarto_summary/figures_index.tsv"):
        assert p in lines


def test_generated_files_not_tracked():
    if not (ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    out = subprocess.run(["git", "ls-files", "docs/quarto_summary/figures", "docs/quarto_summary/tables",
                          "docs/quarto_summary/figures_index.tsv"], cwd=ROOT,
                         capture_output=True, text=True, check=True).stdout
    assert out.strip() == ""


def test_repo_manifest_parses():
    m = json.loads((ROOT / "docs" / "quarto_summary" / "figures_release.json").read_text())
    assert KEYS <= set(m)
    assert m["tag"] == m["asset"].split(".tar")[0]
    assert len(m["sha256"]) == 64
    assert m["n_png"] >= 200
