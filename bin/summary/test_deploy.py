"""Tests for the GitHub Pages deployment: workflow, CNAME and the local render."""
import os
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
SITE = ROOT / "docs" / "quarto_summary"
WORKFLOW = ROOT / ".github" / "workflows" / "pages.yml"
DOMAIN = "bfd.fungalgenomes.org"


def _wf():
    return yaml.safe_load(WORKFLOW.read_text())


def _on(wf):
    # PyYAML parses the bare key `on` as boolean True.
    return wf.get("on", wf.get(True))


def _steps(job):
    return _wf()["jobs"][job]["steps"]


def _uses(job):
    return [s["uses"] for s in _steps(job) if "uses" in s]


def _runs(job):
    return "\n".join(s["run"] for s in _steps(job) if "run" in s)


def test_workflow_triggers():
    on = _on(_wf())
    assert on["push"]["branches"] == ["main"]
    assert set(on["push"]["paths"]) == {"docs/quarto_summary/**", ".github/workflows/pages.yml"}
    assert "workflow_dispatch" in on
    pr = on["pull_request"]
    assert set(pr["paths"]) == set(on["push"]["paths"])


def test_deploy_skipped_on_pull_request():
    wf = _wf()
    assert wf["jobs"]["deploy"]["if"] == "github.event_name != 'pull_request'"
    assert "if" not in wf["jobs"]["build"]
    assert wf["jobs"]["build"]["timeout-minutes"] == 30


def test_requirements_pinned_and_used():
    req = (SITE / "requirements.txt").read_text().split()
    pins = {r.split("==")[0].lower(): r.split("==")[1] for r in req}
    for pkg in ("pandas", "matplotlib", "nbformat", "nbclient", "ipykernel", "pyyaml", "tabulate", "numpy"):
        assert pkg in pins and pins[pkg][0].isdigit(), pkg
    assert all("==" in r for r in req)
    py = next(s for s in _steps("build") if s.get("uses", "").startswith("actions/setup-python"))
    assert py["with"]["cache"] == "pip"
    assert py["with"]["cache-dependency-path"] == "docs/quarto_summary/requirements.txt"


def test_workflow_permissions_and_concurrency():
    wf = _wf()
    assert wf["permissions"] == {"contents": "read", "pages": "write", "id-token": "write"}
    assert wf["concurrency"] == {"group": "pages", "cancel-in-progress": False}


def test_build_job_steps():
    uses = _uses("build")
    for u in ("actions/checkout@v4", "actions/setup-python@v5",
              "quarto-dev/quarto-actions/setup@v2", "actions/upload-pages-artifact@v3"):
        assert u in uses
    setup = next(s for s in _steps("build") if s.get("uses", "").startswith("quarto-dev"))
    assert setup["with"]["version"] == "1.7.31"
    py = next(s for s in _steps("build") if s.get("uses", "").startswith("actions/setup-python"))
    assert str(py["with"]["python-version"]) == "3.12"
    runs = _runs("build")
    assert "pip install -r docs/quarto_summary/requirements.txt" in runs
    assert "quarto render docs/quarto_summary" in runs
    assert "CNAME" in runs and "12" in runs  # assertion step
    assert "du -sm" in runs and "900" in runs  # size check


def test_upload_path_matches_render_output():
    runs = _runs("build")
    up = next(s for s in _steps("build") if s.get("uses", "").startswith("actions/upload-pages-artifact"))
    assert up["with"]["path"] == "docs/quarto_summary/_site"
    cfg = yaml.safe_load((SITE / "_quarto.yml").read_text())
    assert cfg["project"]["output-dir"] == "_site"
    assert "--output-dir" not in runs


def test_deploy_job():
    job = _wf()["jobs"]["deploy"]
    assert job["needs"] == "build"
    assert job["environment"]["name"] == "github-pages"
    assert any(s.get("uses") == "actions/deploy-pages@v4" for s in job["steps"])


def test_cname_exact():
    assert (SITE / "CNAME").read_bytes() == (DOMAIN + "\n").encode()


def test_quarto_lists_cname_resource():
    cfg = yaml.safe_load((SITE / "_quarto.yml").read_text())
    assert "CNAME" in cfg["project"]["resources"]


def test_figures_tracked():
    out = subprocess.run(["git", "ls-files", "docs/quarto_summary/figures"], cwd=ROOT,
                         capture_output=True, text=True, check=True).stdout.split()
    assert len(out) > 0


def test_local_render_places_cname(tmp_path):
    quarto = shutil.which("quarto")
    if not quarto or not os.environ.get("QUARTO_PYTHON"):
        pytest.skip("quarto or QUARTO_PYTHON not available; the controller renders manually")
    src = tmp_path / "src" / "quarto_summary"
    shutil.copytree(SITE, src, ignore=shutil.ignore_patterns("_site", ".quarto", "site_libs", "__pycache__"))
    out = tmp_path / "site"
    subprocess.run([quarto, "render", str(src), "--output-dir", str(out)], check=True)
    assert (out / "CNAME").read_text().strip() == DOMAIN
