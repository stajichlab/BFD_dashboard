"""Tests for export_clade_json: the JSON for the Explore page."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from export_clade_json import export_clade_json  # noqa: E402

LEVELS = ("phylum", "subphylum", "class", "order")


def _records(j):
    return [dict(zip(j["fields"], r)) for r in j["rows"]]


def _frames(n=200, seed=0):
    rng = np.random.default_rng(seed)
    ph = rng.choice(["Ascomycota", "Basidiomycota", "Mucoromycota"], n)
    cl = [f"{p[:3]}_class{rng.integers(0, 4)}" for p in ph]
    core = pd.DataFrame({
        "species_prefix": [f"sp_{i}" for i in range(n)],
        "PHYLUM": ph, "CLASS": cl, "ORDER": [f"{c}_ord{rng.integers(0, 3)}" for c in cl],
        "FAMILY": "F", "GENUS": [f"G{rng.integers(0, 30)}" for _ in range(n)], "SPECIES": "s",
        "is_representative": rng.random(n) < 0.3,
        "total_length_bp": rng.uniform(1e7, 9e7, n), "gc_pct": rng.uniform(30, 60, n),
        "busco_complete": rng.uniform(50, 100, n), "N50_bp": rng.uniform(1e4, 1e6, n),
        "n_genes": rng.integers(5000, 15000, n).astype(float), "genes_per_mb": rng.uniform(100, 500, n),
        "coding_fraction": rng.uniform(0.3, 0.7, n), "median_intron_len_bp": rng.uniform(40, 120, n),
        "median_intergenic_bp": rng.uniform(300, 2000, n),
    })
    core.loc[0, ["PHYLUM", "CLASS", "ORDER", "GENUS"]] = np.nan      # NULL taxonomy
    core.loc[1, ["gc_pct", "n_genes"]] = np.nan                       # NaN numbers
    st = pd.DataFrame({"species_prefix": core.species_prefix[2:], "median_protein_len_aa": rng.uniform(300, 500, n - 2)})
    return core, st


def _run(tmp_path, n=200):
    core, st = _frames(n)
    out = export_clade_json(core, st, tmp_path / "genomes.json")
    return core, out, json.loads(out.read_text())


def test_parses_and_counts(tmp_path):
    core, out, j = _run(tmp_path)
    assert len(_records(j)) == len(core)
    assert len({r["id"] for r in _records(j)}) == len(core)


def test_unknown_taxonomy_and_nulls(tmp_path):
    _, out, j = _run(tmp_path)
    r = _records(j)[0]
    assert (r["phylum"], r["class"], r["order"], r["genus"]) == ("Unknown",) * 4
    assert _records(j)[1]["gc"] is None and _records(j)[1]["n_genes"] is None
    assert _records(j)[0]["median_protein_len_aa"] is None      # absent from structure
    assert _records(j)[5]["median_protein_len_aa"] is not None


def test_no_nan_literal(tmp_path):
    _, out, _ = _run(tmp_path)
    t = out.read_text()
    for bad in ("NaN", "Infinity"):
        assert bad not in t


def test_rounding(tmp_path):
    _, out, j = _run(tmp_path)
    for r in _records(j):
        for k, v in r.items():
            if isinstance(v, float):
                assert round(v, 3) == v, (k, v)


def test_colors_cover_every_group(tmp_path):
    _, out, j = _run(tmp_path)
    assert set(j["colors"]) == set(LEVELS)
    for lv in LEVELS:
        for r in _records(j):
            assert r[lv] in j["colors"][lv], (lv, r[lv])
        for c in j["colors"][lv].values():
            assert c.startswith("#") and len(c) == 7


def test_booleans_and_meta(tmp_path):
    core, out, j = _run(tmp_path)
    assert all(isinstance(r["is_representative"], bool) for r in _records(j))
    assert "genome_mb" in j["columns"] and "label" in j["columns"]["genome_mb"] and "unit" in j["columns"]["genome_mb"]
    assert set(j["levels"]) == set(LEVELS)


def test_missing_is_representative_is_null(tmp_path):
    core, st = _frames(20)
    core["is_representative"] = core["is_representative"].astype(object)
    core.loc[3, "is_representative"] = None
    j = json.loads(export_clade_json(core, st, tmp_path / "g.json").read_text())
    assert _records(j)[3]["is_representative"] is None


def test_size_under_5mb_for_full_table(tmp_path):
    _, out, _ = _run(tmp_path, n=23680)
    assert out.stat().st_size < 5_000_000


def test_meta_fields_and_missing_columns(tmp_path):
    core, out, j = _run(tmp_path)
    assert j["meta"]["n_genomes"] == len(core) and "generated_from" in j["meta"]
    assert j["meta"]["rows_layout"] == "array of arrays; column order = fields"
    assert not j["meta"]["generated_from"].startswith("/")
    import pytest
    core, st = _frames(10)
    with pytest.raises(ValueError, match="gc_pct.*N50_bp|N50_bp.*gc_pct"):
        export_clade_json(core.drop(columns=["gc_pct", "N50_bp"]), st, tmp_path / "x.json")


def test_meta_matches_real_table():
    import pytest
    p = Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_tables/genome_all.parquet")
    g = Path(__file__).resolve().parents[2] / "docs" / "quarto_summary" / "data" / "genomes.json"
    if not p.exists():
        pytest.skip("real tables not available")
    assert json.loads(g.read_text())["meta"]["n_genomes"] == len(pd.read_parquet(p))
