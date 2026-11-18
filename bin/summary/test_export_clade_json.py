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
    meta = json.loads(g.read_text())["meta"]
    # genome_all.parquet keeps the suppressed genomes (suppression is applied at read time), the JSON drops them.
    assert meta["n_genomes"] == len(pd.read_parquet(p)) - meta["n_suppressed"]


# ---- species and strain labels ----
def _strains(core, skip=(3,)):
    ids = [x for i, x in enumerate(core.species_prefix) if i not in skip]
    return pd.DataFrame({"species_prefix": ids, "STRAIN": [f"ST{i}" for i in range(len(ids))]})


def test_species_and_strain_fields(tmp_path):
    core, st = _frames()
    core["SPECIES"] = [f"G{i} sp{i}" for i in range(len(core))]
    core.loc[2, "SPECIES"] = np.nan
    core.loc[3, "SPECIES"] = ""
    core.loc[4, ["SPECIES", "GENUS"]] = np.nan
    j = json.loads(export_clade_json(core, st, tmp_path / "g.json", strains=_strains(core)).read_text())
    assert "species" in j["fields"] and "strain" in j["fields"] and "id" in j["fields"]
    r = _records(j)
    assert r[5]["species"] == "G5 sp5"
    assert r[2]["species"] == f"{core.GENUS[2]} sp."
    assert r[3]["species"] == f"{core.GENUS[3]} sp."
    assert r[4]["species"] == "Unknown"
    assert r[0]["species"] != ""
    assert r[5]["strain"].startswith("ST")
    assert r[3]["strain"] == ""
    assert len({x["id"] for x in r}) == len(r)
    assert all(x["species"] for x in r)


def test_strains_none_gives_empty_strain(tmp_path):
    _, _, j = _run(tmp_path)
    assert all(x["strain"] == "" for x in _records(j))


def test_strain_null_is_empty_string_not_null(tmp_path):
    core, st = _frames()
    sd = _strains(core)
    sd.loc[0, "STRAIN"] = None
    j = json.loads(export_clade_json(core, st, tmp_path / "g.json", strains=sd).read_text())
    assert all(isinstance(x["strain"], str) for x in _records(j))


def test_label_fields_not_numeric_columns(tmp_path):
    _, _, j = _run(tmp_path)
    assert "species" not in j["columns"] and "strain" not in j["columns"]
    assert j["label_columns"]["species"]["label"] and j["label_columns"]["strain"]["label"]


def test_no_nan_with_species(tmp_path):
    core, st = _frames()
    core.loc[2, "SPECIES"] = np.nan
    t = export_clade_json(core, st, tmp_path / "g.json", strains=_strains(core)).read_text()
    assert "NaN" not in t and "None" not in t


def test_size_under_5mb_with_species(tmp_path):
    core, st = _frames(23680)
    core["SPECIES"] = [f"Genus{i % 3000} species{i % 7000}" for i in range(len(core))]
    out = export_clade_json(core, st, tmp_path / "g.json", strains=_strains(core))
    assert out.stat().st_size < 5_000_000


def _cli_tables(tmp_path):
    core, st = _frames(20)
    t = tmp_path / "t"
    t.mkdir()
    core.to_parquet(t / "genome_all.parquet")
    st.to_parquet(t / "gene_structure_genome.parquet")
    return t


def test_cli_species_table_none_and_missing(tmp_path):
    import subprocess
    t = _cli_tables(tmp_path)
    script = Path(__file__).resolve().parent / "export_clade_json.py"
    for val in ("none", str(tmp_path / "nope.parquet")):
        out = tmp_path / "o.json"
        r = subprocess.run([sys.executable, str(script), str(t), str(out), "--species-table", val],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        j = json.loads(out.read_text())
        assert all(x["strain"] == "" for x in _records(j))
        if val != "none":
            assert "nope.parquet" in r.stderr


def test_cli_species_table_reads_strains(tmp_path):
    import subprocess
    t = _cli_tables(tmp_path)
    core = pd.read_parquet(t / "genome_all.parquet")
    pd.DataFrame({"LOCUSTAG": core.species_prefix, "SPECIES": "x", "STRAIN": "S1"}).to_parquet(tmp_path / "sp.parquet")
    script = Path(__file__).resolve().parent / "export_clade_json.py"
    out = tmp_path / "o.json"
    r = subprocess.run([sys.executable, str(script), str(t), str(out), "--species-table", str(tmp_path / "sp.parquet")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert all(x["strain"] == "S1" for x in _records(json.loads(out.read_text())))
