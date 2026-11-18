"""Tests for the TF landscape figures (figs_tf.py) and the TF data helpers (tf_data.py). Fixtures are synthetic.

Run on the HPCC:  python3.12 -m pytest -q -W error::UserWarning bin/summary/test_figs_tf.py
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_figures as mf  # noqa: E402
import tf_data  # noqa: E402
from figs_tf import FIGURES, RARE_ACCESSIONS, page_tf_landscape  # noqa: E402
from test_make_figures import synth  # noqa: E402

NAMES = ["T01_tf_fraction", "T02_tf_composition", "T03_tf_family_heatmap", "T04_tf_rare_presence", "T05_tf_shannon",
         "T06_tf_vs_proteome"]
FAMS = ["Zn(II)2Cys6 zinc cluster", "C2H2 zinc finger", "Basic leucine zipper", "Homeodomain", "GATA zinc finger",
        "Myb-like DNA-binding domain", "WRKY", "E2F_TDP", "Tub"]
RARE = {"WRKY", "E2F_TDP", "Tub"}


def tf_synth(root, seed=1):
    """TF tables for the genomes of genome_all: a few empty proteomes, zero-TF genomes and one-family genomes."""
    rng = np.random.default_rng(seed)
    ga = pd.read_parquet(root / "genome_all.parquet")
    ids = ga["species_prefix"].tolist()
    rows, summ = [], []
    for k, g in enumerate(ids):
        n_total = int(ga.loc[ga.species_prefix == g, "n_genes"].iloc[0])
        if k % 50 == 7:   # empty proteome
            summ.append((g, 0, 0, np.nan, 0, np.nan))
            continue
        if k % 50 == 11:  # proteome, no TF
            summ.append((g, n_total, 0, 0.0, 0, np.nan))
            continue
        fams = FAMS[:1] if k % 50 == 13 else [f for f in FAMS if f not in RARE or rng.random() < 0.4]
        counts = {f: int(rng.integers(1, 120)) for f in fams}
        rows += [(g, f, c) for f, c in counts.items()]
        tot = sum(counts.values())
        p = np.array(list(counts.values())) / tot
        shannon = 0.0 if len(counts) == 1 else float(-(p * np.log(p)).sum())
        summ.append((g, n_total, tot, tot / n_total, len(counts), shannon))
    fam = pd.DataFrame(rows, columns=["species_prefix", "family", "n_proteins"])
    summ = pd.DataFrame(summ, columns=["species_prefix", "n_proteins_total", "n_tf_proteins", "tf_fraction",
                                       "n_families", "shannon"])
    fam.to_parquet(root / "tf_family_genome.parquet")
    summ.to_parquet(root / "tf_summary_genome.parquet")
    return fam, summ


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    root = tmp_path_factory.mktemp("data")
    synth(root)
    fam, summ = tf_synth(root)
    return root, fam, summ


def _fig(out, level, name, ext="png"):
    return Path(out) / "figures" / level / f"{name}.{ext}"


NO_SUPPRESS = pd.DataFrame({"species_prefix": []})


def test_all_tf_figures_written(data, tmp_path):
    root, _, _ = data
    made = mf.main(root, tmp_path / "out", ("phylum", "order"), suppress=NO_SUPPRESS)
    for level in ("phylum", "order"):
        for name in NAMES:
            for ext in ("png", "svg"):
                assert _fig(tmp_path / "out", level, name, ext).stat().st_size > 1000, (level, name, ext)
        by = {m: m for m in made[level]}
        assert {by[n].topic for n in NAMES} == {"tf"}
        assert by["T06_tf_vs_proteome"].scope == "all_genomes"
        assert by["T01_tf_fraction"].scope == "representatives"
        assert not any(getattr(m, "kind", "") == "skipped" for m in made[level])
    idx = pd.read_csv(tmp_path / "out" / "figures_index.tsv", sep="\t")
    assert set(NAMES) <= set(idx["name"])
    assert set(idx.loc[idx["name"].isin(NAMES), "topic"]) == {"tf"}
    assert FIGURES == NAMES
    for p in (tmp_path / "out" / "figures").rglob("*.*"):
        assert re.fullmatch(r"[A-Za-z0-9_]+\.(png|svg|json)", p.name), p


def test_main_skips_without_tf_tables(tmp_path, capsys):
    root = tmp_path / "data"
    root.mkdir()
    synth(root)   # no TF tables
    made = mf.main(root, tmp_path / "out", ("phylum", "order"), suppress=NO_SUPPRESS)
    out = capsys.readouterr().out
    assert out.count("TF landscape skipped") == 1
    assert "F01_genomes_by_group" in made["phylum"]
    assert not any(m in NAMES or m == "T_tf_groups" for lv in ("phylum", "order") for m in made[lv])


def _ga(rows):
    return pd.DataFrame(rows, columns=["species_prefix", "GENUS", "PHYLUM", "CLASS", "ORDER", "busco_complete", "N50_bp"])


def _tabs(ids):
    fam = pd.DataFrame([(i, f, k + 1) for i in ids for k, f in enumerate(["A fam", "B fam"])],
                       columns=["species_prefix", "family", "n_proteins"])
    summ = pd.DataFrame({"species_prefix": ids, "n_proteins_total": 1000, "n_tf_proteins": 10, "tf_fraction": 0.01,
                         "n_families": 2, "shannon": 0.6})
    return fam, summ


def test_genus_table_one_row_per_genus_and_best_busco():
    ga = _ga([("A1", "G1", "P", "C", "O", 90.0, 100), ("A2", "G1", "P", "C", "O", 95.0, 10),   # best busco wins
              ("B1", "G2", "P", "C", "O", 90.0, 100), ("B2", "G2", "P", "C", "O", 90.0, 200),   # tie busco: N50
              ("C1", "G3", "P", "C", "O", 90.0, 100), ("C0", "G3", "P", "C", "O", 90.0, 100),   # tie: lowest id
              ("N1", None, "P", "C", "O", 99.0, 5),
              ("D1", "G4", "P", "C", "O", 99.0, 5)])    # best genome of G4 has no TF data
    fam, summ = _tabs(["A1", "A2", "B1", "B2", "C1", "C0", "N1"])
    t = tf_data.genus_table(ga, fam, summ)
    assert list(t["genus"]) == ["G1", "G2", "G3"]
    assert dict(zip(t["genus"], t["species_prefix"])) == {"G1": "A2", "G2": "B2", "G3": "C0"}
    assert t.attrs["n_no_genus"] == 1
    for c in ("genus", "species_prefix", "PHYLUM", "CLASS", "ORDER", "busco_complete", "n_proteins_total",
              "n_tf_proteins", "tf_fraction", "shannon", "fam_A fam", "fam_B fam"):
        assert c in t.columns, c
    assert (t["fam_B fam"] == 2).all() and (t["fam_A fam"] == 1).all()


def test_genus_table_zero_filled_and_complete_data():
    ga = _ga([("A1", "G1", "P", "C", "O", 99.0, 100), ("A2", "G1", "P", "C", "O", 80.0, 100)])
    fam = pd.DataFrame({"species_prefix": ["A2"], "family": ["B fam"], "n_proteins": [4]})
    summ = pd.DataFrame({"species_prefix": ["A1", "A2"], "n_proteins_total": [0, 900], "n_tf_proteins": [0, 4],
                         "tf_fraction": [np.nan, 4 / 900], "n_families": [0, 1], "shannon": [np.nan, 0.0]})
    t = tf_data.genus_table(ga, fam, summ)
    assert list(t["species_prefix"]) == ["A2"]   # A1 has an empty proteome: not complete data
    assert t.loc[0, "fam_B fam"] == 4


def test_genus_table_busco_filter():
    ga = _ga([("A1", "G1", "P", "C", "O", 60.0, 100), ("B1", "G2", "P", "C", "O", 90.0, 100),
              ("B2", "G2", "P", "C", "O", 70.0, 100), ("X1", None, "P", "C", "O", 99.0, 1)])
    fam, summ = _tabs(["A1", "B1", "B2", "X1"])
    t = tf_data.genus_table(ga, fam, summ, busco_min=80)
    assert list(t["genus"]) == ["G2"] and list(t["species_prefix"]) == ["B1"]
    assert t.attrs["n_no_genus"] == 1
    t0 = tf_data.genus_table(ga, fam, summ)
    assert list(t0["genus"]) == ["G1", "G2"]


def test_load_tf_applies_suppress(data):
    root, fam, summ = data
    sup = pd.DataFrame({"species_prefix": [summ["species_prefix"].iloc[3]]})
    f, s = tf_data.load_tf(root, sup)
    assert sup["species_prefix"].iloc[0] not in set(f["species_prefix"]) | set(s["species_prefix"])
    assert len(s) == len(summ) - 1
    f2, s2 = tf_data.load_tf(root, NO_SUPPRESS)
    assert len(s2) == len(summ)


LONG1 = "Fungal-specific transcription factor domain (C-terminal to Zn_clus)"
LONG2 = "Fungal-specific transcription factor domain 2 (Fungal_trans_2)"


def _handmade(tmp_path, n_single=True, shannon_nan=False):
    """Direct call of the page: one phylum with a one-genome group, a zero-TF and an empty-proteome genome."""
    ids = [f"S{i:02d}" for i in range(13)]
    grp = ["Big"] * 10 + ["Single", "Unknown", "Unknown"]
    core = pd.DataFrame({"species_prefix": ids, "grp": grp, "PHYLUM": "Ascomycota", "gc_pct": 50.0})
    levels = ["Big", "Single", "Unknown"]
    colors = {"Big": "#0072B2", "Single": "#0072B2", "Unknown": "#d4d4d4"}
    rng = np.random.default_rng(3)
    fam = pd.DataFrame([(i, f, int(rng.integers(1, 40))) for i in ids[:11] for f in ["Zn(II)2Cys6 zinc cluster", "WRKY", LONG1, LONG2]],
                       columns=["species_prefix", "family", "n_proteins"])
    summ = pd.DataFrame({"species_prefix": ids, "n_proteins_total": [900] * 11 + [0, 500],
                         "n_tf_proteins": [20] * 11 + [0, 0], "tf_fraction": [0.02] * 11 + [np.nan, 0.0],
                         "n_families": [2] * 11 + [0, 0],
                         "shannon": ([np.nan] * 13 if shannon_nan else [0.5] * 11 + [np.nan, np.nan])})
    ga = pd.DataFrame({"species_prefix": ids, "PHYLUM": "Ascomycota"})
    o = mf.Out(tmp_path, "phylum")
    o.set_groups(core, levels)
    page_tf_landscape(core, levels, colors, fam, summ, o, genome_all=ga)
    return o


def test_nan_and_empty_groups_do_not_crash(tmp_path):
    o = _handmade(tmp_path)
    kinds = {m: m.kind for m in o.made}
    for n in NAMES:
        assert kinds[n] == "figure", (n, kinds)
    # a rare family absent everywhere (E2F_TDP) is a zero row, not an error
    assert (tmp_path / "figures" / "phylum" / "T04_tf_rare_presence.png").stat().st_size > 1000


def test_all_nan_column_is_skipped_not_raised(tmp_path):
    o = _handmade(tmp_path, shannon_nan=True)
    kinds = {m: m.kind for m in o.made}
    assert kinds["T05_tf_shannon"] == "skipped"
    assert kinds["T01_tf_fraction"] == "figure"
    assert not (tmp_path / "figures" / "phylum" / "T05_tf_shannon.png").exists()


def test_zero_tf_genome_is_nan_fraction_safe(tmp_path):
    o = _handmade(tmp_path)
    t = pd.read_csv(tmp_path / "tables" / "phylum" / "T_tf_groups.tsv", sep="\t", keep_default_na=False, na_values=[""])
    u = t.set_index("group").loc["Unknown"]
    assert u["n_genomes"] == 2 and u["n_with_proteome"] == 1 and u["n_with_tf"] == 0
    assert u["median_tf_fraction"] == 0.0           # the empty proteome (NaN) is left out, the zero-TF genome counts
    assert np.isnan(u["median_shannon"])            # no family at all: NaN, never 0
    assert t.set_index("group").loc["Single", "n_genomes"] == 1


def test_unknown_group_and_group_names_not_in_file_names(tmp_path):
    _handmade(tmp_path)
    names = {p.name for p in (tmp_path / "figures" / "phylum").iterdir()}
    assert all(re.fullmatch(r"[A-Za-z0-9_]+\.(png|svg|json)", n) for n in names)
    assert not any("Unknown" in n or "Single" in n for n in names)


def test_rare_families_resolved_from_map_by_accession():
    labels = tf_data.rare_family_labels(RARE_ACCESSIONS)
    assert "E2F_TDP" in labels and "Tub" in labels and "STAT" in labels and "WRKY" in labels
    assert len(labels) == len(set(labels)) == 6


def test_suppressed_genome_absent_from_tf_outputs(data, tmp_path):
    root, _, summ = data
    gone = pd.read_parquet(root / "genome_core.parquet")["species_prefix"].iloc[5]
    sup = pd.DataFrame({"species_prefix": [gone]})
    mf.main(root, tmp_path / "out", ("phylum",), suppress=sup)
    mf.main(root, tmp_path / "base", ("phylum",), suppress=NO_SUPPRESS)
    t = pd.read_csv(tmp_path / "out" / "tables" / "phylum" / "T_tf_groups.tsv", sep="\t")
    b = pd.read_csv(tmp_path / "base" / "tables" / "phylum" / "T_tf_groups.tsv", sep="\t")
    assert t["n_genomes"].sum() == b["n_genomes"].sum() - 1
    # the all-genome scatter drops it as well
    s = pd.read_json(tmp_path / "out" / "figures" / "phylum" / "T06_tf_vs_proteome.json", typ="series")
    sb = pd.read_json(tmp_path / "base" / "figures" / "phylum" / "T06_tf_vs_proteome.json", typ="series")
    assert s["name"] == sb["name"]


def test_deterministic_output(data, tmp_path):
    root, _, _ = data
    mf.main(root, tmp_path / "a", ("order",), suppress=NO_SUPPRESS)
    mf.main(root, tmp_path / "b", ("order",), suppress=NO_SUPPRESS)
    for n in NAMES:
        a = _fig(tmp_path / "a", "order", n).read_bytes()
        b = _fig(tmp_path / "b", "order", n).read_bytes()
        assert a == b, n
    ta = (tmp_path / "a" / "tables" / "order" / "T_tf_groups.tsv").read_bytes()
    assert ta == (tmp_path / "b" / "tables" / "order" / "T_tf_groups.tsv").read_bytes()


def test_short_labels_stay_unique():
    from figs_tf import _short_unique
    a, b = "X" * 60 + " one", "X" * 60 + " two"
    sm = _short_unique([a, b, "WRKY"])
    assert len(set(sm.values())) == 3 and sm["WRKY"] == "WRKY"
    sm = _short_unique([LONG1, LONG2])
    assert len(set(sm.values())) == 2 and all(len(v) <= 44 for v in sm.values())


def test_family_z_values():
    from figs_tf import family_z
    d = pd.DataFrame({"species_prefix": ["a", "b", "c", "d", "e"], "grp": ["G1", "G1", "G2", "G2", "G3"],
                      "n_proteins_total": [100, 100, 100, 100, 0]})
    fam = pd.DataFrame({"species_prefix": ["a", "b", "c", "d", "a", "b", "c", "d"],
                        "family": ["X"] * 4 + ["Y"] * 4, "n_proteins": [10, 10, 30, 30, 5, 5, 5, 5]})
    z = family_z(d, fam, ["G1", "G2", "G3"])
    assert list(z.index) == ["X", "Y"]
    assert z.loc["X", "G1"] < 0 < z.loc["X", "G2"] and np.isnan(z.loc["X", "G3"])   # real contrast, empty group NaN
    assert (z.loc["Y", ["G1", "G2"]] == 0).all()                                      # same in every group: 0
    assert z.loc["X"].abs().max() > 0.5
