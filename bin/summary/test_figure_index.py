"""Tests for figure_index.build_index and the sidecars written by Out.save."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_figures as mf  # noqa: E402
from figure_index import build_index, COLUMNS  # noqa: E402


def _png(root, level, name, side=None):
    d = Path(root) / "figures" / level
    d.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots()
    ax.plot([1, 2])
    fig.savefig(d / f"{name}.png")
    plt.close(fig)
    if side is not None:
        (d / f"{name}.json").write_text(json.dumps(side))


def _anova(root, level, rows):
    d = Path(root) / "tables" / level
    d.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["stat", "effect_size"]).to_csv(d / "T_anova.tsv", sep="\t", index=False)


def _side(**kw):
    base = {"topic": "structure", "stat": None, "n_genomes": 100, "n_groups": 8, "pooled_share": 0.0, "scope": "representatives", "kind": "figure"}
    base.update(kw)
    return base


@pytest.fixture
def root(tmp_path):
    _png(tmp_path, "phylum", "S01_gene_len", _side(stat="median_gene_len_bp", pooled_share=0.2))
    _png(tmp_path, "phylum", "S02_gene_density", _side(stat="genes_per_mb", pooled_share=0.0))
    _png(tmp_path, "phylum", "S08_heatmap", _side())                       # no stat: no score
    _png(tmp_path, "phylum", "F01_old")                                    # no sidecar
    _png(tmp_path, "class", "S01_gene_len", _side(stat="median_gene_len_bp", pooled_share=0.5))   # no T_anova for class
    _png(tmp_path, "all", "A01_x", _side(topic="scatter", scope="all_genomes", n_genomes=300, n_groups=None))
    _png(tmp_path, "ordination", "O01_codon", _side(topic="ordination", scope="genus_reps", n_genomes=50, n_groups=None))
    _anova(tmp_path, "phylum", [("median_gene_len_bp", 0.5), ("genes_per_mb", 0.25)])
    return tmp_path


def test_one_row_per_png_and_file_written(root):
    df = build_index(root)
    n = len(list((root / "figures").glob("*/*.png")))
    assert len(df) == n == 7
    assert list(df.columns) == COLUMNS
    assert not df.duplicated(["level", "name"]).any()
    t = pd.read_csv(root / "figures_index.tsv", sep="\t")
    assert len(t) == n


def test_pooled_share_in_unit_interval(root):
    df = build_index(root)
    ps = df["pooled_share"].dropna()
    assert len(ps) > 0 and ((ps >= 0) & (ps <= 1)).all()


def test_scores_exact(root):
    df = build_index(root).set_index(["level", "name"])
    assert df.loc[("phylum", "S01_gene_len"), "score"] == pytest.approx(0.5 * (1 - 0.2))
    assert df.loc[("phylum", "S02_gene_density"), "score"] == pytest.approx(0.25)
    assert df.loc[("phylum", "S01_gene_len"), "effect_size"] == pytest.approx(0.5)


def test_no_anova_row_means_empty_score_and_sorted_last(root):
    df = build_index(root)
    s = df.set_index(["level", "name"])
    assert pd.isna(s.loc[("phylum", "S08_heatmap"), "score"])
    assert pd.isna(s.loc[("class", "S01_gene_len"), "score"])   # missing T_anova.tsv
    scored = df["score"].notna().to_numpy()
    assert scored.tolist() == sorted(scored.tolist(), reverse=True)
    assert df["score"].dropna().is_monotonic_decreasing
    t = pd.read_csv(root / "figures_index.tsv", sep="\t")
    assert t["score"].isna().sum() == 5


def test_all_and_ordination_levels_are_indexed(root):
    s = build_index(root).set_index(["level", "name"])
    assert s.loc[("all", "A01_x"), "scope"] == "all_genomes"
    assert s.loc[("ordination", "O01_codon"), "scope"] == "genus_reps"


def test_png_without_sidecar_gets_row(root):
    s = build_index(root).set_index(["level", "name"])
    r = s.loc[("phylum", "F01_old")]
    assert pd.isna(r["topic"]) and pd.isna(r["score"]) and pd.isna(r["n_genomes"])


def test_non_finite_effect_size_has_empty_score(tmp_path):
    _png(tmp_path, "phylum", "S01", _side(stat="a"))
    _anova(tmp_path, "phylum", [("a", np.nan)])
    df = build_index(tmp_path)
    assert pd.isna(df["score"].iloc[0]) and pd.isna(df["effect_size"].iloc[0])


def test_corrupt_sidecar_does_not_crash(tmp_path):
    _png(tmp_path, "phylum", "S01")
    (tmp_path / "figures" / "phylum" / "S01.json").write_text("{not json")
    assert len(build_index(tmp_path)) == 1


def test_empty_root_gives_empty_index(tmp_path):
    df = build_index(tmp_path)
    assert len(df) == 0 and list(df.columns) == COLUMNS


def test_out_save_writes_sidecar(tmp_path):
    core = pd.DataFrame({"grp": ["A"] * 6 + [mf.OTHER] * 4})
    o = mf.Out(tmp_path, "phylum")
    o.set_groups(core, ["A", mf.OTHER])
    o.stat_of["S01_x"] = "genes_per_mb"
    fig, ax = plt.subplots()
    ax.plot([1])
    o.save(fig, "S01_x", topic="structure")
    side = json.loads((tmp_path / "figures" / "phylum" / "S01_x.json").read_text())
    assert side["topic"] == "structure" and side["stat"] == "genes_per_mb" and side["scope"] == "representatives"
    assert side["n_genomes"] == 10 and side["n_groups"] == 2 and side["pooled_share"] == pytest.approx(0.4)


def test_out_save_all_level_uses_scatter_n(tmp_path):
    o = mf.Out(tmp_path, "all")
    o.meta["A01"] = {"n": 77}
    fig, ax = plt.subplots()
    ax.plot([1])
    o.save(fig, "A01", topic="scatter")
    side = json.loads((tmp_path / "figures" / "all" / "A01.json").read_text())
    assert side["n_genomes"] == 77 and side["pooled_share"] == 0 and side["scope"] == "all_genomes"


# ---------------------------------------------------------------- signal columns
def test_signal_columns_and_sort(tmp_path):
    _png(tmp_path, "phylum", "S01", _side(stat="a", pooled_share=0.0))
    _png(tmp_path, "all", "A01", _side(topic="scatter", scope="all_genomes", signal=0.9, signal_kind="spearman_r"))
    _png(tmp_path, "all", "A02", _side(topic="scatter", scope="all_genomes", signal=0.3, signal_kind="spearman_r"))
    _png(tmp_path, "ordination", "O01", _side(topic="ordination", scope="genus_reps", signal=0.5, signal_kind="permanova_R2"))
    _png(tmp_path, "phylum", "F01")
    _anova(tmp_path, "phylum", [("a", 0.1)])
    df = build_index(tmp_path)
    assert list(df.columns) == COLUMNS and COLUMNS[-2:] == ["signal", "signal_kind"] and COLUMNS[:10][-1] == "score"
    assert df["name"].tolist() == ["S01", "A01", "O01", "A02", "F01"]   # scored first, then signal desc, then unsignalled
    s = df.set_index("name")
    assert s.loc["S01", "signal"] == pytest.approx(0.1) and s.loc["S01", "signal_kind"] == "epsilon2"
    assert s.loc["O01", "signal_kind"] == "permanova_R2"
    assert pd.isna(s.loc["F01", "signal"]) and pd.isna(s.loc["F01", "signal_kind"])
    assert pd.isna(s.loc["A01", "score"])


def test_out_save_scatter_signal(tmp_path):
    o = mf.Out(tmp_path, "all")
    o.meta["A01"] = {"n": 5, "spearman_r": -0.8}
    fig, ax = plt.subplots()
    ax.plot([1])
    o.save(fig, "A01", topic="scatter")
    side = json.loads((tmp_path / "figures" / "all" / "A01.json").read_text())
    assert side["signal"] == pytest.approx(0.8) and side["signal_kind"] == "spearman_r"


def test_sidecar_write_failure_is_isolated(tmp_path, capsys):
    o = mf.Out(tmp_path, "phylum")
    (o.fig / "S01.json").mkdir()          # the sidecar path is a directory: write raises OSError
    fig, ax = plt.subplots()
    ax.plot([1])
    o.save(fig, "S01", topic="structure")
    assert [(str(m), m.kind) for m in o.made] == [("S01", "figure")]
    assert "sidecar" in capsys.readouterr().err
    df = build_index(tmp_path)
    assert len(df) == 1 and pd.isna(df["topic"].iloc[0])


def test_out_save_coerces_numpy(tmp_path):
    o = mf.Out(tmp_path, "phylum")
    o.stat_of["S01"] = "x"
    o.info = {"n_genomes": np.int64(5), "n_groups": np.int64(2), "pooled_share": np.float64(0.2)}
    fig, ax = plt.subplots()
    ax.plot([1])
    o.save(fig, "S01", topic="t")
    assert json.loads((o.fig / "S01.json").read_text())["n_genomes"] == 5


def test_old_pages_get_a_topic(tmp_path):
    o = mf.Out(tmp_path, "phylum")
    o.default_topic = "gene_content"
    fig, ax = plt.subplots()
    ax.plot([1])
    o.save(fig, "F01")
    assert json.loads((o.fig / "F01.json").read_text())["topic"] == "gene_content"


# ---------------------------------------------------------------- end to end
def test_end_to_end_index(tmp_path):
    from test_make_figures import synth
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    ga = pd.read_parquet(data / "genome_all.parquet")
    ga["GENUS"] = [f"G{i % 60}" for i in range(len(ga))]
    ga.to_parquet(data / "genome_all.parquet")
    out = tmp_path / "out"
    mf.main(data, out, ("phylum",))
    side = [json.loads(p.read_text()) for p in (out / "figures").glob("*/*.json")]
    anova = pd.read_csv(out / "tables" / "phylum" / "T_anova.tsv", sep="\t")
    assert {s["stat"] for s in side if s.get("stat")} & set(anova["stat"])          # (a)
    idx = pd.read_csv(out / "figures_index.tsv", sep="\t")
    sc = idx[idx["score"].notna()]
    assert sc["name"].str.startswith("S").any() and sc["name"].str.startswith("F").any()   # (b)
    for pre in ("A", "R"):                                                          # (c)
        r = idx[idx["name"].str.match(rf"{pre}\d") & idx["signal"].notna()]
        assert len(r) > 0 and (r["signal_kind"] == "spearman_r").all()
    assert (idx[idx["name"].str.startswith("A")]["scope"] == "all_genomes").all()
    o = idx[idx["level"] == "ordination"]                                           # (d)
    assert o["name"].str.startswith("O01").any() and o["name"].str.startswith("O02").any()
    assert (o["scope"] == "genus_reps").all() and o["n_genomes"].notna().all() and o["n_groups"].notna().all()
    assert o["signal"].notna().all() and (o["signal_kind"] == "permanova_R2").all()
    assert o["signal"].between(0, 1).all()
    assert idx[idx["level"] == "phylum"]["topic"].notna().all()
