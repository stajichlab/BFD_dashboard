"""Behaviour tests for figs_common / figs_structure / figs_composition (known answers and failure isolation)."""
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
import figs_common  # noqa: E402
import figs_composition as fc  # noqa: E402
import figs_structure as fs  # noqa: E402

ALA = ["GCT", "GCC", "GCA", "GCG"]


def test_rscu_known_answer():
    f = pd.DataFrame({"codon_GCT": [.4], "codon_GCC": [.2], "codon_GCA": [.2], "codon_GCG": [.2],
                      "codon_ATG": [.05], "codon_TGG": [.03]})
    r = fc._rscu(f)
    assert r[[f"codon_{c}" for c in ALA]].iloc[0].tolist() == pytest.approx([1.6, .8, .8, .8])
    assert r["codon_ATG"].iloc[0] == pytest.approx(1.0) and r["codon_TGG"].iloc[0] == pytest.approx(1.0)


def test_gc3_known_answer():
    # sense codons: AAA (A3), AAG (G3), GCC (C3), GCT (T3); stop TAA must be ignored
    codon = pd.DataFrame({"species_prefix": ["x", "y"], "codon_AAA": [.1, 0.], "codon_AAG": [.2, 0.], "codon_GCC": [.3, 1.0],
                          "codon_GCT": [.4, 0.], "codon_TAA": [.5, .5]})
    g = fc._gc3(codon)
    assert g.tolist() == pytest.approx([0.5, 1.0])   # (0.2+0.3)/1.0 ; 1.0/1.0


def test_codon_cols_and_group_mean():
    codon = pd.DataFrame({"species_prefix": ["a", "b", "c"], "codon_TTT": [.1, .3, .9], "codon_TAA": [0, 0, 0], "codon_AAA": [.9, .7, .1]})
    assert fc._codon_cols(codon) == ["codon_TTT", "codon_AAA"] or fc._codon_cols(codon) == ["codon_AAA", "codon_TTT"]
    assert "codon_TAA" not in fc._codon_cols(codon)
    core = pd.DataFrame({"species_prefix": ["a", "b", "c"], "grp": ["G1", "G1", "G2"], "gc_pct": 50.0})
    d = figs_common.join_groups(core, codon)
    m = d.groupby("grp")[["codon_TTT"]].mean()
    assert m.loc["G1", "codon_TTT"] == pytest.approx(0.2) and m.loc["G2", "codon_TTT"] == pytest.approx(0.9)


class _FakeOut(mf.Out):
    pass


def _fig():
    fig, ax = plt.subplots()
    ax.plot([1, 2])
    return fig


def test_safe_failure_is_skipped_and_later_figure_saves(tmp_path, capsys):
    o = mf.Out(tmp_path, "phylum")

    def boom():
        raise ValueError("bad")
    assert figs_common.safe(o, "X01_bad", "t", boom) is False
    assert figs_common.safe(o, "X02_ok", "t", _fig) is True
    kinds = {str(m): m.kind for m in o.made}
    assert kinds == {"X01_bad": "skipped", "X02_ok": "figure"}
    assert (tmp_path / "figures" / "phylum" / "X02_ok.png").exists()


def test_safe_save_failure_is_skipped(tmp_path):
    o = mf.Out(tmp_path, "phylum")
    real = o.save
    calls = []

    def flaky(fig, name, topic=None):
        calls.append(name)
        if name == "X01_nosave":
            raise OSError("disk full")
        return real(fig, name, topic=topic)
    o.save = flaky
    assert figs_common.safe(o, "X01_nosave", "t", _fig) is False
    assert figs_common.safe(o, "X02_ok", "t", _fig) is True
    assert {str(m): m.kind for m in o.made} == {"X01_nosave": "skipped", "X02_ok": "figure"}


def _edge_core():
    n = 14
    ph = ["Big"] * 10 + ["Solo"] * 1 + ["Big"] * 3
    core = pd.DataFrame({"species_prefix": [f"E{i}" for i in range(n)], "PHYLUM": ph,
                         "CLASS": ["c"] * n, "ORDER": [None] * 5 + ["o1"] * 9,
                         "gc_pct": np.linspace(35, 60, n)})
    return core


def _edge_tables(core):
    rng = np.random.default_rng(3)
    n = len(core)
    ids = core["species_prefix"]
    s = pd.DataFrame({"species_prefix": ids, "n_genes": 8000, "median_gene_len_bp": rng.uniform(1200, 2000, n),
                      "median_cds_len_bp": 1500.0, "median_protein_len_aa": 500.0,   # constant columns
                      "genes_per_mb": rng.uniform(200, 400, n), "coding_fraction": rng.uniform(.4, .6, n),
                      "median_exons_per_gene": 2.0, "intron_frac_genes_with_intron": rng.uniform(.2, .8, n),
                      "median_intron_len_bp": np.nan, "median_intergenic_bp": rng.uniform(300, 2000, n)})
    s.loc[:4, "median_intron_len_bp"] = [60.0, 80, 90, 70, 100]       # the rest, incl. the Solo genome, are NaN
    cods = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG" if a + b + c not in ("TAA", "TAG", "TGA")]
    f = rng.dirichlet(np.ones(61) * 8, n)
    codon = pd.DataFrame(f, columns=[f"codon_{c}" for c in cods]).assign(species_prefix=ids.values)
    aa = pd.DataFrame(rng.dirichlet(np.ones(20) * 20, n), columns=[f"aa_{a}" for a in "ACDEFGHIKLMNPQRSTVWY"]).assign(species_prefix=ids.values)
    return s, codon, aa


def test_edge_cases_pages_do_not_raise(tmp_path):
    core0 = _edge_core()
    s, codon, aa = _edge_tables(core0)
    # phylum level: Solo has one genome; order level: ORDER None -> Unknown group; a tiny phylum may be pooled to Other
    for level in ("phylum", "order"):
        core, levels, colors = mf.assign_groups(core0, level)
        o = mf.Out(tmp_path, level)
        fs.page_gene_structure(core, levels, colors, s, o)
        fc.page_composition(core, levels, colors, codon, aa, o)
        figdir = tmp_path / "figures" / level
        have = {p.stem for p in figdir.glob("*.png")}
        required = {"S01_gene_len", "S02_gene_density", "S03_coding_fraction", "S04_exons", "S06_intron_frac", "S07_intergenic",
                    "C01_codon_heatmap", "C02_rscu_boxes", "C03_aa_composition", "C04_gc3", "C05_aa_vs_gc"}
        assert required <= have, (level, required - have)
        assert "S05_intron_len" in have           # some genomes have introns, so it is drawn
        skipped = {str(m) for m in o.made if m.kind == "skipped"}
        assert not skipped, (level, skipped)
        if level == "phylum":
            assert "Solo" in levels or mf.OTHER in levels


def test_all_nan_column_is_skipped_not_raised(tmp_path):
    core0 = _edge_core()
    s, codon, aa = _edge_tables(core0)
    s["median_intron_len_bp"] = np.nan
    core, levels, colors = mf.assign_groups(core0, "phylum")
    o = mf.Out(tmp_path, "phylum")
    fs.page_gene_structure(core, levels, colors, s, o)
    assert {str(m): m.kind for m in o.made}["S05_intron_len"] == "skipped"
    assert (tmp_path / "figures" / "phylum" / "S06_intron_frac.png").exists()
