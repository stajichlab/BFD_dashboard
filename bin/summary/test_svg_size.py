"""SVG size: point clouds are rasterized, so an SVG with thousands of genome points stays small."""
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_figures as mf  # noqa: E402
import figs_composition as fc  # noqa: E402
import figs_structure as fs  # noqa: E402
import summary_ordination as so  # noqa: E402
from test_make_figures import synth  # noqa: E402

N = 3000
MAX_BYTES = 400_000     # a vector scatter of 3,000 points per panel is several MB; text and axes are < 100 KB
MAX_POINT_ELEMENTS = 200
# C05 has 20 panels, each with its own embedded raster (about 70 KB) and its own ticks: measured 0.8 MB at 3,000
# points once rasterized, against about 6.5 MB as vectors.
C05_MAX_BYTES = 1_500_000
C05_MAX_ELEMENTS = 500


def _frames(n=N, seed=3):
    rng = np.random.default_rng(seed)
    levels = [f"G{i}" for i in range(6)]
    cmap = {g: plt.get_cmap("tab10")(i) for i, g in enumerate(levels)}
    colors = {g: matplotlib.colors.to_hex(c) for g, c in cmap.items()}
    ids = [f"S{i}" for i in range(n)]
    core = pd.DataFrame({"species_prefix": ids, "grp": rng.choice(levels, n), "gc_pct": rng.uniform(30, 60, n),
                         "PHYLUM": "P"})
    core.attrs["parent_colors"] = {}
    return rng, levels, colors, ids, core


def _svg_check(fig, tmp_path, name):
    p = tmp_path / f"{name}.svg"
    fig.savefig(p, format="svg")
    plt.close(fig)
    text = p.read_text()
    n_el = len(re.findall(r"<(?:path|use|circle)\b", text))
    return p.stat().st_size, n_el


def _assert_small(size, n_el, name, max_bytes=MAX_BYTES, max_elements=MAX_POINT_ELEMENTS):
    assert size < max_bytes, (name, size)
    assert n_el < max_elements, (name, n_el)


def test_c04b_gc3_vs_gc_svg_small(tmp_path):
    rng, levels, colors, ids, core = _frames()
    g3 = pd.Series(rng.uniform(0.2, 0.8, len(ids)), index=pd.Index(ids, name="species_prefix"), name="gc3")
    _assert_small(*_svg_check(fc._gc3_vs_gc(g3, core, levels, colors), tmp_path, "c04b"), "C04b")


def test_c05_aa_vs_gc_svg_small(tmp_path):
    rng, levels, colors, ids, core = _frames()
    g3 = pd.Series(rng.uniform(0.2, 0.8, len(ids)), index=pd.Index(ids, name="species_prefix"), name="gc3")
    aa = pd.DataFrame(rng.uniform(0.01, 0.1, (len(ids), 20)), columns=[f"aa_{a}" for a in fc.AA_ORDER])
    aa.insert(0, "species_prefix", ids)
    _assert_small(*_svg_check(fc._aa_vs_gc(aa, g3, core, levels, colors), tmp_path, "c05"), "C05", C05_MAX_BYTES, C05_MAX_ELEMENTS)


def test_s09_size_fig_svg_small(tmp_path):
    rng, levels, colors, ids, core = _frames()
    d = core.copy()
    d["median_gene_len_bp"] = rng.uniform(1000, 3000, len(d))
    d["median_intergenic_bp"] = rng.uniform(100, 5000, len(d))
    _assert_small(*_svg_check(fs._size_fig(d, core, levels, colors), tmp_path, "s09"), "S09")


def test_ordination_scatter_svg_small(tmp_path):
    rng, levels, colors, ids, core = _frames()
    coords = pd.DataFrame({"PCo1": rng.normal(size=len(ids)), "PCo2": rng.normal(size=len(ids))}, index=ids)
    groups = core.set_index("species_prefix")["grp"]
    _assert_small(*_svg_check(so._scatter(coords, groups, colors, [10.0, 5.0], core, "t"), tmp_path, "pcoa"), "PCoA")


def test_box_points_are_rasterized(tmp_path):
    rng, levels, colors, ids, core = _frames()
    fig, ax = plt.subplots()
    mf.box_by_group(ax, core.assign(v=rng.normal(size=len(core))), "v", levels, colors)
    _assert_small(*_svg_check(fig, tmp_path, "box"), "box")


def test_every_scatter_in_a_full_run_is_rasterized(tmp_path, monkeypatch):
    calls = []
    orig = matplotlib.axes.Axes.scatter

    def spy(self, *a, **k):
        calls.append(k.get("rasterized"))
        return orig(self, *a, **k)
    monkeypatch.setattr(matplotlib.axes.Axes, "scatter", spy)
    data = tmp_path / "data"
    data.mkdir()
    synth(data)
    mf.main(data, tmp_path / "out", ("phylum",))
    assert calls and all(c is True for c in calls), [i for i, c in enumerate(calls) if c is not True][:10]
