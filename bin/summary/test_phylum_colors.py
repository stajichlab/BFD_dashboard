"""One phylum -> colour map for the whole site (ranked by genome count in genome_all), and the A scatters."""
import json
import re
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
from make_figures import PALETTE, GREY, phylum_order, phylum_color_map, assign_groups  # noqa: E402
from figs_scatter import scatter_all  # noqa: E402

COUNTS = {"Asco": 500, "Basidio": 300, "Muco": 50, "Kick": 40, "Micro": 30, "Mort": 20, "Chytrid": 15, "Glom": 12,
          "Rozello": 8, "Tiny": 3}


def _all(n_extra_unknown=5, seed=0):
    rng = np.random.default_rng(seed)
    ph = [p for p, c in COUNTS.items() for _ in range(c)] + [None] * n_extra_unknown
    n = len(ph)
    size = rng.integers(10_000_000, 60_000_000, n)
    return pd.DataFrame({"PHYLUM": ph, "CLASS": [f"{(p or 'U')[:3]}_c{i % 3}" for i, p in enumerate(ph)],
                         "ORDER": [f"{(p or 'U')[:3]}_o{i % 4}" for i, p in enumerate(ph)],
                         "total_length_bp": size, "n_genes": (size / 2500 + rng.normal(0, 500, n)).astype(int),
                         "species_prefix": [f"s{i}" for i in range(n)]})


def test_map_ranks_by_count_and_greys_beyond_palette():
    df = _all()
    order = phylum_order(df)
    assert order[:3] == ["Asco", "Basidio", "Muco"] and order[-1] == "Unknown"
    cmap = phylum_color_map(order)
    assert [cmap[p] for p in order[:8]] == PALETTE
    assert "Unknown" not in cmap and "Rozello" not in cmap and "Tiny" not in cmap


def test_assign_groups_uses_the_global_colours_at_every_level():
    df = _all()
    order = phylum_order(df)
    cmap = phylum_color_map(order)
    reps = df.iloc[::7].reset_index(drop=True)       # a subset with a different rank order
    for level in ("phylum", "subphylum", "class", "order"):
        _, _, _ = assign_groups(reps, level, phylum_order=order)
        d, _, colors = assign_groups(reps, level, phylum_order=order)
        for p, c in d.attrs["parent_colors"].items():
            if p in cmap:
                assert c == cmap[p], (level, p)
    d, levels, colors = assign_groups(reps, "phylum", phylum_order=order)
    for p in cmap:
        if p in colors:
            assert colors[p] == cmap[p]


def test_assign_groups_default_unchanged_without_order():
    df = _all()
    d, levels, colors = assign_groups(df, "phylum")
    assert colors["Asco"] == PALETTE[0] and colors["Basidio"] == PALETTE[1]


def test_scatter_all_colours_match_global_map_and_draw_every_genome(tmp_path, monkeypatch):
    df = _all()
    cmap = phylum_color_map(phylum_order(df))
    calls = []
    orig = matplotlib.axes.Axes.scatter

    def spy(self, x, y, *a, **k):
        calls.append({"n": len(x), "color": k.get("color"), "alpha": k.get("alpha"), "s": k.get("s"),
                      "rasterized": k.get("rasterized")})
        return orig(self, x, y, *a, **k)
    monkeypatch.setattr(matplotlib.axes.Axes, "scatter", spy)
    hexbin_calls = []
    monkeypatch.setattr(matplotlib.axes.Axes, "hexbin", lambda *a, **k: hexbin_calls.append(1))
    o = mf.Out(tmp_path, "all")
    assert scatter_all(df, "total_length_bp", "n_genes", cmap, o, "A07_x", logx=True, logy=True)
    n_ok = [m for m in o.made if m == "A07_x"][0].stats["n"]
    assert n_ok == len(df) and sum(c["n"] for c in calls) == n_ok and not hexbin_calls
    assert all(c["rasterized"] is True and c["alpha"] <= 0.35 and c["s"] <= 4 for c in calls)
    assert {c["color"] for c in calls} == set(cmap.values()) | {GREY}
    assert calls[0]["color"] == GREY                       # grey points first, coloured phyla stay visible
    sizes = [c["n"] for c in calls if c["color"] != GREY]
    assert sizes == sorted(sizes, reverse=True)             # small phyla last


def test_scatter_all_legend_lists_phyla_with_n_and_other(tmp_path, monkeypatch):
    df = _all()
    cmap = phylum_color_map(phylum_order(df))
    seen = {}
    orig = plt.Figure.savefig

    def grab(self, *a, **k):
        seen["labels"] = [t.get_text() for ax in self.axes for t in (ax.get_legend().get_texts() if ax.get_legend() else [])]
        return orig(self, *a, **k)
    monkeypatch.setattr(plt.Figure, "savefig", grab)
    o = mf.Out(tmp_path, "all")
    scatter_all(df, "total_length_bp", "n_genes", cmap, o, "A07_x", logx=True, logy=True)
    labs = seen["labels"]
    assert "Asco (n=500)" in labs and "Glom (n=12)" in labs and "Other phyla" in labs


def test_json_phylum_colours_match_global_map(tmp_path):
    from export_clade_json import export_clade_json
    from test_export_clade_json import _frames
    core, st = _frames(300)
    j = json.loads(export_clade_json(core, st, tmp_path / "g.json").read_text())
    cmap = phylum_color_map(phylum_order(core))
    for p, c in j["colors"]["phylum"].items():
        if p in cmap:
            assert c == cmap[p]
    assert j["colors"]["phylum"]["Ascomycota"] == cmap["Ascomycota"]


def test_a_scatter_svg_small_on_full_size_frame(tmp_path):
    n = 23680
    rng = np.random.default_rng(2)
    ph = rng.choice(list(COUNTS), n, p=np.array(list(COUNTS.values())) / sum(COUNTS.values()))
    df = pd.DataFrame({"PHYLUM": ph, "total_length_bp": rng.uniform(1e7, 9e7, n), "busco_complete": rng.uniform(50, 100, n)})
    o = mf.Out(tmp_path, "all")
    assert scatter_all(df, "total_length_bp", "busco_complete", phylum_color_map(phylum_order(df)), o, "A01_x", logx=True)
    svg = (tmp_path / "figures" / "all" / "A01_x.svg").stat().st_size
    png = (tmp_path / "figures" / "all" / "A01_x.png").stat().st_size
    print("A01 on 23,680 points: svg", svg, "png", png)
    assert svg < 600_000
    assert len(re.findall(r"<(?:path|use|circle)\b", (tmp_path / "figures" / "all" / "A01_x.svg").read_text())) < 300
