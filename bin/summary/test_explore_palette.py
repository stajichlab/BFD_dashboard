"""Tests for the pure JS functions in 12_explore.qmd (distinctPalette, parseExploreParams), run with node."""
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

QMD = Path(__file__).resolve().parents[2] / "docs" / "quarto_summary" / "12_explore.qmd"
NODE = shutil.which("node")
pytestmark = pytest.mark.skipif(NODE is None, reason="node is not on PATH")


def _src(name):
    t = QMD.read_text()
    m = re.search(rf"// BEGIN {name}\n(.*?)// END {name}", t, re.S)
    assert m, f"marker {name} not found in 12_explore.qmd"
    return m.group(1)


def _palette():
    t = QMD.read_text()
    m = re.search(r"distinctColours = (\[.*?\])\n```", t, re.S)
    assert m
    return json.loads(re.sub(r",\s*]", "]", m.group(1)))


def _node(fn, call):
    js = f"{_src(fn)}\n{call}"
    r = subprocess.run([NODE, "-e", f"const out = (function(){{ {js} }})(); console.log(JSON.stringify(out));"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.loads(r.stdout)


def _dp(counts, palette=None):
    pal = palette or _palette()
    call = (f"const c = new Map(Object.entries({json.dumps(counts)}));"
            f"return Array.from(distinctPalette(Object.keys({json.dumps(counts)}).reverse(), c, {json.dumps(pal)}));")
    return _node("distinctPalette", call)


def test_palette_has_24_distinct_colours():
    pal = _palette()
    assert len(pal) == 24 and len(set(pal)) == 24
    assert "#9a9a9a" not in pal


def test_ten_groups_in_count_order():
    counts = {f"g{i}": 100 - i for i in range(10)}
    res = _dp(counts)
    assert [g for g, _ in res] == [f"g{i}" for i in range(10)]
    assert len({c for _, c in res}) == 10
    assert [c for _, c in res] == _palette()[:10]


def test_thirty_groups_others_grey():
    counts = {f"g{i:02d}": 1000 - i for i in range(30)}
    res = _dp(counts)
    cols = [c for _, c in res]
    assert len(set(cols[:24])) == 24 and cols[:24] == _palette()
    assert cols[24:] == ["#9a9a9a"] * 6
    assert sum(1 for c in cols if c == "#9a9a9a") == 6     # legend text: "Other (6 groups)"


def test_ties_broken_by_name_and_deterministic():
    counts = {"b": 5, "a": 5, "c": 5, "d": 9}
    res = _dp(counts)
    assert [g for g, _ in res] == ["d", "a", "b", "c"]
    assert res == _dp({"c": 5, "d": 9, "a": 5, "b": 5})


def test_top_colours_stable_when_small_group_removed():
    counts = {f"g{i:02d}": 500 - i * 3 for i in range(26)}
    before = dict(_dp(counts))
    smaller = dict(counts)
    del smaller["g25"]
    after = dict(_dp(smaller))
    assert all(before[f"g{i:02d}"] == after[f"g{i:02d}"] for i in range(24))


CTX = '{numFields: ["genome_mb", "busco_complete"], levels: ["phylum", "class", "order"], phyla: ["All", "Ascomycota"]}'


def _parse(q):
    return _node("parseExploreParams", f"return parseExploreParams({json.dumps(q)}, {CTX});")


def test_parse_valid():
    r = _parse("?x=genome_mb&y=busco_complete&level=order&phylum=Ascomycota&reps=1&logx=0&logy=1&palette=site")
    assert r == {"x": "genome_mb", "y": "busco_complete", "level": "order", "phylum": "Ascomycota",
                 "reps": True, "logx": False, "logy": True, "palette": "site"}


def test_parse_invalid_ignored():
    assert _parse("?x=nope&level=family&phylum=Foo&reps=yes&logx=2&palette=rainbow&zzz=1") == {}
    assert _parse("") == {}


def test_page_reads_parameter_names():
    t = QMD.read_text()
    for k in ("x", "y", "level", "phylum", "reps", "logx", "logy", "palette"):
        assert f'q.get("{k}")' in t or f'["{k}"' in t, k
    assert "window.location.search" in t and "Distinct colours" in t and "Site hues" in t and '"Auto"' in t
