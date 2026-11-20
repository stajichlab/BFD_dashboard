"""Fixture tests for t3_lib: the stratified permutation must find a planted within-genus shift, must not find one when the
shift is only between genera, and must ignore genera that are not eligible. BH comes from the already-tested stats_util."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import t3_lib as L
from stats_util import bh_qvalues


def _frame(shift_within: float, shift_between: float, seed: int = 1):
    rng = np.random.default_rng(seed)
    rows = []
    for gi in range(8):
        base = shift_between * (1 if gi % 2 == 0 else -1)
        for k in range(10):
            lab = k < 3
            rows.append(("G%d" % gi, lab, base + (shift_within if lab else 0.0) + rng.normal(0, 0.5)))
    return pd.DataFrame(rows, columns=["GENUS", "label", "M"])


def test_planted_within_genus_shift_is_found():
    r = L.stratified_permutation(_frame(2.0, 0.0), "M", n_perm=2000)
    assert r["n_genera"] == 8 and r["obs"] > 1.5 and r["p"] < 0.01


def test_between_genus_difference_is_not_a_within_genus_signal():
    r = L.stratified_permutation(_frame(0.0, 5.0), "M", n_perm=2000)
    assert abs(r["obs"]) < 0.6 and r["p"] > 0.05


def test_ineligible_genera_are_ignored():
    df = pd.DataFrame({"GENUS": ["A"] * 4 + ["B"] * 12, "label": [True, False, False, False] + [True] * 3 + [False] * 9,
                       "M": list(range(4)) + list(range(12))})
    t = L.genus_shift(df, "M")
    assert list(t.GENUS) == ["B"]                         # A has only one labelled species
    assert L.stratified_permutation(df[df.GENUS == "A"], "M", n_perm=50)["n_genera"] == 0


def test_dry_epithet_is_the_th1_pattern():
    names = pd.Series(["Aspergillus desertorum", "Aspergillus niger", "Cryomyces antarcticus", "Penicillium salinarum"])
    assert list(names.str.contains(L.DRY_EPITHET, regex=True)) == [True, False, True, True]


def test_bh_is_the_fixed_implementation():
    q = bh_qvalues([0.04, 0.001, 0.03])
    assert np.allclose(q, [0.04, 0.003, 0.04])
