"""Fixture tests for 04_gc_adjusted.py (uses the synthetic cube of test_deeploc_idp.py)."""
import importlib.util
from pathlib import Path

import numpy as np

here = Path(__file__).parent


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, here / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


T = load("t", "test_deeploc_idp.py")
G = load("g", "04_gc_adjusted.py")


def test_contrast_independent_of_gc_keeps_its_size():
    """GC raises IDR in every category alike, so the within-genome contrast must not depend on GC."""
    reps, cube, _ = T.synth()
    nuc = cube.loc_ == "nucleus"
    c = G.A.matched_contrast(cube, cube.tf & nuc, ~cube.tf & nuc)
    d, slope, rho, n = G.adjust_for_gc(c, reps)
    assert abs(slope) < 0.002
    assert abs(d.adjusted.median() - d.contrast.median()) < 0.01


def test_gc_dependent_contrast_is_removed_by_adjustment():
    reps, cube, _ = T.synth()
    c = G.A.matched_contrast(cube, cube.sec, ~cube.sec)
    gc = reps.set_index("species_prefix").gc_pct
    c2 = c + 0.01 * (gc.reindex(c.index) - 47)          # inject a GC-driven part of the contrast
    d, slope, rho, n = G.adjust_for_gc(c2, reps)
    assert slope > 0.005
    assert abs(np.polyfit(d.gc_dm, d.adjusted, 1)[0]) < 1e-9   # nothing GC-linear remains
