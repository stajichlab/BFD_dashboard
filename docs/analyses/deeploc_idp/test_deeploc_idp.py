"""Fixture tests for 02_analyze.py (run on the HPCC: python3.12 -m pytest -q test_deeploc_idp.py).

A synthetic cube with known effects checks the direction and verdict logic: nuclear TFs are made more disordered than other nuclear
proteins, secreted proteins less, unannotated more, and Microsporidia less than the rest; GC is tied to IDR by construction.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("ana", Path(__file__).with_name("02_analyze.py"))
ana = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ana)

LBINS = ["<100", "100-299", "300-999", ">=1000"]
LOCS = ["cytoplasm", "nucleus", "extracellular", "mitochondrion", "cell_membrane"]


def synth(n_per_class=12, seed=0):
    rng = np.random.default_rng(seed)
    classes = {f"C{i}": ("Ascomycota" if i < 8 else "Microsporidia") for i in range(9)}
    reps, rows, ssp = [], [], []
    k = 0
    for cl, ph in classes.items():
        for _ in range(n_per_class if ph != "Microsporidia" else 10):
            sp = f"G{k:04d}"
            k += 1
            gc = rng.uniform(35, 60)
            base = 0.2 + 0.004 * (gc - 35) + rng.normal(0, 0.01) - (0.1 if ph == "Microsporidia" else 0)
            reps.append({"species_prefix": sp, "PHYLUM": ph, "CLASS": cl, "GENUS": f"Gen{k}", "gc_pct": gc})
            for lb in LBINS:
                for loc in LOCS:
                    for sec in (False, True):
                        for tf in (False, True):
                            if tf and loc != "nucleus":
                                continue
                            for unann in (False, True):
                                n = int(rng.integers(20, 60))
                                eff = base + (0.15 if tf else 0) - (0.1 if sec else 0) + (0.1 if unann else 0)
                                eff += 0.05 if loc == "nucleus" else (-0.05 if loc == "mitochondrion" else 0)
                                rows.append({"species_prefix": sp, "lbin": lb, "loc_": loc, "sec": sec, "tf": tf,
                                             "has_pfam": not unann, "unann": unann,
                                             "mito_carrier": loc == "mitochondrion" and not sec and not tf,
                                             "effector_fam": sec and not unann and rng.random() < 0.5,
                                             "n": n, "n_row": n, "sum_idr": n * max(eff, 0.0)})
            ssp.append({"species_prefix": sp, "n_ssp": int(rng.integers(5, 30)), "n_proteins": 8000})
    return pd.DataFrame(reps), pd.DataFrame(rows), pd.DataFrame(ssp)


def test_bh_matches_known_values():
    q = ana.bh([0.01, 0.04, 0.03, 0.005])
    assert np.allclose(q, [0.02, 0.04, 0.04, 0.02])


def test_matched_contrast_direction_and_weights():
    reps, cube, _ = synth()
    c = ana.matched_contrast(cube, cube.tf & (cube.loc_ == "nucleus"), ~cube.tf & (cube.loc_ == "nucleus"))
    assert (c > 0.1).all()                                    # built-in effect +0.15
    c2 = ana.matched_contrast(cube, cube.sec, ~cube.sec)
    assert (c2 < -0.05).all()


def test_sign_test_and_verdicts():
    frac, n, p = ana.sign_test([0.1, 0.2, 0.3, -0.1], "+")
    assert (frac, n) == (0.75, 4) and 0.3 < p < 0.4
    assert ana.verdict_sign(0.95, 0.90, 0.60) == "supported"
    assert ana.verdict_sign(0.40, 0.90, 0.60) == "killed"
    assert ana.verdict_sign(0.75, 0.90, 0.60) == "mixed"
    assert ana.verdict_sign(0.95, 0.90, 0.60, pooled_ok=False) == "mixed"   # effect size below the pooled threshold


def test_full_run(tmp_path):
    reps, cube, ssp = synth()
    ft = pd.DataFrame({"GENUS": reps.GENUS, "primary_lifestyle": ["plant_pathogen" if i % 2 else "soil_saprotroph"
                                                                    for i in range(len(reps))]})
    ctl, res = ana.run(cube, reps, ssp, ft, tmp_path)
    assert ctl.set_index("control").loc["C5", "passed"]       # nucleus above mitochondrion by construction
    r = res.set_index("hypothesis")
    assert r.loc["H1", "verdict"] == "supported" and r.loc["H2", "verdict"] == "supported"
    assert r.loc["H3", "verdict"] == "supported" and r.loc["H6", "verdict"] == "supported"
    assert r.loc["H4", "value"] > 0.3 and r.loc["H4", "verdict"] == "supported"   # IDR tied to GC by construction
    assert (tmp_path / "verdicts.txt").exists() and (tmp_path / "hypotheses.tsv").exists()
