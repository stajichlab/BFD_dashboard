"""Fixture tests for 02_analyze.py (run on the HPCC: python3.12 -m pytest -q test_t4t6.py)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("t4t6", Path(__file__).with_name("02_analyze.py"))
t4t6 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t4t6)

FAM = [f"PF9{i:04d}" for i in range(60)]


def synth(seed=0, focal_ath_share=0.9, nth_drop=0.3):
    rng = np.random.default_rng(seed)
    reps, copies = [], {}

    def add(n, **tax):
        for _ in range(n):
            sp = f"G{len(reps):04d}"
            reps.append({"species_prefix": sp, "GENUS": f"Gen{len(reps)}", "CLASS": tax.get("CLASS", "Other"),
                         "ORDER_": tax.get("ORDER_", ""), "PHYLUM": tax.get("PHYLUM", "Ascomycota"),
                         "n_proteins": int(tax.get("n_proteins", 10000)), "primary_lifestyle": tax.get("life", "soil_saprotroph")})
            row = {f: rng.poisson(3.0 * tax.get("scale", 1.0)) + 1 for f in FAM}
            row[t4t6.NTH] = rng.poisson(3.0 * tax.get("nth_scale", 1.0)) + 1
            row[t4t6.TPS] = rng.poisson(3.0 + tax.get("tps_extra", 0))
            row[t4t6.ATH] = int(rng.random() < tax.get("ath_p", 0.5)) * (1 + rng.poisson(1))
            row[t4t6.S8] = rng.poisson(8.0 + tax.get("prot_extra", 0))
            row[t4t6.GH18] = rng.poisson(8.0 + tax.get("prot_extra", 0))
            copies[sp] = row
    add(20, CLASS="Sordariomycetes", life="animal_parasite", tps_extra=2.0, ath_p=focal_ath_share, prot_extra=4.0)
    add(60, CLASS="Sordariomycetes", life="plant_pathogen", ath_p=0.5)
    add(10, CLASS="Zygomycetes", ORDER_="Kickxellales", scale=0.8, nth_scale=nth_drop, n_proteins=6000)
    add(12, CLASS="Zygomycetes", ORDER_="Mucorales")
    add(10, CLASS="Microsporidea", PHYLUM="Microsporidia", scale=0.8, nth_scale=nth_drop, n_proteins=3000)
    add(10, CLASS="Chytridiomycetes", PHYLUM="Chytridiomycota")
    add(10, CLASS="Saccharomycetes", scale=0.8, nth_scale=nth_drop, n_proteins=5500)
    add(30, CLASS="Eurotiomycetes")      # Pezizomycotina reference, also control K1 side
    add(8, CLASS="Malasseziomycetes", scale=0.8, nth_scale=nth_drop, n_proteins=4500)
    add(10, CLASS="Ustilaginomycetes")
    add(10, CLASS="Dothideomycetes")
    add(10, CLASS="Leotiomycetes")
    reps = pd.DataFrame(reps)
    wide = pd.DataFrame.from_dict(copies, orient="index").fillna(0).astype(int)
    prev = {f: 0.99 for f in FAM}
    prev.update({t4t6.NTH: 0.99, t4t6.TPS: 0.95, t4t6.ATH: 0.5, t4t6.S8: 0.97, t4t6.GH18: 0.97})
    names = pd.DataFrame({"acc": list(prev), "pfam_id": ["Ribosomal_x" if i < 5 else "Fam" for i in range(len(prev))],
                          "prevalence": list(prev.values())})
    reps["subphylum"] = reps.CLASS.map(t4t6.CLASS_TO_SUBPHYLUM)
    return reps, wide, names


def test_cliffs_delta_and_bootstrap():
    assert t4t6.cliffs_delta([3, 4, 5], [0, 1, 2]) == 1.0
    assert t4t6.cliffs_delta([0, 1, 2], [3, 4, 5]) == -1.0
    lo, hi = t4t6.boot_median_diff([5, 6, 7, 8] * 5, [1, 2, 3, 4] * 5)
    assert lo >= 2.5 and hi <= 5.5     # true median difference 4; discrete data give interval ends such as 3.0 and 5.0


def test_t4_verdict_rules():
    good = {"diff": 1.0, "ci_lo": 0.5, "ci_hi": 1.5, "cliffs_delta": 0.6}
    assert t4t6.t4_verdict(good, good, 0.9) == ("supported", "supported", "supported")
    assert t4t6.t4_verdict(good, good, 0.3)[1:] == ("killed", "killed")          # ATH absent in most focal genera
    weak = {"diff": 0.1, "ci_lo": -0.2, "ci_hi": 0.4, "cliffs_delta": 0.1}
    assert t4t6.t4_verdict(weak, good, 0.9)[0] == "mixed"
    assert t4t6.t4_verdict({"diff": -1, "ci_lo": -2, "ci_hi": 0, "cliffs_delta": -0.5}, good, 0.9)[0] == "killed"


def test_t6_verdict_rules():
    assert t4t6.t6_verdict([0.05, 0.08, 0.02, 0.5], [0.5, 0.6]) == "supported"
    assert t4t6.t6_verdict([0.05, 0.08, 0.02, 0.5], [0.05, 0.04]) == "killed"        # same pattern without reduction
    assert t4t6.t6_verdict([0.6, 0.7, 0.4, 0.9], [0.5, 0.6]) == "killed"             # not dropped
    assert t4t6.t6_verdict([0.05, 0.2, 0.15, 0.5], [0.5, 0.6]) == "mixed"
    assert t4t6.t6_verdict([0.05, None, None, 0.1], [0.5, 0.6]).startswith("not evaluable")


def test_retention_percentile_picks_the_dropped_family():
    wide = pd.DataFrame({"a": [1, 1, 4, 4], "b": [1, 1, 4, 4], "c": [1, 1, 4, 4], "x": [0, 0, 4, 4]}, index=list("pqrs"))
    pct, ratio = t4t6.retention_percentile(wide, ["p", "q"], ["r", "s"], ["a", "b", "c", "x"], "x")
    assert pct == 0.0 and ratio < 0.05            # x is lost most, below every other family


def test_full_run_detects_built_in_effects(tmp_path):
    reps, wide, names = synth()
    t4, t6, lines = t4t6.run(reps, wide, names, tmp_path)
    text = "\n".join(lines)
    assert "T4 overall: supported" in text and "positive control" in text and "passes" in text
    assert "T6 overall: supported" in text        # NTH dropped in all four reduced groups, not in the controls
    killed, _, _ = synth(focal_ath_share=0.2)
    t4k, _, lk = t4t6.run(*synth(focal_ath_share=0.2), tmp_path / "k")
    assert "T4 overall: killed" in "\n".join(lk)
    _, _, ln = t4t6.run(*synth(nth_drop=1.0), tmp_path / "n")
    assert "T6 overall: supported" not in "\n".join(ln)
