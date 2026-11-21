#!/usr/bin/env python3
"""POST HOC GC-adjusted versions of H1-H3 (not pre-registered; do not replace the pre-registered verdicts).

Usage: 04_gc_adjusted.py CUBE_DIR OUT_DIR
Each H1-H3 contrast is already computed within genome and within length bin, so a genome-level GC difference cannot create it
directly. What this checks is whether the size of the contrast depends on genome GC, which would show composition driving it:
  (a) linear adjustment: per class, genome GC is centred; the contrast is regressed on the centred GC across all genomes
      (pooled OLS slope, Spearman rho); the adjusted contrast = contrast - slope * centred GC. Class medians, the pooled
      median and the sign test are repeated on the adjusted contrast, with the same thresholds as README.md.
  (b) GC strata: within each class, genomes are split at the class median GC; the pooled median contrast is reported for the
      low and high halves.
Limit: this adjusts for genome GC, not for the amino-acid composition of each protein (the cube has no per-protein composition).
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("analyze", HERE / "02_analyze.py")
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)


def adjust_for_gc(contrast, reps):
    """Return (frame with contrast, gc_dm, adjusted; slope; rho; n). Class-centred GC, classes with >= MIN_GENOMES genomes."""
    d = contrast.rename("contrast").to_frame().join(reps.set_index("species_prefix")[["CLASS", "gc_pct"]], how="inner").dropna()
    d = d[d.groupby("CLASS")["contrast"].transform("size") >= A.MIN_GENOMES].copy()
    d["gc_dm"] = d.gc_pct - d.groupby("CLASS").gc_pct.transform("mean")
    slope = float(np.polyfit(d.gc_dm, d.contrast, 1)[0]) if len(d) > 2 else float("nan")
    rho = float(stats.spearmanr(d.gc_dm, d.contrast)[0]) if len(d) > 2 else float("nan")
    d["adjusted"] = d.contrast - slope * d.gc_dm
    return d, slope, rho, len(d)


def gc_strata(d):
    d = d.copy()
    d["half"] = np.where(d.gc_pct >= d.groupby("CLASS").gc_pct.transform("median"), "high_GC", "low_GC")
    return d.groupby("half").contrast.agg(n_genomes="size", median="median").reset_index()


def run(cube, reps, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    nuc = cube.loc_ == "nucleus"
    specs = [("H1", matched_pair(cube, cube.tf & nuc, ~cube.tf & nuc), "+", 0.90, 0.60, 0.05),
             ("H2", matched_pair(cube, cube.unann, ~cube.unann), "+", 0.90, 0.60, None),
             ("H3", matched_pair(cube, cube.sec, ~cube.sec), "-", 0.90, 0.60, None)]
    rows, strata = [], []
    for name, c, direction, supp, kill, pooled_min in specs:
        d, slope, rho, n = adjust_for_gc(c, reps)
        for label, col in (("as_preregistered", "contrast"), ("gc_adjusted", "adjusted")):
            cm = A.class_medians(d[col].rename("contrast"), reps)
            frac, ncl, p = A.sign_test(cm["median"], direction)
            pooled = float(d[col].median())
            ok = True if pooled_min is None else pooled >= pooled_min
            rows.append({"hypothesis": name, "version": label, "n_genomes": n, "n_classes": ncl, "share_classes_in_direction": frac,
                         "pooled_median": pooled, "p_sign": p, "verdict_by_original_rule": A.verdict_sign(frac, supp, kill, ok),
                         "gc_slope": slope, "gc_spearman_rho": rho})
        s = gc_strata(d)
        s.insert(0, "hypothesis", name)
        strata.append(s)
    r = pd.DataFrame(rows)
    r.to_csv(out / "gc_adjusted_H1_H3.tsv", sep="\t", index=False)
    s = pd.concat(strata)
    s.to_csv(out / "gc_strata_H1_H3.tsv", sep="\t", index=False)
    print(r.to_string(index=False))
    print(s.to_string(index=False))
    return r, s


def matched_pair(cube, m1, m2):
    return A.matched_contrast(cube, m1, m2)


def main(cube_dir, out):
    cube = pd.read_parquet(Path(cube_dir) / "cube.parquet").rename(columns={"loc": "loc_"})
    reps = pd.read_parquet(Path(cube_dir) / "reps.parquet")
    run(cube, reps, out)


if __name__ == "__main__":
    main(*sys.argv[1:3])
