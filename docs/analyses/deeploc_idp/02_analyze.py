#!/usr/bin/env python3
"""Controls C1-C5 and hypotheses H1-H7 of the DeepLoc/IDR analysis, by the rules in README.md.

Usage: 02_analyze.py CUBE_DIR FUNGALTRAITS_CSV OUT_DIR
Reads reps/cube/ssp_genome parquet files written by 01_build_cube.py.

Matched contrast of two protein categories in one genome: in each length bin where both categories have at least MIN_N proteins,
the difference of IDR_all0 (sum of IDP_fraction / number of proteins; no row counts as 0) is taken, and the genome contrast is the
weighted mean over bins with weights min(n1, n2); a genome needs total weight >= MIN_W. Class summaries use classes with at least
MIN_GENOMES genomes with a valid contrast.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

MIN_N, MIN_W, MIN_GENOMES, MIN_GROUP = 5, 50, 8, 6
LOC_CLASSES_NUCLEAR = "nucleus"


# ---------------------------------------------------------------- generic helpers
def bh(p):
    """Benjamini-Hochberg q-values (monotone), same order as p."""
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, idx in zip(range(n, 0, -1), order[::-1]):
        prev = min(prev, p[idx] * n / rank)
        q[idx] = prev
    return q


def matched_contrast(cube, m1, m2, min_n=MIN_N, min_w=MIN_W):
    """Per-genome weighted mean over length bins of IDR_all0(category 1) - IDR_all0(category 2)."""
    a = cube[m1].groupby(["species_prefix", "lbin"]).agg(n1=("n", "sum"), s1=("sum_idr", "sum"))
    b = cube[m2].groupby(["species_prefix", "lbin"]).agg(n2=("n", "sum"), s2=("sum_idr", "sum"))
    j = a.join(b, how="inner")
    j = j[(j.n1 >= min_n) & (j.n2 >= min_n)]
    j["w"] = np.minimum(j.n1, j.n2)
    j["d"] = (j.s1 / j.n1 - j.s2 / j.n2) * j.w
    g = j.groupby("species_prefix").agg(d=("d", "sum"), w=("w", "sum"))
    g = g[g.w >= min_w]
    return (g.d / g.w).rename("contrast")


def class_medians(contrast, reps, min_genomes=MIN_GENOMES):
    d = contrast.to_frame().join(reps.set_index("species_prefix")["CLASS"], how="inner")
    t = d.groupby("CLASS")["contrast"].agg(n_genomes="size", median="median")
    return t[t.n_genomes >= min_genomes].reset_index()


def sign_test(medians, direction):
    """One-sided sign test on class medians in the predicted direction; returns (fraction in direction, n classes, p)."""
    x = np.asarray(medians, float)
    k = int((x > 0).sum() if direction == "+" else (x < 0).sum())
    n = len(x)
    if n == 0:
        return float("nan"), 0, float("nan")
    return k / n, n, float(stats.binomtest(k, n, 0.5, alternative="greater").pvalue)


def verdict_sign(frac, supp, kill, pooled_ok=True):
    if frac != frac:
        return "not evaluable"
    if frac >= supp and pooled_ok:
        return "supported"
    if frac < kill:
        return "killed"
    return "mixed"


def demeaned_spearman(df, x, y, min_genomes=MIN_GENOMES):
    """Spearman rho of class-demeaned x and y across genomes in classes with at least min_genomes genomes."""
    d = df.copy()
    keep = d.groupby("CLASS")[x].transform("size") >= min_genomes
    d = d[keep]
    for c in (x, y):
        d[c + "_dm"] = d[c] - d.groupby("CLASS")[c].transform("mean")
    rho, p = stats.spearmanr(d[x + "_dm"], d[y + "_dm"])
    p_one = p / 2 if rho > 0 else 1 - p / 2
    return float(rho), float(p_one), len(d)


def boot_ci(x, n=2000, seed=1):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    m = [np.median(rng.choice(x, len(x))) for _ in range(n)]
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


# ---------------------------------------------------------------- analysis
def genome_table(cube, reps):
    g = cube.groupby("species_prefix").agg(n=("n", "sum"), s=("sum_idr", "sum")).join(
        cube[cube.loc_ == "nucleus"].groupby("species_prefix")["n"].sum().rename("n_nuc"), how="left")
    g["n_nuc"] = g["n_nuc"].fillna(0)
    g["idr_all0"] = g.s / g.n
    g["nuclear_share"] = g.n_nuc / g.n
    g["log_proteins"] = np.log10(g.n)
    return reps.set_index("species_prefix").join(g, how="inner").reset_index()


def controls(cube, reps):
    rows = []
    cls = reps.set_index("species_prefix")["CLASS"]
    # C1: TFs are nuclear (class median of the per-genome share)
    tf = cube[cube.tf].groupby("species_prefix").agg(n=("n", "sum"))
    tfn = cube[cube.tf & (cube.loc_ == "nucleus")].groupby("species_prefix")["n"].sum()
    share = (tfn / tf.n).fillna(0).rename("share").to_frame().join(cls)
    ok = share.groupby("CLASS")["share"].agg(["size", "median"])
    ok = ok[ok["size"] >= MIN_GENOMES]
    frac = float((ok["median"] >= 0.90).mean()) if len(ok) else float("nan")
    rows.append(("C1", "TF proteins DeepLoc nuclear: share of classes with median >= 0.90", frac, 0.90, frac >= 0.90))
    # C2
    m = cube[cube.mito_carrier]
    v = float(m[m.loc_ == "mitochondrion"].n.sum() / m.n.sum()) if m.n.sum() else float("nan")
    rows.append(("C2", "PF00153 proteins DeepLoc mitochondrion", v, 0.70, v >= 0.70))
    # C3
    s = cube[cube.sec]
    v = float(s[s.loc_.isin(["extracellular", "cell_membrane"])].n.sum() / s.n.sum()) if s.n.sum() else float("nan")
    rows.append(("C3", "strict-secretome proteins DeepLoc extracellular or membrane", v, 0.60, v >= 0.60))
    # C4
    e = cube[cube.effector_fam]
    v = float(e[e.sec].n.sum() / e.n.sum()) if e.n.sum() else float("nan")
    rows.append(("C4", "effector-like families recovered by the secretion filter", v, 0.70, v >= 0.70))
    # C5
    c = matched_contrast(cube, cube.loc_ == "nucleus", cube.loc_ == "mitochondrion")
    cm = class_medians(c, reps)
    f, n, p = sign_test(cm["median"], "+")
    rows.append(("C5", "nuclear minus mitochondrial IDR_all0 positive: share of classes", f, 0.90, f >= 0.90))
    return pd.DataFrame(rows, columns=["control", "what", "observed", "criterion", "passed"])


def run(cube, reps, ssp, lifestyle, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    gen = genome_table(cube, reps)
    ctl = controls(cube, reps)
    ctl.to_csv(out / "controls.tsv", sep="\t", index=False)
    res = []

    def sign_h(name, text, c, direction, supp, kill, pooled_rule=None):
        cm = class_medians(c, reps)
        cm.to_csv(out / f"{name}_class_contrasts.tsv", sep="\t", index=False)
        frac, n, p = sign_test(cm["median"], direction)
        pooled = float(c.median()) if len(c) else float("nan")
        pooled_ok = True if pooled_rule is None else pooled_rule(pooled)
        res.append({"hypothesis": name, "what": text, "n": n, "statistic": "share of classes in predicted direction",
                    "value": frac, "pooled_median": pooled, "p": p, "verdict": verdict_sign(frac, supp, kill, pooled_ok)})

    nuc = cube.loc_ == "nucleus"
    sign_h("H1", "TF minus other nuclear IDR_all0", matched_contrast(cube, cube.tf & nuc, ~cube.tf & nuc), "+", 0.90, 0.60,
           lambda m: m >= 0.05)
    sign_h("H2", "unannotated minus annotated IDR_all0", matched_contrast(cube, cube.unann, ~cube.unann), "+", 0.90, 0.60)
    sign_h("H3", "strict secretome minus other IDR_all0", matched_contrast(cube, cube.sec, ~cube.sec), "-", 0.90, 0.60)

    rho, p, n = demeaned_spearman(gen, "idr_all0", "gc_pct")
    res.append({"hypothesis": "H4", "what": "genome IDR_all0 vs GC, class-demeaned", "n": n, "statistic": "Spearman rho",
                "value": rho, "pooled_median": np.nan, "p": p,
                "verdict": "supported" if rho >= 0.30 else ("killed" if abs(rho) < 0.15 else "mixed")})
    rho, p, n = demeaned_spearman(gen, "nuclear_share", "log_proteins")
    res.append({"hypothesis": "H5", "what": "nuclear share vs log proteome size, class-demeaned", "n": n, "statistic": "Spearman rho",
                "value": rho, "pooled_median": np.nan, "p": p,
                "verdict": "supported" if rho >= 0.20 else ("killed" if rho <= 0.05 else "mixed")})

    # H6: Microsporidia against the pooled rest, matched by length bin
    micro = set(reps[reps.PHYLUM == "Microsporidia"].species_prefix)
    bins = cube.groupby(["species_prefix", "lbin"]).agg(n=("n", "sum"), s=("sum_idr", "sum")).reset_index()
    rest = bins[~bins.species_prefix.isin(micro)].groupby("lbin").agg(n=("n", "sum"), s=("s", "sum"))
    ref = (rest.s / rest.n).rename("ref")
    mb = bins[bins.species_prefix.isin(micro)].join(ref, on="lbin")
    mb["d"] = (mb.s / mb.n - mb.ref) * mb.n
    c6 = (mb.groupby("species_prefix").d.sum() / mb.groupby("species_prefix").n.sum())
    if len(c6) >= 3:
        lo, hi = boot_ci(c6)
        med = float(c6.median())
        w = stats.wilcoxon(c6, alternative="less").pvalue if len(c6) >= 5 else np.nan
        v6 = "supported" if (med < 0 and hi < 0) else ("killed" if med >= 0 else "mixed")
        res.append({"hypothesis": "H6", "what": f"Microsporidia minus rest IDR_all0 (CI {lo:.3f} to {hi:.3f})", "n": len(c6),
                    "statistic": "median genome contrast", "value": med, "pooled_median": med, "p": float(w), "verdict": v6})
        c6.rename("contrast").to_frame().to_csv(out / "H6_microsporidia_contrasts.tsv", sep="\t")
    else:
        res.append({"hypothesis": "H6", "what": "Microsporidia minus rest", "n": len(c6), "statistic": "median genome contrast",
                    "value": np.nan, "pooled_median": np.nan, "p": np.nan, "verdict": "not evaluable"})

    # H7: effector-like candidates per 1,000 proteins, plant pathogen minus saprotroph, within class
    d = reps.merge(ssp, on="species_prefix").merge(lifestyle, on="GENUS", how="left")
    d["rate"] = 1000 * d.n_ssp / d.n_proteins
    d["group"] = np.where(d.primary_lifestyle == "plant_pathogen", "plant_pathogen",
                          np.where(d.primary_lifestyle.fillna("").str.endswith("_saprotroph"), "saprotroph", "other"))
    rows7 = []
    for cl, sub in d.groupby("CLASS"):
        pp, sp = sub[sub.group == "plant_pathogen"].rate, sub[sub.group == "saprotroph"].rate
        if len(pp) >= MIN_GROUP and len(sp) >= MIN_GROUP:
            rows7.append({"CLASS": cl, "n_plant_pathogen": len(pp), "n_saprotroph": len(sp),
                          "median_pp": pp.median(), "median_sap": sp.median(), "median": pp.median() - sp.median()})
    t7 = pd.DataFrame(rows7)
    t7.to_csv(out / "H7_class_contrasts.tsv", sep="\t", index=False)
    if len(t7) >= 3:
        frac, n, p = sign_test(t7["median"], "+")
        v7 = "supported" if frac >= 0.75 else ("killed" if frac < 0.50 else "mixed")
        res.append({"hypothesis": "H7", "what": "effector-like candidates per 1,000 proteins, plant pathogen minus saprotroph",
                    "n": n, "statistic": "share of classes in predicted direction", "value": frac,
                    "pooled_median": float(t7["median"].median()), "p": p, "verdict": v7})
    else:
        res.append({"hypothesis": "H7", "what": "effector-like candidates, plant pathogen minus saprotroph", "n": len(t7),
                    "statistic": "share of classes in predicted direction", "value": np.nan, "pooled_median": np.nan,
                    "p": np.nan, "verdict": "not evaluable (fewer than 3 evaluable classes)"})
    r = pd.DataFrame(res)
    pv = r["p"].astype(float)
    ok = pv.notna()
    r["q_BH"] = np.nan
    r.loc[ok, "q_BH"] = bh(pv[ok].values)
    r.to_csv(out / "hypotheses.tsv", sep="\t", index=False)
    (d.groupby(["CLASS", "group"]).agg(genera=("rate", "size"), median_rate=("rate", "median"))
        .reset_index().to_csv(out / "H7_rate_by_class_and_lifestyle.tsv", sep="\t", index=False))
    lines = ["CONTROLS"] + [f"{c.control}: {c.what}: {c.observed:.3f} (criterion {c.criterion}) -> "
                           f"{'pass' if c.passed else 'FAIL'}" for c in ctl.itertuples()]
    lines += ["", "HYPOTHESES"] + [f"{x.hypothesis}: {x.what}: {x.statistic} {x.value:.3f} (n={x.n}, p={x.p:.2g}, q={x.q_BH:.2g}) -> "
                                   f"{x.verdict}" for x in r.itertuples()]
    (out / "verdicts.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return ctl, r


def main(cube_dir, ft_csv, out):
    cube_dir = Path(cube_dir)
    cube = pd.read_parquet(cube_dir / "cube.parquet").rename(columns={"loc": "loc_"})  # 'loc' clashes with DataFrame.loc
    reps = pd.read_parquet(cube_dir / "reps.parquet")
    ssp = pd.read_parquet(cube_dir / "ssp_genome.parquet")
    ft = pd.read_csv(ft_csv, encoding="utf-8")
    ft.columns = [c.strip() for c in ft.columns]
    lifestyle = ft.drop_duplicates("GENUS")[["GENUS", "primary_lifestyle"]]
    run(cube, reps, ssp, lifestyle, out)


if __name__ == "__main__":
    main(*sys.argv[1:4])
