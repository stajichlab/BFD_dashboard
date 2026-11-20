#!/usr/bin/env python3
"""T4 and T6 by the rules in README.md.

Usage: 02_analyze.py TABLE_DIR FUNGALTRAITS_CSV OUT_DIR
TABLE_DIR holds reps.parquet, copies.parquet, names.parquet from 01_build_tables.py.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3] / "bin" / "summary"))   # bin/summary/taxonomy.py (subphylum derived from class)
from taxonomy import CLASS_TO_SUBPHYLUM  # noqa: E402

TPS, TPP, NTH, ATH = "PF00982", "PF02358", "PF01204", "PF03632"
GSY, GPH = "PF05693", "PF00343"
S8, GH18 = "PF00082", "PF00704"
MIN_SIDE = 5
N_BOOT = 2000


def cliffs_delta(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    gt = sum((x > b).sum() for x in a)
    lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (len(a) * len(b))


def boot_median_diff(a, b, seed=1):
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(N_BOOT)]
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def contrast(a, b):
    lo, hi = boot_median_diff(a, b)
    return {"n_focal": len(a), "n_ref": len(b), "median_focal": float(np.median(a)), "median_ref": float(np.median(b)),
            "diff": float(np.median(a) - np.median(b)), "ci_lo": lo, "ci_hi": hi, "cliffs_delta": float(cliffs_delta(a, b)),
            "p": float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)}


def t4_verdict(tps, ath, ath_share_focal):
    """Pre-registered rules. tps, ath: contrast dicts; ath_share_focal: share of focal genera with >= 1 ATH copy."""
    a = "supported" if (tps["diff"] > 0 and tps["ci_lo"] > 0 and tps["cliffs_delta"] >= 0.33) else (
        "killed" if tps["cliffs_delta"] <= 0 else "mixed")
    if ath_share_focal < 0.50:
        b = "killed"
    elif ath_share_focal >= 0.75 and ath["diff"] > 0 and ath["cliffs_delta"] >= 0.33:
        b = "supported"
    else:
        b = "mixed"
    overall = "supported" if (a == "supported" and b == "supported") else ("killed" if b == "killed" else "mixed")
    return a, b, overall


def retention_percentile(copies, in_red, in_ref, fam_ids, acc, eps=0.01):
    """Share of the other near-universal families whose retention ratio is below that of `acc` (low = dropped more)."""
    m_red = copies.loc[in_red, fam_ids].mean()
    m_ref = copies.loc[in_ref, fam_ids].mean()
    r = (m_red + eps) / (m_ref + eps)
    others = r.drop(acc, errors="ignore")
    return float((others < r[acc]).mean()), float(r[acc])


def t6_verdict(red_pcts, ctl_pcts):
    """red_pcts, ctl_pcts: NTH percentiles of the evaluable reduction / control contrasts (None for not evaluable)."""
    red = [p for p in red_pcts if p is not None]
    ctl = [p for p in ctl_pcts if p is not None]
    if len(red) < 3 or len(ctl) < 2:
        return "not evaluable (fewer than 3 reduction or 2 control contrasts)"
    low = sum(p <= 0.10 for p in red)
    high = sum(p > 0.25 for p in red)
    if low >= 3 and all(p > 0.25 for p in ctl):
        return "supported"
    if high >= 3 or all(p <= 0.10 for p in ctl):
        return "killed"
    return "mixed"


def load(table_dir, ft_csv):
    table_dir = Path(table_dir)
    reps = pd.read_parquet(table_dir / "reps.parquet")
    cp = pd.read_parquet(table_dir / "copies.parquet")
    names = pd.read_parquet(table_dir / "names.parquet")
    wide = cp.pivot_table(index="species_prefix", columns="acc", values="n", aggfunc="sum", fill_value=0)
    wide = wide.reindex(reps.species_prefix, fill_value=0)
    ft = pd.read_csv(ft_csv, encoding="utf-8")
    ft.columns = [c.strip() for c in ft.columns]
    ft = ft.drop_duplicates("GENUS")[["GENUS", "primary_lifestyle"]]
    reps = reps.merge(ft, on="GENUS", how="left")
    reps["subphylum"] = reps.CLASS.map(CLASS_TO_SUBPHYLUM)
    return reps, wide, names


def run(reps, wide, names, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    reps = reps.set_index("species_prefix")
    wide.index.name = "species_prefix"
    per1000 = wide.div(reps.loc[wide.index, "n_proteins"], axis=0) * 1000
    lines = []

    # ---------------- T4
    sor = reps[(reps.CLASS == "Sordariomycetes") & reps.primary_lifestyle.notna()]
    focal = sor[sor.primary_lifestyle == "animal_parasite"].index
    ref = sor[sor.primary_lifestyle != "animal_parasite"].index
    rows = []
    for name, acc in [("TPS", TPS), ("ATH", ATH), ("TPP", TPP), ("NTH", NTH), ("GSY", GSY), ("GPH", GPH), ("S8 protease", S8), ("GH18 chitinase", GH18)]:
        if acc in per1000.columns and len(focal) >= MIN_SIDE and len(ref) >= MIN_SIDE:
            rows.append({"family": name, "pfam": acc, **contrast(per1000.loc[focal, acc], per1000.loc[ref, acc])})
    t4 = pd.DataFrame(rows)
    t4.to_csv(out / "t4_contrasts.tsv", sep="\t", index=False)
    if len(t4):
        c = t4.set_index("family")
        ath_share = float((wide.loc[focal, ATH] >= 1).mean()) if ATH in wide.columns else 0.0
        a, b, overall = t4_verdict(c.loc["TPS"], c.loc["ATH"], ath_share)
        ctrl = all(c.loc[f, "cliffs_delta"] > 0 and c.loc[f, "ci_lo"] > 0 for f in ("S8 protease", "GH18 chitinase"))
        lines += [f"T4 focal set: {len(focal)} Sordariomycetes genera with primary_lifestyle animal_parasite; reference {len(ref)} other labelled genera",
                  f"T4a TPS per 1,000: diff {c.loc['TPS','diff']:.3f} (CI {c.loc['TPS','ci_lo']:.3f} to {c.loc['TPS','ci_hi']:.3f}), delta {c.loc['TPS','cliffs_delta']:.2f} -> {a}",
                  f"T4b ATH: {ath_share:.0%} of focal genera carry >= 1 copy; per 1,000 diff {c.loc['ATH','diff']:.3f}, delta {c.loc['ATH','cliffs_delta']:.2f} -> {b}",
                  "T4c trehalose transport: not tested (no specific family in the tables)",
                  f"positive control (S8 protease and GH18 chitinase higher in focal): {'passes' if ctrl else 'FAILS'}",
                  f"T4 overall: {overall}"]
    else:
        lines.append("T4: not evaluable (focal or reference set smaller than %d genera)" % MIN_SIDE)

    # ---------------- T6
    prev = names.set_index("acc")["prevalence"]
    univ = [a for a in wide.columns if prev.get(a, 0) >= 0.90]
    sacc = reps.subphylum.fillna("")
    groups = {
        "R1 Kickxellales": (reps.ORDER_ == "Kickxellales", reps.ORDER_ == "Mucorales"),
        "R2 Microsporidia": (reps.PHYLUM == "Microsporidia", reps.PHYLUM == "Chytridiomycota"),
        "R3 Saccharomycotina": (sacc == "Saccharomycotina", sacc == "Pezizomycotina"),
        "R4 Malasseziomycetes": (reps.CLASS == "Malasseziomycetes", reps.CLASS.isin(["Ustilaginomycetes", "Exobasidiomycetes"])),
        "K1 Eurotiomycetes": (reps.CLASS == "Eurotiomycetes", reps.CLASS == "Sordariomycetes"),
        "K2 Dothideomycetes": (reps.CLASS == "Dothideomycetes", reps.CLASS == "Leotiomycetes"),
    }
    hk = [a for a in univ if str(names.set_index("acc").loc[a, "pfam_id"]).startswith("Ribosomal")]
    t6, red_p, ctl_p = [], [], []
    for name, (m_red, m_ref) in groups.items():
        ir, ig = reps.index[m_red.values], reps.index[m_ref.values]
        ir, ig = [x for x in ir if x in wide.index], [x for x in ig if x in wide.index]
        is_ctl = name.startswith("K")
        if len(ir) < MIN_SIDE or len(ig) < MIN_SIDE or NTH not in wide.columns:
            (ctl_p if is_ctl else red_p).append(None)
            t6.append({"contrast": name, "n_reduced": len(ir), "n_reference": len(ig), "note": "not evaluable"})
            continue
        ratio = float(reps.loc[ir, "n_proteins"].median() / reps.loc[ig, "n_proteins"].median())
        fam_ids = univ if NTH in univ else univ + [NTH]
        row = {"contrast": name, "n_reduced": len(ir), "n_reference": len(ig), "proteome_ratio": ratio}
        for label, acc in [("NTH", NTH), ("TPS", TPS), ("GSY", GSY), ("GPH", GPH)]:
            if acc in wide.columns:
                fids = fam_ids if acc in fam_ids else fam_ids + [acc]
                pct, r = retention_percentile(wide, ir, ig, fids, acc)
                row[f"{label}_percentile"], row[f"{label}_ratio"] = pct, r
        if hk:
            m_red, m_ref = wide.loc[ir, hk].mean(), wide.loc[ig, hk].mean()
            row["housekeeping_median_ratio"] = float(((m_red + 0.01) / (m_ref + 0.01)).median())
        (ctl_p if is_ctl else red_p).append(row.get("NTH_percentile"))
        t6.append(row)
    t6 = pd.DataFrame(t6)
    t6.to_csv(out / "t6_retention.tsv", sep="\t", index=False)
    lines += ["", f"T6: {len(univ)} near-universal families"]
    for r in t6.itertuples():
        if hasattr(r, "NTH_percentile") and r.NTH_percentile == r.NTH_percentile:
            lines.append(f"  {r.contrast}: NTH percentile {r.NTH_percentile:.2f} (ratio {r.NTH_ratio:.2f}); proteome ratio {r.proteome_ratio:.2f}"
                         f"; TPS percentile {getattr(r, 'TPS_percentile', float('nan')):.2f}")
        else:
            lines.append(f"  {r.contrast}: not evaluable")
    lines.append(f"T6 overall: {t6_verdict(red_p, ctl_p)}")
    (out / "verdicts.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return t4, t6, lines


def main(table_dir, ft_csv, out):
    reps, wide, names = load(table_dir, ft_csv)
    run(reps, wide, names, out)


if __name__ == "__main__":
    main(*sys.argv[1:4])
