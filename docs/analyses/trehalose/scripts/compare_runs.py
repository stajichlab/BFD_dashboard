#!/usr/bin/env python3
"""Compare the first run (data/, MotherDuck, laptop) with the J13 rebuild (data_j13/): prints the old-versus-new values quoted in README section 0.

Usage: compare_runs.py <trehalose analysis dir>   (reads small result tables only)
"""
import csv
import sys
from pathlib import Path

base = Path(sys.argv[1])


def rd(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


for tag in ("data", "data_j13"):
    d = base / tag
    print("=====", tag)
    ms = {r["module"]: r for r in rd(d / "th_module_summary.tsv")}
    print("presence %", {k: round(100 * float(ms[k]["frac_genomes_present"]), 1) for k in ("TPS", "TPP", "NTH", "ATH", "GSY", "GPH")})
    pairs = {(r["module_a"], r["module_b"]): round(float(r["rho_class_demeaned"]), 2) for r in rd(d / "th_module_pairs.tsv")}
    print("class-demeaned rho", {k: v for k, v in pairs.items() if k in [("TPP", "TPS"), ("GPD", "GSY"), ("GPH", "GSY"), ("GSY", "UGP")]})
    al = rd(d / "th_allometry_modules.tsv")
    print("allometry percentile", {r["module"]: round(float(r["percentile_in_null"]) * 100) for r in al
                                   if r["module"] in ("TPS", "TPP", "NTH", "GPD", "ATH", "PGI (control)", "HSP70 (control)")})
    a = rd(d / "th_predictive_auc.tsv")
    print("AUC table rows", len(a), "q<0.1:", sum(float(x["q_BH"]) < 0.1 for x in a), "p<0.05:", sum(float(x["p"]) < 0.05 for x in a))
    for x in a:
        if x["marker"] in ("NTH_k", "GPH_k", "TPS_over_GSY") and x["focal_set"] in ("xerophile genera", "entomopathogen genera"):
            print("   ", x["focal_set"], x["marker"], "AUC", round(float(x["auc"]), 3), "p", round(float(x["p"]), 4), "q", round(float(x["q_BH"]), 3))
    c = rd(d / "th_coloss_vs_TPS.tsv")
    print("co-absence rows", len(c), "q==0:", sum(float(x["q_BH"]) == 0 for x in c), "q<0.1:", sum(float(x["q_BH"]) < 0.1 for x in c),
          "q<0.01:", sum(float(x["q_BH"]) < 0.01 for x in c))
    if "q_BH_class" in c[0]:
        cc = [x for x in c if x["q_BH_class"] not in ("", "nan")]
        print("   class-stratified null: rows", len(cc), "q<0.1:", sum(float(x["q_BH_class"]) < 0.1 for x in cc),
              "q<0.01:", sum(float(x["q_BH_class"]) < 0.01 for x in cc))
    for x in c:
        if x["target"] in ("PF02358", "PF01204", "PF05693", "PF03632"):
            print("   ", x["target"], "absent TPS", x["n_focal_absent"], "obs", x["obs_coabsent"], "exp", x["expected"], "obs/exp", x["obs_exp"],
                  "q", x["q_BH"], "| class null exp", x.get("expected_class", ""), "obs/exp", x.get("obs_exp_class", ""), "q", x.get("q_BH_class", ""))
    tc = rd(d / "th_trait_contrasts.tsv")
    print("trait contrasts rows", len(tc), "(q_BH column present)" if "q_BH" in tc[0] else "(no q column)")
