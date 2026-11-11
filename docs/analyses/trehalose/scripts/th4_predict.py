#!/usr/bin/env python3
"""How predictive is the trehalose module of a dry-substrate lifestyle?

Each candidate marker is converted to a within-class percentile (so that the
dominant class-level allometry, the biggest single effect in analysis B, cannot
masquerade as a lifestyle signal), then compared between focal and control
genera with the Mann-Whitney AUC. Focal sets are deliberately small and curated
from names, so this is screening, not inference.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"

XERO = {"Wallemia", "Xeromyces", "Baudoinia", "Knufia", "Cryomyces", "Antarctomyces",
        "Cryoendolithus", "Friedmanniomyces", "Hortaea", "Exophiala", "Aureobasidium",
        "Zygosaccharomyces", "Debaryomyces", "Coniochaeta", "Gymnoascus"}
ENTOMO = {"Beauveria", "Metarhizium", "Ophiocordyceps", "Cordyceps", "Lecanicillium",
          "Purpureocillium", "Torrubiella", "Cordyceps", "Ophiocordyceps"}
MARKERS = ["TPS_k", "TPP_k", "NTH_k", "ATH_k", "GSY_k", "GPH_k", "UGP_k", "PGM_k", "GPD_k",
           "MIP_k", "MPD_k", "sHSP_k", "HSP12_k", "HYDRO_I_k", "HYDRO_II_k", "DIT1_k",
           "TPS_over_GSY", "TPS_over_NTH", "n_prot", "mb", "gc_pct", "introns_per_transcript"]


def main():
    G = pd.read_csv(DATA / "th_genus_medians_enriched.tsv", sep="\t")
    X = pd.read_csv(DATA / "th_genome_matrix.tsv", sep="\t")
    G["dry_epithet"] = G.GENUS.isin(set(X.loc[X.SPECIES.fillna("").str.contains(
        r"(?i)antarct|arctic|desert|endolith|xero|halophil|salin|cavernicol|arid"), "GENUS"].dropna()))
    G["xero"] = G.GENUS.isin(XERO)
    G["entomo"] = G.GENUS.isin(ENTOMO)
    G["lichen"] = G.CLASS.isin(["Lecanoromycetes", "Lecideomycetes", "Candelariomycetes",
                                "Coniocybomycetes", "Arthopyreniomycetes"])
    # within-class percentile for every marker
    G = G.set_index("GENUS")
    P = G.groupby("CLASS")[MARKERS].transform(lambda s: s.rank(pct=True))
    P = P.join(G[["CLASS", "xero", "dry_epithet", "entomo", "lichen"]])
    P = P.dropna(subset=MARKERS, how="all")

    rows = []
    for name, flag in [("xerophile genera", "xero"), ("dry species epithet", "dry_epithet"),
                       ("entomopathogen genera", "entomo"), ("lichen-forming classes", "lichen")]:
        for m in MARKERS:
            a = P.loc[P[flag] & P[m].notna(), m]
            b = P.loc[~P[flag] & P[m].notna(), m]
            if len(a) < 3 or len(b) < 10:
                continue
            u = stats.mannwhitneyu(a, b, alternative="two-sided")
            rows.append(dict(focal_set=name, marker=m, n_focal=len(a), n_ctrl=len(b),
                             auc=u.statistic / (len(a) * len(b)), p=u.pvalue,
                             med_focal_classpct=a.median(), med_ctrl_classpct=b.median()))
    R = pd.DataFrame(rows)
    k = len(R)
    R["q_BH"] = np.sort(R.p.values) * k / np.arange(1, k + 1)
    R["q_BH"] = np.minimum.accumulate(R.q_BH.values[::-1])[::-1]
    R = R.sort_values("auc", ascending=False)
    R.to_csv(DATA / "th_predictive_auc.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    for name in ["xerophile genera", "dry species epithet", "entomopathogen genera", "lichen-forming classes"]:
        s = R[R.focal_set == name].copy()
        s["|AUC-0.5|"] = (s.auc - .5).abs()
        print(f"\n== {name}: markers ranked by AUC distance from chance (within-class percentiles) ==")
        print(s.sort_values("|AUC-0.5|", ascending=False).head(8)[
            ["marker", "n_focal", "n_ctrl", "auc", "p", "q_BH", "med_focal_classpct"]].round(4).to_string(index=False))

    print("\n== how many markers reach q<0.1 for each focal set, and best single-marker AUC ==")
    for name in R.focal_set.unique():
        s = R[R.focal_set == name]
        best = s.loc[(s.auc - .5).abs().idxmax()]
        print(f"  {name:26s} markers q<0.1: {int((s.q_BH < 0.1).sum()):2d}/{len(s):2d}   "
              f"best {best.marker:22s} AUC {best.auc:.3f} (n focal {int(best.n_focal)})")


if __name__ == "__main__":
    main()
