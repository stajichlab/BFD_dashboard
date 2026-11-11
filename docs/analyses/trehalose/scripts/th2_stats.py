#!/usr/bin/env python3
"""Second pass: allometry against the kingdom-wide null, pathway balance ratios,
and the desiccation / lichenization / insect-pathogen contrasts.

Inputs: th_genus_medians.tsv (from th1_survey.py), the kingdom genus x Pfam
table exported as th_loose_genus.csv, and the atlas genus structural table
docs/analyses/kickxellomycota/data/genus_medians_all.tsv (proteome size per genus).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
FIGS = HERE.parent / "figs"
RNG = np.random.default_rng(7)


def cliffs(a, b):
    a = np.asarray(a, float)[:, None]
    b = np.asarray(b, float)[None, :]
    return float(np.mean(np.sign(a - b)))


def boot(a, b, n=5000):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = np.array([np.median(RNG.choice(a, a.size)) - np.median(RNG.choice(b, b.size)) for _ in range(n)])
    return float(np.median(a) - np.median(b)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def contrast(focal, sisters, name, out, extra=""):
    if len(focal) < 3 or len(sisters) < 5:
        return
    d, lo, hi = boot(focal, sisters)
    out.append(dict(contrast=name, focal_group=extra, n_focal=len(focal), n_sisters=len(sisters),
                    med_focal=float(np.median(focal)), med_sisters=float(np.median(sisters)),
                    diff=d, ci_lo=lo, ci_hi=hi, cliffs_delta=cliffs(focal, sisters),
                    mannwhitney_p=float(stats.mannwhitneyu(focal, sisters).pvalue)))


def main():
    G = pd.read_csv(DATA / "th_genus_medians.tsv", sep="\t")
    X = pd.read_csv(DATA / "th_genome_matrix.tsv", sep="\t")
    L = pd.read_csv(DATA / "th_loose_genus.csv")           # genus x PF00982/PF02358 hit counts
    strc = pd.read_csv("/Users/jstajich/projects/BFD_dashboard/docs/analyses/kickxellomycota/data/"
                       "genus_medians_all.tsv", sep="\t")  # 1,441 genera: structural metrics
    G = G.merge(strc[["GENUS", "introns_per_transcript", "gene_density", "median_protein_len",
                      "masked_pct", "median_intergenic"]], on="GENUS", how="left")

    # --- pathway balance ratios ------------------------------------------------
    for a, b, nm in [("TPS", "GSY", "TPS_over_GSY"), ("TPS", "NTH", "TPS_over_NTH"),
                     ("TPS", "GPH", "TPS_over_GPH"), ("ATH", "NTH", "ATH_over_NTH")]:
        G[nm] = (G[a + "_k"] + 0.05) / (G[b + "_k"] + 0.05)
    G.to_csv(DATA / "th_genus_medians_enriched.tsv", sep="\t", index=False)

    # --- 1. allometry against the kingdom-wide family null --------------------
    # kw_genus_pfam mean_copies per genus for all 15,052 families vs genus proteome size.
    # th_allobin.csv: mean copies per family within deciles of proteome size
    # (built in MotherDuck by scripts/md_allobin.sql; 10 bins x 15k families)
    kp = pd.read_csv(DATA / "th_allobin.csv")
    slopes = {}
    for acc, s in kp.groupby("pfam_acc"):
        if s.n_bin_genus.min() >= 30:
            rho, p = stats.spearmanr(s.mean_copies, np.log10(s.mean_prot))
            slopes[acc] = rho
    SL = pd.Series(slopes, name="rho_copies_vs_log_proteome")
    SL.to_csv(DATA / "th_allometry_null.tsv", sep="\t")
    print(f"== allometry null: {len(SL)} families tested (>=300 genera), median rho {SL.median():+.3f}, "
          f"IQR {SL.quantile(.25):+.3f} to {SL.quantile(.75):+.3f}")
    fam = {"PF00982": "TPS", "PF02358": "TPP", "PF01204": "NTH", "PF03632": "ATH",
           "PF05693": "GSY", "PF00343": "GPH", "PF01704": "UGP", "PF07479": "GPD",
           "PF00011": "sHSP", "PF04119": "HSP12", "PF01185": "HYDRO_I", "PF05141": "DIT1",
           "PF00012": "HSP70 (control)", "PF00342": "PGI (control)", "PF01232": "MPD"}
    rows = []
    for acc, nm in fam.items():
        if acc in SL.index:
            rows.append((nm, acc, SL[acc], float((SL < SL[acc]).mean())))
    A = pd.DataFrame(rows, columns=["module", "pfam", "rho_copies_vs_proteome", "percentile_in_null"])
    print(A.sort_values("rho_copies_vs_proteome").to_string(index=False))
    A.to_csv(DATA / "th_allometry_modules.tsv", sep="\t", index=False)

    # --- 2. desiccation sets, class-matched -----------------------------------
    dry = {"Lecanoromycetes": ["xero_genus_set", "lichen"],
           "Dothideomycetes": ["xero_genus_set"], "Eurotiomycetes": ["xero_genus_set"]}
    XERO = {"Wallemia", "Xeromyces", "Baudoinia", "Knufia", "Cryomyces", "Antarctomyces",
            "Cryoendolithus", "Friedmanniomyces", "Hortaea", "Exophiala", "Aureobasidium",
            "Zygosaccharomyces", "Debaryomyces", "Coniochaeta", "Gymnoascus"}
    DRY_EP = r"(?i)antarct|arctic|desert|endolith|xero|halophil|salin|cavernicol|arid"
    G["dry_epithet"] = G.GENUS.map(X.groupby("GENUS").SPECIES.apply(
        lambda s: bool(s.fillna("").str.contains(DRY_EP).any()))).fillna(False)
    G["xero"] = G.GENUS.isin(XERO)
    res = []
    metrics = ["TPS_k", "TPP_k", "NTH_k", "ATH_k", "GSY_k", "GPH_k", "HSP12_k", "sHSP_k",
               "HYDRO_I_k", "DIT1_k", "TPS_over_GSY", "TPS_over_NTH", "ATH_over_NTH", "n_prot", "mb"]
    for label, mask in [("xerophile genera (curated)", G.xero),
                        ("dry-habitet species epithet", G.dry_epithet),
                        ("lichen-forming class Lecanoromycetes", G.CLASS == "Lecanoromycetes"),
                        ("osmotolerant yeasts (Zygosacch./Debaryomyces)",
                         G.GENUS.isin({"Zygosaccharomyces", "Debaryomyces", "Zygotorulaspora", "Nakaseomyces"})),
                        ("entomopathogen genera", G.GENUS.isin(
                            {"Beauveria", "Metarhizium", "Ophiocordyceps", "Cordyceps",
                             "Lecanicillium", "Purpureocillium", "Cordyceps", "Torrubiella"}))]:
        focal = G[mask]
        if focal.empty:
            continue
        for m in metrics:
            for cls, fsub in focal.groupby("CLASS"):
                sib = G[(G.CLASS == cls) & (~mask)]
                contrast(fsub[m].dropna(), sib[m].dropna(), m, res, extra=f"{label} | {cls}")
    R = pd.DataFrame(res)
    R.to_csv(DATA / "th_trait_contrasts.tsv", sep="\t", index=False)
    print("\n== trait contrasts (genus medians, class-matched sisters; |delta|>=0.4 and CI excluding 0) ==")
    keep = R[(R.cliffs_delta.abs() >= 0.4) & ((R.ci_lo > 0) | (R.ci_hi < 0))].sort_values("cliffs_delta")
    print(keep.to_string(index=False))

    # --- 3. credible losses: who, and are they clustered ----------------------
    print("\n== credible TPS/TPP loss (qcov>=0.985, BUSCO>=85, no Pfam or SwissProt hit) ==")
    cl = X[X.credible_TPS_loss].groupby(["PHYLUM", "CLASS"]).agg(n_genomes=("genome", "size"),
                                                                 n_genera=("GENUS", "nunique"))
    tot = X[X.keep].groupby(["PHYLUM", "CLASS"]).agg(n_gen_annot=("genome", "size"), n_genus=("GENUS", "nunique"))
    print(cl.join(tot).assign(freq=lambda d: (d.n_genomes / d.n_gen_annot).round(3)).to_string())
    X[X.credible_TPS_loss][["GENUS", "CLASS", "ord", "SPECIES", "qcov", "complete_pct",
                            "n_prot", "TPS", "TPP", "n_tps", "n_tpp", "NTH", "ATH"]].sort_values(
        ["CLASS", "GENUS"]).to_csv(DATA / "th_credible_loss_genomes.tsv", sep="\t", index=False)
    # co-loss of TPS and TPP, and of neighbours
    sub = X[X.keep]
    print("\nco-absence (credible) TPS&TPP:", int((sub.credible_TPS_loss & sub.credible_TPP_loss).sum()),
          "of", int(sub.credible_TPS_loss.sum()), "TPS-loss genomes")
    print("NTH still present in TPS-loss genomes:",
          float((sub[sub.credible_TPS_loss].NTH > 0).mean()))
    print("GSY still present in TPS-loss genomes:",
          float((sub[sub.credible_TPS_loss].GSY > 0).mean()))

    # --- 4. GH65 acid trehalase distribution ----------------------------------
    print("\n== ATH (GH65 acid trehalase) presence by class (n>=8 genera) ==")
    at = G.groupby("CLASS").agg(n_gen=("GENUS", "size"), ath_present=("ATH", lambda s: float((s > 0).mean())),
                                ath_med=("ATH_k", "median"), nth_med=("NTH_k", "median"))
    print(at[at.n_gen >= 8].sort_values("ath_present", ascending=False).to_string())
    at.to_csv(DATA / "th_ath_by_class.tsv", sep="\t")

    # --- 5. trehalose module and reduction, Kickxellomycota focus -------------
    print("\n== trehalose module in reduced genomes: quartiles of proteome size (all genera) ==")
    G["q"] = pd.qcut(G.n_prot, 4, labels=["Q1 smallest", "Q2", "Q3", "Q4 largest"])
    print(G.groupby("q", observed=True)[["n_prot", "TPS_k", "TPP_k", "NTH_k", "GSY_k", "GPH_k",
                                         "sHSP_k", "HSP12_k", "TPS_over_GSY", "introns_per_transcript"]]
          .median().round(3).to_string())
    for grp, sel in [("Saccharomycotina", G.yeast_origin == "Saccharomycotina"),
                     ("Kickxellomycetes", G.CLASS == "Kickxellomycetes"),
                     ("Mucoromycetes", G.CLASS == "Mucoromycetes"),
                     ("Microsporidia", G.PHYLUM == "Microsporidia")]:
        s = G[sel]
        if len(s) == 0:
            continue
        print(f"  {grp:18s} n_gen={len(s):3d} prot={s.n_prot.median():7.0f} TPS_k={s.TPS_k.median():.3f} "
              f"TPP_k={s.TPP_k.median():.3f} NTH_k={s.NTH_k.median():.3f} GSY_k={s.GSY_k.median():.3f} "
              f"TPS/GSY={s.TPS_over_GSY.median():.2f} introns/tx={s.introns_per_transcript.median():.2f}")

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
    ax[0].hist(SL.values, bins=60, color=".7")
    shown = [(SL[acc], nm) for acc, nm in fam.items() if acc in SL.index]
    for i, (v, nm) in enumerate(sorted(shown)):
        ax[0].axvline(v, lw=1.1, color="C0")
        ax[0].text(v + 0.008, 20 + 42 * (i % 7), nm, fontsize=7, rotation=90, va="center", color="C0")
    ax[0].set_xlabel("rho(copies, log proteome size) across genera")
    ax[0].set_title("Allometry vs the kingdom-wide family null")
    XERO = {"Wallemia", "Xeromyces", "Baudoinia", "Knufia", "Cryomyces", "Antarctomyces",
            "Cryoendolithus", "Friedmanniomyces", "Hortaea", "Exophiala", "Aureobasidium",
            "Zygosaccharomyces", "Debaryomyces", "Coniochaeta", "Gymnoascus"}
    lo, hi = 0.02, 30
    ax[1].scatter(G.GSY_k + .05, G.TPS_k + .05, s=8, label="genus median")
    xp = G[G.GENUS.isin(XERO)]
    ax[1].scatter(xp.GSY_k + .05, xp.TPS_k + .05, s=40, facecolors="none", edgecolors="crimson",
                  label="curated xerophiles")
    tg = G[G.TPS == 0]
    ax[1].scatter(tg.GSY_k + .05, np.full(len(tg), .045), marker="x", c="darkorange", s=28,
                  label="no TPS detected")
    ax[1].plot([lo, hi], [lo, hi], ls=":", c=".4", lw=1)
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlim(lo, hi); ax[1].set_ylim(.04, hi)
    ax[1].set_xlabel("GSY per 1,000 proteins (glycogen)")
    ax[1].set_ylabel("TPS per 1,000 proteins (trehalose)")
    ax[1].set_title("Trehalose-vs-glycogen branch balance")
    ax[1].legend(fontsize=7, loc="upper left")
    b = G.groupby("q", observed=True).TPS_k.median()
    ax[2].bar(range(4), b.values, tick_label=[str(i) for i in b.index])
    ax[2].set_ylabel("TPS per 1,000 proteins")
    ax[2].set_title("Trehalose density by proteome-size quartile")
    fig.tight_layout(); fig.savefig(FIGS / "f4_allometry_and_balance.png", dpi=150)


if __name__ == "__main__":
    main()
