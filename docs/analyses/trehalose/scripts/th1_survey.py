#!/usr/bin/env python3
"""Trehalose-pathway survey across the BFD kingdom, plus dessication case-controls.

Inputs are the CSVs in ../data/, exported from MotherDuck by
scripts/export_queries.sql (scratch tables bfd_scratch.th_*).

Conventions follow docs/summary_atlas/index.md: genus medians are the unit of
replication for cross-clade tests, effect sizes come first (median difference
with a 5,000-resample bootstrap CI, Cliff's delta), Mann-Whitney is screening
only, and copy numbers are normalised to proteome size (per 1,000 proteins).

Outputs: ../data/th_*.tsv and ../figs/*.png.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
FIGS = HERE.parent / "figs"
FIGS.mkdir(exist_ok=True)
ATLAS = Path("/Users/jstajich/projects/BFD_dashboard/docs/summary_atlas/tables")
KXDIR = Path("/Users/jstajich/projects/BFD_dashboard/docs/analyses/kickxellomycota/data")
RNG = np.random.default_rng(1)
N_BOOT = 5000

# --- module definition -------------------------------------------------------
# One representative Pfam accession per activity, chosen from the SwissProt
# annotation of characterised fungal enzymes (see README table 1). Summing
# several sub-domains of the same protein would overcount, so a single
# diagnostic domain is used per activity.
MODULES = {
    "TPS":      ("PF00982", "trehalose-6-P synthase (Glyco_transf_20)"),
    "TPP":      ("PF02358", "trehalose-6-P phosphatase (Trehalose_PPase)"),
    "NTH":      ("PF01204", "neutral trehalase (Trehalase, GH20-like)"),
    "NTH_Ca":   ("PF07492", "neutral trehalase Ca-binding subdomain"),
    "ATH":      ("PF03632", "acid trehalase (Glyco_hydro_65m, GH65)"),
    "ATH_N":    ("PF03636", "acid trehalase N-domain"),
    "ATH_C":    ("PF03633", "acid trehalase C-domain"),
    "GSY":      ("PF05693", "glycogen synthase"),
    "GPH":      ("PF00343", "glycogen phosphorylase"),
    "GDB":      ("PF06202", "glycogen debranching enzyme"),
    "GH13":     ("PF00128", "alpha-amylase / GH13 (incl. Glc3 branching, TreS)"),
    "UGP":      ("PF01704", "UTP-glucose-1-P uridylyltransferase (UDP-glucose)"),
    "PGM":      ("PF02878", "phosphoglucomutase (G1P<->G6P)"),
    "PGI":      ("PF00342", "glucose-6-P isomerase"),
    "HXK":      ("PF00349", "hexokinase"),
    "GPD":      ("PF07479", "glycerol-3-P dehydrogenase"),
    "MIP":      ("PF00230", "aquaglyceroporin (Fps1-type)"),
    "MPD":      ("PF01232", "mannitol-1-P dehydrogenase"),
    "sHSP":     ("PF00011", "small HSP (HSP20)"),
    "HSP70":    ("PF00012", "HSP70"),
    "HSP12":    ("PF04119", "hydrophilin HSP9/HSP12"),
    "HPF1":     ("PF11765", "haze protective factor / Hyphal_reg_CWP"),
    "CAT":      ("PF00199", "catalase"),
    "SOD_Cu":   ("PF00080", "Cu-Zn superoxide dismutase"),
    "HYDRO_I":  ("PF01185", "class I hydrophobin"),
    "HYDRO_II": ("PF06766", "class II hydrophobin"),
    "DIT1":     ("PF05141", "ascospore wall maturation DIT1/PvcA"),
    "SPO7":     ("PF03907", "Spo7 (spore cortex, meiosis exit)"),
}
N_UNIV = 1409  # size of the near-universal family set in scr.th_univ

# --- dessication-associated genus sets (curated, see README section 6) -------
XERO_GENERA = {  # free-living xerophiles / osmo-halotolerant / rock-inhabiting
    "Wallemia", "Xeromyces", "Baudoinia", "Knufia", "Cryomyces", "Antarctomyces",
    "Cryoendolithus", "Friedmanniomyces", "Hortaea", "Exophiala", "Aureobasidium",
    "Zygosaccharomyces", "Debaryomyces", "Coniochaeta", "Gymnoascus",
}
DRY_EPITHET = r"(?i)antarct|arctic|desert|endolith|xero|halophil|salin|cavernicol|arid|mumm"


def cliffs_delta(a, b):
    a = np.asarray(a, float)[:, None]
    b = np.asarray(b, float)[None, :]
    return float(np.mean(np.sign(a - b)))


def boot_diff(a, b, n=N_BOOT):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    d = np.empty(n)
    for i in range(n):
        d[i] = np.median(RNG.choice(a, a.size)) - np.median(RNG.choice(b, b.size))
    return float(np.median(a) - np.median(b)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def load():
    g = pd.read_csv(DATA / "th_genomes.csv")
    d = pd.read_csv(DATA / "th_dom.csv")
    wide = d.pivot_table(index="genome", columns="acc", values="n_prot", aggfunc="sum")
    X = g.set_index("genome").join(wide)
    for name, (acc, _) in MODULES.items():
        X[name] = X[acc].fillna(0) if acc in X.columns else 0
    sp = pd.read_csv(DATA / "th_sp.csv").set_index("genome")
    X = X.join(sp[["n_tps", "n_tpp", "n_nth", "n_ath", "n_gsy", "n_gph", "n_gpd", "n_mfn"]])
    uc = pd.read_csv(DATA / "th_univcov.csv").set_index("genome")
    X = X.join(uc[["n_univ_found"]])
    # share of the genome's own proteins that carry at least one Pfam hit: catches
    # genomes whose domain scan (not their annotation) is incomplete, which BUSCO misses
    pf = pd.read_csv(DATA / "th_pfprot.csv").set_index("genome")
    X = X.join(pf)
    X["n_prot_hit"] = X.n_prot_hit.fillna(0)
    X["frac_pfam"] = X.n_prot_hit / X.n_prot
    X["prot_per_mb"] = X.n_prot / (X.total_length_bp / 1e6)
    X["qcov"] = X.n_univ_found / N_UNIV
    # atlas Q-busco-rel: BUSCO complete at or above the 25th percentile of its own
    # class, so that real gene loss in reduced lineages is not treated as poor data
    X["busco_class_p25"] = X.groupby("CLASS").complete_pct.transform(lambda s: s.quantile(0.25))
    lab = pd.read_csv(ATLAS / "V31_genus_medians.tsv", sep="\t")[
        ["GENUS", "growth_form", "primary_lifestyle", "yeast_origin"]]
    X = X.reset_index()  # merge() would otherwise drop the genome index
    X = X.merge(lab, on="GENUS", how="left")
    X["mb"] = X.total_length_bp / 1e6
    X["k_prot"] = 1000.0 / X.n_prot
    for m in MODULES:
        X[m + "_k"] = X[m] * X.k_prot
    X["dry_epithet"] = X.SPECIES.fillna("").str.contains(DRY_EPITHET, regex=True)
    X["xero_genus"] = X.GENUS.isin(XERO_GENERA)
    # credible absence: neither a Pfam record at any e-value nor a SwissProt best
    # hit to the named enzyme, in a genome whose annotation is essentially complete
    X["tps_any"] = X.TPS > 0
    X["credible_TPS_loss"] = (X.qcov >= 0.985) & (X.complete_pct >= 85) & (X.TPS == 0) & (X.n_tps.fillna(0) == 0)
    X["credible_TPP_loss"] = (X.qcov >= 0.985) & (X.complete_pct >= 85) & (X.TPP == 0) & (X.n_tpp.fillna(0) == 0)
    X["keep"] = (X.qcov >= 0.98) & (X.complete_pct >= X.busco_class_p25)
    X["keep_strict"] = X.keep & (X.complete_pct >= 90)
    return X.reset_index()


def genus_table(X):
    K = X[X.keep].copy()
    num = [m + "_k" for m in MODULES] + list(MODULES) + ["mb", "n_prot", "med_len", "gc_pct",
                                                          "masked_pct", "frac_pfam", "prot_per_mb"]
    agg = {c: "median" for c in num}
    agg.update({"PHYLUM": "first", "CLASS": "first", "ord": "first",
                "growth_form": "first", "primary_lifestyle": "first", "yeast_origin": "first"})
    agg["genome"] = "size"
    G = K.groupby("GENUS").agg(agg).rename(columns={"genome": "n_sp", "ord": "ORDER_"})
    pres = K.groupby("GENUS").agg(**{f"f_{m.lower()}": (m, lambda s: float((s > 0).mean())) for m in MODULES})
    return G.join(pres).reset_index()


def main():
    cache = DATA / ".cache.pkl"
    if cache.exists():
        X, G = pd.read_pickle(cache)
    else:
        X = load()
        G = genus_table(X)
        pd.to_pickle((X, G), cache)
    X.to_csv(DATA / "th_genome_matrix.tsv", sep="\t", index=False)
    print(f"genomes {len(X)}  quality-kept {int(X.keep.sum())}  genera {X.GENUS.nunique()}")

    print("\n== 1. false absences track annotation completeness, not phylogeny ==")
    miss = X[X.TPS == 0]
    print(f"TPS absent in {len(miss)} species reps ({miss.GENUS.nunique()} genera); "
          f"present genomes qcov median {X[X.TPS > 0].qcov.median():.4f}, absent {miss.qcov.median():.4f}")
    print(f"TPS absent AND SwissProt-anchored absent: {int(((X.TPS==0)&(X.n_tps.fillna(0)==0)).sum())} genomes")
    print(f"  of those, credible (qcov>=0.985 & BUSCO>=85): {int(X.credible_TPS_loss.sum())} genomes, "
          f"{X.loc[X.credible_TPS_loss,'GENUS'].nunique()} genera")
    cl = X[X.credible_TPS_loss].groupby(["PHYLUM", "CLASS"]).GENUS.nunique().sort_values(ascending=False)
    print(cl.head(10).to_string())
    X[X.credible_TPS_loss].groupby("GENUS").agg(n=("genome", "size"), qcov=("qcov", "median"),
        busco=("complete_pct", "median"), n_tpp=("n_tpp", "median"),
        cls=("CLASS", "first")).sort_values("n", ascending=False).to_csv(DATA / "th_credible_TPS_loss.tsv", sep="\t")

    G.to_csv(DATA / "th_genus_medians.tsv", sep="\t", index=False)
    print(f"\n== 2. genus medians: {len(G)} genera from {int(X.keep.sum())} quality-kept genomes ==")

    # kingdom summary by class
    rows = []
    for m, (acc, desc) in MODULES.items():
        s = X[X.keep][m]
        rows.append((m, acc, desc, float((s > 0).mean()), float(s.median()),
                     float((s * X.loc[s.index, "k_prot"]).median()), float(s.max())))
    S = pd.DataFrame(rows, columns=["module", "pfam", "description", "frac_genomes_present",
                                    "median_copies", "median_per_1000_prot", "max_copies"])
    S.to_csv(DATA / "th_module_summary.tsv", sep="\t", index=False)
    print(S.sort_values("frac_genomes_present").head(12).to_string(index=False))

    print("\n== 3. trehalose module vs genome reduction (genus level, Spearman) ==")
    for m in ["TPS", "TPP", "NTH", "ATH", "GSY", "GPH", "sHSP", "HSP12", "HYDRO_I", "DIT1"]:
        r, p = stats.spearmanr(G[m + "_k"], G.n_prot)
        r2, p2 = stats.spearmanr(G[m], G.n_prot)
        print(f"  {m:8s} per-1000 rho={r:+.3f} (p={p:.1e})   raw copies rho={r2:+.3f} (p={p2:.1e})")
    print("  within-class (Spearman, pooled residual-free) Saccharomycotina and Mucoromycota:")
    for grp, sel in [("Saccharomycotina", G.yeast_origin == "Saccharomycotina"),
                     ("Mucoromycota", G.PHYLUM.isin(["Mucoromycota", "Kickxellomycota", "Zoopagomycota"]))]:
        s = G[sel]
        if len(s) < 5:
            continue
        print(f"   {grp} n={len(s)}: " + "  ".join(
            f"{m} rho={stats.spearmanr(s[m+'_k'], s.n_prot)[0]:+.2f}" for m in ["TPS", "TPP", "NTH", "GSY"]))

    print("\n== 4. module co-occurrence across genera (class-demeaned Spearman) ==")
    cols = ["TPS", "TPP", "NTH", "ATH", "GSY", "GPH", "GDB", "GH13", "UGP", "PGM",
            "GPD", "MIP", "MPD", "sHSP", "HSP12", "CAT", "SOD_Cu", "HYDRO_I", "HYDRO_II", "DIT1"]
    Dm = G[["GENUS", "CLASS"] + [c + "_k" for c in cols]].dropna()
    DM = Dm.groupby("CLASS")[ [c + "_k" for c in cols] ].transform(lambda s: s - s.median())
    C = DM.corr(method="spearman")
    C.to_csv(DATA / "th_module_correlation.tsv", sep="\t")
    pairs = [(a, b, C.loc[a + "_k", b + "_k"], int((Dm[[a + "_k", b + "_k"]] != 0).all(1).sum()))
             for a in cols for b in cols if a < b]
    P = pd.DataFrame(pairs, columns=["module_a", "module_b", "rho_class_demeaned", "n_genera_both"])
    P["p"] = [stats.spearmanr(DM[a + "_k"], DM[b + "_k"])[1] for a, b, *_ in P.itertuples(index=False)]
    P = P.sort_values("rho_class_demeaned", ascending=False)
    P.to_csv(DATA / "th_module_pairs.tsv", sep="\t", index=False)
    print("top positive:\n" + P.head(12).to_string(index=False))
    print("most negative:\n" + P.tail(6).to_string(index=False))

    print("\n== 5. TPS-TPP-NTH-ATH pathway states ==")
    st = X[X.keep].assign(state=lambda d: np.select(
        [(d.TPS > 0) & (d.TPP > 0) & (d.ATH > 0), (d.TPS > 0) & (d.TPP > 0) & (d.ATH == 0)],
        ["TPS+TPP+ATH", "TPS+TPP, no ATH"], default="other/incomplete"))
    print(pd.crosstab(st.CLASS, st.state, normalize="index").round(3)[
        ["TPS+TPP+ATH", "TPS+TPP, no ATH", "other/incomplete"]].sort_values(
        "other/incomplete", ascending=False).head(12).to_string())

    print("\n== 6. dessication case-controls: focal genera vs their own class ==")
    out = []
    for label, focal_genera in [
            ("xerophile_genus_set", sorted(set(G.loc[G.GENUS.isin(XERO_GENERA), "GENUS"]))),
            ("dry_epithet_genus_set", sorted(set(G.loc[G.GENUS.isin(
                X.loc[X.dry_epithet, "GENUS"].dropna()), "GENUS"])))]:
        for gen in focal_genera:
            row = G[G.GENUS == gen].iloc[0]
            sib = G[(G.CLASS == row.CLASS) & (G.GENUS != gen)]
            if len(sib) < 5:
                sib = G[(G["ORDER_"] == row["ORDER_"]) & (G.GENUS != gen)]
            if len(sib) < 5:
                continue
            for m in ["TPS", "TPP", "NTH", "ATH", "GSY", "GPH", "sHSP", "HSP12", "HYDRO_I", "CAT"]:
                v, lo, hi = boot_diff([row[m + "_k"]], sib[m + "_k"])
                out.append(dict(set=label, GENUS=gen, CLASS=row.CLASS, metric=m,
                                focal=row[m + "_k"], n_sisters=len(sib),
                                sib_median=sib[m + "_k"].median(), diff=v, ci_lo=lo, ci_hi=hi,
                                # percentile rank of the focal genus among its sisters
                                pct=float((sib[m + "_k"] <= row[m + "_k"]).mean())))
    CC = pd.DataFrame(out)
    CC.to_csv(DATA / "th_xerophile_casecontrols.tsv", sep="\t", index=False)
    piv = CC.pivot_table(index=["set", "GENUS", "CLASS"], columns="metric", values="pct")
    print(piv[["TPS", "TPP", "NTH", "ATH", "GSY", "GPH", "sHSP", "HSP12", "HYDRO_I", "CAT"]].round(2).to_string())

    print("\n== 7. lichen-forming (Lecanoromycetes) vs free-living Dothideomycetes/Eurotiomycetes ==")
    lic = G[G.CLASS == "Lecanoromycetes"]
    non = G[G.CLASS.isin(["Dothideomycetes", "Eurotiomycetes"])]
    for m in ["TPS", "TPP", "NTH", "ATH", "sHSP", "HSP12", "HYDRO_I", "DIT1", "n_prot", "mb"]:
        if m in ("n_prot", "mb"):
            a, b = lic[m], non[m]
        else:
            a, b = lic[m + "_k"], non[m + "_k"]
        d, lo, hi = boot_diff(a, b)
        print(f"  {m:9s} lichen med {np.median(a):8.2f} vs other {np.median(b):8.2f}  diff {d:+8.2f} "
              f"[{lo:+.2f},{hi:+.2f}] delta {cliffs_delta(a,b):+.2f} p={stats.mannwhitneyu(a,b).pvalue:.1e}")

    print("\n== 8. Kickxellomycota and their relatives ==")
    kx = G[G.PHYLUM.isin(["Kickxellomycota", "Mucoromycota", "Zoopagomycota"])].copy()
    print(kx[["GENUS", "CLASS", "n_sp", "n_prot", "mb", "TPS_k", "TPP_k", "NTH_k", "ATH_k",
              "GSY_k", "GPH_k", "HSP12_k", "sHSP_k", "f_tps", "f_tpp", "f_ath"]].sort_values(
        ["CLASS", "GENUS"]).to_string(index=False))
    kx.to_csv(DATA / "th_genus_mucoromycota_kickx.tsv", sep="\t", index=False)

    make_figs(X, G, C, cols)


def make_figs(X, G, C, cols):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    present, absent = X[X.TPS > 0], X[X.TPS == 0]
    ax[0].scatter(100 * (1 - present.qcov), present.TPS, s=4, alpha=.25, label="TPS present")
    ax[0].scatter(100 * (1 - absent.qcov), np.full(len(absent), .8), s=16, c="crimson",
                  label="TPS absent")
    ax[0].set_xscale("symlog", linthresh=0.05)
    ax[0].set_xlabel("% of the 1,409 near-universal Pfam families NOT found (annotation gaps)")
    ax[0].set_yscale("log")
    ax[0].set_ylabel("TPS copies per genome")
    ax[0].set_title("TPS absence vs annotation completeness")
    ax[0].legend(fontsize=8)
    ax[1].scatter(G.n_prot, G.TPS_k, s=10)
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("proteome size (genus median)")
    ax[1].set_ylabel("TPS per 1,000 proteins")
    ax[1].set_title("TPS copy number vs genome reduction")
    fig.tight_layout()
    fig.savefig(FIGS / "f1_qc_and_allometry.png", dpi=150)

    fig, ax = plt.subplots(figsize=(8, 6.5))
    im = ax.imshow(C.values, vmin=-.4, vmax=.4, cmap="RdBu_r")
    ax.set_xticks(range(len(C.columns)), [c.replace("_k", "") for c in C.columns], rotation=90, fontsize=8)
    ax.set_yticks(range(len(C.index)), [c.replace("_k", "") for c in C.index], fontsize=8)
    ax.set_title("Class-demeaned Spearman correlation of module copy number (per 1,000 proteins)")
    fig.colorbar(im, shrink=.8)
    fig.tight_layout()
    fig.savefig(FIGS / "f2_module_correlation.png", dpi=150)

    cls = G.groupby("CLASS").agg(n=("GENUS", "size"), TPS=("TPS_k", "median"), TPP=("TPP_k", "median"),
                                 NTH=("NTH_k", "median"), ATH=("ATH_k", "median")).query("n >= 8")
    cls = cls.sort_values("TPS", ascending=False)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    (cls[["TPS", "TPP", "NTH", "ATH"]] / 10).plot.bar(ax=ax, width=.8)
    ax.set_ylabel("copies per 100 proteins")
    ax.set_title("Trehalose module by class (genus medians, n>=8 genera)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "f3_class_module_medians.png", dpi=150)


if __name__ == "__main__":
    main()
