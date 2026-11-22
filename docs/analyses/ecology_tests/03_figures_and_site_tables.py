#!/usr/bin/env python3
"""Figures ET1-ET3 and the slim tables for the Ecology page, from the result tables of 02_tests.py.

Usage: 03_figures_and_site_tables.py RESULTS_DIR OUT_DIR
Writes OUT_DIR/figures/ecology/ET*.png|svg|json and OUT_DIR/tables/ecology/*.tsv. Reads only RESULTS_DIR (no database).
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

LABEL = {"total_length_bp": "genome length", "n_genes": "gene count", "genes_per_mb": "genes per Mb", "coding_fraction": "coding fraction",
         "gc_pct": "GC %", "median_intron_len_bp": "median intron length", "median_exons_per_gene": "median exons per gene",
         "intron_frac_genes_with_intron": "genes with an intron", "median_protein_len_aa": "median protein length",
         "cazy_GH_per1000": "CAZy GH", "cazy_GT_per1000": "CAZy GT", "cazy_PL_per1000": "CAZy PL", "cazy_CE_per1000": "CAZy CE",
         "cazy_AA_per1000": "CAZy AA", "cazy_CBM_per1000": "CAZy CBM", "secretome_per1000": "strict secretome",
         "merops_per1000": "MEROPS proteases", "tf_per1000": "transcription factors"}
COL = {"taxonomy": "#4c78a8", "shared": "#b8b8b8", "ecology": "#e45756", "residual": "#f0f0f0"}


def save(fig, out, name, meta):
    d = out / "figures" / "ecology"
    d.mkdir(parents=True, exist_ok=True)
    fig.savefig(d / f"{name}.png", dpi=150, bbox_inches="tight")
    fig.savefig(d / f"{name}.svg", bbox_inches="tight")
    (d / f"{name}.json").write_text(json.dumps(meta, indent=1))
    plt.close(fig)


def et1(res, out):
    a = pd.read_csv(res / "A_variance_partition.tsv", sep="\t")
    a = a[a.model == "order + coarse lifestyle"].copy()
    a["tax_only"] = a.unique_taxonomy.clip(lower=0)
    a["eco_only"] = a.unique_ecology.clip(lower=0)
    a["shared_c"] = a.shared.clip(lower=0)
    a["resid"] = (1 - (a.tax_only + a.eco_only + a.shared_c)).clip(lower=0)
    a = a.sort_values("eco_only")
    fig, ax = plt.subplots(figsize=(8, 6.2))
    y = np.arange(len(a))
    left = np.zeros(len(a))
    for col, key, lab in (("tax_only", "taxonomy", "taxonomy only (order)"), ("shared_c", "shared", "shared taxonomy and ecology"),
                          ("eco_only", "ecology", "ecology only, after order"), ("resid", "residual", "unexplained")):
        ax.barh(y, a[col], left=left, color=COL[key], label=lab, edgecolor="white", linewidth=0.5)
        left += a[col].to_numpy()
    for i, (q, ok) in enumerate(zip(a.q_BH, a.ecology_explains)):
        if ok:
            ax.text(1.01, i, "*", va="center", ha="left", fontsize=14, color=COL["ecology"])
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[r] for r in a.response])
    ax.set_xlim(0, 1.04)
    ax.set_xlabel("share of variance (adjusted R2)")
    ax.set_title("What explains genome properties and protein content: order vs coarse lifestyle")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.28), ncol=2, frameon=False, fontsize=8)
    fig.text(0.01, -0.02, "* ecology after order: unique part >= 0.02 and q < 0.05 (permutation within order)", fontsize=8)
    save(fig, out, "ET1_variance_partition", {"n_genomes": int(a.n.max()), "n_groups": 6, "kind": "variance partition"})


def et2(res, out):
    e = pd.read_csv(res / "B_per_clade_effects_order.tsv", sep="\t")
    e = e[e.contrast == "mycorrhizal vs saprotroph"]
    b = pd.read_csv(res / "B_matched_contrasts_order.tsv", sep="\t")
    b = b[b.contrast == "mycorrhizal vs saprotroph"].set_index("response")
    resp = ["cazy_GH_per1000", "cazy_AA_per1000", "secretome_per1000", "merops_per1000"]
    orders = sorted(e.clade.unique())
    fig, axes = plt.subplots(1, 4, figsize=(11, 3.4), sharey=True)
    for ax, r in zip(axes, resp):
        d = e[e.response == r].set_index("clade").reindex(orders)
        ax.axvline(0, color="#777", lw=0.8)
        ax.scatter(d.effect, range(len(orders)), s=50, color=COL["ecology"], zorder=3)
        ax.axvline(b.loc[r, "pooled_effect"], color="#333", lw=1.2, ls="--")
        ax.set_title(LABEL[r], fontsize=10)
        ax.set_xlabel("per 1,000 genes")
        ax.grid(axis="y", color="#eee")
    axes[0].set_yticks(range(len(orders)))
    axes[0].set_yticklabels(orders)
    fig.suptitle("Mycorrhizal minus saprotrophic genera, within each order (dashed line: pooled effect)", fontsize=10, y=1.04)
    save(fig, out, "ET2_mycorrhizal_vs_saprotroph", {"n_genomes": int(b.n_a.max() + b.n_b.max()), "n_groups": 2, "orders": len(orders)})


def et3(res, out):
    c = pd.read_csv(res / "C_domain_scan_order.tsv", sep="\t")
    c = c[(c.contrast == "mycorrhizal vs saprotroph") & c.feature.str.startswith("cazy:") & (c.q_BH < 0.05) & (c.share_clades_same_sign >= 0.75)]
    c = c.assign(a=c.pooled_effect_log10.abs()).nlargest(20, "a").sort_values("pooled_effect_log10")   # the 20 largest effects
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.barh(range(len(c)), c.pooled_effect_log10, color=[COL["ecology"] if v < 0 else COL["taxonomy"] for v in c.pooled_effect_log10])
    ax.set_yticks(range(len(c)))
    ax.set_yticklabels([f.split(":")[1] for f in c.feature], fontsize=8)
    ax.axvline(0, color="#777", lw=0.8)
    ax.set_xlabel("mycorrhizal minus saprotroph: difference of log10(1 + copies per 1,000 genes)")
    ax.set_title("The 20 CAZy families with the largest difference (q < 0.05, same sign in at least 75% of the 5 orders)", fontsize=8)
    save(fig, out, "ET3_mycorrhizal_cazy_families", {"n_genomes": 214, "n_groups": 2, "n_families_shown": int(len(c))})


def site_tables(res, out):
    d = out / "tables" / "ecology"
    d.mkdir(parents=True, exist_ok=True)
    a = pd.read_csv(res / "A_variance_partition.tsv", sep="\t")
    a = a[a.model.isin(["order + coarse lifestyle", "family + coarse lifestyle"])]
    a = a[["model", "response", "kind", "n", "R2adj_taxonomy", "R2adj_ecology", "unique_ecology", "shared", "q_BH", "ecology_explains"]]
    a["response"] = a.response.map(LABEL)
    a.to_csv(d / "ecology_variance_partition.tsv", sep="\t", index=False)
    p = pd.read_csv(res / "A_permanova_profiles.tsv", sep="\t")
    p[["profile", "n_features", "trait", "test", "n_genera", "n_groups", "R2_adj", "p"]].to_csv(d / "ecology_permanova_profiles.tsv", sep="\t", index=False)
    rows = []
    for lvl in ("family", "order"):
        b = pd.read_csv(res / f"B_matched_contrasts_{lvl}.tsv", sep="\t")
        for ct, g in b.groupby("contrast", sort=False):
            vc = g.verdict.value_counts()
            rows.append({"clade_level": lvl, "contrast": ct, "n_clades": int(g.n_clades.max()), "genera_A": int(g.n_a.max()), "genera_B": int(g.n_b.max()),
                         "supported": int(vc.get("supported", 0)), "mixed": int(vc.get("mixed", 0)), "no_evidence": int(vc.get("no evidence", 0)),
                         "not_evaluable": int(vc.get("not evaluable", 0))})
    pd.DataFrame(rows).to_csv(d / "ecology_matched_summary.tsv", sep="\t", index=False)
    hits = []
    for lvl in ("family", "order"):
        b = pd.read_csv(res / f"B_matched_contrasts_{lvl}.tsv", sep="\t")
        b = b[b.verdict.isin(["supported", "mixed"])].copy()
        b["response"] = b.response.map(LABEL)
        hits.append(b[["clade_level", "contrast", "response", "n_clades", "pooled_effect", "q_BH", "share_clades_same_sign", "verdict"]])
    pd.concat(hits).to_csv(d / "ecology_matched_hits.tsv", sep="\t", index=False)
    s = []
    for lvl in ("family", "order"):
        c = pd.read_csv(res / f"C_domain_scan_{lvl}.tsv", sep="\t")
        for ct, g in c.groupby("contrast", sort=False):
            hit = g[g.q_BH < 0.05]
            s.append({"clade_level": lvl, "contrast": ct, "n_features_tested": len(g), "n_q_lt_0.05": len(hit),
                      "n_same_sign_ge_75pct": int((hit.share_clades_same_sign >= 0.75).sum())})
    pd.DataFrame(s).to_csv(d / "ecology_domain_scan_summary.tsv", sep="\t", index=False)


def main(res, out):
    res, out = Path(res), Path(out)
    et1(res, out)
    et2(res, out)
    et3(res, out)
    site_tables(res, out)
    print("done", sorted(p.name for p in (out / "figures" / "ecology").glob("*.png")))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
