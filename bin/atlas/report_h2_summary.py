#!/usr/bin/env python3
"""Summary figures for the H2 report (yeast vs filamentous protein length).

Reads only the small TSV tables that render_views.py / render_p2.py wrote to
docs/summary_atlas/tables/, so it runs in seconds on any node.

Writes to docs/summary_atlas/figures/:
  H2R_forest.{png,svg}       percent length difference by proxy, filter set, yeast origin
  H2R_family_scatter.{png,svg} yeast vs filamentous median length per P-A family

Percent difference = (2 ** median_log2_ratio - 1) * 100. Negative = yeast shorter.

Usage:
    atlas_env/venv/bin/python bin/atlas/report_h2_summary.py --atlas docs/summary_atlas
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

YEAST = "#c2410c"
FILA = "#1d4ed8"
GREY = "#6b7280"


def pct(x):
    return (np.power(2.0, np.asarray(x, dtype=float)) - 1.0) * 100.0


def forest(ax, labels, est, lo, hi, ns, title, colors=None, notes=None):
    y = np.arange(len(labels))[::-1]
    colors = colors or [YEAST] * len(labels)
    for yi, e, l, h, c in zip(y, est, lo, hi, colors):
        if np.isfinite(l) and np.isfinite(h):
            ax.plot([l, h], [yi, yi], color=c, lw=2)
        ax.plot(e, yi, "o", color=c, ms=6)
    ax.axvline(0, color=GREY, lw=1, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{lab}  (n={n:,})" for lab, n in zip(labels, ns)], fontsize=8)
    ax.set_title(title, fontsize=10, loc="left")
    ax.set_xlabel("yeast vs filamentous length (%)", fontsize=8)
    ax.tick_params(axis="x", labelsize=8)
    if notes:
        for yi, note in zip(y, notes):
            if note:
                ax.annotate(note, (1.0, yi), xycoords=("axes fraction", "data"),
                            xytext=(4, 0), textcoords="offset points",
                            fontsize=7, color=GREY, va="center")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", type=Path, required=True)
    a = ap.parse_args()
    t, f = a.atlas / "tables", a.atlas / "figures"

    # ---- panel A: proxies ---------------------------------------------------
    v10 = pd.read_csv(t / "V10_anchor_contrast.tsv", sep="\t")
    v11 = pd.read_csv(t / "V11_arch_summary.tsv", sep="\t") if (t / "V11_arch_summary.tsv").exists() else None
    rows = []
    pa = v10[v10.proxy.str.startswith("P-A")].iloc[0]
    rows.append(("P-A shared SwissProt best hit", pa))
    names = {"YEAST": "S. cerevisiae", "NEUCR": "N. crassa", "SCHPO": "S. pombe", "EMENI": "A. nidulans"}
    for code, nm in names.items():
        r = v10[v10.proxy == f"P-B anchor {code}"].iloc[0]
        rows.append((f"P-B anchor: {nm}", r))
    lab = [r[0] for r in rows]
    est = pct([r[1].median_log2_ratio for r in rows])
    lo = pct([r[1].ci_lo for r in rows])
    hi = pct([r[1].ci_hi for r in rows])
    ns = [int(r[1].n_families) for r in rows]
    # P-C numbers are taken from the V11 table on the H2 page (n=5,026; log2 -0.029, CI -0.032 to -0.026)
    if v11 is not None:
        r = v11.iloc[0]
        lab.append("P-C Pfam architecture")
        est = np.append(est, pct(r.median_log2_ratio)); lo = np.append(lo, pct(r.ci_lo)); hi = np.append(hi, pct(r.ci_hi))
        ns.append(int(r.n_families))
    else:
        lab.append("P-C Pfam architecture")
        est = np.append(est, pct(-0.029)); lo = np.append(lo, pct(-0.032)); hi = np.append(hi, pct(-0.026))
        ns.append(5026)

    # ---- panel B: filter grid -------------------------------------------------
    v12 = pd.read_csv(t / "V12_filter_grid.tsv", sep="\t")
    lab_b = [s.split(" ", 1)[0] + " " + s.split(" ", 1)[1][:38] for s in v12.filter_set]

    # ---- panel C: yeast origin -------------------------------------------------
    v09 = pd.read_csv(t / "V09_offset_by_origin.tsv", sep="\t")
    v09 = v09[v09.yeast_origin != "filamentous"]
    notes_c = ["only 2 genera" if n < 3 else "" for n in v09.n_genera]
    col_c = [YEAST if "dimorphic" not in o else GREY for o in v09.yeast_origin]

    fig, axes = plt.subplots(3, 1, figsize=(8.5, 10.5),
                             gridspec_kw={"height_ratios": [len(lab), len(v12), len(v09)]})
    forest(axes[0], lab, est, lo, hi, ns, "A. Within-family length, by ortholog proxy (n = families)")
    forest(axes[1], lab_b, pct(v12.median_log2_ratio), pct(v12.ci_lo), pct(v12.ci_hi),
           v12.n_families.astype(int), "B. P-A result under 7 quality-filter sets (n = families)")
    forest(axes[2], list(v09.yeast_origin), pct(v09.median_offset_log2), pct(v09.ci_lo),
           pct(v09.ci_hi), v09.n_genera.astype(int),
           "C. Genome offset by yeast origin (n = genera; grey = dimorphic)",
           colors=col_c, notes=notes_c)
    lim = max(abs(np.nanmin(pct(v09.ci_lo))), abs(np.nanmax(pct(v09.ci_hi))), 4.5)
    for ax in axes:
        ax.set_xlim(-lim - 0.5, lim + 0.5)
    fig.suptitle("H2: within ortholog families, yeast proteins are 0.6-2 % shorter; only Saccharomycotina",
                 fontsize=11)
    fig.text(0.01, 0.005, "Points = median; bars = bootstrap 95 % CI. Negative = yeast shorter. "
             "Percent = (2^median log2 ratio - 1) x 100.", fontsize=7, color=GREY)
    fig.tight_layout(rect=(0, 0.015, 1, 0.97))
    for ext in ("png", "svg"):
        fig.savefig(f / f"H2R_forest.{ext}", dpi=150)
    plt.close(fig)

    # ---- family scatter: absolute lengths --------------------------------------
    v08 = pd.read_csv(t / "V08_family_ratios.tsv", sep="\t", comment="#")
    x = v08["median_filamentous_mycelium"].astype(float)
    yv = v08["median_yeast"].astype(float)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.6))
    ax1.scatter(x, yv, s=3, alpha=0.35, color=YEAST, rasterized=True, lw=0)
    m = max(x.max(), yv.max())
    ax1.plot([30, m], [30, m], color=GREY, lw=1, ls="--", label="equal length")
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("filamentous median length (aa)"); ax1.set_ylabel("yeast median length (aa)")
    ax1.set_title(f"A. {len(v08):,} P-A families: yeast vs filamentous", fontsize=10, loc="left")
    ax1.legend(fontsize=8, loc="upper left")
    diff = yv - x
    bins = np.arange(-150, 151, 5)
    ax2.hist(diff.clip(-150, 150), bins=bins, color=YEAST, alpha=0.8)
    ax2.axvline(0, color=GREY, lw=1, ls="--")
    med = float(np.median(diff))
    ax2.axvline(med, color="black", lw=1.2)
    ax2.set_xlabel("yeast minus filamentous median length (aa); clipped at +/-150")
    ax2.set_ylabel("families")
    ax2.set_title(f"B. Difference per family; median {med:+.1f} aa", fontsize=10, loc="left")
    fig.text(0.01, 0.005, "Point layer rasterized in the SVG. One point = one SwissProt best-hit family "
             "(>= 10 genera per side, strict gates).", fontsize=7, color=GREY)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    for ext in ("png", "svg"):
        fig.savefig(f / f"H2R_family_scatter.{ext}", dpi=150)
    plt.close(fig)

    # numbers the report quotes
    print(f"P-A families: {len(v08)}; median aa difference {med:+.2f}; "
          f"median filamentous {np.median(x):.1f}, yeast {np.median(yv):.1f}; "
          f"yeast shorter {(diff < 0).sum()}, longer {(diff > 0).sum()}, tied {(diff == 0).sum()}")
    for l, e, n in zip(lab, est, ns):
        print(f"{l}: {e:+.2f} % (n={n})")


if __name__ == "__main__":
    main()
