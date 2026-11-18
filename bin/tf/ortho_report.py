"""Ortholog report: presence, copy number and architecture by phylum and order; tables and figures.

    python bin/tf/ortho_report.py report --stage DIR [--summary-dir DIR]

Writes, through the Out class (level 'tf') under DIR/tables/tf and DIR/figures/tf:
    T20_ortho_<slug>.tsv            slug wc2, zap1, bzip5; first line '# thresholds'; rows for phyla and orders
    T20_ortho_<slug>_presence       figure: orders by prevalence, with median copies
    T20_ortho_<slug>_tree           figure: gene tree (when the tree exists)
    T20_ortho_seeds.tsv             seed list, source, status, best BFD match and its architecture, sanity-check warning
    T20_ortho_thresholds.tsv        thresholds, rules, counts per target and tree settings
A genus representative without an ortholog that passes the thresholds and the architecture rule counts as absent.
Similarity and tree position are evidence for orthology, not proof. No function is inferred.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "bin" / "summary"))
from ortho_search import SLUG, write_commented  # noqa: E402

RUNS = Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs")
SUMMARY_COLS = ["target", "lineage", "parent", "n_genera", "n_with_ortholog", "prevalence", "median_copies", "architecture_counts"]
PHYLA_SHOWN = ["Mucoromycota", "Mortierellomycota", "Glomeromycota", "Chytridiomycota", "Ascomycota", "Basidiomycota"]


def _arch_counts(h, top=5):
    if not len(h):
        return ""
    vc = h["architecture"].value_counts()
    return "; ".join(f"{a}:{n}" for a, n in vc.head(top).items())


def summarise(hits: pd.DataFrame, genus_reps: pd.DataFrame, target=None) -> pd.DataFrame:
    """Presence by phylum and order. `hits`: ortholog proteins (species_prefix, protein_id, architecture) that passed the
    thresholds and the architecture rule. `genus_reps`: one row per genus representative (species_prefix, PHYLUM, ORDER).
    Columns: target, lineage, parent (the phylum of an order; empty for a phylum), n_genera (all representatives of the
    lineage), n_with_ortholog, prevalence (n_with/n_genera), median_copies (median number of ortholog proteins among the
    representatives that have one; NaN when none), architecture_counts (top architectures among the ortholog proteins)."""
    target = target if target is not None else (hits["target"].iloc[0] if "target" in hits.columns and len(hits) else "")
    copies = hits.groupby("species_prefix").size() if len(hits) else pd.Series(dtype=int)
    reps = genus_reps.copy()
    reps["PHYLUM"] = reps["PHYLUM"].fillna("Unknown")
    reps["ORDER"] = reps["ORDER"].fillna("Unknown")
    reps["copies"] = reps["species_prefix"].map(copies).fillna(0).astype(int)
    rows = []

    def one(lineage, parent, g):
        have = g[g["copies"] > 0]
        hh = hits[hits["species_prefix"].isin(set(g["species_prefix"]))] if len(hits) else hits
        rows.append({"target": target, "lineage": lineage, "parent": parent, "n_genera": len(g), "n_with_ortholog": len(have),
                     "prevalence": len(have) / len(g) if len(g) else np.nan,
                     "median_copies": float(have["copies"].median()) if len(have) else np.nan,
                     "architecture_counts": _arch_counts(hh)})
    for ph, g in reps.groupby("PHYLUM"):
        one(ph, "", g)
    for (ph, od), g in reps.groupby(["PHYLUM", "ORDER"]):
        one(od, ph, g)
    return pd.DataFrame(rows, columns=SUMMARY_COLS)


# ---------------------------------------------------------------- figures
def fig_presence(summ: pd.DataFrame, colors: dict, porder, target, min_genera=3):
    """Orders (at least `min_genera` representatives) grouped by phylum: prevalence and median copies."""
    import matplotlib.pyplot as plt
    od = summ[summ["parent"] != ""].copy()
    od = od[od["n_genera"] >= min_genera]
    if not len(od):
        return None
    rank = {p: i for i, p in enumerate(porder)}
    od["_r"] = od["parent"].map(lambda p: rank.get(p, len(rank)))
    od = od.sort_values(["_r", "parent", "prevalence", "lineage"], ascending=[True, True, False, True])
    n = len(od)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.4, min(max(3.0, n * 0.17 + 1.2), 40)), sharey=True,
                                 gridspec_kw={"width_ratios": [3, 1.2], "wspace": 0.05})
    y = np.arange(n)
    a1.imshow(od[["prevalence"]].to_numpy(dtype=float), aspect="auto", cmap="viridis", vmin=0, vmax=1)
    for i, r in enumerate(od.itertuples()):
        v = r.prevalence
        a1.text(0, i, f"{r.n_with_ortholog}/{r.n_genera}", ha="center", va="center", fontsize=6,
                color="white" if v < 0.6 else "black")
    a1.set_xticks([0])
    a1.set_xticklabels(["genera with an ortholog"], fontsize=7)
    a1.set_yticks(y)
    a1.set_yticklabels([f"{r.lineage}" for r in od.itertuples()], fontsize=6)
    for tl, p in zip(a1.get_yticklabels(), od["parent"]):
        tl.set_color(colors.get(p, "#555555"))
    mc = od["median_copies"].to_numpy(dtype=float)
    a2.imshow(np.nan_to_num(mc, nan=0.0)[:, None], aspect="auto", cmap="Greys", vmin=0, vmax=max(2.0, np.nanmax(mc) if np.isfinite(mc).any() else 2.0))
    for i, v in enumerate(mc):
        a2.text(0, i, "" if np.isnan(v) else f"{v:g}", ha="center", va="center", fontsize=6)
    a2.set_xticks([0])
    a2.set_xticklabels(["median copies"], fontsize=7)
    a2.tick_params(left=False)
    for ax in (a1, a2):
        ax.spines[:].set_visible(False)
    from matplotlib.lines import Line2D
    ph = [p for p in porder if p in set(od["parent"]) and p in colors]
    a1.legend(handles=[Line2D([], [], marker="s", ls="", color=colors[p], label=p, markersize=5) for p in ph],
              fontsize=6, loc="upper left", bbox_to_anchor=(0, -0.02 if n > 20 else -0.15), frameon=False, ncol=2)
    a1.set_title(f"{target}: orders with >= {min_genera} genus representatives", fontsize=8)
    return fig


# ---------------------------------------------------------------- pipeline
def phylum_colors(summary_dir):
    from make_figures import phylum_color_map, phylum_order
    from suppress import drop_suppressed, load_suppress
    ga = pd.read_parquet(Path(summary_dir) / "genome_all.parquet")
    ga, _ = drop_suppressed(ga, load_suppress(None))
    order = phylum_order(ga)
    return order, phylum_color_map(order)


def build_seed_table(stage: Path, counts: pd.DataFrame):
    s = pd.read_csv(stage / "seeds.tsv", sep="\t", keep_default_na=False)
    s["length"] = s["sequence"].str.len()
    thr = counts.set_index("target")
    s["min_pident"] = s["target"].map(thr["min_pident"])
    s["min_cov"] = s["target"].map(thr["min_cov"])
    s["architecture_rule"] = s["target"].map(thr["architecture_rule"])
    cols = ["target", "label", "source", "key", "seed_id", "status", "length", "bfd_match", "bfd_pident", "bfd_min_cov",
            "bfd_architecture_names", "warning", "min_pident", "min_cov", "architecture_rule"]
    out = s[cols].rename(columns={"label": "seed", "seed_id": "accession"})
    return out.rename(columns={"bfd_match": "bfd_match_FBF4BF29", "bfd_min_cov": "bfd_coverage",
                               "bfd_architecture_names": "bfd_architecture"})


def step_report(stage: Path, summary_dir: Path, out_root=None):
    from figs_common import safe
    from make_figures import Out
    import ortho_tree as ort
    out_root = Path(out_root or stage)
    o = Out(out_root, "tf")
    o.scope = "genus_reps"
    porder, colors = phylum_colors(summary_dir)
    reps = pd.read_csv(stage / "genus_reps.tsv", sep="\t")
    counts = pd.read_csv(stage / "counts.tsv", sep="\t")
    summaries = {}
    for target, slug in SLUG.items():
        ph = pd.read_csv(stage / f"pass_{slug}.tsv", sep="\t", comment="#")
        summ = summarise(ph, reps, target=target)
        summaries[slug] = summ
        c = counts[counts["target"] == target].iloc[0]
        params = {"min_pident": c["min_pident"], "min_cov": c["min_cov"], "cov_mode": 0,
                  "architecture_rule": c["architecture_rule"], "n_genera": c["n_reps"]}
        write_commented(summ, o.tab / f"T20_ortho_{slug}.tsv", params)
        o.made.append(_mk(f"T20_ortho_{slug}", o))
        ng = int(c["n_arch_pass_genera"])
        o.info = {"n_genomes": int(c["n_reps"]), "n_groups": int((summ["parent"] != "").sum()), "pooled_share": 0.0}
        o.meta[f"T20_ortho_{slug}_presence"] = {"signal": ng, "signal_kind": "n_genera_with_ortholog"}
        safe(o, f"T20_ortho_{slug}_presence", "tf", lambda s=summ, t=target: fig_presence(s, colors, porder, t))
        nwk = stage / f"tree_{slug}.nwk"
        tips = stage / f"tree_{slug}_tips.tsv"
        if nwk.exists() and tips.exists():
            meta = pd.read_csv(tips, sep="\t")
            o.meta[f"T20_ortho_{slug}_tree"] = {"signal": len(meta), "signal_kind": "n_tree_tips", "n": len(meta)}
            safe(o, f"T20_ortho_{slug}_tree", "tf",
                 lambda n=nwk, m=meta, t=target: ort.tree_figure(n, m, colors, title=f"{t}: {len(m)} tips (best hit per genus plus seeds)"))
        else:
            print(f"tree for {slug} not found: figure skipped", flush=True)
            o.skip(f"T20_ortho_{slug}_tree", "tf")
    seeds = build_seed_table(stage, counts)
    seeds.to_csv(o.tab / "T20_ortho_seeds.tsv", sep="\t", index=False)
    thr = counts.copy()
    thr["tree"] = "mafft --auto; columns with >50% gaps removed; FastTree -lg; one best hit per genus plus seeds"
    thr["search"] = "MMseqs2 easy-search seeds vs candidates; identity and coverage as given; coverage of both sequences"
    thr["note"] = "Candidates carry the target Pfam only. The architecture rule is applied after the search. Evidence for orthology, not proof."
    thr.to_csv(o.tab / "T20_ortho_thresholds.tsv", sep="\t", index=False)
    return o.made


def _mk(name, o):
    from make_figures import Made
    return Made(name, "tables", "table", o.scope)


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["report"])
    ap.add_argument("--stage", required=True)
    ap.add_argument("--summary-dir", default=str(RUNS / "analysis/summary_tables"))
    ap.add_argument("--out", default=None, help="output root (default: the stage dir)")
    a = ap.parse_args(argv)
    made = step_report(Path(a.stage), Path(a.summary_dir), a.out)
    for m in made:
        print(f"{m.kind:8s} {m}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
