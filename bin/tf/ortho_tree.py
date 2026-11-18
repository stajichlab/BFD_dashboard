"""Gene trees for the ortholog candidates: MAFFT, gap trimming, FastTree, and a Bio.Phylo/matplotlib drawing.

    build_tree(fasta, out_prefix, scratch) -> Path to the Newick file, or None (with a message) when mafft or FastTree is missing
    draw_tree(newick, meta, colors, out_png)   tip points coloured by phylum, anchor species marked
    python bin/tf/ortho_tree.py trees --stage DIR --scratch DIR   builds the three trees of a run (skips existing)

A tree is a picture of one alignment and one model (FastTree, LG). Position in it is evidence for orthology, not proof.
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "bin" / "summary"))

MAX_GAP = 0.5
MAX_TIPS = 2500
SLUG = {"WC-2": "wc2", "Zap1-like": "zap1", "bZip5": "bzip5"}


def read_alignment(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name is not None:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def trim_gappy(aln: dict, max_gap=MAX_GAP) -> dict:
    """Drop the alignment columns in which more than `max_gap` of the sequences have a gap."""
    names = list(aln)
    if not names:
        return {}
    arr = np.array([list(aln[n]) for n in names])
    keep = (arr == "-").mean(axis=0) <= max_gap
    return {n: "".join(arr[i, keep]) for i, n in enumerate(names)}


def build_tree(fasta: Path, out_prefix: Path, scratch: Path):
    """mafft --auto, trim columns with more than 50% gaps, FastTree -lg. Writes <out_prefix>.aln.fa (trimmed alignment)
    and <out_prefix>.nwk, returns the Newick path. Returns None and prints a message when mafft or FastTree is not on PATH."""
    mafft, ft = shutil.which("mafft"), shutil.which("FastTree") or shutil.which("fasttree")
    if mafft is None or ft is None:
        miss = " and ".join(n for n, e in (("mafft", mafft), ("FastTree", ft)) if e is None)
        print(f"{miss} not on PATH: tree skipped for {fasta} (module load mafft/7.505 fasttree/2.1.11)", flush=True)
        return None
    scratch = Path(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    out_prefix = Path(out_prefix)
    raw = scratch / (out_prefix.name + ".mafft.fa")
    with open(raw, "w") as o:
        subprocess.run([mafft, "--auto", "--thread", str(_threads()), "--quiet", str(fasta)], stdout=o, check=True)
    trimmed = trim_gappy(read_alignment(raw), MAX_GAP)
    aln = Path(str(out_prefix) + ".aln.fa")
    aln.write_text("".join(f">{k}\n{v}\n" for k, v in trimmed.items()))
    nwk = Path(str(out_prefix) + ".nwk")
    with open(nwk, "w") as o:
        subprocess.run([ft, "-lg", "-quiet", str(aln)], stdout=o, check=True)
    raw.unlink(missing_ok=True)
    return nwk


def _threads():
    import os
    return int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))


def tree_figure(newick: Path, meta: pd.DataFrame, colors: dict, title=None):
    """Figure of a rectangular tree. meta columns: tip, phylum, anchor (bool, optional), is_seed (bool, optional).
    Tips are points coloured by phylum (grey when the phylum is not in `colors`); anchor tips get a black ring;
    seeds get a star. The tree is midpoint-rooted for display."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D
    from Bio import Phylo
    tree = Phylo.read(str(newick), "newick")
    try:
        tree.root_at_midpoint()
    except Exception:
        pass
    tips = tree.get_terminals()
    ypos = {t: i for i, t in enumerate(tips)}
    xpos = tree.depths()
    if not any(xpos.values()):
        xpos = tree.depths(unit_branch_lengths=True)

    for c in reversed(list(tree.find_clades(order="level"))):      # children before parents, no recursion
        if not c.is_terminal():
            ys = [ypos[k] for k in c.clades]
            ypos[c] = (min(ys) + max(ys)) / 2
    segs = []
    for c in tree.find_clades(order="level"):
        for k in c.clades:
            segs.append([(xpos[c], ypos[c]), (xpos[c], ypos[k])])    # vertical
            segs.append([(xpos[c], ypos[k]), (xpos[k], ypos[k])])    # horizontal
    n = len(tips)
    fig, ax = plt.subplots(figsize=(7.5, min(max(4.0, n * 0.012 + 2.5), 18)))
    ax.add_collection(LineCollection(segs, colors="#555555", linewidths=0.4))
    m = meta.set_index("tip")
    x = np.array([xpos[t] for t in tips])
    y = np.array([ypos[t] for t in tips])
    names = [t.name for t in tips]
    phy = [m.loc[nm, "phylum"] if nm in m.index else "Unknown" for nm in names]
    cols = [colors.get(p, "#9a9a9a") for p in phy]
    anchor = np.array([bool(m.loc[nm, "anchor"]) if nm in m.index and "anchor" in m.columns else False for nm in names])
    seed = np.array([bool(m.loc[nm, "is_seed"]) if nm in m.index and "is_seed" in m.columns else False for nm in names])
    ax.scatter(x[~seed], y[~seed], s=6, c=[c for c, s in zip(cols, seed) if not s], linewidths=0, zorder=3)
    if anchor.any():
        ax.scatter(x[anchor], y[anchor], s=26, facecolors="none", edgecolors="black", linewidths=0.8, zorder=4)
    if seed.any():
        ax.scatter(x[seed], y[seed], s=48, marker="*", c=[c for c, s in zip(cols, seed) if s], edgecolors="black",
                   linewidths=0.5, zorder=5)
    ax.set_ylim(n, -1)
    ax.set_xlim(0, max(x.max() * 1.05, 1e-6))
    ax.set_yticks([])
    ax.set_xlabel("substitutions per site")
    ax.spines["left"].set_visible(False)
    present = [p for p in colors if p in set(phy)]
    handles = [Line2D([], [], marker="o", ls="", color=colors[p], label=p, markersize=5) for p in present]
    if "Unknown" in phy or any(p not in colors for p in phy):
        handles.append(Line2D([], [], marker="o", ls="", color="#9a9a9a", label="other phyla", markersize=5))
    if anchor.any():
        handles.append(Line2D([], [], marker="o", ls="", mfc="none", mec="black", label="anchor species", markersize=7))
    if seed.any():
        handles.append(Line2D([], [], marker="*", ls="", color="#777777", mec="black", label="seed", markersize=9))
    ax.legend(handles=handles, fontsize=7, loc="lower right", frameon=False)
    ax.set_title(title or f"{n} tips")
    return fig


def draw_tree(newick: Path, meta: pd.DataFrame, colors: dict, out_png: Path):
    import matplotlib.pyplot as plt
    fig = tree_figure(newick, meta, colors)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def tree_inputs(stage: Path, slug: str, target: str, max_tips=MAX_TIPS):
    """Tree FASTA and metadata: the best hit of each genus (highest bit score) that passed the thresholds and the
    architecture rule, plus the seeds with a sequence. More than max_tips tips: keep the highest bit scores.
    Returns (fasta text, meta DataFrame)."""
    import ortho_search as ors
    ph = pd.read_csv(stage / f"pass_{slug}.tsv", sep="\t", comment="#")
    reps = pd.read_csv(stage / "genus_reps.tsv", sep="\t")
    ph = ph.sort_values("bits", ascending=False).drop_duplicates("species_prefix").head(max_tips)
    seeds = pd.read_csv(stage / "seeds.tsv", sep="\t", keep_default_na=False)
    seeds = seeds[(seeds["target"] == target) & (seeds["status"] == "ok")]
    prot = {}
    for r in _read_fa(stage / f"cand_{slug}.faa.gz").items():
        prot[r[0]] = r[1]
    ph = ph[ph["protein_id"].isin(prot)]
    phy = reps.set_index("species_prefix")["PHYLUM"].to_dict()
    rows, fa = [], []
    for r in ph.itertuples():
        fa.append(f">{r.protein_id}\n{prot[r.protein_id]}\n")
        rows.append({"tip": r.protein_id, "phylum": phy.get(r.species_prefix, "Unknown"),
                     "anchor": r.species_prefix in ors.ANCHORS, "is_seed": False})
    have = {x["tip"] for x in rows}
    for r in seeds.itertuples():
        if r.seed_id in have or not r.sequence:
            continue
        fa.append(f">{r.seed_id}\n{r.sequence}\n")
        have.add(r.seed_id)
        rows.append({"tip": r.seed_id, "phylum": r.phylum or "Unknown", "anchor": bool(r.genome in ors.ANCHORS), "is_seed": True})
    return "".join(fa), pd.DataFrame(rows)


def _read_fa(path):
    import jgi_map as jm
    seqs, name = {}, None
    with jm._open_text(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def step_trees(stage: Path, scratch: Path):
    for target, slug in SLUG.items():
        nwk = stage / f"tree_{slug}.nwk"
        if nwk.exists():
            print(f"tree {slug}: exists, skipped", flush=True)
            continue
        text, meta = tree_inputs(stage, slug, target)
        fa = Path(scratch) / f"tree_{slug}.faa"
        Path(scratch).mkdir(parents=True, exist_ok=True)
        fa.write_text(text)
        meta.to_csv(stage / f"tree_{slug}_tips.tsv", sep="\t", index=False)
        print(f"tree {slug}: {len(meta)} tips ({int(meta['is_seed'].sum()) if len(meta) else 0} seeds)", flush=True)
        if len(meta) < 3:
            print(f"tree {slug}: fewer than 3 tips, no tree", flush=True)
            continue
        out = build_tree(fa, stage / f"tree_{slug}", Path(scratch) / slug)
        if out is None:
            print(f"tree {slug}: skipped", flush=True)


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["trees"])
    ap.add_argument("--stage", required=True)
    ap.add_argument("--scratch", required=True)
    a = ap.parse_args(argv)
    step_trees(Path(a.stage), Path(a.scratch))


if __name__ == "__main__":
    main(sys.argv[1:])
