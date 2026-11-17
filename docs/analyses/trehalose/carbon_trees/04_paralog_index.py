#!/usr/bin/env python3
"""Paralog index (PI) of the carbon-storage family trees and the pre-registered TP1 verdict (README.md).

Usage: 04_paralog_index.py BENCH_DIR OUT_DIR
Reads BENCH/tips.tsv and BENCH/trees/<FAM>.fasttree.nwk (and <FAM>.treefile for the IQ-TREE families when present, reported separately).
PI of a class = median over genomes with >= 2 copies of
    (median patristic distance among the genome's own copies) / (median patristic distance between copies of different genera of the class),
the reference being MAX_PAIRS random pairs per class (fixed seed). PI < 1 means copies of one genome are closer than copies of different
genera (lineage-specific duplication); PI near or above 1 means the copies diverged before the genera did.
TP1: supported if PI >= 0.8 for TPS and TPP in both Mucoromycetes and Mortierellomycetes; killed if any of the four is < 0.5; else mixed;
not evaluable if a class has fewer than MIN_GENOMES genomes with >= 2 copies.
"""
import csv
import itertools
import statistics
import sys
import threading
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "kx2_ku_pathway" / "family_trees"))   # docs/analyses/kx2_ku_pathway/family_trees
import fam_lib  # noqa: E402

MAX_PAIRS, MIN_GENOMES, SEED, N_BOOT = 2000, 5, 20261116, 2000
TARGET = ("Mucoromycetes", "Mortierellomycetes")
TEST_FAMILIES = ("TPS", "TPP")
BASELINE = ("Sordariomycetes", "Agaricomycetes", "Eurotiomycetes")


def leaf_map(root):
    return {n.name: n for n in fam_lib.leaves(root)}


def paralog_index(root, tips, rng):
    """tips: rows (protein_id, locustag, genus, class). Returns {class: dict(n_genomes, n_multi, median_within, ref, pi, ci)}."""
    lm = leaf_map(root)
    by_class = defaultdict(lambda: defaultdict(list))
    genus_of = {}
    for pid, loc, genus, cls in tips:
        key = f"{pid}|{loc}"
        if key in lm:
            by_class[cls][loc].append(lm[key])
            genus_of[loc] = genus
    out = {}
    for cls, genomes in by_class.items():
        within, closest = {}, {}
        for loc, nodes in genomes.items():
            if len(nodes) >= 2:
                d = [fam_lib.path_distance(a, b) for a, b in itertools.combinations(nodes, 2)]
                within[loc] = statistics.median(d)
                closest[loc] = min(d)   # post hoc supplement: the closest pair of copies in the genome
        locs = [loc for loc in genomes]
        ref_d = []
        tries = 0
        while len(ref_d) < MAX_PAIRS and tries < 20 * MAX_PAIRS and len(locs) > 1:
            tries += 1
            a, b = rng.choice(len(locs), 2, replace=False)
            la, lb = locs[a], locs[b]
            if genus_of[la] == genus_of[lb]:
                continue
            na, nb = genomes[la][rng.integers(len(genomes[la]))], genomes[lb][rng.integers(len(genomes[lb]))]
            ref_d.append(fam_lib.path_distance(na, nb))
        ref = statistics.median(ref_d) if ref_d else float("nan")
        pis = [w / ref for w in within.values()] if ref and ref == ref and ref > 0 else []
        # post hoc supplement (not part of the TP1 rule): share of genomes whose closest pair is closer than half the reference
        recent = float(np.mean([c / ref < 0.5 for c in closest.values()])) if closest and ref and ref > 0 else float("nan")
        if pis:
            boot = [float(np.median(rng.choice(pis, len(pis)))) for _ in range(N_BOOT)]
            ci = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))
        else:
            ci = (float("nan"), float("nan"))
        out[cls] = {"n_genomes": len(genomes), "n_multi": len(within), "median_within": statistics.median(within.values()) if within else float("nan"),
                    "ref": ref, "pi": float(np.median(pis)) if pis else float("nan"), "ci_lo": ci[0], "ci_hi": ci[1],
                    "recent_pair_frac": recent}
    return out


def tp1_verdict(table):
    """table: {(family, class): row}. Pre-registered rule."""
    vals = []
    for fam in TEST_FAMILIES:
        for cls in TARGET:
            r = table.get((fam, cls))
            if r is None or r["n_multi"] < MIN_GENOMES or r["pi"] != r["pi"]:
                return "not evaluable (a test class has fewer than %d genomes with >= 2 copies)" % MIN_GENOMES
            vals.append(r["pi"])
    if any(v < 0.5 for v in vals):
        return "killed (recent duplications: PI < 0.5 in at least one test cell)"
    if all(v >= 0.8 for v in vals):
        return "supported (ancient paralogs: PI >= 0.8 in all four test cells)"
    return "mixed (PI between 0.5 and 0.8 in at least one test cell, none below 0.5)"


def run(bench, out):
    bench, out = Path(bench), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    tips = defaultdict(list)
    for r in csv.DictReader(open(bench / "tips.tsv"), delimiter="\t"):
        tips[r["family"]].append((r["protein_id"], r["locustag"], r["genus"], r["class"]))
    rng = np.random.default_rng(SEED)
    rows, table = [], {}
    for fam in tips:
        for kind, fn in (("fasttree", f"{fam}.fasttree.nwk"), ("iqtree", f"{fam}.treefile")):
            path = bench / "trees" / fn
            if not path.exists():
                continue
            root = fam_lib.parse_newick(path.read_text())
            res = paralog_index(root, tips[fam], rng)
            for cls, r in res.items():
                if r["n_genomes"] >= 8:
                    rows.append({"family": fam, "tree": kind, "class": cls, **r})
                    if kind == "fasttree":
                        table[(fam, cls)] = r
    cols = ["family", "tree", "class", "n_genomes", "n_multi", "median_within", "ref", "pi", "ci_lo", "ci_hi", "recent_pair_frac"]
    with open(out / "paralog_index.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    verdict = tp1_verdict(table)
    lines = ["TP1: " + verdict, ""]
    for fam in TEST_FAMILIES:
        for cls in TARGET + BASELINE:
            r = table.get((fam, cls))
            if r:
                lines.append(f"{fam} {cls}: PI {r['pi']:.2f} (95% CI {r['ci_lo']:.2f} to {r['ci_hi']:.2f}); {r['n_multi']} of {r['n_genomes']} genomes with >= 2 copies; "
                         f"genomes with a recent pair (closest pair < 0.5 x reference, post hoc): {r['recent_pair_frac']:.0%}")
    (out / "tp1_verdict.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return verdict


def main(bench, out):
    sys.setrecursionlimit(1_000_000)
    threading.stack_size(512 * 1024 * 1024)
    result = {}
    t = threading.Thread(target=lambda: result.update(v=run(bench, out)))
    t.start()
    t.join()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
