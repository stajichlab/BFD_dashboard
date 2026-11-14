#!/usr/bin/env python3
"""Co-absence test with an annotation-propensity null.

Question: in genomes where trehalose-6-P synthase (PF00982) has no Pfam record,
is trehalose-6-P phosphatase (PF02358) also absent more often than expected?
The naive 2x2 Fisher test is confounded because a genome that misses one family
tends to miss others (annotation gaps). Expected co-absence is therefore
computed Mantel-Haenszel style, stratifying genomes by their total number of
missing near-universal families, following the "per-family loss-rate null"
called for in docs/FUNCTIONAL_ANALYSES.md (agenda item C4 / A2).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from stats_util import bh_qvalues, expected_coabsence  # noqa: E402

DATA = Path(os.environ.get("TH_DATA", HERE.parent / "data"))


def main():
    m = pd.read_csv(DATA / "th_uvmiss.csv")
    genomes = np.sort(m.genome.unique())
    gi = {g: i for i, g in enumerate(genomes)}
    fams = np.sort(m.pfam_acc.unique())
    fi = {f: i for i, f in enumerate(fams)}
    M = np.zeros((len(genomes), len(fams)), bool)
    M[m.genome.map(gi).values, m.pfam_acc.map(fi).values] = True
    n_missing = M.sum(1)
    print(f"genomes {len(genomes)}  families with >=1 absence {len(fams)}  "
          f"median missing/genome {np.median(n_missing):.0f}  max {n_missing.max()}")

    # the trehalose genes are themselves not in the near-universal set (that is the
    # point), so focal presence/absence comes from the per-genome module matrix
    key = pd.read_csv(DATA / "th_genome_matrix.tsv", sep="\t").set_index("genome")
    q = [g for g in genomes if key.loc[g, "keep"]]
    qi = np.array([gi[g] for g in q])
    M, n_missing = M[qi], n_missing[qi]
    genomes = np.array(q)
    focal_extra = {"PF00982": "TPS", "PF02358": "TPP", "PF01204": "NTH", "PF03632": "ATH",
                   "PF05693": "GSY", "PF00343": "GPH"}
    abs_extra = {k: (key.loc[genomes, v].values == 0) for k, v in focal_extra.items()}

    cls = np.array([key.loc[g, "CLASS"] for g in genomes])

    def stratified(focal_acc, target_acc):
        """Both nulls: n_missing bins only (first version) and n_missing bins x class (added after the Opus review)."""
        f = abs_extra[focal_acc] if focal_acc in abs_extra else M[:, fi[focal_acc]]
        t = abs_extra[target_acc] if target_acc in abs_extra else M[:, fi[target_acc]]
        a = expected_coabsence(f, t, n_missing)
        b = expected_coabsence(f, t, n_missing, group=cls)
        if a is None:
            return None
        row = dict(focal=focal_acc, target=target_acc, n_focal_absent=int(f.sum()),
                   obs_coabsent=a["obs"], expected=round(a["exp"], 2),
                   obs_exp=round(a["obs"] / a["exp"], 1), p_two=round(a["p_two"], 6))
        row.update(expected_class=round(b["exp"], 2) if b else np.nan,
                   obs_exp_class=round(b["obs"] / b["exp"], 1) if b else np.nan,
                   p_two_class=round(b["p_two"], 6) if b else np.nan)
        return row

    # 1. focal = TPS, every other family as target (near-universal set + pathway genes)
    targets = sorted((set(fams) | set(focal_extra)) - {"PF00982"})
    rows = [r for r in (stratified("PF00982", a) for a in targets) if r]
    R = pd.DataFrame(rows).sort_values("obs_exp", ascending=False)
    R["q_BH"] = bh_qvalues(R.p_two.values)                # BH on p in input order (was applied to effect-size order: all q = 0)
    R["q_BH_class"] = bh_qvalues(R.p_two_class.values)
    R["rank"] = np.arange(1, len(R) + 1)
    R.to_csv(DATA / "th_coloss_vs_TPS.tsv", sep="\t", index=False)
    ids = {"PF02358": "TPP", "PF01204": "NTH", "PF05693": "GSY", "PF00343": "GPH", "PF01704": "UGP",
           "PF03632": "ATH", "PF07479": "GPD", "PF02878": "PGM", "PF00342": "PGI", "PF00128": "GH13",
           "PF00011": "sHSP", "PF04119": "HSP12"}
    print("\n== co-absence with TPS, ranked by obs/expected (annotation-propensity-stratified) ==")
    print(R.head(12).assign(module=lambda d: d.target.map(ids).fillna("")).to_string(index=False))
    print("\nwhere the other trehalose/glycogen genes rank:")
    print(R[R.target.isin(ids)].assign(module=lambda d: d.target.map(ids))[
        ["target", "module", "n_focal_absent", "obs_coabsent", "expected", "obs_exp", "p_two", "q_BH",
         "expected_class", "obs_exp_class", "p_two_class", "q_BH_class", "rank"]
    ].to_string(index=False))

    # 2. symmetric: focal = TPP
    print("\n== focal = TPP: co-absence of TPS ==")
    print(stratified("PF02358", "PF00982"))

    # 3. who are the TPS-absent genomes
    f = abs_extra["PF00982"]
    tab = pd.DataFrame({"GENUS": [key.loc[g, "GENUS"] for g in genomes[f]],
                        "CLASS": [key.loc[g, "CLASS"] for g in genomes[f]],
                        "n_missing": n_missing[f],
                        "qcov": [key.loc[g, "qcov"] for g in genomes[f]],
                        "BUSCO": [key.loc[g, "complete_pct"] for g in genomes[f]]})
    print(f"\n{len(tab)} quality-kept genomes (of {len(genomes)}) with no gated PF00982 record; "
          f"they sit in {tab.GENUS.nunique()} genera and {tab.CLASS.nunique()} classes")
    print(tab.CLASS.value_counts().head(8).to_string())
    print("\nTPS-absent genomes missing >=5 universal families (i.e. gap-prone):",
          int((tab.n_missing >= 5).sum()), "of", len(tab))
    tab.to_csv(DATA / "th_tpsabsent_genomes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
