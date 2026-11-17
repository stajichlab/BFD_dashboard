"""T3 helpers: within-genus contrasts of trehalose-synthesis dose and balance between labelled and unlabelled species.

Everything here is a pure function of small tables so it can be tested on fixtures (test_t3_lib.py).

Labels. The only species-level label available in the repo or BFD is `dry_epithet`: the species name matches DRY_EPITHET (the
same pattern as th1_survey.py). It is a name-based flag, not a measured water activity or a curated xerophily call, and it is
reported as such. Genus-level XERO_GENERA cannot be used for within-genus contrasts.

Metrics (fixed before any outcome was computed):
  dose     D = n_tps * 1000 / n_prot       TPS proteins per 1,000 proteins (SwissProt best-hit counts, th_sp.csv)
  arm      A = (n_gsy + n_gph) * 1000 / n_prot
  balance  B = log2((n_tps + 0.5) / (n_gsy + n_gph + 0.5))
Statistic: for each eligible genus (>= MIN_FLAG labelled and >= MIN_NON unlabelled kept species, one genome per species),
mean(metric | labelled) - mean(metric | unlabelled); the combined statistic is the mean over genera (each genus counts once).
Null: the labelled/unlabelled flag is permuted within each genus (stratified null), so genus and everything shared inside a genus
(clade, genome size, annotation style) is held fixed.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

DRY_EPITHET = r"(?i)antarct|arctic|desert|endolith|xero|halophil|salin|cavernicol|arid|mumm"  # identical to th1_survey.py
MIN_FLAG = 2
MIN_NON = 3
N_PERM = 20000


def load_species(data: Path) -> pd.DataFrame:
    """One row per species representative with the T3 metrics and the `keep` quality flag (same rule as th1_survey.load)."""
    g = pd.read_csv(data / "th_genomes.csv")
    sp = pd.read_csv(data / "th_sp.csv")
    uc = pd.read_csv(data / "th_univcov.csv")
    n_univ = len(pd.read_csv(data / "th_univ.csv"))
    X = g.merge(sp, on="genome", how="left").merge(uc[["genome", "n_univ_found"]], on="genome", how="left")
    for c in ["n_tps", "n_gsy", "n_gph"]:
        X[c] = X[c].fillna(0)
    X["qcov"] = X.n_univ_found.fillna(0) / n_univ
    X["busco_class_p25"] = X.groupby("CLASS").complete_pct.transform(lambda s: s.quantile(0.25))
    X["keep"] = (X.qcov >= 0.98) & (X.complete_pct >= X.busco_class_p25)
    X["label"] = X.SPECIES.fillna("").str.contains(DRY_EPITHET, regex=True)
    X["D"] = X.n_tps * 1000.0 / X.n_prot
    X["A"] = (X.n_gsy + X.n_gph) * 1000.0 / X.n_prot
    X["B"] = np.log2((X.n_tps + 0.5) / (X.n_gsy + X.n_gph + 0.5))
    return X


def design_table(X: pd.DataFrame) -> pd.DataFrame:
    """Per genus: kept species, how many carry the label. No outcome metric is touched."""
    K = X[X.keep]
    t = K.groupby("GENUS").agg(CLASS=("CLASS", "first"), n_species=("genome", "size"), n_label=("label", "sum")).reset_index()
    t["n_unlabelled"] = t.n_species - t.n_label
    t["eligible"] = (t.n_label >= MIN_FLAG) & (t.n_unlabelled >= MIN_NON)
    return t.sort_values(["eligible", "n_label"], ascending=False)


def genus_shift(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Per-genus difference of means (labelled minus unlabelled) for the eligible genera in `df` (columns GENUS, label, metric)."""
    rows = []
    for gname, s in df.groupby("GENUS"):
        a, b = s.loc[s.label, metric], s.loc[~s.label, metric]
        if len(a) >= MIN_FLAG and len(b) >= MIN_NON:
            rows.append((gname, len(a), len(b), float(a.mean() - b.mean())))
    return pd.DataFrame(rows, columns=["GENUS", "n_label", "n_unlabelled", "shift"])


def stratified_permutation(df: pd.DataFrame, metric: str, n_perm: int = N_PERM, seed: int = 20261115) -> dict:
    """Mean over eligible genera of the labelled-minus-unlabelled shift, with a two-sided permutation p (flag shuffled within genus)."""
    rng = np.random.default_rng(seed)
    groups = []
    for _, s in df.groupby("GENUS"):
        lab = s.label.to_numpy()
        if lab.sum() >= MIN_FLAG and (~lab).sum() >= MIN_NON:
            groups.append((s[metric].to_numpy(float), int(lab.sum())))
    if not groups:
        return {"n_genera": 0, "obs": float("nan"), "p": float("nan")}
    def stat(perm: bool) -> float:
        out = []
        for v, k in groups:
            if perm:
                idx = rng.permutation(v.size)
                a, b = v[idx[:k]], v[idx[k:]]
            else:
                a, b = v[:k], v[k:]
            out.append(a.mean() - b.mean())
        return float(np.mean(out))
    # observed uses the real labels: reorder each genus so labelled values come first
    groups_real = []
    for _, s in df.groupby("GENUS"):
        lab = s.label.to_numpy()
        if lab.sum() >= MIN_FLAG and (~lab).sum() >= MIN_NON:
            v = s[metric].to_numpy(float)
            groups_real.append((np.concatenate([v[lab], v[~lab]]), int(lab.sum())))
    obs = float(np.mean([v[:k].mean() - v[k:].mean() for v, k in groups_real]))
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = stat(True)
    p = float((1 + np.sum(np.abs(null) >= abs(obs))) / (n_perm + 1))
    return {"n_genera": len(groups), "obs": obs, "p": p, "null_sd": float(null.std())}
