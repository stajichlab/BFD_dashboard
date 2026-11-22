#!/usr/bin/env python3
"""Tests A, B and C of README.md.

Usage: 02_tests.py INPUT_DIR OUT_DIR [N_PERM_SCAN]
INPUT_DIR holds units.parquet, features_long.parquet, feature_names.parquet from 01_build_inputs.py.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "bin" / "summary"))

MIN_FINE, N_PERM_A, N_PERM_B = 8, 999, 19999
GENOME_PROPS = ["total_length_bp", "n_genes", "genes_per_mb", "coding_fraction", "gc_pct", "median_intron_len_bp",
                "median_exons_per_gene", "intron_frac_genes_with_intron", "median_protein_len_aa"]
CONTENT = ["cazy_GH_per1000", "cazy_GT_per1000", "cazy_PL_per1000", "cazy_CE_per1000", "cazy_AA_per1000", "cazy_CBM_per1000",
           "secretome_per1000", "merops_per1000", "tf_per1000"]
CONTRASTS = [("plant_pathogen", "saprotroph", "lifestyle_coarse"), ("animal_parasite", "saprotroph", "lifestyle_coarse"),
             ("mycorrhizal", "saprotroph", "lifestyle_coarse"), ("lichenized", "saprotroph", "lifestyle_coarse"),
             ("mycoparasite", "saprotroph", "lifestyle_coarse"), ("yeast", "filamentous_mycelium", "growth_form")]


def bh(p):
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    v = p[ok]
    m = len(v)
    if m == 0:
        return q
    o = np.argsort(v)
    adj = np.minimum.accumulate((v[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


# ------------------------------------------------------------------ Test A
def _dummies(lab, drop_first):
    d = pd.get_dummies(pd.Series(lab).astype(str), drop_first=drop_first).to_numpy(float)
    return d


def _r2(y, X):
    """R2 and rank of the model y ~ [1, X]."""
    Xc = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(Xc, y, rcond=None)
    res = y - Xc @ beta
    tss = float(((y - y.mean()) ** 2).sum())
    return 1 - float((res ** 2).sum()) / tss, int(np.linalg.matrix_rank(Xc)) - 1


def _adj(r2, n, p):
    return 1 - (1 - r2) * (n - 1) / (n - p - 1) if n - p - 1 > 0 else np.nan


def partition(y, tax, eco, n_perm=N_PERM_A, seed=1):
    """Variance partition of y by taxonomy (strata) and ecology. Returns a dict. p: permutation of eco labels within tax."""
    y = np.asarray(y, float)
    tax = np.asarray(tax)
    eco = np.asarray(eco)
    n = len(y)
    Xt, Xe = _dummies(tax, True), _dummies(eco, True)
    r_t, p_t = _r2(y, Xt)
    r_e, p_e = _r2(y, Xe)
    r_b, p_b = _r2(y, np.column_stack([Xt, Xe]))
    a_t, a_e, a_b = _adj(r_t, n, p_t), _adj(r_e, n, p_e), _adj(r_b, n, p_b)
    uniq_e, uniq_t = a_b - a_t, a_b - a_e
    # permutation F of eco after tax (Frisch-Waugh)
    Xt1 = np.column_stack([np.ones(n), Xt])
    Q, _ = np.linalg.qr(Xt1)
    rank_t = int(np.linalg.matrix_rank(Xt1))
    yr = y - Q @ (Q.T @ y)

    def fstat(lab):
        E = _dummies(lab, True)
        Er = E - Q @ (Q.T @ E)
        Qe, Re = np.linalg.qr(Er)
        keep = np.abs(np.diag(Re)) > 1e-9
        Qe = Qe[:, keep]
        q = Qe.shape[1]
        if q == 0:
            return np.nan
        ss = float(((Qe.T @ yr) ** 2).sum())
        rss = float((yr ** 2).sum()) - ss
        df2 = n - rank_t - q
        return (ss / q) / (rss / df2) if df2 > 0 and rss > 0 else np.nan
    f0 = fstat(eco)
    rng = np.random.default_rng(seed)
    idx_by = [np.where(tax == s)[0] for s in np.unique(tax)]
    ge = 0
    for _ in range(n_perm):
        lab = eco.copy()
        for idx in idx_by:
            if len(idx) > 1:
                lab[idx] = eco[rng.permutation(idx)]
        f = fstat(lab)
        if np.isfinite(f) and f >= f0:
            ge += 1
    return {"n": n, "n_tax": len(idx_by), "n_eco": len(np.unique(eco)), "R2adj_taxonomy": a_t, "R2adj_ecology": a_e,
            "R2adj_both": a_b, "unique_ecology": uniq_e, "unique_taxonomy": uniq_t, "shared": a_e - uniq_e,
            "F_ecology_after_taxonomy": f0, "p": (1 + ge) / (1 + n_perm)}


def run_test_a(units, out):
    rows = []
    variants = [("ORDER_", "lifestyle", "order + lifestyle"), ("ORDER_", "lifestyle_coarse", "order + coarse lifestyle"),
                ("FAMILY", "lifestyle_coarse", "family + coarse lifestyle")]
    for tax_col, eco_col, name in variants:
        for resp in GENOME_PROPS + CONTENT:
            d = units[[tax_col, eco_col, resp]].dropna()
            d = d[d[eco_col] != "other"]
            if eco_col == "lifestyle":
                cnt = d[eco_col].value_counts()
                d = d[d[eco_col].isin(cnt[cnt >= MIN_FINE].index)]
            if d[eco_col].nunique() < 2 or len(d) < 50:
                continue
            r = partition(d[resp].to_numpy(), d[tax_col].to_numpy(), d[eco_col].to_numpy())
            rows.append({"model": name, "response": resp, "kind": "genome" if resp in GENOME_PROPS else "protein_content", **r})
    t = pd.DataFrame(rows)
    t["q_BH"] = t.groupby("model")["p"].transform(lambda s: pd.Series(bh(s.to_numpy()), index=s.index))
    t["ecology_explains"] = (t.unique_ecology >= 0.02) & (t.q_BH < 0.05)
    t.to_csv(out / "A_variance_partition.tsv", sep="\t", index=False)
    return t


def run_permanova(units, feats, out):
    import summary_ecology as se
    from scipy.spatial.distance import pdist, squareform
    rows = []
    sets = {"pfam_profile": "pfam:", "cazy_family_profile": "cazy:", "tf_family_profile": "tf:"}
    u = units.set_index("species_prefix")
    for name, pre in sets.items():
        f = feats[feats.feature.str.startswith(pre)]
        m = f.pivot_table(index="species_prefix", columns="feature", values="rate", fill_value=0.0)
        m = m.reindex(u.index[u.has_function]).fillna(0.0).dropna(how="all")
        present = (m > 0).mean()
        m = m.loc[:, present >= 0.05]
        m = m[m.sum(axis=1) > 0]                      # Bray-Curtis is undefined for an all-zero profile
        dist = pd.DataFrame(squareform(pdist(m.to_numpy(), "braycurtis")), index=m.index, columns=m.index)
        for trait in ("lifestyle", "lifestyle_coarse"):
            lab0 = u.loc[m.index, trait].where(lambda s: s != "other")
            for test, st in (("unrestricted", None), ("within_order", u.loc[m.index, "ORDER_"])):
                lab, n_rare = se.mask_rare(lab0, st, MIN_FINE)
                r = se.stratified_permanova(dist, lab, st)
                rows.append({"profile": name, "n_features": m.shape[1], "trait": trait, "test": test, "n_genera": r["n"],
                             "n_groups": r["n_groups"], "pseudo_F": r["F"], "R2_adj": r["R2_adj"], "p": r["p"], "caveat": r["caveat"]})
    t = pd.DataFrame(rows)
    t.to_csv(out / "A_permanova_profiles.tsv", sep="\t", index=False)
    return t


# ------------------------------------------------------------------ Tests B and C
def clade_design(df, grp_col, a, b, clade_col="FAMILY", min_per=3):
    """Rows in groups a or b, restricted to clades with >= min_per of each. Returns (sub, clade array, is_a array)."""
    s = df[df[grp_col].isin([a, b]) & df[clade_col].notna()]
    cnt = s.groupby([clade_col, grp_col]).size().unstack(fill_value=0)
    ok = cnt.index[(cnt.get(a, 0) >= min_per) & (cnt.get(b, 0) >= min_per)]
    s = s[s[clade_col].isin(ok)]
    return s, s[clade_col].to_numpy(), (s[grp_col] == a).to_numpy()


def label_vector(clade, is_a):
    """Weights v with stat = v @ y = sum_c w_c (mean_A,c - mean_B,c) / sum w_c, w_c = min(nA_c, nB_c)."""
    v = np.zeros(len(clade))
    W = 0.0
    for c in np.unique(clade):
        m = clade == c
        na, nb = int(is_a[m].sum()), int((~is_a[m]).sum())
        w = min(na, nb)
        W += w
        v[m] = np.where(is_a[m], w / na, -w / nb)
    return v / W


def permute_vectors(v, clade, n_perm, rng):
    """n_perm x n matrix of label vectors with the labels shuffled within clade (v values move with the labels)."""
    V = np.empty((n_perm, len(v)))
    groups = [np.where(clade == c)[0] for c in np.unique(clade)]
    for k in range(n_perm):
        p = np.arange(len(v))
        for idx in groups:
            p[idx] = idx[rng.permutation(len(idx))]
        V[k] = v[p]
    return V


def matched_stat_test(Y, clade, is_a, n_perm=N_PERM_B, seed=1, block=400):
    """Y: n x D (no NaN). Returns pooled effect (D,), two-sided p (D,), per-clade effects (K x D), sign agreement (D,)."""
    rng = np.random.default_rng(seed)
    v = label_vector(clade, is_a)
    obs = v @ Y
    V = permute_vectors(v, clade, n_perm, rng)
    ge = np.zeros(Y.shape[1])
    for s in range(0, Y.shape[1], block):
        P = V @ Y[:, s:s + block]
        ge[s:s + block] = (np.abs(P) >= np.abs(obs[s:s + block]) - 1e-12).sum(0)
    p = (1 + ge) / (1 + n_perm)
    clades = np.unique(clade)
    eff = np.vstack([Y[(clade == c) & is_a].mean(0) - Y[(clade == c) & ~is_a].mean(0) for c in clades])
    agree = (np.sign(eff) == np.sign(obs)).mean(0)
    return obs, p, eff, agree, clades


def verdict(q, agree, n_clades):
    if n_clades < 4:
        return "not evaluable"
    if q < 0.05 and agree >= 0.75:
        return "supported"
    if q < 0.05:
        return "mixed"
    return "no evidence"


def run_test_b(units, out, clade_col="FAMILY", tag="family"):
    rows, clade_rows = [], []
    for a, b, col in CONTRASTS:
        for resp in GENOME_PROPS + CONTENT:
            d = units.dropna(subset=[resp])
            sub, clade, is_a = clade_design(d, col, a, b, clade_col)
            ncl = len(np.unique(clade))
            if ncl < 4:
                rows.append({"contrast": f"{a} vs {b}", "response": resp, "n_clades": ncl, "n_a": int(is_a.sum()), "n_b": int((~is_a).sum()),
                             "verdict": "not evaluable"})
                continue
            Y = sub[[resp]].to_numpy(float)
            obs, p, eff, agree, clades = matched_stat_test(Y, clade, is_a)
            rows.append({"contrast": f"{a} vs {b}", "response": resp, "n_clades": ncl, "n_a": int(is_a.sum()), "n_b": int((~is_a).sum()),
                         "pooled_effect": float(obs[0]), "p": float(p[0]), "share_clades_same_sign": float(agree[0])})
            for c, e in zip(clades, eff[:, 0]):
                m = clade == c
                clade_rows.append({"contrast": f"{a} vs {b}", "response": resp, "clade": c, "n_a": int((is_a & m).sum()),
                                   "n_b": int((~is_a & m).sum()), "effect": float(e)})
    t = pd.DataFrame(rows)
    t["q_BH"] = np.nan
    for c, g in t.groupby("contrast"):
        t.loc[g.index, "q_BH"] = bh(g["p"].to_numpy(float))
    t["verdict"] = [verdict(q, ag, n) if pd.notna(q) else "not evaluable"
                    for q, ag, n in zip(t.q_BH, t.share_clades_same_sign, t.n_clades)]
    t.insert(0, "clade_level", tag)
    t.to_csv(out / f"B_matched_contrasts_{tag}.tsv", sep="\t", index=False)
    pd.DataFrame(clade_rows).to_csv(out / f"B_per_clade_effects_{tag}.tsv", sep="\t", index=False)
    return t


def run_test_c(units, feats, names, out, n_perm, clade_col="FAMILY", tag="family"):
    rows = []
    u = units[units.has_function].copy()
    fw = feats.pivot_table(index="species_prefix", columns="feature", values="rate", fill_value=0.0)
    nm = dict(zip(names.feature, names.name))
    for a, b, col in CONTRASTS:
        sub, clade, is_a = clade_design(u, col, a, b, clade_col)
        ncl = len(np.unique(clade))
        if ncl < 4:
            continue
        M = fw.reindex(sub.species_prefix).fillna(0.0)
        M = M.loc[:, (M > 0).mean() >= 0.10]
        Y = np.log10(1 + M.to_numpy())
        obs, p, eff, agree, clades = matched_stat_test(Y, clade, is_a, n_perm=n_perm)
        q = bh(p)
        for j, f in enumerate(M.columns):
            rows.append({"clade_level": tag, "contrast": f"{a} vs {b}", "feature": f, "name": nm.get(f, ""), "n_clades": ncl, "n_a": int(is_a.sum()),
                         "n_b": int((~is_a).sum()), "pooled_effect_log10": float(obs[j]), "p": float(p[j]), "q_BH": float(q[j]),
                         "share_clades_same_sign": float(agree[j])})
    t = pd.DataFrame(rows)
    t.to_csv(out / f"C_domain_scan_{tag}.tsv", sep="\t", index=False)
    return t


def main(inp, out, n_perm=N_PERM_B, parts="a,perm,b,c"):
    inp, out = Path(inp), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    parts = set(parts.split(","))
    units = pd.read_parquet(inp / "units.parquet")
    feats = pd.read_parquet(inp / "features_long.parquet")
    names = pd.read_parquet(inp / "feature_names.parquet")
    units = units.dropna(subset=["n_genes"])                 # 15 genomes without gene statistics
    if "a" in parts:
        a = run_test_a(units, out)
        print(a[["model", "response", "n", "R2adj_taxonomy", "R2adj_ecology", "unique_ecology", "shared", "p", "q_BH", "ecology_explains"]]
              .round(3).to_string(index=False), flush=True)
    if "perm" in parts:
        print(run_permanova(units, feats, out).round(4).to_string(index=False), flush=True)
    for clade_col, tag in (("FAMILY", "family"), ("ORDER_", "order")):
        if "b" in parts:
            b = run_test_b(units, out, clade_col, tag)
            print(tag, b.groupby("contrast").verdict.value_counts().to_string(), flush=True)
        if "c" in parts:
            c = run_test_c(units, feats, names, out, n_perm, clade_col, tag)
            if len(c):
                s = c.groupby("contrast").agg(n_features=("feature", "size"), n_q05=("q_BH", lambda x: int((x < 0.05).sum())))
                print(tag, s.to_string(), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else N_PERM_B, sys.argv[4] if len(sys.argv) > 4 else "a,perm,b,c")
