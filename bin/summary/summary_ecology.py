"""Ecology (FungalTraits) in the summary ordinations: a separate step that runs after make_figures.

Reads the same per-genome tables as make_figures read-only, joins lifestyle and growth form from the DB view
v_genome_labels, and writes figures and tables under its own output directory. It changes no shared file.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import colorsys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import summary_ordination as so  # noqa: E402
from figs_common import safe  # noqa: E402
from make_figures import Out  # noqa: E402


def permute_within(labels: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Return a copy of labels with the values shuffled inside each stratum only."""
    out = np.array(labels, copy=True)
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]
        if len(idx) > 1:
            out[idx] = out[rng.permutation(idx)]
    return out


def _na(n, k, caveat, n_strata=0, n_drop=0):
    return {"F": np.nan, "R2": np.nan, "R2_adj": np.nan, "p": np.nan, "n": int(n), "n_groups": int(k),
            "n_strata_used": int(n_strata), "n_dropped_no_stratum": int(n_drop), "n_permutations_run": 0, "caveat": caveat}


def stratified_permanova(dist: pd.DataFrame, groups: pd.Series, strata: pd.Series | None = None,
                         n_permutations: int = so.N_PERM, seed: int = so.SEED) -> dict:
    """PERMANOVA of `groups`. With `strata`, labels are permuted only inside a stratum. Rows with a NULL group
    are left out; with strata, rows with a NULL stratum are also left out and counted. Strata with one row do
    not move. Degenerate designs give NaN and a caveat."""
    g = groups.reindex(dist.index)
    keep = g.notna()
    n_drop = 0
    st = None
    if strata is not None:
        st = strata.reindex(dist.index)
        n_drop = int((keep & st.isna()).sum())
        keep = keep & st.notna()
    n = int(keep.sum())
    lab = g[keep].to_numpy()
    k = len(np.unique(lab))
    if n < 3 or k < 2:
        return _na(n, k, "not tested: fewer than 3 genera or fewer than 2 groups", n_drop=n_drop)
    sv = None
    n_used = 0
    if st is not None:
        sv = st[keep].to_numpy()
        sizes = pd.Series(sv).value_counts()
        n_used = int((sizes >= 2).sum())
        if n_used < 2:
            return _na(n, k, "not tested: fewer than 2 strata with at least 2 genera", n_used, n_drop)
    d = dist.loc[keep, keep].to_numpy()
    f = float(so.lib_ordination._pseudo_f(d, lab))
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_permutations):
        perm = permute_within(lab, sv, rng) if sv is not None else lab[rng.permutation(n)]
        if so.lib_ordination._pseudo_f(d, perm) >= f:
            ge += 1
    p = (1 + ge) / (1 + n_permutations)
    r2 = f * (k - 1) / (f * (k - 1) + (n - k)) if np.isfinite(f) and n > k else np.nan
    caveat = ""
    if np.isinf(f):
        caveat = "pseudo-F is infinite: no within-group variation"
    return {"F": f, "R2": float(r2), "R2_adj": so._r2_adj(r2, n, k), "p": float(p), "n": n, "n_groups": k,
            "n_strata_used": n_used, "n_dropped_no_stratum": n_drop, "n_permutations_run": n_permutations, "caveat": caveat}


def bh_adjust(p: pd.Series) -> pd.Series:
    """Benjamini-Hochberg q values. NaN stays NaN and is not counted."""
    q = pd.Series(np.nan, index=p.index, dtype=float)
    ok = p.notna()
    v = p[ok].to_numpy(dtype=float)
    m = len(v)
    if m == 0:
        return q
    order = np.argsort(v)
    adj = v[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def pool_groups(labels: pd.Series, min_group: int, other: str = "other") -> pd.Series:
    """Groups with fewer than min_group members become 'other (k groups)'. NULL stays NULL."""
    counts = labels.value_counts()
    small = set(counts[counts < min_group].index)
    if not small:
        return labels.copy()
    name = f"{other} ({len(small)} groups)"
    return labels.where(~labels.isin(small), name)


LEVEL_COLUMN = {"phylum": "PHYLUM", "class": "CLASS", "order": "ORDER"}
STRATUM_COLUMN = {"phylum": "ORDER", "class": "ORDER", "order": "FAMILY"}


def load_traits(db_path: str) -> pd.DataFrame:
    """species_prefix, primary_lifestyle, growth_form from v_genome_labels, read-only."""
    import duckdb
    con = duckdb.connect(db_path, read_only=True)
    try:
        have = {r[0] for r in con.sql("SELECT column_name FROM information_schema.columns WHERE table_name = 'v_genome_labels'").fetchall()}
        for col in ("primary_lifestyle", "growth_form"):
            if col not in have:
                raise RuntimeError(f"{db_path}: v_genome_labels has no column {col}; rebuild the database with workflows/build_bfd_db.sh")
        return con.sql("SELECT LOCUSTAG AS species_prefix, primary_lifestyle, growth_form FROM v_genome_labels").df()
    finally:
        con.close()


@dataclass
class Clade:
    level: str
    name: str
    ids: list


def select_clades(reps: pd.DataFrame, labels: pd.Series, levels: tuple = ("phylum", "class", "order"),
                  min_genera: int = 30, min_group: int = 5):
    """Clades with at least min_genera labeled genera and at least 2 labels with at least min_group genera each.
    Only labeled genera are counted. Clade.ids holds every genus of the clade, labeled or not."""
    clades, skipped = [], []
    for level in levels:
        col = LEVEL_COLUMN[level]
        for name, sub in reps[reps[col].notna()].groupby(col, sort=True):
            lab = labels.reindex(sub["species_prefix"]).dropna()
            n_lab = len(lab)
            n_ge = int((lab.value_counts() >= min_group).sum())
            if n_lab < min_genera:
                reason = f"fewer than {min_genera} labeled genera"
            elif n_ge < 2:
                reason = f"fewer than 2 lifestyles with at least {min_group} genera"
            else:
                clades.append(Clade(level, str(name), sub["species_prefix"].tolist()))
                continue
            skipped.append({"level": level, "clade": str(name), "n_labeled": n_lab, "n_groups_ge_min": n_ge, "reason": reason})
    cols = ["level", "clade", "n_labeled", "n_groups_ge_min", "reason"]
    return clades, pd.DataFrame(skipped, columns=cols)


NO_LABEL = "#d0d0d0"
OTHER_GREY = "#9a9a9a"
TOPIC = "ecology"
TRAITS = (("primary_lifestyle", "lifestyle", "primary lifestyle"), ("growth_form", "growthform", "growth form"))


def lifestyle_palette(groups: list) -> dict:
    """One distinct, non-grey colour per group. The two reserved greys are never used for a real group."""
    cand = []
    for name in ("tab20", "tab20b", "tab20c", "Set1", "Dark2"):
        for c in plt.get_cmap(name).colors:
            r, g, b = matplotlib.colors.to_rgb(c)
            if colorsys.rgb_to_hsv(r, g, b)[1] >= 0.25:   # skip greys
                hexc = matplotlib.colors.to_hex((r, g, b))
                if hexc not in cand:
                    cand.append(hexc)
    n = len(groups)
    while len(cand) < n:   # many groups: add evenly spaced hues
        k = len(cand)
        cand.append(matplotlib.colors.to_hex(colorsys.hsv_to_rgb((k * 0.618034) % 1, 0.65, 0.8)))
    return {g: cand[i] for i, g in enumerate(groups)}


def scatter_labels(coords: pd.DataFrame, labels: pd.Series, palette: dict, pct: list, title: str):
    """PCo1 against PCo2 coloured by labels (indexed like coords). NULL labels draw in NO_LABEL and are listed as
    'no label (n)'; groups named 'other (...)' draw in OTHER_GREY. Legend right of the axes, 2 columns, 6 pt."""
    lab = labels.reindex(coords.index) if labels.index.equals(coords.index) is False and set(coords.index) <= set(labels.index) else labels
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    y = coords["PCo2"] if "PCo2" in coords else 0
    un = lab.isna().to_numpy()
    if un.any():
        ax.scatter(coords["PCo1"][un], (y[un] if "PCo2" in coords else 0), s=9, c=NO_LABEL, alpha=0.85, linewidths=0, rasterized=True)
    counts = lab.value_counts()
    handles = []
    for g in counts.index:
        col = OTHER_GREY if str(g).startswith("other (") else palette.get(g, OTHER_GREY)
        m = (lab == g).to_numpy()
        ax.scatter(coords["PCo1"][m], (y[m] if "PCo2" in coords else 0), s=9, c=col, alpha=0.85, linewidths=0, rasterized=True)
        handles.append(Line2D([], [], marker="o", ls="", color=col, markersize=4, label=f"{g} ({int(counts[g])})"))
    if un.any():
        handles.append(Line2D([], [], marker="o", ls="", color=NO_LABEL, markersize=4, label=f"no label ({int(un.sum())})"))
    ax.set_xlabel(f"PCo1 ({pct[0]:.1f}%)")
    ax.set_ylabel(f"PCo2 ({pct[1]:.1f}%)" if len(pct) > 1 else "PCo2")
    ax.set_title(title)
    leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), ncol=2, fontsize=6,
                    title=labels.name or None, title_fontsize=6.5, frameon=False, handletextpad=0.3, columnspacing=1.0)
    fig.canvas.draw()
    w_in = leg.get_window_extent().width / fig.dpi
    right = 1 - (w_in + 0.05 * fig.get_figwidth() / 8.4 + 0.1) / fig.get_figwidth()
    fig.subplots_adjust(right=max(0.5, right))
    return fig


def _test_rows(name, trait, dist, labels, strata, min_group):
    """Unrestricted and within-stratum PERMANOVA rows for one matrix and trait. Labels with fewer than min_group
    genera in the tested set are dropped from the test and counted."""
    lab = labels.reindex(dist.index)
    counts = lab.value_counts()
    rare = lab.isin(counts[counts < min_group].index)
    lab = lab.mask(rare)
    rows = []
    for test, st in (("unrestricted", None), ("within_order", strata.reindex(dist.index))):
        r = stratified_permanova(dist, lab, st)
        rows.append({"matrix": name, "trait": trait, "test": test, "n_genera": r["n"], "n_groups": r["n_groups"],
                     "n_dropped_rare_group": int(rare.sum()), "n_dropped_no_stratum": r["n_dropped_no_stratum"],
                     "n_strata_used": r["n_strata_used"], "pseudo_F": r["F"], "R2": r["R2"], "R2_adj": r["R2_adj"],
                     "p_value": r["p"], "caveat": r["caveat"]})
    return rows


def run_global(reps: pd.DataFrame, mats: dict, traits: pd.DataFrame, out: Path, min_group: int = 5) -> dict:
    """Global PCoA of the genus representatives (same matrices, distances and PCoA as so.run_ordination), coloured
    by primary lifestyle and growth form, with lifestyle and growth-form PERMANOVA. Returns {matrix: coords}."""
    o = Out(out, TOPIC)
    o.scope = "genus_reps"
    tr = traits.drop_duplicates("species_prefix").set_index("species_prefix")
    ids = reps["species_prefix"]
    strata = reps.set_index("species_prefix")["ORDER"]
    cov, pooled, palettes = [], {}, {}
    for col, _, _ in TRAITS:
        lab = tr[col].reindex(ids.to_numpy()).rename(col)
        n_lab = int(lab.notna().sum())
        cov.append({"trait": col, "n_genera": len(lab), "n_labeled": n_lab, "n_unlabeled": len(lab) - n_lab})
        pooled[col] = se_pool = pool_groups(lab, min_group)
        palettes[col] = lifestyle_palette([g for g in se_pool.dropna().unique() if not str(g).startswith("other (")])
    o.table(pd.DataFrame(cov), "ecology_traits_coverage")
    coords_by, rows = {}, []
    for name in ("codon", "aa", "structure"):
        if name not in mats:
            continue
        mat, _, _ = so._feature_matrix(name, mats[name], reps)
        if len(mat) < 3:
            print(f"  {name}: {len(mat)} genera, too few for PCoA", flush=True)
            continue
        dist = so.lib_ordination.distance_matrix(mat, so.METRIC[name])
        ord_ = so.lib_ordination.pcoa(dist, n_axes=so.N_AXES)
        coords = ord_.coords.copy()
        coords.index = mat.index
        coords.index.name = "species_prefix"
        coords_by[name] = coords
        pct = [float(x) for x in ord_.pct_variance[:coords.shape[1]]]
        for col, short, nice in TRAITS:
            trows = _test_rows(name, col, dist, tr[col].reindex(mat.index).rename(col), strata, min_group)
            rows += trows
            lab = pooled[col].reindex(mat.index)
            stem = f"E0{1 if col == 'primary_lifestyle' else 2}_{name}_pcoa_{short}"
            r2 = trows[0]["R2"]
            o.meta[stem] = {"n": int(len(mat)), "n_groups": int(lab.nunique()), "signal_kind": "permanova_R2",
                            "signal": r2 if np.isfinite(r2) else None}
            title = f"{name} PCoA of {len(mat):,} genus representatives, coloured by {nice}"
            safe(o, stem, TOPIC, lambda lab=lab, col=col, title=title, coords=coords, pct=pct:
                 scatter_labels(coords, lab.rename(col), palettes[col], pct, title))
    o.table(pd.DataFrame(rows, columns=["matrix", "trait", "test", "n_genera", "n_groups", "n_dropped_rare_group",
                                        "n_dropped_no_stratum", "n_strata_used", "pseudo_F", "R2", "R2_adj", "p_value", "caveat"]),
            "ecology_permanova_global")
    print(f"ecology global: {len(coords_by)} matrices, {len(o.made)} outputs", flush=True)
    return coords_by


def clade_slug(level: str, name: str) -> str:
    """Lower-case file-name stem: level, then the name with non-alphanumeric runs turned into one '_'."""
    import re
    import unicodedata
    ascii_name = unicodedata.normalize("NFKD", str(name)).encode("ascii", "ignore").decode()
    return f"{level}_{re.sub(r'[^a-z0-9]+', '_', ascii_name.lower()).strip('_')}"


def clade_matrix(name: str, table: pd.DataFrame, reps: pd.DataFrame) -> pd.DataFrame:
    """Feature matrix of one clade's genera. Structure z-scoring is done on these genera only."""
    return so._feature_matrix(name, table, reps)[0]


CLADE_COLUMNS = ["level", "clade", "matrix", "test", "n_genera", "n_labeled", "n_groups", "n_dropped_rare_group",
                 "n_dropped_no_stratum", "n_strata_used", "pseudo_F", "R2", "R2_adj", "p_value", "q_BH", "caveat"]


def run_clades(reps: pd.DataFrame, mats: dict, traits: pd.DataFrame, clades: list, skipped: pd.DataFrame,
               out: Path, min_group: int = 5) -> pd.DataFrame:
    """For each clade and matrix: PCoA on the clade's genera (reps are not re-picked), a scatter coloured by
    lifestyle (E03_<slug>_<matrix>_pcoa_lifestyle), and unrestricted and stratified lifestyle PERMANOVA.
    BH q values are computed within (level, matrix, test). Returns the clade table."""
    o = Out(out, TOPIC)
    o.scope = "genus_reps"
    tr = traits.drop_duplicates("species_prefix").set_index("species_prefix")["primary_lifestyle"]
    by_id = reps.set_index("species_prefix")
    rows, skip = [], skipped.to_dict("records")
    for c in clades:
        sub = reps[reps["species_prefix"].isin(c.ids)]
        strata = by_id[STRATUM_COLUMN[c.level]]
        slug = clade_slug(c.level, c.name)
        for name in ("codon", "aa", "structure"):
            if name not in mats:
                continue
            mat = clade_matrix(name, mats[name], sub)
            lab = tr.reindex(mat.index).rename("primary_lifestyle")
            n_lab = int(lab.notna().sum())
            if len(mat) < 3:
                skip.append({"level": c.level, "clade": c.name, "n_labeled": n_lab, "n_groups_ge_min": 0,
                             "reason": f"{name}: fewer than 3 genera with data"})
                continue
            dist = so.lib_ordination.distance_matrix(mat, so.METRIC[name])
            ord_ = so.lib_ordination.pcoa(dist, n_axes=so.N_AXES)
            coords = ord_.coords.copy()
            coords.index = mat.index
            pct = [float(x) for x in ord_.pct_variance[:coords.shape[1]]]
            counts = lab.value_counts()
            rare = lab.isin(counts[counts < min_group].index)
            tested = lab.mask(rare)
            first = None
            for test, st in (("unrestricted", None), ("stratified", strata.reindex(mat.index))):
                r = stratified_permanova(dist, tested, st)
                first = first or r
                rows.append({"level": c.level, "clade": c.name, "matrix": name, "test": test, "n_genera": r["n"],
                             "n_labeled": n_lab, "n_groups": r["n_groups"], "n_dropped_rare_group": int(rare.sum()),
                             "n_dropped_no_stratum": r["n_dropped_no_stratum"], "n_strata_used": r["n_strata_used"],
                             "pseudo_F": r["F"], "R2": r["R2"], "R2_adj": r["R2_adj"], "p_value": r["p"], "q_BH": np.nan,
                             "caveat": r["caveat"]})
            pooled = pool_groups(lab, min_group)
            pal = lifestyle_palette([g for g in pooled.dropna().unique() if not str(g).startswith("other (")])
            stem = f"E03_{slug}_{name}_pcoa_lifestyle"
            o.meta[stem] = {"n": int(len(mat)), "n_groups": int(pooled.nunique()), "signal_kind": "permanova_R2",
                            "signal": first["R2"] if np.isfinite(first["R2"]) else None}
            title = f"{c.name} ({c.level}): {name} PCoA of {len(mat):,} genera, coloured by primary lifestyle"
            safe(o, stem, TOPIC, lambda coords=coords, pooled=pooled, pal=pal, pct=pct, title=title:
                 scatter_labels(coords, pooled, pal, pct, title))
    t = pd.DataFrame(rows, columns=CLADE_COLUMNS)
    for _, idx in t.groupby(["level", "matrix", "test"]).groups.items():
        t.loc[idx, "q_BH"] = bh_adjust(t.loc[idx, "p_value"]).to_numpy()
    o.table(t, "ecology_permanova_clades")
    o.table(pd.DataFrame(skip, columns=["level", "clade", "n_labeled", "n_groups_ge_min", "reason"]), "ecology_clades_skipped")
    print(f"ecology clades: {len(clades)} clades, {len(t)} test rows, {len(skip)} skipped", flush=True)
    return t
