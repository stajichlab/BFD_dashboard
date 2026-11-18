"""Calibration of BFD's Pfam-based TF counts against the paper's own InterPro-based counts (Lax et al. 2025, MOESM3), and the
paper's claims tested on the paper's species table.

The xlsx stays in resources/. Only derived TSVs are committed: config/paper_tf_counts.tsv (species x paper family),
config/paper_tf_totals.tsv (total TFs and total genes per species) and config/paper_species_map.tsv (paper species -> BFD genome).
Unit: one value per paper species (62 species), matched to the best-BUSCO BFD genome of the same genus and epithet.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

import tf_claims as tc

CONFIG = Path(__file__).resolve().parents[2] / "config"
SOURCE_LINE = ("# Source: Lax et al. 2025, Nature Communications 16:10277, Supplementary Data (MOESM3), sheet 'TF count' "
               "(counts) and sheet 'Genomes' (full species names, paired by column order). Whitespace in names stripped.")
# Paper species names with a typo or an old name, to the BFD name. Marked match = alias in the species map.
ALIAS_REASONS = {
    "Malasseizia sympodialis": "typo",
    "Ustilago maydis": "renamed taxon",
    "Saksenea vasiformis": "typo",
    "Mortierella elongata": "renamed taxon",
    "Synchephalis fuscata": "typo",
}
ALIASES = {
    "Malasseizia sympodialis": "Malassezia sympodialis",
    "Ustilago maydis": "Mycosarcoma maydis",
    "Saksenea vasiformis": "Saksenaea vasiformis",
    "Mortierella elongata": "Linnemannia elongata",
    "Synchephalis fuscata": "Syncephalis fuscata",
}
BACTERIAL = "Bacterial regulatory protein"


def convert_paper_counts(xlsx, out_tsv, totals_tsv=None):
    """Read sheets 'TF count' and 'Genomes' of MOESM3. Writes out_tsv (columns species, species_full, paper_family, count; a
    header comment gives the source) and totals_tsv (species, species_full, total_tfs, total_genes; default
    paper_tf_totals.tsv next to out_tsv). Returns the counts DataFrame. species = the paper's abbreviated name,
    species_full = the name in the 'Genomes' sheet (same column order; the first letters are checked)."""
    import openpyxl
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    rows = [list(r) for r in wb["TF count"].iter_rows(values_only=True)]
    full = [str(r[0]).strip() for r in wb["Genomes"].iter_rows(values_only=True, min_row=2) if r[0] is not None]
    wb.close()
    short = [str(s).strip() for s in rows[0][1:] if s is not None]
    if len(full) < len(short):
        raise ValueError(f"Genomes sheet has {len(full)} species, TF count has {len(short)}")
    for s, f in zip(short, full):
        if s[0].lower() != f[0].lower():
            raise ValueError(f"species order differs between sheets: '{s}' and '{f}'")
    fulls = dict(zip(short, full))
    recs, totals = [], {}
    for r in rows[1:]:
        name = None if r[0] is None else str(r[0]).strip()
        if not name:
            continue
        vals = r[1:1 + len(short)]
        if name in ("Total TFs", "Total Genes"):
            totals[name] = vals
            continue
        for s, v in zip(short, vals):
            recs.append((s, fulls[s], name, 0 if v is None else v))
    df = pd.DataFrame(recs, columns=["species", "species_full", "paper_family", "count"])
    df["count"] = pd.to_numeric(df["count"]).astype(int)
    tdf = pd.DataFrame({"species": short, "species_full": [fulls[s] for s in short],
                        "total_tfs": [int(v) for v in totals["Total TFs"][:len(short)]],
                        "total_genes": [int(v) for v in totals["Total Genes"][:len(short)]]})
    out_tsv = Path(out_tsv)
    totals_tsv = Path(totals_tsv) if totals_tsv else out_tsv.parent / "paper_tf_totals.tsv"
    for d, p in ((df, out_tsv), (tdf, totals_tsv)):
        with open(p, "w") as fh:
            fh.write(SOURCE_LINE + "\n")
            d.to_csv(fh, sep="\t", index=False)
    return df


def load_paper_counts(counts_tsv=None, totals_tsv=None):
    """(counts, totals) DataFrames from the derived TSVs (default: the files in config/)."""
    counts_tsv = Path(counts_tsv) if counts_tsv else CONFIG / "paper_tf_counts.tsv"
    totals_tsv = Path(totals_tsv) if totals_tsv else counts_tsv.parent / "paper_tf_totals.tsv"
    return (pd.read_csv(counts_tsv, sep="\t", comment="#", keep_default_na=False),
            pd.read_csv(totals_tsv, sep="\t", comment="#", keep_default_na=False))


def _genus_epithet(name):
    """(list of genera, epithet) of a paper name: 'Ashbya(Eremothecium) gossypii' -> (['Ashbya', 'Eremothecium'], 'gossypii').
    A one-word name has no epithet (None)."""
    s = str(name).strip()
    parts = s.split()
    if len(parts) < 2:
        return [], None
    head = " ".join(parts[:-1])
    gen = [g.strip() for g in head.replace(")", "").split("(") if g.strip()]
    return gen, parts[-1]


def match_species(paper_species, genome_all, aliases=None):
    """For each paper species name (full name from the 'Genomes' sheet) the best BFD genome with the same genus and epithet:
    highest busco_complete, then highest N50_bp, then lowest species_prefix. A name with a genus in brackets is tried with
    both genera. When nothing matches, `aliases` (paper name -> BFD 'Genus epithet') is tried (match = alias). Columns:
    paper_species, species_prefix, strain, match (exact, alias, none), n_candidates, bfd_species, busco_complete,
    alias_reason (typo, renamed taxon or synonym for an alias; empty otherwise)."""
    aliases = ALIASES if aliases is None else aliases
    ga = genome_all.copy()
    sp = ga["SPECIES"].fillna("").astype(str).str.split()
    ga["_epi"] = sp.str[1]
    n50 = "N50_bp" if "N50_bp" in ga.columns else None
    recs = []
    for name in paper_species:
        gen, epi = _genus_epithet(name)
        how, cand = "none", ga.iloc[0:0]
        if epi:
            cand = ga[ga["GENUS"].isin(gen) & (ga["_epi"] == epi)]
            how = "exact" if len(cand) else "none"
        if not len(cand) and str(name).strip() in aliases:
            g2, e2 = _genus_epithet(aliases[str(name).strip()])
            cand = ga[ga["GENUS"].isin(g2) & (ga["_epi"] == e2)]
            how = "alias" if len(cand) else "none"
        if not len(cand):
            recs.append((name, None, "", "none", 0, "", np.nan, ""))
            continue
        by = ["busco_complete"] + ([n50] if n50 else []) + ["species_prefix"]
        best = cand.sort_values(by, ascending=[False] + [False] * (len(by) - 2) + [True], na_position="last").iloc[0]
        strain = str(best["STRAIN"]) if "STRAIN" in cand.columns and pd.notna(best.get("STRAIN")) else ""
        reason = ALIAS_REASONS.get(str(name).strip(), "synonym") if how == "alias" else ""
        recs.append((name, best["species_prefix"], strain, how, len(cand), best["SPECIES"], best["busco_complete"], reason))
    return pd.DataFrame(recs, columns=["paper_species", "species_prefix", "strain", "match", "n_candidates", "bfd_species",
                                       "busco_complete", "alias_reason"])


def _bridge(fmap):
    fm = fmap[(fmap["paper_family"] != "") & (fmap["family"] != "")][["family", "paper_family"]].drop_duplicates()
    return fm


def calibrate(bfd_family_genome, paper_counts, species_map, fmap, bfd_paper_family=None, genus_nonzero=None):
    """Compare BFD and paper counts per paper family over the matched species. Returns (summary, per_species).
    summary: paper_family, n_species, n_nonzero_paper, n_nonzero_bfd (matched species with a count above 0), spearman,
    pearson_log (Pearson on log10(count + 1)), median_ratio (BFD / paper over species with a paper count above 0), note.
    One row for each paper family that has a Pfam in the family map and a row in the paper table. A family with no BFD rows
    at all gets n_species 0 and a note. A family is "not comparable" when it is non-zero in none of the matched species or
    (genus_nonzero: {paper family: number of genera with a protein}) in fewer than MIN_NONZERO_GENERA genera: spearman,
    pearson_log and median_ratio are NaN and the note says so.
    BFD counts: bfd_paper_family (tf_paper_family_genome: species_prefix, paper_family, n_proteins; distinct proteins) when
    given; otherwise the counts of the BFD families of one paper family are added (a note says double counting is possible).
    per_species: paper_species, species, species_prefix, paper_family, paper_count, bfd_count, ratio. A species or family
    absent from the BFD table counts as 0."""
    br = _bridge(fmap)
    nlab = br.groupby("paper_family")["family"].nunique()
    sm = species_map[(species_map["match"] != "none") & species_map["species_prefix"].notna() &
                     (species_map["species_prefix"].astype(str) != "")][["paper_species", "species_prefix"]]
    if bfd_paper_family is not None:
        fam = bfd_paper_family[bfd_paper_family["paper_family"].isin(set(br["paper_family"]))]
    else:
        fam = bfd_family_genome.merge(br, on="family", how="inner")
    wide = fam.groupby(["species_prefix", "paper_family"])["n_proteins"].sum()
    pc = paper_counts[paper_counts["paper_family"].isin(set(br["paper_family"]))]
    per = pc.merge(sm, left_on="species_full", right_on="paper_species", how="inner")
    per["bfd_count"] = [int(wide.get((s, f), 0)) for s, f in zip(per["species_prefix"], per["paper_family"])]
    per = per.rename(columns={"count": "paper_count"})
    per["ratio"] = np.where(per["paper_count"] > 0, per["bfd_count"] / per["paper_count"].where(per["paper_count"] > 0), np.nan)
    per = per[["paper_species", "species", "species_prefix", "paper_family", "paper_count", "bfd_count", "ratio"]]
    has_bfd = set(fam["paper_family"])
    rows = []
    for f in sorted(set(pc["paper_family"])):
        d = per[per["paper_family"] == f]
        npp, nb = int((d["paper_count"] > 0).sum()), int((d["bfd_count"] > 0).sum())
        sumnote = tc.SUM_NOTE if bfd_paper_family is None and nlab.get(f, 1) > 1 else ""
        if f not in has_bfd:
            rows.append((f, 0, npp, 0, np.nan, np.nan, np.nan, "no BFD rows for this family (family absent from the BFD table)"))
            continue
        x, y = d["paper_count"].to_numpy(float), d["bfd_count"].to_numpy(float)
        k = genus_nonzero.get(f) if genus_nonzero else None
        if nb == 0 or (k is not None and k < tc.MIN_NONZERO_GENERA):
            note = tc.not_comparable_note(k) if k is not None else \
                f"BFD Pfam definition does not detect this family (non-zero in {nb} matched species); not comparable with the paper"
            rows.append((f, len(d), npp, nb, np.nan, np.nan, np.nan, tc._join(note, sumnote)))
            continue
        sp = pl = np.nan
        if len(d) >= 3 and np.ptp(x) > 0 and np.ptp(y) > 0:
            sp = float(st.spearmanr(x, y).statistic)
            pl = float(st.pearsonr(np.log10(x + 1), np.log10(y + 1)).statistic)
        rt = d["ratio"].dropna()
        rows.append((f, len(d), npp, nb, sp, pl, float(rt.median()) if len(rt) else np.nan, sumnote))
    summ = pd.DataFrame(rows, columns=["paper_family", "n_species", "n_nonzero_paper", "n_nonzero_bfd", "spearman", "pearson_log",
                                       "median_ratio", "note"])
    return summ, per.reset_index(drop=True)


def evaluate_claims_on_paper(paper_counts, species_phylum, claims, totals, n_boot=2000, seed=1, fmap=None, not_comparable=None):
    """The claims on the paper's own table: one value per species (not per genus). fraction = count / Total Genes. Phylum
    comes from species_phylum (columns species_full, PHYLUM: the BFD taxonomy of the matched genome); species without a phylum
    are left out. Shannon index over the paper families that have a Pfam in the family map (all non-bacterial families when
    fmap is None). Same columns as tf_claims.evaluate_claims, source = 'paper', n_genera_* = number of species. Claims on a
    family in `not_comparable` (a set) keep their result and get the note 'paper definition; not reproduced by BFD Pfam'."""
    wide = paper_counts.pivot_table(index="species_full", columns="paper_family", values="count", aggfunc="sum", fill_value=0)
    tot = totals.drop_duplicates("species_full").set_index("species_full")["total_genes"].astype(float)
    ph = species_phylum.dropna(subset=["PHYLUM"]).drop_duplicates("species_full").set_index("species_full")["PHYLUM"]
    idx = [s for s in wide.index if s in ph.index]
    wide, tot, ph = wide.loc[idx], tot.reindex(idx), ph.reindex(idx)
    df = pd.DataFrame({"genus": idx, "PHYLUM": ph.to_numpy()})
    if fmap is not None:
        keep = set(_bridge(fmap)["paper_family"])
    else:
        keep = {c for c in wide.columns if not c.startswith(BACTERIAL)}
    cols = [c for c in wide.columns if c in keep]
    p = wide[cols].to_numpy(float)
    s = p.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        q = p / s[:, None]
        sh = -np.nansum(np.where(q > 0, q * np.log(q), 0.0), axis=1)
    sh = np.where(s > 0, sh, np.nan)

    def value_fn(family, metric):
        if metric == "shannon":
            return sh, ""
        if family not in wide.columns:
            return None, f"family '{family}' is not a row of the paper table"
        v = wide[family].to_numpy(float)
        return (v / tot.to_numpy(float) if metric == "fraction" else v), ""

    nc = set(not_comparable or ())
    flag_fn = (lambda f: (False, "paper definition; not reproduced by BFD Pfam") if f in nc else None)
    return tc.evaluate_table(df, claims, value_fn, n_boot, seed, "paper", "species", flag_fn, "n/a (paper table)")


def _calibration_fig(summ, per):
    import matplotlib.pyplot as plt
    fams = [f for f in summ.sort_values("n_species", ascending=False)["paper_family"]
            if (per["paper_family"] == f).sum() >= 1]
    if not fams:
        return None
    nc = 5
    nr = int(np.ceil(len(fams) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(2.4 * nc, 2.3 * nr + 0.8), layout="constrained", squeeze=False)
    s = summ.set_index("paper_family")
    for ax, f in zip(axs.ravel(), fams):
        d = per[per["paper_family"] == f]
        x, y = np.log10(d["paper_count"] + 1), np.log10(d["bfd_count"] + 1)
        ax.scatter(x, y, s=10, color="#0072B2", alpha=0.8, edgecolor="none")
        m = float(max(x.max(), y.max(), 0.5)) * 1.05
        ax.plot([0, m], [0, m], color="black", lw=0.7)
        ax.set_xlim(-0.05, m)
        ax.set_ylim(-0.05, m)
        sp = s.loc[f, "spearman"]
        ax.set_title(f"{f[:32]}\nn={int(s.loc[f, 'n_species'])}, rho={'NA' if not np.isfinite(sp) else f'{sp:.2f}'}", fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in axs.ravel()[len(fams):]:
        ax.set_visible(False)
    fig.supxlabel("Paper count (InterPro), log10(count + 1)", fontsize=8)
    fig.supylabel("BFD count (Pfam), log10(count + 1)", fontsize=8)
    fig.suptitle("Calibration: BFD Pfam counts against the paper's counts, one point per matched species (line = identity)",
                 fontsize=8)
    return fig


def page_tf_calibration(bfd_family_genome, paper_counts, species_map, fmap, o, bfd_paper_family=None, genus_nonzero=None):
    """T09_tf_calibration.tsv, T09_tf_calibration_species.tsv and the figure T09_tf_calibration. o: an Out for 'tf'.
    Returns the summary table."""
    from figs_common import safe
    summ, per = calibrate(bfd_family_genome, paper_counts, species_map, fmap, bfd_paper_family, genus_nonzero)
    o.table(summ, "T09_tf_calibration")
    o.table(per, "T09_tf_calibration_species")
    safe(o, "T09_tf_calibration", "tf", lambda: _calibration_fig(summ, per))
    return summ


def run_tf_claims(genome_all, fam, summ, out, busco_min=80, n_boot=2000, config_dir=None, paper_fam=None):
    """The claims page and the calibration page under tables/tf and figures/tf. BFD claims use genus representatives with
    busco_complete >= busco_min; T10_tf_claims_nofilter.tsv uses all genera. When the derived paper tables exist in
    config_dir (default config/), the claims table also has the paper's own rows (source = paper) and the calibration is
    written; otherwise only the BFD claims are. paper_fam: tf_paper_family_genome (distinct proteins per paper family); when
    None the counts of the BFD families are added with a printed warning and a note on each row. Returns the Made list."""
    from make_figures import Out
    from tf_data import genus_table
    cfg = Path(config_dir) if config_dir else CONFIG
    o = Out(out, "tf")
    o.scope = "genus_reps"
    claims = tc.load_claims(cfg / "tf_claims.tsv")
    fmap = tc.load_fmap(cfg / "tf_family_map.tsv")
    if paper_fam is None:
        print("WARNING: tf_paper_family_genome absent: paper-family counts are summed over Pfam families "
              "(double counting possible)", flush=True)
    g80 = genus_table(genome_all, fam, summ, busco_min=busco_min, paper_fam=paper_fam)
    g0 = genus_table(genome_all, fam, summ, paper_fam=paper_fam)
    nz = tc.nonzero_genera(g0, fmap)
    not_comp = {f for f, k in nz.items() if k < tc.MIN_NONZERO_GENERA}
    o.info = {"n_genomes": len(g80), "n_groups": None, "pooled_share": 0.0}
    paper_rows, calib = None, None
    cf, tf_, mf_ = cfg / "paper_tf_counts.tsv", cfg / "paper_tf_totals.tsv", cfg / "paper_species_map.tsv"
    if cf.exists() and tf_.exists() and mf_.exists():
        counts, totals = load_paper_counts(cf, tf_)
        smap = pd.read_csv(mf_, sep="\t", comment="#", keep_default_na=False, dtype=str)
        smap["species_prefix"] = smap["species_prefix"].replace("", np.nan)
        ph = smap.merge(genome_all[["species_prefix", "PHYLUM"]], on="species_prefix", how="inner")
        ph = ph.rename(columns={"paper_species": "species_full"})[["species_full", "PHYLUM"]]
        paper_rows = evaluate_claims_on_paper(counts, ph, claims, totals, n_boot=n_boot, fmap=fmap, not_comparable=not_comp)
        smap_ok = smap[smap["species_prefix"].isin(set(genome_all["species_prefix"]))]   # drops suppressed genomes
        calib = (fam, counts, smap_ok)
    tc.page_tf_claims(g80, claims, o, paper_rows=paper_rows, genus_df_nofilter=g0, fmap=fmap, n_boot=n_boot)
    if calib is not None:
        page_tf_calibration(calib[0], calib[1], calib[2], fmap, o, bfd_paper_family=paper_fam, genus_nonzero=nz)
    return o.made
