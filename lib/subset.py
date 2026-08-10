"""Resolve a taxon subset + subsampling plan into a concrete LOCUSTAG list.

The whole point of a generic BFD dashboard is to run on any subset of a large BFD
database (23k+ genomes), so this is where a build decides WHICH genomes get views.
Two independent mechanisms, composed in order:

1. Taxonomy filters -- any combination of phylum / subphylum / class / order /
   family / genus / species (case-insensitive substring match on the species table,
   same convention as bfd_data.resolve_taxa_filter), plus an explicit locustag
   allowlist and a "only_with_functional" flag that restricts to genomes present in
   at least one functional-annotation table.

2. Subsampling -- ``max_per_species`` (at most N strains per species; 1 = one
   representative per species) and/or ``max_genomes`` (a seeded random cap on the
   total). Subsampling is deterministic for a given seed so a rebuild reproduces the
   same set. Representative picks prefer the best-assembled strain (highest gene
   count, then N50) so a "one per species" panel shows the most complete genomes,
   not whichever row duckdb happened to return first.

Also reports what was dropped/why, so a build tells the user e.g. "genus filter kept
412 genomes, only 38 of which have functional annotation".
"""
from __future__ import annotations

import sys

import pandas as pd

import bfd_data as bd
from config import SubsetConfig

_RANK_COLUMNS = [
    ("phylum", "PHYLUM"),
    ("subphylum", "SUBPHYLUM"),
    ("klass", "CLASS"),
    ("order", "ORDER"),
    ("family", "FAMILY"),
    ("genus", "GENUS"),
    ("species", "SPECIES"),
]


def resolve_subset(con, cfg: SubsetConfig) -> tuple[list[str], dict]:
    """Return (locustags_in_subset_order, report) for a subset config.

    The report dict is human-readable accounting (n_before, filters applied, n
    dropped by subsampling, functional coverage) for the build's stderr/log.
    """
    report: dict = {"n_before": None, "steps": [], "n_after": None}

    sp = con.execute(
        'SELECT LOCUSTAG, PHYLUM, SUBPHYLUM, CLASS, "ORDER", FAMILY, GENUS, SPECIES '
        'FROM species ORDER BY LOCUSTAG'
    ).fetchdf()
    report["n_before"] = len(sp)

    mask = pd.Series(True, index=sp.index)
    for field, column in _RANK_COLUMNS:
        value = getattr(cfg, field)
        if value:
            before = int(mask.sum())
            mask &= sp[column].str.contains(value, case=False, na=False)
            report["steps"].append(
                f"{field} ~ '{value}': {before} -> {int(mask.sum())}"
            )
    if cfg.locustags:
        wanted = set(cfg.locustags)
        before = int(mask.sum())
        mask &= sp["LOCUSTAG"].isin(wanted)
        report["steps"].append(f"locustags allowlist: {before} -> {int(mask.sum())}")

    if cfg.only_with_functional:
        func = bd.functional_species(con)
        before = int(mask.sum())
        mask &= sp["LOCUSTAG"].isin(func)
        report["steps"].append(
            f"only_with_functional: {before} -> {int(mask.sum())} "
            f"(of {len(func)} functional genomes in db)"
        )

    sub = sp.loc[mask]
    locustags = sub["LOCUSTAG"].tolist()

    # ---- subsampling ----------------------------------------------------------
    if cfg.max_per_species:
        locustags = _representatives_per_species(
            con, locustags, max_per_species=cfg.max_per_species
        )
        report["steps"].append(
            f"max_per_species={cfg.max_per_species}: -> {len(locustags)}"
        )

    if cfg.max_genomes:
        before = len(locustags)
        locustags = _seeded_cap(locustags, cfg.max_genomes, cfg.seed)
        report["steps"].append(f"max_genomes={cfg.max_genomes}: {before} -> {len(locustags)}")

    report["n_after"] = len(locustags)
    return locustags, report


def _representatives_per_species(
    con, locustags: list[str], max_per_species: int
) -> list[str]:
    """Keep at most max_per_species strains per species. 'Best' strain = highest
    gene_count, then N50 (from get_species_table, which falls back to asm_stats when
    the annotation summary is sparse) -- stable picks that don't depend on table
    order. Non-picked strains within a species are dropped deterministically."""
    if max_per_species < 1:
        raise ValueError("max_per_species must be >= 1")

    sp = con.execute(
        "SELECT LOCUSTAG, SPECIES, GENUS FROM species ORDER BY LOCUSTAG"
    ).fetchdf()
    sp = sp.set_index("LOCUSTAG").reindex(locustags)

    stats = bd.get_species_table(con, rows=locustags)

    chosen: list[str] = []
    for name, group in sp.groupby("SPECIES", dropna=False):
        group = group.index.tolist()
        if len(group) <= max_per_species:
            chosen.extend(group)
            continue
        scored = []
        for lt in group:
            if lt in stats.index:
                row = stats.loc[lt]
                gc = _num(row["gene_count"]) if "gene_count" in stats else None
                nc = _num(row["N50_bp"]) if "N50_bp" in stats else None
            else:
                gc = nc = None
            scored.append((_best_key(gc, nc), lt))
        # Stable within the species (tie-break by locustag), not RNG-dependent.
        scored.sort(key=lambda t: (t[0], t[1]), reverse=True)
        chosen.extend(lt for _, lt in scored[:max_per_species])
    return chosen


def _best_key(gene_count, n50):
    return (gene_count if gene_count is not None else -1,
            n50 if n50 is not None else -1)


def _num(v) -> float | None:
    try:
        if v is None or (isinstance(v, float) and v != v):  # NaN
            return None
        return float(v)
    except (TypeError, ValueError):
        return None


def _seeded_cap(locustags: list[str], max_genomes: int, seed: int) -> list[str]:
    import random

    if len(locustags) <= max_genomes:
        return locustags
    rng = random.Random(seed)
    chosen = rng.sample(locustags, max_genomes)
    return sorted(chosen)  # stable, reproducible order


def report_subset(rows: list[str], report: dict, out=sys.stderr) -> None:
    print(f"  subset: {report['n_before']} genomes in db -> {report['n_after']} in build", file=out)
    for step in report["steps"]:
        print(f"    {step}", file=out)
