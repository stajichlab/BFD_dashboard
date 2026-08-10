"""Per-dataset configuration for a BFD_dashboard build.

Every knob a generic dashboard needs lives here instead of being hardcoded in a
module or in the HTML template: which duckdb to read, where BUSCO summaries /
phyling results live, the explicit treefile or the glob that finds it, the
outgroup genus used to root the tree, the displayed title/subtitle, which feature
sets to include, and the optional taxon subset (taxonomy ranks + subsampling).

A build is driven by a ``DashboardConfig`` dataclass. It can be constructed
directly in code, read from a YAML file (see config/BFD_dashboard.yaml), or
assembled by bin/make_dashboard.py, where explicit CLI flags override YAML values
and YAML overrides built-in defaults. ``None`` means "unset" so layering is
unambiguous.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml


@dataclass
class SubsetConfig:
    """Which genomes to include, resolved against the species table.

    Rank filters are case-insensitive substring matches against the matching
    ``species`` column (the same behaviour as bfd_data.resolve_taxa_filter);
    combine any number of ranks. ``only_with_functional`` restricts to genomes
    present in at least one functional-annotation table (pfam/cazy/merops/
    signalp/tmhmm/targetp/idp/wolfpsort/predgpi) -- important against a large
    BFD database where functional tables may only be populated for some genomes.
    ``max_per_species`` subsamples to at most N strains per species (1 = one
    representative per species); ``max_genomes`` caps the total (random, seeded).
    """

    phylum: str | None = None
    subphylum: str | None = None
    klass: str | None = None  # 'class' is a Python keyword; column is CLASS
    order: str | None = None
    family: str | None = None
    genus: str | None = None
    species: str | None = None
    locustags: list[str] | None = None
    only_with_functional: bool = False
    max_per_species: int | None = None
    max_genomes: int | None = None
    seed: int = 42

    def is_empty(self) -> bool:
        """True when no filter/subsample would change the default (all genomes)."""
        return not any(
            [
                self.phylum, self.subphylum, self.klass, self.order, self.family,
                self.genus, self.species, self.locustags, self.only_with_functional,
                self.max_per_species, self.max_genomes,
            ]
        )

    def to_dict(self) -> dict:
        return {
            "phylum": self.phylum, "subphylum": self.subphylum, "class": self.klass,
            "order": self.order, "family": self.family, "genus": self.genus,
            "species": self.species, "locustags": self.locustags,
            "only_with_functional": self.only_with_functional,
            "max_per_species": self.max_per_species,
            "max_genomes": self.max_genomes, "seed": self.seed,
        }


@dataclass
class DashboardConfig:
    """Everything bin/make_dashboard.py needs for one build."""

    # --- data sources ---------------------------------------------------------
    db_path: Path
    busco_dir: Path | None = None        # glob for */*.BUSCO_summary.*.txt
    results_dir: Path | None = None      # base dir for the phyling tree glob
    tree_path: Path | None = None        # explicit treefile; else glob under results_dir
    merops_families_path: Path | None = None  # optional; merops skipped when absent

    # --- analysis -------------------------------------------------------------
    outgroup_genus: str | None = None    # used to root the tree; None = no reroot
    feature_sets: list[str] | None = None  # None = all registered feature sets
    alt_ordinations: bool = False        # add NMDS / t-SNE / UMAP per feature set
    per_species_permanova: bool | None = None  # None = auto (on when n <= cap)

    # --- display --------------------------------------------------------------
    title: str = "BFD Comparative Explorer"
    subtitle: str | None = None          # shown under the title; "N genomes" appended

    # --- subsetting -----------------------------------------------------------
    subset: SubsetConfig = field(default_factory=SubsetConfig)

    @classmethod
    def from_yaml(cls, path: Path | str) -> "DashboardConfig":
        """Load a config file. Keys mirror the dataclass fields; unknown keys are
        ignored so one config file can serve this tool and be the seed for a
        build-subset-db sibling."""
        raw = yaml.safe_load(Path(path).read_text()) or {}
        return cls.from_dict(raw)

    @classmethod
    def from_dict(cls, raw: dict) -> "DashboardConfig":
        sub = raw.get("subset") or {}
        subset = SubsetConfig(
            phylum=sub.get("phylum"),
            subphylum=sub.get("subphylum"),
            klass=sub.get("class") or sub.get("klass"),
            order=sub.get("order"),
            family=sub.get("family"),
            genus=sub.get("genus"),
            species=sub.get("species"),
            locustags=sub.get("locustags"),
            only_with_functional=bool(sub.get("only_with_functional", False)),
            max_per_species=sub.get("max_per_species"),
            max_genomes=sub.get("max_genomes"),
            seed=int(sub.get("seed", 42)),
        )
        db_path = _opt_path(raw.get("db_path") or raw.get("db"))
        if db_path is None:
            raise ValueError("config requires a 'db' / 'db_path' key")
        return cls(
            db_path=db_path,
            busco_dir=_opt_path(raw.get("busco_dir")),
            results_dir=_opt_path(raw.get("results_dir")),
            tree_path=_opt_path(raw.get("tree") or raw.get("tree_path")),
            merops_families_path=_opt_path(raw.get("merops_families_path")),
            outgroup_genus=raw.get("outgroup_genus"),
            feature_sets=raw.get("feature_sets"),
            alt_ordinations=bool(raw.get("alt_ordinations", False)),
            per_species_permanova=raw.get("per_species_permanova"),
            title=str(raw.get("title", cls.title)),
            subtitle=raw.get("subtitle"),
            subset=subset,
        )

    def to_dict(self) -> dict:
        return {
            "db_path": str(self.db_path),
            "busco_dir": str(self.busco_dir) if self.busco_dir else None,
            "results_dir": str(self.results_dir) if self.results_dir else None,
            "tree_path": str(self.tree_path) if self.tree_path else None,
            "merops_families_path": str(self.merops_families_path) if self.merops_families_path else None,
            "outgroup_genus": self.outgroup_genus,
            "feature_sets": self.feature_sets,
            "alt_ordinations": self.alt_ordinations,
            "per_species_permanova": self.per_species_permanova,
            "title": self.title,
            "subtitle": self.subtitle,
            "subset": {
                "phylum": self.subset.phylum,
                "subphylum": self.subset.subphylum,
                "class": self.subset.klass,
                "order": self.subset.order,
                "family": self.subset.family,
                "genus": self.subset.genus,
                "species": self.subset.species,
                "locustags": self.subset.locustags,
                "only_with_functional": self.subset.only_with_functional,
                "max_per_species": self.subset.max_per_species,
                "max_genomes": self.subset.max_genomes,
                "seed": self.subset.seed,
            },
        }


def _opt_path(v) -> Path | None:
    return Path(str(v)) if v else None
