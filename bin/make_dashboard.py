#!/usr/bin/env python3
"""Build the BFD comparative-genomics explorer dashboard for any subset of a BFD duckdb.

Two explicit steps (compute is not free -- PERMANOVA runs per feature set per focus
species):

    make_dashboard.py --compute   # duckdb -> cache/payload.json
    make_dashboard.py --render    # cache/payload.json -> explorer.html

Running with neither flag does both, for convenience on a first build.

Subsetting: pass --config config/BFD_dashboard.yaml (or the individual --genus /
--order / --family / --phylum / --species / --only-with-functional /
--max-per-species / --max-genomes flags) to restrict the build to a taxon subset,
including subsampling to e.g. one representative per species. The resolved subset is
printed to stderr before any heavy analysis runs.

The phylogeny is optional and per-dataset: --tree gives an explicit phyling treefile,
otherwise --results-dir is globbed for one (see bfd_data.find_tree_path). When a tree
is found it is pruned to the subset and used for species ordering; without one the
build falls back to taxonomic order.
"""
import argparse
import dataclasses
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

DASHBOARD_DIR = Path(__file__).resolve().parent.parent
for _sub in ("lib", "templates"):
    sys.path.insert(0, str(DASHBOARD_DIR / _sub))

import bfd_data as bd  # noqa: E402
import ordination as ordi  # noqa: E402
import figures as fig  # noqa: E402
import subset as submod  # noqa: E402
from config import DashboardConfig  # noqa: E402

DIFF_TOP_N = 15
# Above this many genomes, per-species PERMANOVA ("species vs rest" for every focus)
# is O(n) full PERMANOVAs and stops being worth the wall-clock time; it is then
# skipped unless explicitly requested. The genus-level test always runs.
PER_SPECIES_PERMANOVA_CAP = 40


def _available_feature_sets(con, cfg: DashboardConfig) -> dict:
    """Registered feature sets, minus any whose underlying tables/aux files are
    missing for this database (e.g. merops without the MEROPS family-lookup file)."""
    fsets = bd.make_feature_sets(cfg.merops_families_path)
    tables = set(bd.list_tables(con))
    if "merops" in fsets and "merops" not in tables:
        print("  feature set 'merops': table not present, skipped", file=sys.stderr)
        del fsets["merops"]
    if "merops" in fsets:
        try:
            bd._ensure_merops_family_view(con, cfg.merops_families_path)
        except bd.MeropsUnavailableError as e:
            print(f"  feature set 'merops': {e}; skipped", file=sys.stderr)
            del fsets["merops"]
    if cfg.feature_sets:
        wanted = set(cfg.feature_sets)
        missing = sorted(wanted - set(fsets))
        if missing:
            print(f"  requested feature sets not available: {missing}", file=sys.stderr)
        fsets = {k: v for k, v in fsets.items() if k in wanted}
    return fsets


def build_payload(cfg: DashboardConfig) -> dict:
    con = bd.get_connection(cfg.db_path)

    # ---- subset resolution ----------------------------------------------------
    if not cfg.subset.is_empty():
        rows, report = submod.resolve_subset(con, cfg.subset)
        submod.report_subset(rows, report)
    else:
        rows = None
        print(f"  subset: all {len(bd.get_species_table(con))} genomes in {cfg.db_path.name}", file=sys.stderr)

    if rows is not None and len(rows) == 0:
        sys.exit(
            "ERROR: subset resolved to 0 genomes -- check the taxonomy filters / "
            "--only-with-functional against this database (see feature coverage in the "
            "report above). Nothing to build."
        )
    if rows is not None and len(rows) == 1:
        print("  WARNING: subset has a single genome -- ordinations will be empty; "
              "widen the subset for comparative views.", file=sys.stderr)

    feature_sets = _available_feature_sets(con, cfg)
    print(f"  feature sets: {', '.join(feature_sets)}", file=sys.stderr)

    # ---- species block --------------------------------------------------------
    species_df = bd.get_species_table(con, rows=rows)
    busco_df = bd.get_busco_scores(cfg.busco_dir)
    species_df = species_df.join(busco_df)

    phylo_payload = None
    try:
        tree = bd.load_tree(
            con, tree_path=cfg.tree_path, outgroup_genus=cfg.outgroup_genus,
            results_dir=cfg.results_dir, rows=rows,
        )
        phylo_payload = bd.phylo_layout(tree)
        sp_order = phylo_payload["tip_order"]
        print(f"  phylogeny: loaded, using tree tip order for species ordering", file=sys.stderr)
    except Exception as e:
        print(f"  phylogeny: unavailable ({type(e).__name__}: {e}); falling back to taxonomic order", file=sys.stderr)
        sp_order = species_df.index.tolist()

    species_df = species_df.reindex(sp_order)
    species_block = {"locustag": sp_order}
    for col in species_df.columns:
        values = species_df[col]
        if values.dtype.kind in "fc":
            species_block[col] = [None if _pd_isna(v) else float(v) for v in values]
        elif values.dtype.kind in "iu":
            species_block[col] = [None if _pd_isna(v) else int(v) for v in values]
        else:
            species_block[col] = [None if _pd_isna(v) else str(v) for v in values]

    genus = species_df["GENUS"]
    genus_counts = genus.value_counts().to_dict()

    # ---- per-feature-set analysis ----------------------------------------------
    feature_sets_payload = {}
    differential_payload = {}
    composition_payload = {}
    domain_counts_payload = {}

    n_here = len(sp_order)
    do_per_species = (
        cfg.per_species_permanova
        if cfg.per_species_permanova is not None
        else n_here <= PER_SPECIES_PERMANOVA_CAP
    )
    if not do_per_species and n_here > PER_SPECIES_PERMANOVA_CAP:
        print(
            f"  per-species PERMANOVA skipped for n={n_here} (cap {PER_SPECIES_PERMANOVA_CAP}); "
            f"pass --per-species-permanova to force",
            file=sys.stderr,
        )

    for name, fcfg in feature_sets.items():
        print(f"  feature set: {name}", file=sys.stderr)
        matrix = fcfg["loader"](con, rows=rows)

        # A large BFD db may have assembly/BUSCO data for every genome while the
        # functional tables (codon/aa/domain/localization) are only populated for a
        # fraction that has finished the pipeline. An empty or all-zero matrix for
        # this subset is "no data", not a degenerate ordination -- skip it cleanly
        # so the dashboard simply doesn't list that view.
        if matrix is None or matrix.shape[1] == 0 or not matrix.any().any():
            print(
                f"    no data for this subset (shape {matrix.shape if matrix is not None else '?'}); skipping",
                file=sys.stderr,
            )
            continue

        matrix = matrix.reindex(sp_order).fillna(0.0)

        if fcfg["method"] == "ca":
            ordination = ordi.correspondence_analysis(matrix)
            dist_for_stats = ordi.distance_matrix(matrix, "euclidean")
        else:
            dist_for_stats = ordi.distance_matrix(matrix, fcfg["metric"])
            ordination = ordi.pcoa(dist_for_stats)

        permanova_genus = dataclasses.asdict(ordi.permanova(dist_for_stats, genus))

        permanova_by_species = {}
        if do_per_species:
            for focus in sp_order:
                grouping = genus.index.to_series().apply(lambda x: "focus" if x == focus else "rest")
                permanova_by_species[focus] = dataclasses.asdict(
                    ordi.permanova(dist_for_stats, grouping)
                )

        coords = ordination.coords.reindex(sp_order)
        ordinations = {
            fcfg["method"]: {
                "pct_variance": [float(x) for x in ordination.pct_variance[:4]],
                "axes": list(coords.columns),
                "coords": {c: [float(v) for v in coords[c]] for c in coords.columns},
                "extra": ordination.extra,
            }
        }

        if cfg.alt_ordinations and fcfg["method"] != "ca":
            for alt_method, alt_fn in [("nmds", ordi.nmds), ("tsne", ordi.tsne)]:
                try:
                    alt_ord = alt_fn(dist_for_stats)
                    alt_coords = alt_ord.coords.reindex(sp_order)
                    ordinations[alt_method] = {
                        "pct_variance": [float(x) for x in alt_ord.pct_variance[:4]],
                        "axes": list(alt_coords.columns),
                        "coords": {c: [float(v) for v in alt_coords[c]] for c in alt_coords.columns},
                        "extra": alt_ord.extra,
                    }
                except Exception as e:
                    print(f"    {alt_method} failed: {e}", file=sys.stderr)
            try:
                umap_ord = ordi.umap_ordination(dist_for_stats)
                umap_coords = umap_ord.coords.reindex(sp_order)
                ordinations["umap"] = {
                    "pct_variance": [float(x) for x in umap_ord.pct_variance[:4]],
                    "axes": list(umap_coords.columns),
                    "coords": {c: [float(v) for v in umap_coords[c]] for c in umap_coords.columns},
                    "extra": umap_ord.extra,
                }
            except ImportError:
                print("    umap: not installed (optional; pip install umap-learn)", file=sys.stderr)
            except Exception as e:
                print(f"    umap failed: {e}", file=sys.stderr)

        feature_sets_payload[name] = {
            "description": fcfg.get("description", ""),
            "metric": fcfg["metric"],
            "primary_method": fcfg["method"],
            "ordinations": ordinations,
            "permanova_genus": permanova_genus,
            "permanova_by_species": permanova_by_species,
        }

        if name in ("pfam", "cazy", "merops"):
            try:
                differential_payload[name] = bd.differential_families(
                    con, name, top_n=DIFF_TOP_N, rows=rows,
                    merops_families_path=cfg.merops_families_path,
                )
            except ValueError as e:
                print(f"    differential ({name}): {e}", file=sys.stderr)

            # Raw per-species counts for exactly the families that surface in the
            # differential table -- not every family in the domain library, which
            # would be thousands of mostly-irrelevant columns. Backs the
            # "counts across all species" popup for a clicked family.
            referenced = {
                row["family"]
                for per_focus in differential_payload.get(name, {}).values()
                for section in ("present_only", "absent_only", "top_fold_change")
                for row in per_focus[section]
            }
            if referenced:
                raw_counts = bd.get_domain_raw_counts(
                    con, name, merops_families_path=cfg.merops_families_path
                ).reindex(sp_order).fillna(0.0)
                domain_counts_payload[name] = {
                    family: {
                        locustag: int(raw_counts.loc[locustag, family]) if family in raw_counts.columns else 0
                        for locustag in sp_order
                    }
                    for family in sorted(referenced)
                }

        if name in ("codon", "amino_acid"):
            raw = fcfg["loader"](con, rows=rows).reindex(sp_order).fillna(0.0)
            composition_payload[name] = {
                "categories": list(raw.columns),
                "values": {locustag: [float(v) for v in raw.loc[locustag]] for locustag in sp_order},
            }

        if name == "localization":
            # Composition wants the human-readable fractions (SignalP secretion, IDP
            # disorder, etc.), not the z-scored version fed to PCoA/PERMANOVA -- plus
            # the raw counts behind each fraction, so bars/tooltips show "412 (6.2%)".
            frac = bd.get_localization_fractions(con).reindex(sp_order).fillna(0.0)
            counts = bd.get_localization_counts(con).reindex(sp_order).fillna(0)
            composition_payload[name] = {
                "categories": list(frac.columns),
                "values": {locustag: [float(v) for v in frac.loc[locustag]] for locustag in sp_order},
                "counts": {locustag: [int(v) for v in counts.loc[locustag]] for locustag in sp_order},
            }

    # ---- assemble --------------------------------------------------------------
    func_n = len(bd.functional_species(con) & set(sp_order))
    notes = bd.dataset_notes(n_here, genus_counts, func_n)
    subtitle = cfg.subtitle or f"{n_here} genomes"
    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "db_path": str(cfg.db_path),
        "title": cfg.title,
        "subtitle": subtitle,
        "notes": notes,
        "subset": cfg.subset.to_dict() if not cfg.subset.is_empty() else None,
        "palette": fig.PALETTE_JSON,
        "species": species_block,
        "feature_sets": feature_sets_payload,
        "differential": differential_payload,
        "composition": composition_payload,
        "domain_counts": domain_counts_payload,
        "phylogeny": phylo_payload,
    }


def _pd_isna(v) -> bool:
    try:
        import math

        return v is None or (isinstance(v, float) and math.isnan(v))
    except Exception:
        return v is None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=DASHBOARD_DIR / "config" / "BFD_dashboard.yaml",
                    help="YAML config (defaults to config/BFD_dashboard.yaml if present; "
                         "CLI flags override its values).")
    ap.add_argument("--db", type=Path, help="BFD.duckdb path (overrides config)")
    ap.add_argument("--busco-dir", type=Path, help="dir globbed for */*.BUSCO_summary.*.txt (overrides config)")
    ap.add_argument("--results-dir", type=Path, help="base dir globbed for the phyling tree (overrides config)")
    ap.add_argument("--tree", type=Path, help="explicit phyling treefile (overrides config)")
    ap.add_argument("--outgroup-genus", type=str, help="GENUS used to root the tree (overrides config)")
    ap.add_argument("--merops-families-path", type=Path, help="MEROPS families.tab (overrides config)")
    ap.add_argument("--title", type=str, help="dashboard title (overrides config)")
    ap.add_argument("--subtitle", type=str, help="dashboard subtitle (overrides config)")
    ap.add_argument("--feature-sets", type=str, nargs="+",
                    help="restrict feature sets to these names (overrides config)")

    g = ap.add_argument_group("subset")
    g.add_argument("--phylum", type=str)
    g.add_argument("--genus", type=str, help="keep genomes whose GENUS matches (substring, case-insensitive)")
    g.add_argument("--order", type=str)
    g.add_argument("--family", type=str)
    g.add_argument("--class", dest="klass", type=str, help="keep genomes whose CLASS matches")
    g.add_argument("--species", type=str)
    g.add_argument("--locustags", type=str, nargs="+", help="explicit LOCUSTAG allowlist")
    g.add_argument("--only-with-functional", action="store_true",
                   help="keep only genomes present in a functional-annotation table")
    g.add_argument("--max-per-species", type=int,
                   help="subsample to at most N strains per species (1 = one representative)")
    g.add_argument("--max-genomes", type=int, help="seeded random cap on total genomes")
    g.add_argument("--seed", type=int, help="RNG seed for subsampling (default 42)")

    ap.add_argument("--cache", type=Path, default=DASHBOARD_DIR / "cache" / "payload.json")
    ap.add_argument("--out", type=Path, default=DASHBOARD_DIR / "explorer.html")
    ap.add_argument("--alt-ordinations", action="store_true",
                    help="include NMDS, t-SNE, UMAP alongside the default ordination")
    ap.add_argument("--per-species-permanova", action="store_true", default=None,
                    help="force per-species PERMANOVA even for large subsets")
    ap.add_argument("--no-per-species-permanova", action="store_false", dest="per_species_permanova",
                    help="skip per-species PERMANOVA even for small subsets")
    ap.add_argument("--compute", action="store_true")
    ap.add_argument("--render", action="store_true")
    args = ap.parse_args()

    # ---- layer config: built-in defaults <- YAML <- CLI ------------------------
    cfg_path = args.config
    cfg = DashboardConfig.from_yaml(cfg_path) if cfg_path.exists() else DashboardConfig.from_dict({})

    # CLI arg name -> DashboardConfig field name (only args the user explicitly set).
    _CLI_TO_FIELD = {
        "db": "db_path", "busco_dir": "busco_dir", "results_dir": "results_dir",
        "tree": "tree_path", "outgroup_genus": "outgroup_genus",
        "merops_families_path": "merops_families_path", "title": "title",
        "subtitle": "subtitle", "feature_sets": "feature_sets",
        "alt_ordinations": "alt_ordinations",
        "per_species_permanova": "per_species_permanova",
    }
    overrides = {}
    for cli_name, field in _CLI_TO_FIELD.items():
        val = getattr(args, cli_name, None)
        if val is not None:
            overrides[field] = val

    subset = cfg.subset
    for field in ("phylum", "klass", "order", "family", "genus", "species",
                  "locustags", "max_per_species", "max_genomes", "seed"):
        cli_val = getattr(args, field, None)
        if cli_val is not None:
            setattr(subset, field, cli_val)
    if args.only_with_functional:
        subset.only_with_functional = True

    cfg = DashboardConfig(
        db_path=overrides.get("db_path", cfg.db_path),
        busco_dir=overrides.get("busco_dir", cfg.busco_dir),
        results_dir=overrides.get("results_dir", cfg.results_dir),
        tree_path=overrides.get("tree_path", cfg.tree_path),
        merops_families_path=overrides.get("merops_families_path", cfg.merops_families_path),
        outgroup_genus=overrides.get("outgroup_genus", cfg.outgroup_genus),
        feature_sets=overrides.get("feature_sets", cfg.feature_sets),
        alt_ordinations=overrides.get("alt_ordinations", cfg.alt_ordinations),
        per_species_permanova=overrides.get("per_species_permanova", cfg.per_species_permanova),
        title=overrides.get("title", cfg.title),
        subtitle=overrides.get("subtitle", cfg.subtitle),
        subset=subset,
    )

    do_compute = args.compute or not (args.compute or args.render)
    do_render = args.render or not (args.compute or args.render)

    if do_compute:
        print(f"Computing payload from {cfg.db_path} ...", file=sys.stderr)
        payload = build_payload(cfg)
        args.cache.parent.mkdir(parents=True, exist_ok=True)
        args.cache.write_text(json.dumps(payload))
        print(f"Wrote {args.cache} ({args.cache.stat().st_size:,} bytes)", file=sys.stderr)

    if do_render:
        payload = json.loads(args.cache.read_text())
        import explorer_template

        html = explorer_template.build_html(payload)
        args.out.write_text(html)
        print(f"Wrote {args.out} ({args.out.stat().st_size:,} bytes)", file=sys.stderr)


if __name__ == "__main__":
    main()
