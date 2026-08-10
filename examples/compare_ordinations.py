"""Compare ordinations across feature sets for a subset of the BFD database.

Demonstrates the generic (subset-aware) API:

  python3 examples/compare_ordinations.py --config config/rhodotorula_reps.yaml
  python3 examples/compare_ordinations.py --db /path/BFD.duckdb --genus Fusarium --max-genomes 30

For each available feature set: build the matrix for the resolved subset, run PCoA,
report variance explained and (for distance ordinations) whether the 2D projection is
trustworthy. Calls the exact same functions the dashboard uses, so this doubles as a
sanity check before a full --compute.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))

import bfd_data as bd
import ordination as ordmod
import subset as submod
from config import DashboardConfig


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", help="dashboard YAML (db_path/subset/tree knobs are reused)")
    ap.add_argument("--db", help="override db path")
    ap.add_argument("--genus", "--family", "--order", dest="taxa", action="append", default=None)
    ap.add_argument("--max-genomes", type=int, default=None)
    ap.add_argument("--only-with-functional", action="store_true")
    args = ap.parse_args()

    cfg = DashboardConfig.from_yaml(args.config) if args.config else DashboardConfig(db_path=args.db)
    if args.db:
        cfg.db_path = Path(args.db)
    if args.taxa:
        for spec in args.taxa:
            rank, _, value = spec.partition("=")
            setattr(cfg.subset, rank, value)
    if args.max_genomes:
        cfg.subset.max_genomes = args.max_genomes
    if args.only_with_functional:
        cfg.subset.only_with_functional = True
    if not cfg.db_path:
        ap.error("need --db (or --config with a db_path)")

    con = bd.get_connection(cfg.db_path)
    rows, report = submod.resolve_subset(con, cfg.subset)
    submod.report_subset(rows, report)
    print(f"subset -> {len(rows)} genomes\n")

    for name, spec in bd.make_feature_sets(cfg.merops_families_path).items():
        loader = spec["loader"]
        try:
            mat = loader(con, rows=rows)
        except bd.MeropsUnavailableError:
            print(f"{name:12s} skipped (merops families file not configured)")
            continue
        except Exception as exc:  # noqa: BLE001 -- report, don't stop the comparison
            print(f"{name:12s} FAILED: {exc}")
            continue
        if mat is None or mat.empty or mat.shape[0] < 3:
            print(f"{name:12s} no data (shape {getattr(mat, 'shape', None)})")
            continue

        distance = spec["metric"]
        method = spec["method"]
        ordination = (ordmod.correspondence_analysis(mat) if method == "ca"
                      else ordmod.pcoa(ordmod.distance_matrix(mat, distance)))
        pct = [float(x) for x in ordination.pct_variance[:2]]
        frac = (pct[0] + pct[1]) / 100.0
        print(f"{name:12s} {mat.shape[0]}x{mat.shape[1]}  "
              f"method={method:4s}  PCo1={pct[0]:6.1f}% PCo2={pct[1]:6.1f}%  "
              f"2D={frac:.0%} of variation"
              + ("" if frac >= 0.4 else "  [low -- projection unreliable]"))

        # Per-species differential families (only meaningful for domain matrices).
        if name in ("pfam", "cazy", "merops") and len(rows) >= 4:
            diff = bd.differential_families(con, name, rows=rows,
                                            merops_families_path=cfg.merops_families_path)
            focus = next(iter(diff))
            top = diff[focus]["present_only"] or diff[focus]["top_fold_change"]
            if top:
                print(f"              example differential for {focus}: "
                      f"{top[0]['family']} ({top[0].get('value_per_1000', '')} per 1000)")


if __name__ == "__main__":
    main()
