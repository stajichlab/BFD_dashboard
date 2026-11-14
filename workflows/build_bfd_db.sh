#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 64G --time 6:00:00 -J build_bfd_db
#SBATCH -o build_bfd_db.%j.log
# Rebuild the local BFD DuckDB (VIEWs over Fungi_BFD_runs/tables/*.parquet).
# This replaces the MotherDuck copy (workflows/md_load_bfd.sh) as the query database.
#
# This job wraps the pipeline's own builder, Fungi_BFD/nextflow/bin/build_BFD_duckDB.sh,
# so the view set, the SwissProt macros and the duplicate-key checks stay identical
# to the pipeline's db/BFD.duckdb. It does not duplicate that logic.
#
# Safety:
#   - It builds in a staging directory and moves the result into place only after
#     the checks pass. A failed job leaves the previous database untouched.
#   - The default OUT is inside BFD_dashboard. It does not overwrite the pipeline's
#     copy at Fungi_BFD_runs/BFD_function/db/BFD.duckdb.
#   - Views read the parquet at query time. Do not submit this while MERGE_* steps are
#     still rewriting tables/. The duplicate-key checks read every large table.
#
# Submit (after the merged tables are regenerated):
#   sbatch workflows/build_bfd_db.sh
# Overrides (environment):
#   TABLES=/path/to/tables  OUT=/path/to/BFD.duckdb  CHECK_UNIQUE=0 (skip key checks)
#   PIPELINE=/path/to/Fungi_BFD
set -euo pipefail

PIPELINE="${PIPELINE:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD}"
TABLES="${TABLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables}"
OUT="${OUT:-/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/db/BFD.duckdb}"
REPO="${REPO:-/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard}"
BUILDER="$PIPELINE/nextflow/bin/build_BFD_duckDB.sh"
SCHEMA="$PIPELINE/sql/table_schema.json"
BIN_DIR="$PIPELINE/nextflow/bin"
PY=/usr/bin/python3.12

for f in "$BUILDER" "$SCHEMA"; do
  [ -e "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done
[ -d "$TABLES" ] || { echo "ERROR: tables dir $TABLES not found" >&2; exit 1; }
TABLES="$(cd "$TABLES" && pwd)"   # the builder writes this absolute path into every view

OUT_DIR="$(dirname "$OUT")"
mkdir -p "$OUT_DIR"
STAGE="$OUT_DIR/.stage.${SLURM_JOB_ID:-manual}"
rm -rf "$STAGE"; mkdir -p "$STAGE"
trap 'rm -rf "$STAGE"' EXIT

# The builder wants to run from a directory where it can create logs/ and db files.
cd "$STAGE"
module load duckdb

echo "Build start $(date -Is) tables=$TABLES stage=$STAGE"
TABLES="$TABLES" DBDIR="$STAGE" BIN_DIR="$BIN_DIR" SCHEMA="$SCHEMA" \
  CHECK_UNIQUE="${CHECK_UNIQUE:-1}" bash "$BUILDER"

# ---- MotherDuck-compatibility and helper views ------------------------------
# sql/bfd_compat_views.sql: MotherDuck table names that the pipeline builder lacks, plus comments.
# sql/bfd_extra_views.sql : helper views and the in_genome() macro (see docs/BFD_VIEW_REVIEW.md).
# Both files carry the tables path on their first line; replace it with the real $TABLES.
# Both are idempotent and read only parquet, so nothing is copied into the database.
SQLDIR="$REPO/sql"
for f in bfd_compat_views.sql bfd_extra_views.sql; do
  [ -e "$SQLDIR/$f" ] || { echo "ERROR: missing $SQLDIR/$f" >&2; exit 1; }
  DEFAULT_T="$(sed -n '1s/^-- TABLES_DIR=//p' "$SQLDIR/$f")"
  sed "s#${DEFAULT_T:-/nonexistent}#$TABLES#g" "$SQLDIR/$f" > "$STAGE/$f"
  echo "Applying $f"
  duckdb "$STAGE/BFD.duckdb" < "$STAGE/$f"
done

# ---- verification before publish -------------------------------------------
# 1. every required parquet has a view and the view row count equals the parquet count
# 2. the view set is not smaller than the one we are replacing
"$PY" - "$STAGE/BFD.duckdb" "$TABLES" "$OUT" <<'PYEOF'
import sys
from pathlib import Path
import duckdb

new_db, tables, old_db = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
con = duckdb.connect(new_db, read_only=True)
views = {r[0] for r in con.sql("SELECT view_name FROM duckdb_views() WHERE NOT internal").fetchall()}
print(f"new db: {len(views)} views")

problems = []
# row counts for the small per-genome tables (cheap); big tables are covered by the key checks
for view, pq in [("species", "species.parquet"), ("asm_stats", "asm_stats.parquet"),
                 ("busco_genome", "busco_genome.parquet")]:
    if view not in views:
        problems.append(f"view {view} missing")
        continue
    nv = con.sql(f"SELECT count(*) FROM {view}").fetchone()[0]
    npq = con.sql(f"SELECT count(*) FROM read_parquet('{tables}/{pq}')").fetchone()[0]
    print(f"  {view:14s} view={nv:,d} parquet={npq:,d}")
    if nv != npq:
        problems.append(f"{view}: view {nv} != parquet {npq}")

# every name the MotherDuck database had must resolve (LIMIT 0 binds the view without scanning it)
MD_NAMES = """aa_freq asm_stats busco_genome cazy_cazymes_hmm cazy_overview codon_freq gene_CDS gene_exons
gene_info gene_intergenic_distances gene_introns gene_proteins gene_transcripts gene_trnas idp idp_summary
infernal_rfam merops pfam predgpi samples signalp_signal_peptide species swissprot_annot swissprot_best
swissprot_subject swissprot_subject_summary targetP telomere_summary telomere_tracts tf_inventory tmhmm
wgd_ks wgd_ks_nodes wgd_ks_summary wolfpsort""".split()
for name in MD_NAMES:
    try:
        con.sql(f'SELECT * FROM "{name}" LIMIT 0')
    except Exception as e:
        problems.append(f"MotherDuck name {name} does not resolve: {str(e)[:120]}")
for name in ["v_genome_labels", "v_genome_secretome", "v_gene_introns_phase"]:
    if name not in views:
        problems.append(f"helper view {name} missing")

if Path(old_db).exists():
    old = duckdb.connect(old_db, read_only=True)
    old_views = {r[0] for r in old.sql("SELECT view_name FROM duckdb_views() WHERE NOT internal").fetchall()}
    lost = sorted(old_views - views)
    print(f"previous db: {len(old_views)} views; views lost: {lost or 'none'}")
    if lost:
        problems.append(f"views present before but not now: {lost}")

if problems:
    sys.exit("VERIFY FAILED:\n  " + "\n  ".join(problems))
print("verify ok")
PYEOF

mv -f "$STAGE/BFD.duckdb" "$OUT.new"
mv -f "$OUT.new" "$OUT"          # same filesystem, so the rename is atomic
echo "Published $OUT ($(du -h "$OUT" | cut -f1)) at $(date -Is)"
