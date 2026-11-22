#!/usr/bin/bash -l
#SBATCH -p epyc -c 4 --mem 16G --time 1:00:00 -J build_labels_db
#SBATCH -o build_labels_db.%j.log
# Build the dashboard's labels file db/BFD_labels.duckdb from a pipeline database (see sql/dashboard_views.sql).
# The pipeline owns the facts; this file holds the analysis labels (v_genome_labels, v_genome_traits, ref_* tables).
# Rebuild it after every pipeline build, or after a change to config/refs/ or sql/dashboard_views.sql.
#
# Safety: it builds in a staging directory, checks the result, and moves it into place only if the checks pass.
# A failed job leaves the previous labels file untouched. The pipeline database is attached READ_ONLY.
#
# Submit from anywhere:   sbatch /bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/workflows/build_labels_db.sh
# Overrides (environment): SRC=/path/to/BFD.duckdb  OUT=/path/to/BFD_labels.duckdb  REPO=/path/to/BFD_dashboard
#   Default SRC is the pipeline copy Fungi_BFD_runs/BFD_function/db/BFD.duckdb (switched from the dashboard's own
#   db/BFD.duckdb on 2026-11-21).
set -euo pipefail

REPO="${REPO:-/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard}"
SRC="${SRC:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/BFD_function/db/BFD.duckdb}"
OUT="${OUT:-$REPO/db/BFD_labels.duckdb}"
REF_DIR="${REF_DIR:-$REPO/config/refs}"
SQLFILE="$REPO/sql/dashboard_views.sql"
PY=/usr/bin/python3.12

[ -s "$SRC" ] || { echo "ERROR: source database $SRC not found or empty" >&2; exit 1; }
[ -e "$SQLFILE" ] || { echo "ERROR: missing $SQLFILE" >&2; exit 1; }
for f in fungaltraits_1.2.utf8.csv growth_form_overrides.tsv class_subphylum.tsv phylum_subphylum.tsv species_subphylum.tsv; do
  [ -s "$REF_DIR/$f" ] || { echo "ERROR: missing or empty $REF_DIR/$f" >&2; exit 1; }
done
REF_DIR="$(cd "$REF_DIR" && pwd)"
SRC="$(cd "$(dirname "$SRC")" && pwd)/$(basename "$SRC")"
case "$SRC$REF_DIR" in *[\#\&\\\'\"]*) echo "ERROR: paths must not contain # & \\ or quotes" >&2; exit 1 ;; esac

OUT_DIR="$(dirname "$OUT")"
mkdir -p "$OUT_DIR"
STAGE="$OUT_DIR/.stage_labels.${SLURM_JOB_ID:-manual}"
rm -rf "$STAGE"; mkdir -p "$STAGE"
trap 'rm -rf "$STAGE"' EXIT
module load duckdb

SRC_MTIME="$(stat -c %Y "$SRC")"
DEFAULT_P="$(sed -n '1s/^-- REF_DIR=//p' "$SQLFILE")"
[ -n "$DEFAULT_P" ] || { echo "ERROR: line 1 of $SQLFILE must be '-- REF_DIR=<path>'" >&2; exit 1; }
sed -e "s#${DEFAULT_P}#${REF_DIR}#g" -e "s#@SOURCE_DB@#${SRC}#g" -e "s#@SOURCE_MTIME@#${SRC_MTIME}#g" "$SQLFILE" > "$STAGE/dashboard_views.sql"

echo "Labels build start $(date -Is) source=$SRC out=$OUT"
{ echo "ATTACH '$SRC' AS bfd (READ_ONLY);"; cat "$STAGE/dashboard_views.sql"; } | duckdb "$STAGE/labels.duckdb"

"$PY" - "$STAGE/labels.duckdb" "$SRC" <<'PYEOF'
import sys
import duckdb

new, src = sys.argv[1], sys.argv[2]
con = duckdb.connect(new, read_only=True)
con.execute(f"ATTACH '{src}' AS bfd (READ_ONLY)")
problems = []
n_sp = con.sql("SELECT count(*) FROM bfd.species").fetchone()[0]
n_lab, n_key = con.sql("SELECT count(*), count(DISTINCT LOCUSTAG) FROM labels.main.v_genome_labels").fetchone()
n_nosub = con.sql("SELECT count(*) FROM labels.main.v_genome_labels WHERE SUBPHYLUM IS NULL").fetchone()[0]
print(f"  v_genome_labels rows={n_lab:,d} distinct LOCUSTAG={n_key:,d} species={n_sp:,d}; SUBPHYLUM NULL for {n_nosub:,d}")
if not (n_lab == n_key == n_sp):
    problems.append(f"v_genome_labels rows {n_lab}, distinct LOCUSTAG {n_key}, species rows {n_sp} must be equal")
for name in ["ref_fungaltraits", "ref_growth_form_override", "ref_class_subphylum", "ref_phylum_subphylum", "ref_species_subphylum"]:
    n = con.sql(f"SELECT count(*) FROM labels.main.{name}").fetchone()[0]
    print(f"  {name:26s} {n:,d} rows")
    if n == 0:
        problems.append(f"table {name} is empty")
if con.sql("SELECT count(*) FROM labels.main.labels_meta").fetchone()[0] != 1:
    problems.append("labels_meta must have one row")
if problems:
    sys.exit("VERIFY FAILED:\n  " + "\n  ".join(problems))
print("verify ok")
PYEOF

mv -f "$STAGE/labels.duckdb" "$OUT.new"
mv -f "$OUT.new" "$OUT"
echo "Published $OUT ($(du -h "$OUT" | cut -f1)) at $(date -Is)"
