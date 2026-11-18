#!/usr/bin/bash
# Copy the figure run output (figures, tables, figures_index.tsv) into the Quarto site folder.
# FOR LOCAL PREVIEWS ONLY. It copies into folders that git ignores. It does not publish the site.
# To publish figures, use workflows/release_site_assets.sh --upload (release asset + manifest).
# Usage: workflows/publish_figures_to_site.sh [--dry-run]
# Override ROOT (repo root) and OUT (figure run output) with env vars.
# Files are merged into the destination. Extra files in the destination are not deleted.
set -euo pipefail
ROOT=${ROOT:-$(git rev-parse --show-toplevel)}
OUT=${OUT:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_figures_real}
SITE="$ROOT/docs/quarto_summary"
DRY=""
if [ "${1:-}" = "--dry-run" ]; then DRY="--dry-run"; fi

if [ ! -f "$OUT/figures_index.tsv" ]; then
  echo "ERROR: $OUT/figures_index.tsv is missing. Run workflows/summary_figures.sh first." >&2
  exit 1
fi
for d in figures tables; do
  if [ ! -d "$OUT/$d" ]; then echo "ERROR: $OUT/$d is missing." >&2; exit 1; fi
done

mkdir -p "$SITE/figures" "$SITE/tables"
n=0; bytes=0
for d in figures tables; do
  list=$(rsync -a $DRY --out-format='%n' "$OUT/$d/" "$SITE/$d/" | grep -v '/$' || true)
  c=$(printf '%s' "$list" | grep -c . || true)
  b=0
  if [ "$c" -gt 0 ]; then
    b=$(printf '%s\n' "$list" | while IFS= read -r f; do stat -c %s "$OUT/$d/$f"; done | awk '{s+=$1} END{print s+0}')
  fi
  echo "$d: $c files changed or new"
  n=$((n + c)); bytes=$((bytes + b))
done
rsync -a $DRY "$OUT/figures_index.tsv" "$SITE/figures_index.tsv"
n=$((n + 1)); bytes=$((bytes + $(stat -c %s "$OUT/figures_index.tsv")))
mb=$(awk -v b="$bytes" 'BEGIN{printf "%.1f", b/1048576}')
if [ -n "$DRY" ]; then echo "dry run: would copy $n files, $mb MB"; else echo "copied $n files, $mb MB"; fi
