#!/usr/bin/bash
# Package the generated site figures as one release asset and update the manifest.
# Usage: workflows/release_site_assets.sh [--out DIR] [--date YYYY-MM-DD] [--upload] [--dry-run]
#   --out DIR   where the archive is written (default: a new temp dir; never inside the repo)
#   --date D    date in the tag and file name (default: today)
#   --upload    run gh release create. Checks gh auth and the tag BEFORE building. If the tag
#               exists it refuses, unless --clobber is given.
#   --clobber   with --upload: replace the asset of an existing tag. Other branches that pin the
#               old sha256 then fail with exit 3. Prefer a new dated tag (--date).
#   --dry-run   print the actions and the gh command, write nothing
# Env: ROOT (repo root, default: git toplevel), OUT (figure run output with figures/,
#      tables/ and figures_index.tsv), REPO (default stajichlab/BFD_dashboard).
# The manifest is written to docs/quarto_summary/figures_release.json. Commit it.
# It is written after a successful upload (with --upload) or after the local build (without).
# Manifest field source.branch_commit is the HEAD at packaging time.
set -euo pipefail
ROOT=${ROOT:-$(git rev-parse --show-toplevel)}
OUT=${OUT:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_figures_real}
REPO=${REPO:-stajichlab/BFD_dashboard}
SITE="$ROOT/docs/quarto_summary"
DEST=""; DATE=$(date +%Y-%m-%d); UPLOAD=0; DRY=0; CLOBBER=0
while [ $# -gt 0 ]; do
  case "$1" in
    --out) DEST=$2; shift 2 ;;
    --date) DATE=$2; shift 2 ;;
    --upload) UPLOAD=1; shift ;;
    --dry-run) DRY=1; shift ;;
    --clobber) CLOBBER=1; shift ;;
    *) echo "ERROR: unknown argument $1" >&2; exit 1 ;;
  esac
done
case "$DATE" in
  [0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]) ;;
  *) echo "ERROR: --date must be YYYY-MM-DD" >&2; exit 1 ;;
esac

if [ ! -f "$OUT/figures_index.tsv" ]; then
  echo "ERROR: $OUT/figures_index.tsv is missing. Run workflows/summary_figures.sh first." >&2; exit 1
fi
for d in figures tables; do
  if [ ! -d "$OUT/$d" ]; then echo "ERROR: $OUT/$d is missing." >&2; exit 1; fi
done

TAG="site-figures-$DATE"
if command -v zstd >/dev/null 2>&1; then EXT="tar.zst"; else EXT="tar.gz"; echo "WARNING: zstd not found, using gzip" >&2; fi
ASSET="$TAG.$EXT"

EXISTS=0
if [ "$UPLOAD" = 1 ]; then
  if ! gh auth status >/dev/null 2>&1; then
    echo "ERROR: gh auth status failed. Run gh auth login first." >&2; exit 1
  fi
  if gh release view "$TAG" --repo "$REPO" >/dev/null 2>&1; then
    EXISTS=1
    if [ "$CLOBBER" = 0 ]; then
      echo "ERROR: release $TAG already exists. Use a new --date, or pass --clobber to replace the asset." >&2
      exit 1
    fi
    echo "WARNING: --clobber replaces the asset of $TAG. Other branches whose manifest pins the old sha256 will fail with exit 3. Prefer a new dated tag (--date)." >&2
  fi
fi

if [ -z "$DEST" ]; then
  if [ "$DRY" = 1 ]; then DEST="<temp dir>"; else DEST=$(mktemp -d "${TMPDIR:-/tmp}/site_assets.XXXXXX"); fi
fi
if [ "$DRY" = 0 ]; then
  mkdir -p "$DEST"
  d_abs=$(cd "$DEST" && pwd -P); r_abs=$(cd "$ROOT" && pwd -P)
  case "$d_abs/" in
    "$r_abs"/*) echo "ERROR: --out $DEST is inside the repo $ROOT. Use a folder outside it." >&2; exit 1 ;;
  esac
fi
ARCHIVE="$DEST/$ASSET"

cmd_create() {
  printf 'gh release create %s %s --repo %s --target main --title %s --notes %s\n' \
    "$TAG" "$ARCHIVE" "$REPO" "$(printf '%q' "Site figures $DATE")" "$(printf '%q' "$1")"
}

if [ "$DRY" = 1 ]; then
  echo "dry run: would build $ARCHIVE from $OUT (figures/, tables/, figures_index.tsv)"
  echo "dry run: would write $SITE/figures_release.json"
  cmd_create "(notes are generated after the build)"
  exit 0
fi

start=$(date +%s)
if [ "$EXT" = "tar.zst" ]; then
  tar -C "$OUT" --sort=name --mtime="$DATE 00:00:00 UTC" --owner=0 --group=0 --numeric-owner \
    -cf - figures tables figures_index.tsv | zstd -T0 -19 -q -f -o "$ARCHIVE"
else
  tar -C "$OUT" --sort=name --mtime="$DATE 00:00:00 UTC" --owner=0 --group=0 --numeric-owner \
    -cf - figures tables figures_index.tsv | gzip -n -9 > "$ARCHIVE"
fi
SHA=$(sha256sum "$ARCHIVE" | cut -d' ' -f1)
BYTES=$(stat -c %s "$ARCHIVE")
N_PNG=$(find "$OUT/figures" -type f -name '*.png' | wc -l)
N_FILES=$(( $(find "$OUT/figures" "$OUT/tables" -type f | wc -l) + 1 ))
SUP=""
for f in "$ROOT/config/summary_suppress.tsv" "$SITE/data/summary_suppress.tsv"; do
  if [ -f "$f" ]; then SUP=$(sha256sum "$f" | cut -d' ' -f1); break; fi
done
COMMIT=$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)

MB=$(awk -v b="$BYTES" 'BEGIN{printf "%.1f", b/1048576}')
NOTES="Generated figures (PNG, SVG, sidecar JSON), tables and figures_index.tsv for the BFD site.
Source: figure run output $OUT, built from the summary tables of the latest build.
Suppress list sha256: $SUP
Files: $N_FILES, PNG: $N_PNG, archive size: $MB MB, sha256: $SHA
Regenerate with workflows/release_site_assets.sh."
echo "archive: $ARCHIVE"
echo "sha256: $SHA"
echo "bytes: $BYTES"
echo "n_png: $N_PNG"
echo "n_files: $N_FILES"
echo "build seconds: $(( $(date +%s) - start ))"

if [ "$UPLOAD" = 1 ]; then
  if [ "$EXISTS" = 1 ]; then
    gh release upload "$TAG" "$ARCHIVE" --repo "$REPO" --clobber
  else
    gh release create "$TAG" "$ARCHIVE" --repo "$REPO" --target main \
      --title "Site figures $DATE" --notes "$NOTES"
  fi
fi

mkdir -p "$SITE"
cat > "$SITE/figures_release.json" <<JSON
{
  "tag": "$TAG",
  "asset": "$ASSET",
  "sha256": "$SHA",
  "bytes": $BYTES,
  "created": "$DATE",
  "n_png": $N_PNG,
  "n_files": $N_FILES,
  "source": {
    "summary_suppress_sha256": "$SUP",
    "branch_commit": "$COMMIT"
  }
}
JSON
echo "manifest: $SITE/figures_release.json"
if [ "$UPLOAD" = 0 ]; then
  echo "NOTE: the manifest points at a release that does not exist yet. Upload before you commit it."
  echo "Not uploading. The command that would run:"
  cmd_create "$NOTES"
fi
