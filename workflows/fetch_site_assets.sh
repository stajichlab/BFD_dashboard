#!/usr/bin/bash
# Download the pinned figure archive and unpack it into docs/quarto_summary/.
# Usage: workflows/fetch_site_assets.sh [--from-archive PATH] [--no-verify]
#   --from-archive PATH  use a local archive instead of downloading (local tests)
#   --no-verify          skip the sha256 check (local experiments only)
# Env: ROOT (repo root, default: git toplevel), REPO (default stajichlab/BFD_dashboard),
#      GH_TOKEN (used by gh when set).
# Exit codes: 2 manifest missing, 3 sha256 mismatch, 4 download or unpack failure.
set -euo pipefail
ROOT=${ROOT:-$(git rev-parse --show-toplevel)}
REPO=${REPO:-stajichlab/BFD_dashboard}
SITE="$ROOT/docs/quarto_summary"
MANIFEST="$SITE/figures_release.json"
LOCAL=""; VERIFY=1
while [ $# -gt 0 ]; do
  case "$1" in
    --from-archive) LOCAL=$2; shift 2 ;;
    --no-verify) VERIFY=0; shift ;;
    *) echo "ERROR: unknown argument $1" >&2; exit 1 ;;
  esac
done
if [ ! -f "$MANIFEST" ]; then
  echo "ERROR: $MANIFEST is missing. Run workflows/release_site_assets.sh and commit the manifest." >&2
  exit 2
fi
field() { python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]])' "$MANIFEST" "$1"; }
TAG=$(field tag); ASSET=$(field asset); WANT=$(field sha256)

TMP=$(mktemp -d "${TMPDIR:-/tmp}/fetch_site.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

if [ -n "$LOCAL" ]; then
  ARCHIVE=$LOCAL; DOWNLOADED=0
  if [ ! -f "$ARCHIVE" ]; then echo "ERROR: $ARCHIVE not found." >&2; exit 4; fi
else
  DOWNLOADED=1; ARCHIVE="$TMP/$ASSET"
  if command -v gh >/dev/null 2>&1; then
    gh release download "$TAG" --repo "$REPO" --pattern "$ASSET" --dir "$TMP" \
      || { echo "ERROR: gh release download failed for $TAG/$ASSET." >&2; exit 4; }
  else
    curl -fL -o "$ARCHIVE" "https://github.com/$REPO/releases/download/$TAG/$ASSET" \
      || { echo "ERROR: curl download failed for $TAG/$ASSET." >&2; exit 4; }
  fi
fi

if [ "$VERIFY" = 1 ]; then
  GOT=$(sha256sum "$ARCHIVE" | cut -d' ' -f1)
  if [ "$GOT" != "$WANT" ]; then
    echo "ERROR: sha256 mismatch for $ASSET: manifest $WANT, file $GOT." >&2
    if [ "$DOWNLOADED" = 1 ]; then rm -f "$ARCHIVE"; fi
    exit 3
  fi
else
  echo "WARNING: sha256 check skipped (--no-verify)." >&2
fi

case "$ASSET" in
  *.tar.zst) UNPACK=(tar -I "zstd -d") ;;
  *.tar.gz) UNPACK=(tar -z) ;;
  *) echo "ERROR: unknown archive type $ASSET." >&2; exit 4 ;;
esac
if "${UNPACK[@]}" -tf "$ARCHIVE" 2>/dev/null | grep -Eq '^/|(^|/)\.\.(/|$)'; then
  echo "ERROR: archive has absolute or parent paths." >&2; exit 4
fi
mkdir "$TMP/x"
if ! "${UNPACK[@]}" -xf "$ARCHIVE" -C "$TMP/x"; then
  echo "ERROR: could not unpack $ARCHIVE." >&2; exit 3
fi
for need in figures tables figures_index.tsv; do
  if [ ! -e "$TMP/x/$need" ]; then echo "ERROR: archive lacks $need." >&2; exit 4; fi
done

rm -rf "$SITE/figures" "$SITE/tables" "$SITE/figures_index.tsv"
mv "$TMP/x/figures" "$TMP/x/tables" "$TMP/x/figures_index.tsv" "$SITE/"
n=$(find "$SITE/figures" -type f -name '*.png' | wc -l)
mb=$(du -sm --apparent-size "$SITE/figures" "$SITE/tables" "$SITE/figures_index.tsv" | awk '{s+=$1} END{print s}')
echo "fetched $ASSET: $n PNG files, $mb MB"
