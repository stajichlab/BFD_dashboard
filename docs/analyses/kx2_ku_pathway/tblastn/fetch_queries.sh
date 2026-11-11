#!/usr/bin/bash -l
# Build queries.faa for the tblastn check. Run on an HPCC login node (light): it downloads 22 public UniProt FASTA entries by
# accession (only accessions are sent) and extracts one BFD protein from its proteome. Header format: >KEY|ACCESSION.
# Every UniProt entry's length is checked against expected_len from queries.tsv (taken from bfd.swissprot_subject.s_len).
# Required: BENCH (output folder, a purgeable one-off folder under ~/projects/BFD_dashboard/bench/).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
: "${BENCH:?set BENCH to the output folder, e.g. ~/projects/BFD_dashboard/bench/kx2_ku_tblastn}"
BENCH=$(readlink -f "$BENCH"); mkdir -p "$BENCH"
MAP=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/locustag_to_pep.tsv   # read-only
OUT="$BENCH/queries.faa"; : > "$OUT.tmp"
tail -n +2 "$HERE/queries.tsv" | while IFS=$'\t' read -r key src acc exp; do
  if [ "$src" = uniprot ]; then
    seq=$(curl -fsS "https://rest.uniprot.org/uniprotkb/${acc}.fasta")
    len=$(printf '%s\n' "$seq" | awk 'NR>1{n+=length($0)} END{print n}')
    if [ "$len" != "$exp" ]; then echo "LENGTH MISMATCH $key $acc: got $len expected $exp" >&2; exit 3; fi
    printf '%s\n' "$seq" | awk -v h=">${key}|${acc}" 'NR==1{print h; next} {print}' >> "$OUT.tmp"
  elif [ "$src" = bfd ]; then
    loc=${acc%%_*}
    pep=$(awk -F'\t' -v l="$loc" '$1==l{print $3}' "$MAP")
    test -n "$pep" || { echo "no proteome for $loc" >&2; exit 4; }
    zcat "$pep" | awk -v id="$acc" -v h=">${key}|${acc}" '/^>/{p=($1==">" id)} p{if(/^>/) print h; else print}' >> "$OUT.tmp"
  fi
done
n=$(grep -c '^>' "$OUT.tmp"); want=$(( $(wc -l < "$HERE/queries.tsv") - 1 ))
test "$n" -eq "$want" || { echo "queries: got $n sequences, expected $want" >&2; exit 5; }
mv "$OUT.tmp" "$OUT"; echo "wrote $OUT with $n sequences"
