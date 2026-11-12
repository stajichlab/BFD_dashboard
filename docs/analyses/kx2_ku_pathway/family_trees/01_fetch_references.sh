#!/usr/bin/bash -l
# Step 1: fetch the reference proteins in references.tsv from UniProt (only accessions are sent) and verify each length against
# expected_len (taken from bfd.swissprot_subject.s_len). Writes $BENCH/refs/<family>.faa with headers >CLADE|ACCESSION|gene|organism.
# Run on an HPCC login node (light). Required: BENCH (a purgeable folder), HERE (this directory).
set -euo pipefail
: "${BENCH:?set BENCH}"; : "${HERE:?set HERE to the folder that holds references.tsv}"
BENCH=$(readlink -f "$BENCH"); mkdir -p "$BENCH/refs"
for fam in lig reca polx; do : > "$BENCH/refs/$fam.faa.tmp"; done
tail -n +2 "$HERE/references.tsv" | while IFS=$'\t' read -r fam clade acc gene org exp; do
  seq=$(curl -fsS "https://rest.uniprot.org/uniprotkb/${acc}.fasta")
  len=$(printf '%s\n' "$seq" | awk 'NR>1{n+=length($0)} END{print n}')
  if [ "$len" != "$exp" ]; then echo "LENGTH MISMATCH $fam $clade $acc: got $len expected $exp" >&2; exit 3; fi
  printf '%s\n' "$seq" | awk -v h=">${clade}|${acc}|${gene}|${org// /_}" 'NR==1{print h; next} {print}' >> "$BENCH/refs/$fam.faa.tmp"
done
for fam in lig reca polx; do
  want=$(awk -F'\t' -v f=$fam 'NR>1 && $1==f' "$HERE/references.tsv" | wc -l); n=$(grep -c '^>' "$BENCH/refs/$fam.faa.tmp")
  test "$n" -eq "$want" || { echo "$fam: got $n of $want" >&2; exit 5; }
  mv "$BENCH/refs/$fam.faa.tmp" "$BENCH/refs/$fam.faa"; echo "$fam: $n references"
done
