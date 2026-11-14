#!/usr/bin/bash -l
# Build queries.faa for the meiotic-toolkit tblastn check (adapted from ../../tblastn/fetch_queries.sh).
# Run on an HPCC login node (light): downloads the public UniProt FASTA entries by accession (only accessions are sent) and extracts
# BFD proteins from their proteomes. Header: >KEY|ACCESSION. UniProt lengths are checked against expected_len (bfd.swissprot_subject.s_len);
# BFD proteins against a 'min-max' length window. Required: BENCH (a purgeable one-off folder under ~/projects/BFD_dashboard/bench/).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
: "${BENCH:?set BENCH, e.g. ~/projects/BFD_dashboard/bench/kx2_meiosis}"
BENCH=$(readlink -f "$BENCH"); mkdir -p "$BENCH"
MAP=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/locustag_to_pep.tsv   # read-only
QTSV="$HERE/queries.tsv"; OUT="$BENCH/queries.faa"; : > "$OUT.tmp"; : > "$OUT.lengths.tmp"
tail -n +2 "$QTSV" | while IFS=$'\t' read -r key src acc exp; do
  if [ "$src" = uniprot ]; then
    seq=$(curl -fsS "https://rest.uniprot.org/uniprotkb/${acc}.fasta")
    len=$(printf '%s\n' "$seq" | awk 'NR>1{n+=length($0)} END{print n}')
    if [ "$len" != "$exp" ]; then echo "LENGTH MISMATCH $key $acc: got $len expected $exp" >&2; exit 3; fi
    printf '%s\n' "$seq" | awk -v h=">${key}|${acc}" 'NR==1{print h; next} {print}' >> "$OUT.tmp"
    echo "$key $acc $len aa" >> "$OUT.lengths.tmp"
  elif [ "$src" = bfd ]; then
    loc=${acc%%_*}
    pep=$(awk -F'\t' -v l="$loc" '$1==l{print $3}' "$MAP")
    test -n "$pep" || { echo "no proteome for $loc" >&2; exit 4; }
    seq=$(zcat "$pep" | awk -v id="$acc" -v h=">${key}|${acc}" '/^>/{p=($1==">" id)} p{if(/^>/) print h; else print}')
    test -n "$seq" || { echo "protein $acc not found in $pep" >&2; exit 6; }
    len=$(printf '%s\n' "$seq" | awk 'NR>1{n+=length($0)} END{print n}')
    if [[ "$exp" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      if [ "$len" -lt "${BASH_REMATCH[1]}" ] || [ "$len" -gt "${BASH_REMATCH[2]}" ]; then
        echo "LENGTH OUT OF RANGE $key $acc: $len not in $exp" >&2; exit 7; fi
    fi
    printf '%s\n' "$seq" >> "$OUT.tmp"
    echo "$key $acc $len aa" >> "$OUT.lengths.tmp"
  fi
done
n=$(grep -c '^>' "$OUT.tmp"); want=$(( $(wc -l < "$QTSV") - 1 ))
test "$n" -eq "$want" || { echo "queries: got $n sequences, expected $want" >&2; exit 5; }
mv "$OUT.tmp" "$OUT"; mv "$OUT.lengths.tmp" "$OUT.lengths"; echo "wrote $OUT with $n sequences"
