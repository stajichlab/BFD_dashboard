#!/usr/bin/bash -l
# Build queries.faa for the T1+T2 tblastn check from queries.tsv (columns key, source, accession, expected_len, start, end).
# Run on an HPCC login node (light): UniProt entries are downloaded by accession only (nothing else is sent) and their length is
# checked against expected_len (bfd.swissprot_subject.s_len); BFD proteins are extracted from their proteome and, when start/end
# are given, trimmed to that residue range (the Pfam domain envelope). Header: >KEY|ACCESSION.
# Required: BENCH (output folder under ~/projects/BFD_dashboard/bench/) with queries.tsv in it.
set -euo pipefail
: "${BENCH:?set BENCH to the bench folder that holds queries.tsv}"
BENCH=$(readlink -f "$BENCH")
MAP=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/locustag_to_pep.tsv   # read-only
OUT="$BENCH/queries.faa"; : > "$OUT.tmp"; : > "$BENCH/queries.lengths.tsv"
tail -n +2 "$BENCH/queries.tsv" | while IFS=$'\t' read -r key src acc exp s e; do
  if [ "$src" = uniprot ]; then
    seq=$(curl -fsS "https://rest.uniprot.org/uniprotkb/${acc}.fasta")
    len=$(printf '%s\n' "$seq" | awk 'NR>1{n+=length($0)} END{print n}')
    if [ "$len" != "$exp" ]; then echo "LENGTH MISMATCH $key $acc: got $len expected $exp" >&2; exit 3; fi
    printf '%s\n' "$seq" | awk -v h=">${key}|${acc}" 'NR==1{print h; next} {print}' >> "$OUT.tmp"
    echo -e "$key\t$acc\t$len" >> "$BENCH/queries.lengths.tsv"
  else
    loc=${acc%%_*}
    pep=$(awk -F'\t' -v l="$loc" '$1==l{print $3}' "$MAP")
    test -n "$pep" || { echo "no proteome for $loc" >&2; exit 4; }
    full=$(zcat "$pep" | awk -v id="$acc" '/^>/{p=($1==">" id); next} p{printf "%s", $0}')
    test -n "$full" || { echo "protein $acc not found in $pep" >&2; exit 6; }
    if [ -n "${s:-}" ] && [ -n "${e:-}" ]; then sub=${full:$((s-1)):$((e-s+1))}; else sub=$full; fi
    printf '>%s|%s\n%s\n' "$key" "$acc" "$sub" | fold -w 80 >> "$OUT.tmp"
    echo -e "$key\t$acc\t${#sub}" >> "$BENCH/queries.lengths.tsv"
  fi
done
n=$(grep -c '^>' "$OUT.tmp"); want=$(( $(wc -l < "$BENCH/queries.tsv") - 1 ))
test "$n" -eq "$want" || { echo "queries: got $n sequences, expected $want" >&2; exit 5; }
mv "$OUT.tmp" "$OUT"; echo "wrote $OUT with $n sequences"; cat "$BENCH/queries.lengths.tsv"
