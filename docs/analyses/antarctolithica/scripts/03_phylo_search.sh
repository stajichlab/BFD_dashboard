#!/bin/bash
# hmmsearch BUSCO fungi_odb10 HMMs against each proteome (NCPU parallel, 1 cpu each; default 14). Run from AO_BASE or set AO_BASE.
cd "${AO_BASE:-.}/phylo"; mkdir -p tbl fa; rm -f tbl/*.tmp
one() { f=$1; b=$(basename $f .proteins.fa.gz); [ -s tbl/$b.tbl ] && return; gzip -dc $f > fa/$b.fa; hmmsearch --cpu 1 -E 1e-10 --tblout tbl/$b.tbl.tmp -o /dev/null db/all.hmm fa/$b.fa && mv tbl/$b.tbl.tmp tbl/$b.tbl; }
export -f one
ls pep/*.fa.gz | xargs -P "${NCPU:-14}" -n 1 bash -c 'one "$0"'
echo SEARCH_DONE
