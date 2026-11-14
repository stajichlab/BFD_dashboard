# kx2_h2: where the megabases went (H2, Kickxellales vs yeast reductions)

All SQL runs in MotherDuck, reading `bfd` and writing only `bfd_scratch.kx2_h2_*`. Nothing was done by hand or in Python; all steps are the SQL in these files, which are the final working versions of the queries run.

## Run order and dependencies

| file | output | depends on |
|---|---|---|
| 01_genomes.sql | kx2_h2_genomes | bfd_scratch.kx_genomes, bfd.asm_stats |
| 02_tx.sql | kx2_h2_tx | 01, bfd.gene_transcripts, bfd.gene_exons |
| 03_ig.sql | kx2_h2_ig | 01, bfd.gene_info, bfd.gene_intergenic_distances |
| 04_acct.sql | kx2_h2_acct | 01, 02, 03 |
| 05_lin.sql | kx2_h2_lin | 04 |
| 06_decomposition.sql | result query only (ordering A and B) | 05 |
| 07_igorient.sql | kx2_h2_igorient (+ result query) | 01, 03 |
| 08_genus.sql | kx2_h2_genus | 04, 07 |
| 09_delta.sql | kx2_h2_delta (+ pivot query) | 08 |

Graph: kx_genomes/asm_stats -> 01 -> {02, 03} -> 04 -> 05 -> 06; 03 -> 07 -> 08 (also from 04) -> 09.

## Notes

- Lineage assignment (`kx2_h2_genomes.lineage`) is a hand-defined CASE list over `kx_genomes.grp`, not derived from taxonomy. KxCore = grp Kx:Coemansia + Kx:Kickxellales-core-other. Saccharomycotina = all `Sacch:*`; Pezizomycotina = all `Pez:*`.
- Filters: Q = busco_complete >= 85 and N50 >= 50 kb (flag `q` in 01); Q70 = busco >= 70 and N50 >= 25 kb (applied in 05); "all" = no filter. No Harpellales or Dimargaritales genome passes Q.
- Baselines: Mucorales for the Zygomycete groups; Pezizomycotina for yeasts and other reduced lineages (crude, ancestral state unknown).
- UTR: not annotated in any genome (exon length equals CDS length everywhere), so the UTR term is 0 by construction. Gene overlap is also absent (distances >= 1 bp apart from a few tandem pairs down to -6).
- Masked bases are not additive with the other terms (they overlap intergenic and intron space); only reported as a memo.
- 06 was run inline from kx2_h2_lin and was not materialised; the file contains the exact queries. Display/check queries in 01, 02, 05, 07, 08, 09 are commented at the foot of each file and were run as given.
- Not reconstructable exactly: the first attempt at kx2_h2_genus/kx2_h2_delta failed on a syntax error (`ref` alias) and was re-run; only the final versions are saved. The 1,919 kx_genomes rows are an existing scratch table, not rebuilt here.

## Verification (run 2026-11-11)

Steps 05, 08, 09 and 06 were re-run into `bfd_scratch.kx2_h2_chk_lin`, `_chk_genus`, `_chk_delta` (then dropped; confirmed none left).

- kx2_h2_chk_lin: 37 rows, identical to kx2_h2_lin on 37/37 (Mb and intergenic Mb).
- kx2_h2_chk_delta: 119 rows, Cliff's delta identical on 119/119.
- kx2_h2_chk_genus: 1,650 rows, same as kx2_h2_genus.
- Decomposition, ordering A, Q: core Kx vs Mucorales 70 / -2 / 9 / 22 / 1 % (gene number / CDS per gene / introns / intergenic / remainder), deficit 15.0 Mb; Saccharomycotina vs Pezizomycotina 72 / 0 / 2 / 23 / 3 %, 27.5 Mb. Matches the report.
- Unfiltered (all): core Kx 77 / -2 / 6 / 20 / -3; Saccharomycotina 72 / -1 / 3 / 22 / 4. Matches the report.
- Ordering B, Q: core Kx 58 / 13 / 32 % (gene number / intron / intergenic); Saccharomycotina 47 / 4 / 43 %. Matches the report.
- The re-run reused the existing kx2_h2_acct, kx2_h2_igorient and upstream tables; steps 01-04 and 07 were not rebuilt from scratch, so they were not re-verified.
