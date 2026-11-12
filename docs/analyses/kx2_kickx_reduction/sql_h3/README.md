# kx2_h3: genome-level codon / tRNA / protein-rate proxies for Kickxellomycota (hypothesis H3)

All SQL runs in MotherDuck (DuckDB) and writes only to `bfd_scratch.kx2_h3_*`; `bfd` is read-only.
Nothing was done by hand or in Python: every number in the H3 report comes from these queries. Statements
were run through the MotherDuck MCP tool and are copied verbatim (table names qualified with `bfd_scratch.`
here; the originals ran with `bfd_scratch` as the default database).

## Run order and dependencies

```
01_codmap ---------------------------+
02_set (kx_genomes, kw_rep_genomes) -+-> 03_codon_enc  -> kx2_h3_cod ------+
                                     +-> 04_trna_tai   -> kx2_h3_trna ----+
                                     +-> 05_aa_idp     -> kx2_h3_aa, _idp +-> 07_genome_genus -> kx2_h3_genome,
                                     +-> 06_protein_rate_proxy             |      kx2_h3_genus, kx2_h3_genus_f
                                            -> kx2_h3_anchor, kx2_h3_prot -+          |
                                                                                      +-> 08_group_medians -> kx2_h3_grpmed
                                                                                      +-> 09_cliff -> kx2_h3_long, kx2_h3_cliff
                                                                                      +-> 10_residuals -> kx2_h3_resid, kx2_h3_resid2
```

Order: 01, 02, then 03-06 (independent of each other), then 07, then 08-10 (independent). Display queries
used for the report tables are included as comments at the end of 08-10.

## Key parameters
- Unit: species representative, collapsed to genus medians; group value = median of genus medians.
- BUSCO>=85 variant: genomes filtered before the genus medians (`kx2_h3_genus_f`, filt='b85').
- ENC: Wright 1990, genome-wide, F approximated by sum p^2 (no per-aa counts); GC3-expected curve
  2+s+29/(s^2+(1-s)^2). Citations to verify: Wright 1990 Gene 87:23; Sun et al. 2013 (frequency-weighted F).
- tAI: dos Reis et al. 2004 wobble weights (0.41, 0.28, 0.9999, 0.68; citation to verify), pooled genome-level.
- Protein proxy: hit qcov and scov >= 80%; anchor accession present in >=30% of genomes; OUTGRP = org_code in
  HUMAN, MOUSE, DROME, CAEEL, ARATH, DANRE, XENLA, BOVIN, RAT, DICDI (106 accessions); SCHPO = 707 accessions.
- Residuals: 5-point GC3 bins (and log10 Mb x 6 size bins), reference = Mucorales, Sacch other, Pezizomycotina,
  Basidiomycota, other; bins need >=8 reference genera.

## Reconstruction notes
- The 106 outgroup anchors were chosen by rule, not by hand: the statement in `06` regenerates them, and the
  list is stored in `bfd_scratch.kx2_h3_anchor` (aset='OUTGRP'). The animal/plant org_code list was fixed in
  advance as above; I did not test other lists.
- The exploratory looks (org_code distribution, tRNA codon semantics check, codon_freq sums, kx_genomes join
  coverage, Hanseniaspora counts in kw_rep_genomes) were one-off read-only queries and are not saved; they
  only informed design choices (the kx_genomes set holds 1,919 genomes, one per genus for most yeasts).
- The size-scatter correlations (r of log10 Mb vs GC3/ENC departure/tAI) are in a comment at the end of `10`.
- `02` reproduces the final set; originally Hanseniaspora was first drawn from kx_genomes and then replaced
  (DELETE+INSERT retained as written).

## Verification (2026-11-11)
Re-ran 07 (genome join, genus, genus_f), 08, 09, 10-resid and an independent copy of 03 (ENC/GC3) into
`bfd_scratch.kx2_h3_chk_*`, using the existing upstream tables (01, 02, 04, 05, 06 were not recomputed).
- Group medians: `EXCEPT` both ways against `kx2_h3_grpmed` = 0 differing rows. Kx core: GC3 67.2, ENC 49.5,
  ENC departure 0.110, tAI 0.287 (matches report).
- Cliff's delta table: 0 differing rows; e.g. Kx core vs Mucorales enc_dep +0.87, gc3 +0.97, tai -0.72, outgroup -1.
- GC3-matched residuals: 0 differing rows. Kx core ENC departure residual +0.054, skew +0.095, size GC3 residual
  +14.3 points; Harpellales -0.017, -0.029, -22.9.
- ENC/GC3 recomputed: max abs difference 3e-14 versus `kx2_h3_cod`.
- Not re-verified by re-running: `kx2_h3_resid2`, the tRNA/tAI step (04), the aa/IDP step (05) and the protein
  proxy (06); they were run once and are the source of the reported values. All chk tables were dropped.
- Everything that was checked reproduced; nothing failed to reproduce.
