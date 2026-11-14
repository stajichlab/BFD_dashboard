# HS_I. Intron position within proteins, genome reduction and regulatory retention (pilot)

Status: **pilot**, run 2026-11-09 on the MotherDuck copy of BFD (`md:bfd`, snapshot 2026-11-06). Queries: `bin/atlas/sql/md_intron_position.sql`. No figures, no bootstrap CIs and no tree-based tests yet; group labels are taxonomic (growth-form labels are not in MotherDuck).

[Back to index](index.md) · extends [H3](H3_gene_structure.md) and [G3](HS_G3_intron_loss_by_origin.md)

## Data note: `gene_introns.codon_position` is wrong in this snapshot

`codon_position` and `codon_frame` hold the CDS offset at the end of the exon *after* the intron, so the last intron of every transcript reports the full CDS length and phase 0 (stored phase 0 share: 60.7 % of 508.6 M introns). In *S. cerevisiae* F0200A73, 0 of 392 introns match the position recomputed from `gene_CDS` (query I0). Fixed in `Fungi_BFD/nextflow/bin/build_genestats_table.py` with a fixture test (`nextflow/tests/test_genestats_intron_position.py`); `gene_introns` needs a rebuild. Everything below recomputes position from `gene_CDS` and does not use the stored columns.

## Hypotheses

| ID | Hypothesis | Supports | Falsifies |
|---|---|---|---|
| I1 | Intron loss removes 3' introns first (reverse-transcript-mediated loss), so surviving introns sit at the 5' end of the CDS | 5' share rises as intron density falls, in each lineage that lost introns | 5' share stays near 0.10 (uniform) in an intron-poor lineage |
| I2 | Intron phase shifts with loss | phase-0 share changes with density within clades | no relation |
| I3 | Retained 5' introns carry regulatory sequence, so they are longer than other introns in compact genomes | 5' / rest length ratio > 1 in intron-poor classes | ratio <= 1 |
| I4 | Retention is concentrated in highly expressed, regulated genes (ribosomal proteins) | RP genes keep introns more often than other genes, and the excess grows as density falls | RP = other at every density |
| I5 | Introns avoid (or mark) protein features: domain interiors, disordered linkers | observed / expected share inside Pfam envelopes or IDRs differs from 1 | ratio = 1 |

Filters: one genome per species (best BUSCO, then N50), BUSCO complete >= class 25th percentile, complete gene models. Unit = genus median. Position = coding bases 5' of the intron / CDS length; uniform expectation for "first tenth of the CDS" is 0.10.

## I1a. Position, phase and length by class (38 classes with >= 3 genera; selection)

| Class | Genera | Genome Mb | Introns / transcript | In first 10 % of CDS | In last 10 % | Within 10 codons of start | Phase 0 | 5' intron bp | Other intron bp |
|---|---|---|---|---|---|---|---|---|---|
| Pichiomycetes | 60 | 12.4 | 0.154 | 0.286 | 0.142 | 0.152 | 0.378 | 75 | 127 |
| Saccharomycetes | 31 | 12.1 | 0.207 | 0.374 | 0.127 | 0.187 | 0.403 | 78 | 119 |
| Dipodascomycetes | 21 | 12.4 | 0.386 | 0.439 | 0.063 | 0.159 | 0.410 | 66 | 63 |
| Microsporidia | 21 | 9.5 | 0.481 | 0.088 | 0.172 | 0.057 | 0.390 | 76 | 67 |
| Ustilaginomycetes | 12 | 18.9 | 0.491 | 0.184 | 0.118 | 0.069 | 0.400 | 99 | 93 |
| Kickxellomycetes | 6 | 21.0 | 0.793 | 0.215 | 0.079 | 0.066 | 0.422 | 118 | 96 |
| Dothideomycetes | 180 | 34.7 | 1.693 | 0.180 | 0.065 | 0.037 | 0.393 | 56 | 55 |
| Sordariomycetes | 240 | 38.9 | 1.803 | 0.175 | 0.068 | 0.035 | 0.389 | 69 | 65 |
| Eurotiomycetes | 52 | 30.7 | 2.126 | 0.183 | 0.065 | 0.036 | 0.399 | 66 | 62 |
| Mucoromycetes | 33 | 34.1 | 3.842 | 0.127 | 0.067 | 0.026 | 0.419 | 66 | 64 |
| Chytridiomycetes | 31 | 31.1 | 3.717 | 0.110 | 0.074 | 0.018 | 0.448 | 67 | 65 |
| Pucciniomycetes | 5 | 166.1 | 3.832 | 0.104 | 0.078 | 0.024 | 0.420 | 82 | 83 |
| Agaricomycetes | 291 | 44.4 | 5.136 | 0.098 | 0.070 | 0.016 | 0.431 | 55 | 55 |
| Tremellomycetes | 25 | 20.8 | 5.439 | 0.106 | 0.081 | 0.019 | 0.414 | 54 | 52 |

- Basidiomycota outside Ustilaginomycotina are at the uniform expectation (0.10 to 0.11), as are most early-diverging classes (0.10 to 0.14; Kickxellomycota 0.17 to 0.22). All Pezizomycotina classes sit at 0.15 to 0.19. Budding-yeast classes reach 0.29 to 0.44.
- Microsporidia: 0.48 introns per transcript but no 5' bias (0.088) and more introns in the last tenth (0.172). Tremellomycetes (20.8 Mb, yeasts) are indistinguishable from Agaricomycetes.

## I1b. Within-clade Spearman (genus medians; genomes with >= 50 CDS introns)

| Group | Genera | Introns / transcript | 5' share (p10 to p90) | rho density vs 5' share | rho genome size vs 5' share | rho genome size vs density | rho density vs phase 0 |
|---|---|---|---|---|---|---|---|
| Saccharomycotina | 123 | 0.18 | 0.332 (0.171 to 0.486) | 0.196 | 0.159 | 0.380 | 0.352 |
| Microsporidia | 15 | 0.76 | 0.094 (0.074 to 0.118) | 0.004 | -0.057 | 0.657 | -0.021 |
| Ustilaginomycotina | 31 | 1.30 | 0.168 (0.098 to 0.192) | -0.848 | -0.370 | 0.294 | 0.626 |
| Taphrinomycotina | 6 | 1.46 | 0.216 (0.122 to 0.265) | -0.886 | 0.200 | -0.486 | 0.771 |
| Pezizomycotina | 745 | 1.87 | 0.173 (0.154 to 0.196) | -0.503 | -0.529 | 0.029 | 0.308 |
| early diverging | 138 | 3.28 | 0.121 (0.100 to 0.170) | -0.590 | -0.222 | 0.190 | 0.165 |
| other Basidiomycota | 342 | 5.14 | 0.099 (0.093 to 0.107) | -0.151 | -0.334 | -0.471 | -0.163 |

## I4. Ribosomal-protein genes keep introns (27 classes with >= 5 genera; selection)

| Class | Genera | Introns / transcript | RP genes with intron | Other genes with intron | RP / other | Share of all introns in RP genes |
|---|---|---|---|---|---|---|
| Microsporidia | 9 | 0.01 | 0.016 | 0.009 | 0.86 | 0.012 |
| Pichiomycetes | 58 | 0.15 | 0.501 | 0.124 | 4.04 | 0.050 |
| Saccharomycetes | 30 | 0.19 | 0.568 | 0.143 | 3.83 | 0.047 |
| Dipodascomycetes | 21 | 0.39 | 0.510 | 0.271 | 1.61 | 0.015 |
| Ustilaginomycetes | 12 | 0.49 | 0.736 | 0.297 | 2.49 | 0.023 |
| Kickxellomycetes | 6 | 0.79 | 0.766 | 0.396 | 1.92 | 0.014 |
| Sordariomycetes | 240 | 1.80 | 1.000 | 0.758 | 1.32 | 0.010 |
| Eurotiomycetes | 52 | 2.13 | 1.000 | 0.792 | 1.26 | 0.010 |
| Mucoromycetes | 33 | 3.84 | 0.964 | 0.862 | 1.10 | 0.011 |
| Agaricomycetes | 288 | 5.14 | 0.985 | 0.921 | 1.07 | 0.004 |

RP genes = best Swiss-Prot hit named as a cytosolic ribosomal subunit protein; genomes with fewer than 40 such genes dropped. Median RP genes found per genome: 48 to 122.

## I5. Position relative to protein features (10 reference genomes, one per lineage)

| Genome | CDS introns | In Pfam obs / exp | In IDR obs / exp | RP with intron | Other with intron | RP intron bp | Other intron bp | RP 5' share | Other 5' share |
|---|---|---|---|---|---|---|---|---|---|
| *E. cuniculi* | 12 | 0.80 | 1.35 | 0.016 | 0.005 | 67 | 47 | 0.000 | 0.182 |
| *S. cerevisiae* | 287 | 0.79 | 1.03 | 0.651 | 0.053 | 408 | 92 | 0.607 | 0.338 |
| *Y. lipolytica* | 1,266 | 0.76 | 1.03 | 0.631 | 0.183 | 332 | 171 | 0.605 | 0.502 |
| *S. pombe* | 3,444 | 1.05 | 0.77 | 0.568 | 0.480 | 82 | 58 | 0.321 | 0.254 |
| *B. dendrobatidis* | 19,139 | 1.12 | 0.71 | 0.831 | 0.610 | 93 | 90 | 0.315 | 0.144 |
| *N. crassa* | 12,355 | 1.13 | 0.72 | 1.000 | 0.739 | 96 | 72 | 0.167 | 0.207 |
| *R. delemar* | 18,753 | 1.10 | 0.75 | 0.912 | 0.766 | 58.5 | 55 | 0.328 | 0.132 |
| *C. cinerea* | 53,442 | 1.19 | 0.64 | 0.985 | 0.883 | 55 | 56 | 0.198 | 0.114 |
| *P. jirovecii* | 5,375 | 1.09 | 0.81 | 1.000 | 0.944 | 45 | 45 | 0.096 | 0.139 |
| *C. neoformans* | 32,052 | 1.19 | 0.60 | 0.986 | 0.953 | 58 | 55 | 0.153 | 0.115 |

Observed / expected = share of introns that fall inside a Pfam envelope (or an IDR) divided by the share of residues covered in the same proteins. Phase-1 share inside vs outside Pfam envelopes differs by <= 0.04 in every genome with > 1,000 introns (no exon-shuffling signal at this resolution).

## Reading (pilot)

- **I1 supported in Ascomycota and Ustilaginomycotina, not a general rule of reduction.** Within Pezizomycotina (745 genera), Ustilaginomycotina (31) and the early-diverging groups (138), genera with fewer introns have a larger 5' share (rho -0.50, -0.85, -0.59). Microsporidia lost introns without any 5' bias (0.094; rho 0.00), and in Pezizomycotina smaller genomes have the larger 5' share (rho -0.53) although genome size and intron density are unrelated there (rho 0.03). Inside Saccharomycotina the relation is weak and has the other sign (rho 0.20): the 5' bias is already high across the subphylum.
- **I2:** phase-0 share rises with intron density inside Ustilaginomycotina (0.63), Saccharomycotina (0.35) and Pezizomycotina (0.31); class medians span only 0.36 to 0.51.
- **I3 not supported as stated; supported for RP genes.** In Pichiomycetes and Saccharomycetes 5' introns are *shorter* than the rest (75 vs 127 bp, 78 vs 119 bp); elsewhere the ratio is 0.94 to 1.16. But RP introns are long in the two intron-poor yeasts checked (408 vs 92 bp in *S. cerevisiae*, 332 vs 171 bp in *Y. lipolytica*) and equal to other introns in intron-rich genomes.
- **I4 supported.** RP genes keep introns 3.8 to 4.0 times as often as other genes in Saccharomycetes and Pichiomycetes, 2.5 times in Ustilaginomycetes, 1.0 to 1.4 times in classes with more than 1.5 introns per transcript. Microsporidia show no RP excess (0.86).
- **I5:** introns are over-represented inside Pfam envelopes (1.05 to 1.19) and under-represented in IDRs (0.60 to 0.81) in the seven intron-rich genomes; the pattern reverses in the three intron-poor genomes (Pfam 0.76 to 0.80, IDR 1.03 to 1.35), as expected if the survivors sit near the N-terminus, ahead of the first domain. Ten genomes only.

## Caveats

- Ab initio gene models without RNA-seq (D8). Short introns close to the start codon can be prediction artefacts: 15 to 19 % of budding-yeast introns are within 10 codons of the start, and those classes are where 5' introns are shorter than the rest. Real yeast RP introns are also start-proximal, so this needs the RP / non-RP split at scale.
- Microsporidia: 21 genera in I1a (0.48 introns per transcript) but the 9 genera with >= 40 RP genes found have 0.01. The high-intron microsporidian annotations are likely poor; treat the Microsporidia rows as provisional.
- CDS introns only (UTR introns are not in `gene_CDS`), complete models only. *S. cerevisiae* has 3,986 complete models of 5,441 with a CDS; partial-model filtering removes about a quarter of genes there.
- No phylogenetic correction; Taphrinomycotina is 6 genera; Saccharomycotina is one origin.

## Next

1. Rebuild `gene_introns` with the fixed `build_genestats_table.py`, then rerun I0.
2. Move I1 to I5 into an HPCC job (`bin/atlas/sql/J11_intron_position.sql`, per-genome parquet in `atlas_data/`), with growth-form and `hs_group` labels, bootstrap CIs, figures and the G3 sister pairs.
3. I5 at scale with interval merging; add signal-peptide cleavage sites (`signalp_signal_peptide`) and intron length against median intergenic distance per genome (regulatory space moved into introns in compact genomes).
4. Same-ortholog test: for Swiss-Prot anchored families, compare intron positions in intron-poor vs intron-rich relatives (loss vs gain, sliding).
