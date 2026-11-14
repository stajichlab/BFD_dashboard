# Intron cross-check: is the low intron count of Coemansia an artifact of the gene models?

Status: **design fixed before any group result was looked at** (2026-11-11); results are added below the line "Results".

Question. Analyses A and E report about 0.43-0.65 introns per transcript in core Kickxellales against 3.8 in Mucorales, from
funannotate gene models. Is the deficit still there when introns are called from spliced protein-to-genome alignments (miniprot)
of **one fixed reference protein set on every assembly**, which do not use the gene finder? (Gene completeness is not the issue:
protein-mode BUSCO equals genome-mode BUSCO in *Coemansia*, median gap -0.2 points.)

## Design (fixed in advance)

| Choice | Value | Why |
|---|---|---|
| Reference proteins | the 1,122 BUSCO `fungi_odb12` ancestral proteins (`.../phyling_db/fungi_odb12/ancestral`, median length 448 aa) | one fixed, taxon-neutral set, so no group is favoured by its own proteins |
| Genomes (83) | the 63 backbone genus representatives, 12 *Coemansia* strains with N50 >= 100 kb and BUSCO >= 60 (within-genus spread), 4 yeast anchors (*Lachancea*, *Petasospora*, *Saccharomycopsis*, *Williopsis*; intron-poor) and 4 Pezizomycotina anchors (*Blumeria*, *Cenococcum*, *Clonostachys*, *Trichophyton*; intron-rich) | the anchors calibrate the method, they are not part of the Kickxellomycota question |
| Assemblies | `Fungi_BFD_runs/input_clean_genomes/<ASMID>.fa.gz` (the `input/dna` links are broken); the scaffold count must equal `bfd.asm_stats.contig_count` | same check as the tblastn work |
| Aligner | miniprot 0.2, `--gff`, best-scoring alignment per reference protein, `Rank=1` only | no gene finder |
| Splice settings | **A** default (`-J 31 -C 11 -G 200k`); **B** permissive (`-J 16 -C 6 -G 5000`); **C** very permissive (`-J 8 -C 3 -G 5000`) | miniprot penalises introns; if that disfavoured introns it would bias every group towards "few introns". The test is whether the *Coemansia* vs Mucorales difference survives lowering the penalty. A is the primary setting because it is the tool default. |
| Qualifying alignment | query coverage >= 0.80 and miniprot `Identity` >= 0.50 (sensitivity: 0.40 and 0.60) | divergent references align less well in diverged lineages and can hide introns, so identity is also used as a stratifier (bins 0.4-0.5, 0.5-0.6, 0.6-0.7, >= 0.7) |
| Intron call | a gap between consecutive CDS blocks of >= 20 nt; gaps of 1-19 nt are counted separately as frameshift or short gaps and not as introns | miniprot splits a CDS at frameshifts |
| Metrics per genome | median introns per qualifying protein; share of proteins with >= 1 intron; introns per kb of aligned CDS; number of qualifying proteins | |
| Unit and effect size | one value per genus (backbone representative); Cliff's delta of each group vs Mucorales (8 genera); the 12 strains only show within-genus spread | genus medians are the atlas replication unit |
| Splice sites | donor/acceptor dinucleotides read from the assembly; shares of GT-AG, GC-AG and other | a high non-canonical share would suggest artefacts |
| Annotation comparison | funannotate introns per transcript (all transcripts) per genome from `gene_introns` and `gene_transcripts`; Spearman correlation across genomes; per-group medians side by side | |
| Robustness | bootstrap over proteins (1,000 resamples per genome) to show the result does not rest on a few proteins | |
| Minimum evidence per genome | a genome with fewer than 50 qualifying proteins at a given setting and identity threshold is excluded from group comparisons there (and listed); in identity bins at least 15 proteins | a genome with a handful of aligned proteins must not drive a group result |
| Success rule for "the deficit is supported independently of the gene models" | in settings A, B and C, and at identity thresholds 0.40, 0.50 and 0.60, core Kickxellales are below Mucorales (Cliff's delta <= -0.5) and still below after stratifying by identity bin | |

## Scripts (numbered, run on the HPCC)

1. `01_prepare_genomes.py` adds assembly ID, contig count, N50, BUSCO to `genomes.tsv` (reads the BFD parquet tables).
2. `02_run_miniprot.sh` (SLURM array, one task per genome): decompress to `$SCRATCH`, check the scaffold count, run miniprot with
   settings A, B, C, parse with `miniprot_introns.py`, write one small table per genome and setting.
3. `miniprot_introns.py` the GFF parser (library and command line); `test_miniprot_introns.py` is the fixture test.
4. `03_annotation_introns.py` funannotate introns per transcript per genome.
5. `04_summarize.py` metrics, group tables, effect sizes, bootstrap, annotation comparison.

Limits stated now: ancestral consensus proteins are diverged from every genome, so alignment sensitivity varies with lineage;
miniprot may miss short exons or very short introns; this tests intron calls on 1,122 conserved single-copy genes, not all genes,
and conserved genes tend to be intron-rich, so absolute counts are higher than genome-wide averages and only the comparison is
meaningful. No literature is cited. Anything about miniprot behaviour beyond what the help text states is a citation to verify.

---
## Results

Run 2026-11-11/12 on the HPCC (miniprot array, jobs 32016275 and 32031682; summary job 32059464, `04_summarize.py`). 83 genomes, 1,122
reference proteins. Small tables are in [`results/`](results/) (`summary.txt` is the full printed output); the per-genome miniprot tables and
GFFs stay in the purgeable HPCC folder `bench/kx2_intron_crosscheck/`. Scripts are numbered above; the parser has a fixture test
(`test_miniprot_introns.py`). The first run used a buggy parser and was discarded (`logs_run1_buggy_parser/` on the HPCC); all numbers here come from the second run.

### Verdict (observed)

**The pre-registered success rule is met.** In all 9 setting × identity-threshold cells, core Kickxellales (5 genera) have Cliff's delta = -1.00
against Mucorales (8 genera) on median introns, share of proteins with at least one intron, and introns per kb; so does every identity bin of setting A
(-1.00 in 0.4-0.5, 0.5-0.6, 0.6-0.7 and 0.7-1.0). In the bootstrap over reference proteins (1,000 resamples, settings A and C) the
share of resamples with delta <= -0.5 is 1.00 for core Kickxellales, Kickxellales all and Harpellales. The low intron count of core Kickxellales is therefore **not
an artifact of the funannotate gene models**, and it does not depend on the aligner's intron penalty. With only 5 and 8 genera per group, a delta of -1.00 means complete
separation of the two sets of genus values, not a large sample.

### Setting A (default), identity >= 0.5, group medians of genus values

| Group | genera | median introns / protein | share with >= 1 intron | introns per kb CDS | funannotate introns / transcript |
|---|---|---|---|---|---|
| Kickxellales core | 5 | 0 | 0.29 | 0.30 | 0.65 |
| Dimargaritales + *Ramicandelaber* | 4 | 0.5 | 0.38 | 0.54 | 0.94 |
| Zoopagales | 8 | 1 | 0.61 | 1.04 | 3.29 |
| Harpellales | 4 | 1 | 0.63 | 1.21 | 1.98 |
| *Mycoemilia* + *Spiromyces* | 2 | 1 | 0.67 | 1.30 | 1.70 |
| Entomophthorales + *Basidiobolus* | 4 | 1.25 | 0.73 | 1.51 | 2.62 |
| Mucorales | 8 | 2 | 0.94 | 2.10 | 3.96 |
| Other Mucoromycota group | 13 | 4 | 0.97 | 3.36 | 3.48 |
| Root (chytrids etc.) | 15 | 3 | 0.90 | 2.46 | 4.35 |
| *Saccharomycetes* anchors (intron-poor) | 4 | 0 | 0.19 | 0.16 | 0.21 |
| Pezizomycotina anchors (intron-rich) | 4 | 1.5 | 0.82 | 1.06 | 1.83 |

(The funannotate column uses all transcripts, the miniprot columns only 1,122 conserved proteins, so compare within a column.)
Core Kickxellales sit at about twice the yeast anchors and one seventh of Mucorales per kb. Within *Coemansia* (13 genomes: the representative and 12 contiguous strains)
the spread is small: 0.21-0.28 introns per kb and 23-28% of proteins with an intron.

### Checks (observed)

1. **The intron count does not rise at higher identity, so a divergence artifact is unlikely.** Mean introns per protein in core Kickxellales are 0.56, 0.40, 0.36 and 0.30 in the
   four identity bins (0.4 to 1.0), against 3.70, 3.27, 2.35 and 2.32 in Mucorales. If poorly aligned divergent proteins hid introns, the well-aligned bins would
   look intron-richer.
2. **Agreement with the gene models.** Spearman rho between miniprot introns per kb and funannotate introns per transcript is 0.86 across all 83 genomes
   (0.72 in the 63 backbone genera; setting C 0.87 and 0.78). The calibration anchors behave as expected: yeasts 0.15-0.19 per kb, Pezizomycotina 0.79-1.24. miniprot gives
   lower values than the annotation in the intron-rich anchors (for example *Trichophyton* 1.24 per kb against 2.02 introns per transcript), so the method is, if anything, conservative.
3. **The penalty dependence was real but is spurious introns.** Lowering the intron penalty raises counts in every group: core Kickxellales introns per kb are 0.30 (A), 0.70 (B) and 4.50 (C); Mucorales 2.10, 3.00 and 6.60. The ratio narrows from 0.14 (A)
   to 0.68 (C). In setting C, the share of called introns with a GT-AG site falls to 0.12 in core Kickxellales and 0.33 in Mucorales (85% and 65% non-canonical), so most introns called there are
   alignment gaps, not splice junctions. The early signal of 5.59 introns per protein in *Coemansia* under setting C was this effect. Counting only canonical (GT-AG or GC-AG) introns,
   which was **added after seeing this and is post hoc**, core Kickxellales have 0.24, 0.38 and 0.67 canonical introns per kb in A, B and C against 1.99, 2.50 and 2.26 in Mucorales (delta -1.00 in every cell).
4. **Splice sites in the primary setting.** Core Kickxellales introns are 79% GT-AG and 18% other (Harpellales 93%, Mucorales 95%, Pezizomycotina anchors 95%, yeast anchors 69%). A higher non-canonical share fits either
   diverged splice sites or a larger share of false calls among few true introns; it cannot be told apart here. The yeast anchors show the same pattern at very low counts.
5. **The gradient is wider than the core clade.** Dimargaritales + *Ramicandelaber* are as intron-poor as core Kickxellales on this measure (0.54 introns per kb; delta against Mucorales -1.00 in every setting-A
   cell), Harpellales, Zoopagales and *Mycoemilia* + *Spiromyces* are about 1.0-1.3, and Entomophthorales + *Basidiobolus* 1.5. Intron loss therefore seems to begin before the core Kickxellales (this is a descriptive reading
   of group medians, not a reconstruction on the tree).

### Caveats

- Group sizes are small (5 core genera, 8 Mucorales genera); delta is a rank statistic and saturates at -1.00. The bootstrap resamples reference proteins only, not genomes, so it does not measure genus sampling.
- Setting C is not a usable intron caller (mostly non-canonical calls, see 3). Its role was to test whether penalty alone produces the deficit; it does not. miniprot aborts (core dump) on *Cryptocalciella* under setting C, also on a re-run (job 32031682, task 25), so that genome has settings A and B only and is excluded from the C rows (12 instead of 13 genera in "other Mucoromycota"). The cause was not investigated.
- The 1,122 ancestral BUSCO proteins are conserved single-copy genes and are diverged from every genome, so alignment sensitivity varies by lineage; whether the lower per-kb values hold for all genes (as the annotation suggests, 0.65 against 3.96 per transcript) is not tested here.
- miniprot's known behaviour beyond the help text (for example how it scores short exons) is not cited; citation to verify.
- The success rule and thresholds were fixed before the runs. Items marked post hoc (canonical-only counts, the splice-site reading) were added afterwards.
