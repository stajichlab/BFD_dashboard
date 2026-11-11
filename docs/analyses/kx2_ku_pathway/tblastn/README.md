# tblastn check of the Ku-pathway absences

Purpose: BFD gene annotation under-calls genes, so "Ku absent" from Pfam is weak evidence. This searches the assemblies directly with
fungal reference proteins. It covers the 8 genes on the Q1 worklist (Ku, DNA-PKcs, Lig4, Pol4, XLF/Nej1, aprataxin, ATR/Mec1,
Sae2/CtIP) plus a Rad51 control, in the 31 worklist genomes (24 Kickxellomycota and sister genomes with absences to check, including
9 *Coemansia* strains, plus positive controls: Ku in *Rhizopus*, *Dimargaris* and *Ramicandelaber*; Rad51 in the stem genera). The
worklist is `../sql_coloss/14_tblastn_worklist.sql`.

| File | Role |
|---|---|
| `genomes.tsv`, `queries.tsv` | inputs (genomes from the worklist; 23 query proteins, 22 from UniProt plus one *Actinomortierella* DNA-PKcs from BFD) |
| `fetch_queries.sh` | builds `queries.faa` (checks each UniProt length against `bfd.swissprot_subject`) |
| `run_tblastn.sh` | SLURM array, one task per genome, `tblastn` against `input_clean_genomes/<ASMID>.fa.gz` (unmasked; masked file as fallback), BLAST+ 2.17.0. The `input/dna/*.scaffolds.fa` links are broken (they point at deleted `genome_annotation/.../predict_results/` copies), which stopped the first submission |
| `summarize_tblastn.py` | per genome and gene: strong / weak / none |

## Settings and thresholds (fixed before the results were seen)

- Search: e-value 1e-3 (permissive), up to 50 targets per query, BLAST defaults for the rest (including low-complexity masking).
- **strong** = best e-value <= 1e-10 and merged HSPs cover >= 30% of the query on one scaffold. **weak** = e-value <= 1e-5 but not
  strong. **none** = otherwise. Coverage merges all HSPs of a query on its best scaffold, so a gene split by introns counts once.
- The best query per gene is the one with the highest coverage x bit score, so one close reference is enough to count.

## How to read it

- A **strong** hit in a genome Pfam called absent means the gene is present in the assembly and was missed by annotation.
- **None** in a contiguous, complete assembly, with the positive controls recovered, supports a true loss. None in a fragmented
  assembly is still weak, because the region could be missing from the assembly. A diverged gene could also fall below the
  thresholds, which is why the **weak** class and the controls matter.
- The positive controls must be strong (Ku in *Rhizopus*, *Dimargaris*, *Ramicandelaber*; Rad51 in the stem genera). If a control
  is not strong, the threshold is too strict for these genomes and the negative calls cannot be trusted.

## Run (HPCC)

    BENCH=~/projects/BFD_dashboard/bench/kx2_ku_tblastn REPO=<checkout> bash fetch_queries.sh      # login node, light
    cd $BENCH && BENCH=$BENCH REPO=<checkout> sbatch $REPO/docs/analyses/kx2_ku_pathway/tblastn/run_tblastn.sh
    python3 summarize_tblastn.py $BENCH genomes.tsv results/tblastn

Hit tables live in `$BENCH/hits/` (purgeable); the summaries are copied into `results/`.
