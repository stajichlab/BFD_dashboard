# RecA/Rad51 family for all 107 Coemansia strains: is Dmc1 absent from every strain?

Date: 2026-11-15. Status: **done for the annotated layer.** The earlier RecA tree ([`../RECA_RESULT.md`](../RECA_RESULT.md)) had only 10 of the 107 *Coemansia* strains.
Here the same pipeline, rules and model are run on all 128 genomes of the tblastn round-3 set (the 107 strains plus 21 other genomes), 593 annotated proteins plus the 31 references (**624 sequences**, 332 alignment sites, Q.PFAM+R8, UFBoot 1,000, SH-aLRT 1,000, seed 20261111).
HPCC folder `bench/kx2_reca_all/` (purgeable). Small tables are in [`results/`](results/).

## Answer

**No *Coemansia* strain has a Dmc1: 0 of 107.** All 301 annotated RecA-fold proteins from the 107 strains are placed in the tree with Rad51, Rad57 or XRCC2 as the nearest reference
(177, 96 and 28 proteins); **none has Dmc1 as its nearest reference** and none is assigned to Dmc1. The same classification finds Dmc1 where it is known to be present: assigned to Dmc1 in
Harpellales (4 proteins, 4 of 4 genomes), Dimargaritales (4 proteins, 3 of 3), Entomophthorales + *Basidiobolus* (4 of 4), Mucoromycota (18 proteins, 18 of 21 genomes) and root taxa (12 of 15), and in *Mycoemilia*.
In the whole core clade (107 *Coemansia* + *Kickxella*, *Linderina*, *Dipsacomyces*, *Martensiomyces*) Dmc1 is not seen in any of 111 genomes. This extends the earlier result from 10 to all 107 strains (`results/reca_by_genome.tsv`, `reca_by_group.tsv`).

## What the tree can and cannot call (leave-one-out, 17 of 31 references)

| Clade | references correct | Use |
|---|---|---|
| Dmc1 | **5 of 5** | usable: this is what the absence call rests on |
| Rad51B, Rad51C | 3 of 3 each | usable |
| Rad51D | 2 of 3 | usable with care |
| RadA, XRCC2 | 2 of 2 each | usable |
| Rad51 | 0 of 7 | **not callable**: the 177 Coemansia proteins nearest to Rad51 are ambiguous, as expected |
| Rad55, Rad57, XRCC3 | 0 of 2 each | not usable (Rad57 passed 2/2 in the earlier 362-sequence tree and now fails: the larger tree is not better for every clade) |

## Tblastn locus count next to the tree call (independent of the annotation; `results/loci_per_genome.tsv`)

| Loci hit by Rad51/Dmc1 queries (e-value <= 1e-25, >= 150 aa, HSPs merged within 20 kb) | *Coemansia* strains |
|---|---|
| 1 | 97 |
| 2 | 6 (identical Rad51-like score at both loci, no Dmc1-like second locus: a duplicated Rad51 locus, as before) |
| 0 | 4 (no hit at all, so the assembly or the annotation is the limit; uninformative) |

In Dmc1-positive genomes the second locus scores higher with the Dmc1 queries than with Rad51 queries (for example *Dimargaris*: Rad51-queries 504 / Dmc1-queries 313 at the first locus, 0 / 120 at the second). No *Coemansia* strain shows that pattern.

## By assembly quality (`results/coemansia_strain_reca.tsv`)

| BUSCO complete | strains | median annotated RecA-fold proteins | tblastn loci |
|---|---|---|---|
| >= 90 | 21 | 4 | 1 in all 21 |
| 80-90 | 33 | 3 | 1: 29, 2: 3, 0: 1 |
| 70-80 | 21 | 2 | 1: 20, 2: 1 |
| < 70 | 32 | 2 | 1: 27, 2: 2, 0: 3 |

Every strain, including all 21 of the best assemblies, has an annotated RecA-fold protein (minimum 1, no strain with zero) and no Dmc1, so the absence does not follow assembly quality.

## A labeling trap in the per-strain table (read before using `n_dmc1_weak`)

`strain_summary.py` reports `n_dmc1_weak` and a `dmc1_call` of "uncertain" for 101 of 107 strains. This counts ambiguous proteins whose ambiguous ancestor clade **contains** Dmc1 among many clades
(`mixed_labels` = Dmc1, Rad51, Rad51B, ... XRCC3). It does not mean the protein is Dmc1-like: `reca_classification.tsv` shows 0 of 301 *Coemansia* proteins with a Dmc1 nearest reference. Use the nearest reference and `reca_by_genome.tsv` (all 107 `not_seen_in_annotation`, written by `07_classify.py`), not the "uncertain" label.
The strain-summary label should be changed to count only proteins whose nearest reference is Dmc1; not done here.

## Scripts and jobs

Steps 02 (candidates), 03 (extract proteins) and 05 (build inputs) are the scripts in `../`, run on the 128 genomes (jobs 32074396, 32074399, 32074400). Here: `06_tree_reca_all.sh` (mafft `--auto`, trimAl `-gt 0.4 -resoverlap 0.5 -seqoverlap 50`, IQ-TREE 3.0.1 with the earlier
model, job 32074401, 1 h 24 min), `loci_per_genome.py` (job 32074402), `../07_classify.py` run as `NO_RESCUE=1 LAYER=annotated` (job 32102316; `NO_RESCUE` was added so that genomes not searched with miniprot get no "absent" call),
`strain_summary.py` (job 32102346), fixture test `test_reca_all.py`. Versions: mafft 7.505, trimAl 1.5.rev1, IQ-TREE 3.0.1.

## Limits

- Annotated layer only; the rescue (miniprot) layer was not run for these strains. Absence in annotation can be a missed gene, but the tblastn locus count agrees (one locus in 97 strains) and so does the earlier result for the 10-strain tree.
- Rad51 cannot be called by the tree, so "one Rad51-like protein per strain" rests on annotation and tblastn, not on a tree call.
- One ML tree; no site-heterogeneous model check; four strains have no tblastn locus and are uninformative.
- The earlier result for *Spiromyces* (BUSCO 44%) and *Ramicandelaber* (one genome) stays weak; higher-quality assemblies were not available.
