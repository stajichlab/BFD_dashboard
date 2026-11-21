# Carbon-storage family trees (one genome per genus)

Status: **first-pass trees and the TP1 test done 2026-11-16; the four IQ-TREE trees added 2026-11-21** (see Results). The design and decision
rules were written, and committed, before any tree was built.
Parent: [`../README.md`](../README.md) (module definitions), [`../t5/`](../t5/) (the GH65 tree, same pipeline and settings).

## Why

The trehalose survey found that copy number, not presence, separates many lineages (for example a median of 7 TPS and 7 TPP copies in
every genus of Mucoromycetes and 5 in Mortierellomycetes, against 3 and 2 in most classes). A count cannot say whether extra copies are
**ancient paralogs kept by every genus** or **recent lineage-specific duplications**. A gene tree can.

## Families and sequences

One representative genome per genus (`is_genus_rep`, annotated; 1,669 genomes). All proteins carrying the family's Pfam domain
(`domain_i_evalue <= 1e-5`, version stripped), trimmed to the best-scoring domain envelope +-10 aa. A protein shorter than 60% of the
family's median trimmed length is dropped.

| Family | Pfam | Sequences (counted before design) | Method |
|---|---|---|---|
| TPS | PF00982 | 5,226 | FastTree first pass |
| TPP | PF02358 | 3,639 | FastTree first pass |
| NTH (neutral trehalase) | PF01204 | 3,885 | FastTree first pass |
| UGP (UDP-glucose pyrophosphorylase) | PF01704 | 3,866 | FastTree first pass |
| GSY (glycogen synthase) | PF05693 | 1,839 | FastTree + IQ-TREE |
| GPH (glycogen phosphorylase) | PF00343 | 1,742 | FastTree + IQ-TREE |
| GDB (glycogen debranching) | PF06202 | 1,897 | FastTree + IQ-TREE |
| GPD (glycerol-3-phosphate dehydrogenase) | PF01210 | 2,262 | FastTree + IQ-TREE |

Not built: **PGM** ([PF02878](https://www.ebi.ac.uk/interpro/entry/pfam/PF02878/), 4,564 sequences, but the domain is only about 75 aa, too short for a
reliable tree) and **GH65** (already built in `../t5/`). Accessions are linked to InterPro on first mention in the parent README.

## Pipeline (fixed before the run, as in T5)

mafft 7.505 `--auto`; trimAl 1.5.1 `-gt 0.4 -resoverlap 0.5 -seqoverlap 50`; FastTree 2.1.11 `-lg -gamma` for all families;
IQ-TREE 3.0.1 `-m Q.PFAM+R8` with 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates for the four families of at most 2,300
sequences. Large families get only the first-pass tree (no bootstrap; FastTree local support values are reported but are not
bootstrap values). Trees are unrooted: no outgroup was added, so root-dependent statements are not made.

## Question and decision rule (one pre-registered test)

**TP1.** In classes with elevated copy number (**Mucoromycetes, Mortierellomycetes**) the extra TPS and TPP copies are **ancient
paralogs retained by every genus** rather than recent duplications.

- Measure: for each genus-representative genome with at least 2 copies, the median patristic distance between its own copies; the
  reference is the median patristic distance between copies drawn from **different genera of the same class** (at most 2,000
  random pairs per class, fixed seed). **Paralog index PI** = median over genomes of (within-genome distance / between-genus
  distance). Computed per family and class.
- **Supported** (ancient paralogs): PI >= 0.8 in both classes and both families. **Killed** (recent duplications): PI < 0.5 in
  either class for either family. Otherwise **mixed**.
- Baseline for context, not part of the rule: PI for Sordariomycetes, Agaricomycetes and Eurotiomycetes (3 or 4 copies).
- A within-genome distance below the between-genus distance means copies are more alike than different genera are, the
  signature of lineage-specific duplication; a ratio near or above 1 means copies diverged before the genera did.

The same index is reported for the other families and classes with at least 8 genera, as description only.

## Limits stated in advance

- First-pass FastTree trees are approximate and unrooted; branch lengths are model-based estimates, so PI carries error that is not
  quantified here (a bootstrap of PI over genomes is reported, not over trees).
- Gene-model errors (split or merged genes) can inflate or deflate within-genome distances.
- One genome per genus does not make genera independent; the classes are not independent either.
- No outgroup, so paralog clades are described by tree topology within the fungal set only; claims about the origin of TPS and TPP
  paralogs (or about function) are citations to verify and none is made here.

## Files

`01_extract_domains.py` (HPCC job 32132586, 10 min), `02_trees.sh` (FastTree array, job 32132838), `03_iqtree.sh` (IQ-TREE array, job 32132839,
all four tasks COMPLETED, the last on 2026-11-19), `04_paralog_index.py`, `test_carbon_trees.py` (5 fixture tests pass on the HPCC), `results/`, `trees/`
(FastTree and IQ-TREE trees, logs; the IQ-TREE `.iqtree` reports and run logs are gzipped). Alignments and the tip table (`tips.tsv`, 22,276 kept sequences from 24,356 domain hits) stay in
`bench/trehalose_carbon_trees/` on the HPCC.

---
## Results (2026-11-16)

Order of events: design and TP1 rule committed first (`93cc3b9`), scripts next (`e8fc462`), then the trees and one analysis run. Two
bugs were fixed before the run (a path to `fam_lib` one directory too high; a too-tight bound in a fixture test). After the TP1
result one **post hoc** statistic was added (closest within-genome pair, below) and is labelled as such.

### Trees built (first pass, FastTree 2.1.11 `-lg -gamma`, unrooted)

| Family | Sequences | Alignment columns (raw to trimmed) |
|---|---|---|
| TPS | 4,942 | 3,536 to 472 |
| TPP | 3,529 | 2,038 to 256 |
| NTH | 2,886 | 3,146 to 576 |
| UGP | 3,718 | 1,617 to 489 |
| GSY | 1,742 | 1,485 to 665 |
| GPH | 1,551 | 1,802 to 736 |
| GDB | 1,675 | 1,354 to 500 |
| GPD | 2,233 | 824 to 190 |

Sequences dropped by the pre-set 60% length rule: TPS 284, TPP 110, NTH 999 (a quarter of the neutral-trehalase domain hits are much
shorter than the median), UGP 148, GSY 97, GPH 191, GDB 222, GPD 29. trimAl kept every sequence of every family. The IQ-TREE trees for GSY, GPH,
GDB and GPD are below.

### IQ-TREE maximum-likelihood trees (added 2026-11-21)

IQ-TREE 3.0.1, `-m Q.PFAM+R8`, 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates, seed 20261116 (job 32132839, four array tasks,
5 h 58 m to 2 d 10 h). Only these four families were run: the pre-set rule was families of at most 2,300 sequences, which excludes
TPS (4,942) and TPP (3,529).

| Family | Sequences | Columns | UFBoot >= 95 | SH-aLRT >= 80 | RF vs FastTree | normalised RF |
|---|---|---|---|---|---|---|
| GSY | 1,742 | 665 | 1,050 / 1,736 (60.5%) | 68.2% | 1,055 | 0.303 |
| GPH | 1,551 | 736 | 1,055 / 1,547 (68.2%) | 70.3% | 817 | 0.264 |
| GDB | 1,675 | 500 | 1,027 / 1,672 (61.4%) | 68.4% | 1,006 | 0.301 |
| GPD | 2,233 | 190 | 1,245 / 2,227 (55.9%) | 60.8% | 1,779 | 0.399 |

Two things to read off this, both of which argue against leaning on any one topology.

**The ML and FastTree topologies disagree substantially.** Between a quarter and two fifths of the splits differ (normalised
Robinson-Foulds 0.26 to 0.40). The first-pass FastTree trees are not a close approximation of the ML trees for these families.

**Neither is well resolved.** Only 56% to 68% of nodes reach UFBoot >= 95 even in the ML trees. GPD is the weakest on both counts,
which its alignment explains: 2,233 sequences over 190 trimmed columns, the least information per sequence of any family here.

**TP1 is unaffected, for two independent reasons.** It is a distance test, not a topology test -- the paralog index is a ratio of
median within-genome copy distance to between-genus distance, computed from branch lengths. And it is a test about TPS and TPP,
neither of which has an ML tree, both being over the 2,300-sequence bound. Nothing in this section revises the TP1 verdict; the
trees stand on their own as first-pass topologies for four carbon-storage families.

### TP1 (pre-registered): supported, by the rule

PI (median within-genome copy distance / between-genus distance) is **at or above 0.8 in all four test cells**: TPS 1.11 (Mucoromycetes) and
1.05 (Mortierellomycetes); TPP 1.20 and 1.33 (all 95% intervals above 1.0). By the rule, the extra copies are not recent
duplications; copies of one genome are as different as copies of different genera.

**That answer is not specific to the high-copy classes.** The baseline classes with 3 to 4 copies give the same values (TPS Sordariomycetes
1.07, Agaricomycetes 1.03, Eurotiomycetes 1.05; TPP Eurotiomycetes 1.08, Sordariomycetes 1.13, Agaricomycetes 2.30). PI of about 1 or
above is expected whenever a genome carries copies from different paralog groups: a within-genome pair is always two different copies,
whereas a random between-genus pair often falls in the same group. So TP1 shows that distinct deep paralog groups exist everywhere; it does not
separate Mucoromycetes from the other classes, and a median over pairs cannot see a few recent duplicates among many ancient paralogs.

### Post hoc supplement (not part of the rule): a recent duplicate pair is common only in the high-copy classes

Share of genomes whose **closest** pair of copies is nearer than half the between-genus reference distance:

| Class | TPS | TPP |
|---|---|---|
| Mucoromycetes (36 genomes) | 89% | 86% |
| Mortierellomycetes (14) | 86% | 93% |
| Sordariomycetes (about 270) | 9% | 4% |
| Agaricomycetes (about 300) | 17% | 18% |
| Eurotiomycetes (71) | 21% | 1% |

Read as suggestive: in the two high-copy classes most genomes have at least one recently duplicated pair on top of the ancient paralog groups,
which would explain the elevated copy number. **Caveats:** the 0.5 cutoff was chosen after seeing the index; a genome with 7 copies has 21 pairs and so is
more likely to contain a small distance by chance than a genome with 2 copies, so the contrast is partly mechanical; a copy-number-matched
comparison was not run; gene-model errors (split or merged genes) can create close pairs.

### Not done

Rooting or outgroups; identification of the number of deep paralog groups per class; a copy-number-matched test of the
closest-pair statistic; trees for PGM; any claim about function or origin (citations to verify, none made).
