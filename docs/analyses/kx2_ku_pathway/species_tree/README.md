# Improved Kickxellomycota species trees: ANI dereplication, site-heterogeneous model, gene-tree summary, rooting

Date: 2026-11-13/14. Status: **partly complete.** Done: ANI dereplication, LG+F+G trees with UFBoot, ASTRAL gene-tree summaries, rooting tests,
comparison with the first-pass FastTree trees. **Not finished: the site-heterogeneous (LG+C60+F+G, PMSF) trees** — see "What is still running".
Scripts: `workflows/phylo/` (`kickx_derep.sh/.py`, `kickx_subset_aln.py`, `kickx_iqtree.sh`, `kickx_iq_submit.sh`, `kickx_tree_tools.py`,
`kickx_tree_report.py`, test `test_kickx_phylo.py`; run on the HPCC, purgeable bench folder `bench/kickx_tree/iq/`).
Inputs are the nf_phyling alignments of the first-pass run (80 `fungi_odb12` markers; the same per-marker alignments, only taxa are removed).

| Folder | Contents |
|---|---|
| [`derep/`](derep/) | fastANI output, clusters at 99, 97 and 95 % ANI, retained taxa table, summary |
| [`coemansia/`](coemansia/) | 76-taxon tree set (62 dereplicated *Coemansia* + 14 other Kickxellomycota genera): report, clade tests, rooted Newick files, figure (PDF/PNG), LG+F+G tree, ASTRAL trees |
| [`backbone/`](backbone/) | 63-taxon backbone, same set of files (outgroup: chytrid-grade fungi) |

## Methods and choices (fixed before the trees were viewed unless marked)

1. **Dereplication.** fastANI 1.34 (default k-mer 16, fragment 3,000 nt), all-vs-all on the 107 cleaned *Coemansia* assemblies
   (`Fungi_BFD_runs/input_clean_genomes`). **Rule: two strains are redundant at ANI >= 99.0 %** (either direction); clusters are connected
   components; the representative has the highest BUSCO complete %, then the largest N50; the 14 non-*Coemansia* genera are always kept. The 99 %
   threshold was chosen to remove near-clonal strains (the first-pass tree had 20 of 29 closest pairs under 0.005 substitutions per site); 97 and
   95 % are reported as sensitivity only.
2. **Alignment.** The per-marker phyling alignments restricted to the chosen taxa (no re-alignment or re-trimming); columns that became all gaps
   are removed; markers with fewer than 8 taxa are dropped. 76 taxa: 79 of 80 markers, 39,744 sites, median gap fraction 31.4 %, highest 91 %
   (*Tieghemiomyces*). 63 taxa (backbone): 80 markers, 51,761 sites, median gap 29.6 %.
3. **Concatenated ML trees.** IQ-TREE 3.0.1. Guide model **LG+F+G** with 1,000 UFBoot and 1,000 SH-aLRT replicates (seed 20261112). Site-heterogeneous
   model **LG+C60+F+G with PMSF on the guide tree** (C60 chosen because memory allowed: about 92 GB for 76 taxa and 110 GB for 63 taxa, not the 250-400 GB
   first assumed), 100 nonparametric bootstrap replicates — **unfinished, see below**.
4. **Gene-tree summary.** One IQ-TREE per marker (ModelFinder among LG, WAG, JTT, Q.PFAM; 1,000 UFBoot), ASTRAL 5.15.5 on all gene trees
   (79 trees for the 76-taxon set, 80 for the backbone; branch support = local posterior probability). Sensitivity: the same with gene-tree branches
   of UFBoot < 50 collapsed (threshold chosen in advance). ASTRAL's full-annotation output is kept as `astral_t2.tre`.
5. **Rooting.** The 76-taxon set has no taxon outside Kickxellomycota, so it is rooted on Dimargaritales (*Dimargaris*, *Dispira*, *Tieghemiomyces*), the
   earliest Kickxellomycota branch in the backbone, which is itself rooted on the 15 chytrid-grade taxa (chytrids, Blastocladiomycota,
   Monoblepharidomycota, Neocallimastigomycota). Sensitivity: alternative outgroups and midpoint rooting (below). Rerooting is done by our own function
   (`kickx_tree_tools.root_on`, tested: same splits, same total branch length, supports move with their branches); Bio.Phylo's own rerooting
   mislabelled supports on these trees and was **not** used for any number below (an intermediate version of this analysis used it; those numbers were discarded).
6. **Comparison.** Unrooted Robinson-Foulds distance on shared tips (normalized by the number of splits in both trees) and clade tests on the rooted trees.

## Dereplication result (observed)

- 1,403 strain pairs received a fastANI value (the others are below the tool's reporting range); highest ANI 99.98 %. Pairs with ANI >= 99 %: 116; >= 97 %: 395; >= 95 %: 407; >= 90 %: 488.
- **At 99 %: 62 clusters from 107 strains** (44 singletons; largest clusters of 12, 8, 5, 4, 4, 4, 3, 3 strains). At 97 %: 38 clusters (largest 24); at 95 %: 36. The 97 % and
  95 % results are close, so the strains fall into about 36-38 groups at a species-like ANI; this is a description, not a species claim (most strains are unnamed "sp.").
- Retained: 62 *Coemansia* representatives + 14 others = **76 taxa** (`derep/derep_taxa.tsv`; membership of every strain in `derep/derep_clusters_99.tsv`).

## Tree results (observed; PMSF trees pending)

Trees compared: LG+F+G (UFBoot), ASTRAL (all gene-tree branches), ASTRAL (UFBoot < 50 collapsed), and the first-pass FastTree. Full tables: `*/..._clades.tsv`, `*_report.txt`.
Support is UFBoot for LG+F+G, local posterior probability for ASTRAL, FastTree local support for the first pass.

| Question | Backbone (63 taxa) | *Coemansia* set (76 taxa) |
|---|---|---|
| Kickxellomycota monophyletic | yes in all four trees (UFBoot 100, ASTRAL 1.0, FastTree 1.0); sister = Entomophthorales (*Conidiobolus*, *Capillidium*, *Neoconidiobolus*) in all four | not tested (no outgroup) |
| Order of the early branches | Dimargaritales, then *Ramicandelaber*, then (Kickxellales + Harpellales), in all four trees | same, in all four trees |
| Kickxellales + Harpellales monophyletic (each order also monophyletic) | yes; UFBoot 100, ASTRAL 1.0, FastTree 1.0 | yes; UFBoot 100, ASTRAL 1.0 (also with collapsed gene trees), first-pass FastTree 1.0 |
| ***Coemansia* monophyletic** | trivial (one genome) | **no, in all four trees: *Kickxella* nests inside it** |
| Position of *Kickxella* | sister to the single *Coemansia* genome | sister to one strain (*Coemansia* sp. Benny D115) inside a clade of strains (UFBoot 67; ASTRAL 0.49) |
| ***Linderina*** | LG+F+G and FastTree: sister to (*Coemansia*, *Kickxella*, *Dipsacomyces*, *Martensiomyces*) (UFBoot 100); **both ASTRAL: sister to (*Dipsacomyces*, *Martensiomyces*)** (PP 0.85, and 0.94 with collapsed gene trees) | LG+F+G (UFBoot 67) and both ASTRAL (PP 0.99): sister to (*Dipsacomyces*, *Martensiomyces*); first-pass FastTree: sister to the other four |
| ***Basidiobolus*** | LG+F+G (UFBoot 98): with Zoopagomycota + Entomophthoromycota + Kickxellomycota; first-pass FastTree (0.43) same; **both ASTRAL (PP 0.86): sister to the Mucoromycota-group clade** | n/a |
| Normalized RF, LG+F+G vs ASTRAL | 0.12 (0.10 collapsed) | 0.23 (0.22 collapsed) |
| Normalized RF, LG+F+G vs first-pass FastTree | 0.08 | 0.15 |
| Internal branches with support >= 95 / >= 70 / < 50 (LG+F+G) | 51 / 56 / 0 of 60 | 58 / 68 / 1 of 73 |

1. **The first-pass findings hold with a second model and a second method:** the Kickxellomycota backbone order, Kickxellales + Harpellales as sisters with *Ramicandelaber* as their
   closest relative and Dimargaritales as the earliest branch, and a non-monophyletic *Coemansia* with *Kickxella* inside.
2. ***Kickxella* inside *Coemansia* is reproduced but its exact position is weak** (UFBoot 67, ASTRAL 0.49). There is only one *Kickxella* genome in BFD, so real paraphyly of
   *Coemansia*, a misidentified or contaminated genome, and a long-branch or missing-data artifact cannot be separated here (the first-pass note on the *Kickxella* genome's 40 % missing sites still applies).
3. ***Linderina* and *Basidiobolus* remain unstable**, and the instability follows the method (concatenation vs ASTRAL), not the model: do not use either placement for a
   reconstruction that depends on it. *Linderina* is a single, flagged genome (24 % duplicated, high gap fraction in the first pass).
4. The dereplicated 76-taxon set resolves the strain-level tree well (58 of 73 internal branches at UFBoot >= 95), with the closest retained pairs still very close
   (0.0014 and 0.0045 substitutions per site), so near-clonal pairs remain below the 99 % ANI cut in places; this does not affect the deeper results.

## Rooting sensitivity (observed)

- **Backbone:** the 15 chytrid-grade taxa form one branch of the unrooted LG+F+G tree; rooting on them, or on them without Blastocladiomycota, gives the same Kickxellomycota result.
  Midpoint rooting puts the root on the Kickxellomycota stem (15 versus 48 tips) because Kickxellomycota branches are long, so **midpoint rooting is not credible here**; it still keeps
  Kickxellales + Harpellales together with *Ramicandelaber* as sister.
- ***Coemansia* set:** the three Dimargaritales genera form one branch of the unrooted tree. Rooting on Dimargaritales + *Ramicandelaber* gives the same clades (the sister of
  Kickxellales + Harpellales is then Dimargaritales + *Ramicandelaber*). Rooting on Harpellales or on *Mycoemilia* + *Spiromyces* breaks Kickxellales + Harpellales (and, for *Mycoemilia* +
  *Spiromyces*, Kickxellales), which contradicts the backbone, which has the proper outgroups. Midpoint rooting falls between (Dimargaritales + *Ramicandelaber*) and the rest, consistent with the
  Dimargaritales rooting. **"*Coemansia* is not monophyletic" holds under every rooting tested.**

## Placement of the Ku loss on the new trees

Ku is absent from all Kickxellales and Harpellales genomes and present in Dimargaritales and *Ramicandelaber* (Ku pathway README; tblastn on 128 genomes). On every tree here
(LG+F+G, both ASTRAL variants, first-pass FastTree; both tree sets) the Ku-negative genera form **one clade, Kickxellales + Harpellales (UFBoot 100, ASTRAL 1.0)**, whose sister is *Ramicandelaber*
(Ku-positive) and then Dimargaritales (Ku-positive). A single loss on the **stem branch of Kickxellales + Harpellales** (between the *Ramicandelaber* split and the crown of the two orders)
therefore explains the pattern parsimoniously on all four topologies, and the placement does not depend on how *Coemansia*, *Linderina* or *Basidiobolus* are resolved, nor on the rooting tested.
This is a placement of one loss on a species tree, not a reconstruction: the gene-presence matrix was not re-run on these trees (Dollo or likelihood on the final tree is still to do), and further
independent losses elsewhere (for example the weak or absent tblastn calls in *Piptocephalis*, *Neoconidiobolus* and *Capillidium*) are outside this test.

## What is still running (finalize recipe)

The site-heterogeneous tree is the one planned step not delivered. With `-b`, IQ-TREE 3.0.1 overwrites `pmsf.treefile` with each bootstrap replicate's tree and keeps no copy of the ML tree
(found the hard way: an intermediate version of this analysis used that file and was discarded), so the ML tree needs its own run. Submitted on the HPCC (all in `bench/kickx_tree/iq/<name>/`):

| Job | What | State when this was written |
|---|---|---|
| 32083430 | backbone PMSF ML tree, no bootstrap (`kickx_iqtree.sh pmsfml`) | running, mixture-model fit iteration 2 |
| 32083431 | 76-taxon PMSF ML tree, no bootstrap | queued |
| 32079503 | backbone PMSF bootstrap (writes `pmsf.boottrees`) | running; 6 of 100 replicates after 3.5 h (about 35 min each) |
| 32080666 | 76-taxon PMSF bootstrap | running; mixture fit still iterating after 3 h; bootstrap not started |

When the ML jobs finish, per dataset: `iqtree3 -sup pmsf_ml.treefile -t pmsf.boottrees --prefix pmsf_ml_sup` (only when enough replicates exist; with a few replicates leave support off), then rerun
`workflows/phylo/kickx_tree_report.py <mode> <iq dir> <taxa.tsv> <first-pass treefile> <outdir>`; it picks up `pmsf_ml.treefile` (or `pmsf_ml_sup.suptree`) automatically as the primary tree, and the tables and
figures here regenerate. The 100-replicate PMSF bootstrap needs about 2-3 days per dataset on the epyc/batch nodes; fewer replicates, or UFBoot on the guide model (already done), are the cheaper support options.

## Limits

- No site-heterogeneous tree is reported (above). LG+F+G is a single-matrix model and can mislead where compositional heterogeneity or long branches matter (the *Kickxella* and *Linderina* placements are the likely cases).
- ASTRAL with 79-80 gene trees from alignments of 12-71 taxa each (median 61 for the 76-taxon set); gene-tree error is large for short alignments, and ASTRAL support is a local posterior probability.
- One genome each of *Kickxella*, *Linderina* and *Basidiobolus*; flagged by high gap fraction or duplication in the first pass.
- Dereplication uses ANI at one pre-stated threshold; strain pairs below fastANI's reporting range are treated as not redundant.
- The first-pass FastTree comparison uses a 121-taxon tree pruned to the 76 shared tips.
- No literature is cited; the contested order of the zygomycete groups is still described as contested without references (citation to verify).
