# Kickxellomycota phylogeny (Analysis E follow-up)

Status: **drafted, not yet run.** Nothing here has been submitted on the HPCC.

Decisions (owner, 2026-11-11): narrow trees only, run on the UCR HPCC; marker set **fungi_odb12**; membership as below
(*Coemansia* at BUSCO >= 60% = 107 genomes; all 15 Kickxellomycota genera kept and flagged `low_busco` when BUSCO < 60%;
backbone outgroups unchanged at 65 taxa). The atlas genus-representative tree is being built on a different machine, so
this pair of trees must be reconciled with it (tip labels are LOCUSTAGs here; confirm the same convention there).

## Why

Every cross-clade result in Analysis E and the Kickxellomycota report is tree-free, so genus medians are treated as
independent and the Mucorales baseline cannot be checked. A tree is the gate for ancestral reconstruction of genome
size and intron density, Dollo or Count gain-loss, CAFE, PGLS, and relative-rate tests of selection vs drift.

## Result so far (2026-11-11, first pass, FastTree and ASTRAL on 80 fungi_odb12 markers; not yet a final tree)

Backbone, 63 taxa, 51,761 sites, median 29.6% gaps. All 63 taxa have markers; the weakest are *Tieghemiomyces* (84% gaps,
55% BUSCO-duplicated, the only genome of its genus), *Linderina*, *Spiromyces* and *Basidiobolus* (58-66%). Checks are
`bin/phylo/check_marker_recovery.py`, `check_tree_clades.py` and `print_kickx_subtree.py` (run on the HPCC with
`python3 - <tree> < script.py`).

- Kickxellomycota is monophyletic in both trees (FastTree 1.00, ASTRAL 0.9998), as are all sampled orders.
- Within it, both trees give: Dimargaritales as the earliest branch; Ramicandelaber sister to (Kickxellales + Harpellales);
  Kickxellales and Harpellales sister to each other (support 1.00 and 1.00). In Kickxellales: (Coemansia, Kickxella),
  (Dipsacomyces, Martensiomyces) and (Mycoemilia, Spiromyces).
- **The one disagreement inside Kickxellomycota is *Linderina***: sister to the (Coemansia, Kickxella) + (Dipsacomyces,
  Martensiomyces) clade in FastTree (0.997), but grouped with (Dipsacomyces, Martensiomyces) in ASTRAL (0.87). It is a flagged
  taxon (24% duplicated, 66% gaps).
- *Basidiobolus* is unstable: with Zoopagomycota + Entomophthoromycota in FastTree at support 0.43, but outside that group
  in ASTRAL. Entomophthorales fall with Kickxellomycota + Zoopagales in both trees. The contested zygomycete backbone is
  therefore not resolved by this run.
- What this means for the reduction analyses is a hypothesis, not a result: the unreduced Harpellales are sister to
  Kickxellales, so Ku loss (absent in both) may be a single event on the stem of that pair, and the intron and size reduction
  happened later, inside Kickxellales. Ancestral-state reconstruction is needed to test it.

### *Coemansia* tree (121 taxa: 107 strains + 14 other Kickxellomycota genera; 40,311 sites; FastTree, first pass)

Scripts: `bin/phylo/check_coemansia_tree.py`, `check_coemansia_redundancy.py`. All 121 taxa are in the concatenated alignment
(median gap 28.9%; strains 26.5%, three strains above 60%). Branch lengths are amino-acid substitutions per site and the tree is
unrooted as printed; the redundancy numbers below do not depend on rooting.

- ***Coemansia* is not monophyletic**: the smallest split holding all 107 strains also holds *Kickxella*. *Kickxella* sits in a
  9-taxon clade of 8 strains + *Kickxella* (support 1.00); its closest relative within that clade is one strain at support
  only 0.53. So *Kickxella* nests inside *Coemansia*, but its exact position is weak. The other 13 genera are outside.
  Possible reasons are real paraphyly of the genus, a mislabelled genome, or a long-branch/coverage effect (the *Kickxella*
  genome has 40% missing sites, 11.4% BUSCO-duplicated); this needs a check against other *Kickxella* genomes and a
  second method before it is stated as a result.
- **The strains are strongly redundant.** Of 29 pairs of strains that are each other's closest relatives, 8 differ by less than
  0.001 substitutions per site and 20 by less than 0.005, so many of the 107 "strains" are near-clonal or very close. This
  supports dereplicating before counting independent lineages.
- **Well resolved:** 66 nested all-*Coemansia* clades have support >= 0.95 (largest 25, 21, 19, 19, 15, 13 strains), and only
  7 of 84 internal nodes among strains are below 0.70.

**Update 2026-11-14:** ANI dereplication, LG+F+G (UFBoot) trees, ASTRAL gene-tree summaries and rooting tests are done for both trees; the
first-pass findings hold (Kickxella still nests inside Coemansia, weakly; Linderina and Basidiobolus still unstable, following concatenation vs
ASTRAL). Results, methods and the Ku-loss placement: [`docs/analyses/kx2_ku_pathway/species_tree/`](../../docs/analyses/kx2_ku_pathway/species_tree/README.md).
The site-heterogeneous (LG+C60+F+G, PMSF) trees are still running. Still not done: a marker-count sensitivity run without the flagged taxa.

### Improved-tree scripts (run on the HPCC; outputs in the purgeable `bench/kickx_tree/iq/`)

| Script | Step |
|---|---|
| `kickx_derep.sh`, `kickx_derep.py` | fastANI all-vs-all of the 107 *Coemansia* assemblies; clusters at >= 99 % ANI (97 and 95 % as sensitivity); retained-taxa table |
| `kickx_subset_aln.py` | restrict the nf_phyling per-marker alignments to a taxa list and rebuild the concatenated alignment |
| `kickx_iqtree.sh`, `kickx_iq_submit.sh` | stages `subset`, `guide` (LG+F+G, UFBoot), `pmsf` (LG+C60+F+G bootstrap), `pmsfml` (PMSF ML tree), `gene` (per-marker trees), `astral` |
| `kickx_tree_tools.py`, `kickx_tree_report.py` | own rerooting, Robinson-Foulds, clade tests, rooting sensitivity, figures |
| `test_kickx_phylo.py` | fixture tests (`/usr/bin/python3.12 -m unittest test_kickx_phylo`; 7 tests pass) |

Two pitfalls found: (1) with `-b`, IQ-TREE 3.0.1 overwrites `<prefix>.treefile` with every bootstrap replicate's tree, so run the ML tree separately (`pmsfml`); (2) Bio.Phylo's
`root_with_outgroup` and `root_at_midpoint` mislabel branch supports on IQ-TREE/FastTree trees (trifurcating root), so use `kickx_tree_tools.root_on`.
The LG+F+G guide fit needs about 2 GB and the C60 mixture 92-110 GB; epyc and stajichlab memory was the scheduling bottleneck, so request what is needed (`--mem 160gb`), not hundreds of GB.

## Lessons from the first runs

- **Genome choice.** Ranking by BUSCO complete favoured duplicated assemblies; genomes with BUSCO duplicated 52-70% recovered
  14-426 of 1,089 markers, and *Massospora* (BUSCO 71.7%) had only 415 proteins in BFD and no markers. The selection now ranks by
  single-copy BUSCO and requires an annotation (see `bin/phylo/kickx_tree_taxa.sql`).
- **Nextflow on a busy scheduler.** A partial `squeue` listing made Nextflow fail healthy SLURM jobs ("terminated for an unknown
  reason", 2026-11-09 and 2026-11-11); the script sets `executor.exitReadTimeout = '12 h'`.
- **Memory.** `PHYLING_ALIGN` needed more than 64 GB for 121 taxa (peak 67 GB); memory now scales with the tree and doubles on retry.
- **Paths.** `~/projects` is a symlink into `/bigdata`; containers need the resolved path (`readlink -f`).

## State of the existing genus tree (checked 2026-11-11, read-only)

The atlas genus tree (1,654 genera, `atlas_data/phylo/`, nf_phyling, `fungi_odb10`, FastTree) has **not finished**.
The Nov 9 run ended with 0 succeeded and 2 failed tasks: `PHYLING_ALIGN` once with exit 137 after 4 h 25 min
(killed; the cause was not established, so memory vs time limit is open), and once at 7 min as "terminated by an
external system" (consistent with the maintenance drain). Resubmitting `run_phyling_atlas.sh` unchanged may hit exit 137
again. Triage that job separately (look at `sacct` for the first task's MaxRSS and state) before re-submitting.

## This draft: two small trees instead

| Tree | Taxa | Question |
|---|---|---|
| `backbone` | 65: 15 Kickxellomycota genera, 15 sister (Zoopagomycota, Entomophthoromycota, Basidiobolomycota), 20 Mucoromycota-group, 15 root | Where on the branches did size, introns and gene content change; how good is the Mucorales baseline |
| `coemansia` | 107 *Coemansia* genomes (BUSCO >= 60%; 32 with N50 >= 50 kb) plus the other Kickxellomycota genera | How many lineages, and relationships among strains (129 of 146 representatives are unnamed "sp.") |

Taxon selection is `bin/phylo/kickx_tree_taxa.sql` (rules in its header). Counts above were checked in MotherDuck.
One ingroup genus has BUSCO of only 44%, so marker coverage per taxon must be inspected before trusting its position.

## Run

    git clone <repo> ~/projects/BFD_dashboard/repo   # or any checkout that has this branch
    REPO=~/projects/BFD_dashboard/repo sbatch workflows/phylo/kickx_tree.sh backbone
    REPO=~/projects/BFD_dashboard/repo sbatch workflows/phylo/kickx_tree.sh coemansia

Outputs go to `~/projects/BFD_dashboard/bench/kickx_tree/` (purgeable). `MARKERSET` defaults to `fungi_odb12`; the
markers are downloaded by phyling into `bench/kickx_tree/phyling_db` (the shared atlas `phyling_db` holds only
`fungi_odb10` and is not written to). The pinned nf_phyling revision (`52b694f356`, the head as of 2026-10-03) defaults to
`fungi_odb12`, and its config says it runs a no-support tree and a bootstrap-support tree per marker set. A comma-separated
list such as `fungi_odb12,mucoromycota_odb12` is supported; the first pass uses `fungi_odb12` only because the root taxa
(chytrids, Blastocladiomycota) are outside Mucoromycota. Taxa with `low_busco = true` in `backbone_taxa.tsv` need their
marker count inspected, and the topology reported with and without them.

## Still to do before this answers the questions

1. **Second method.** Done for LG+F+G and ASTRAL (see the update above); the site-heterogeneous PMSF trees (LG+C60+F+G) are
   still running (`docs/analyses/kx2_ku_pathway/species_tree/README.md`, "What is still running"). The order of the zygomycete groups is
   contested (citation to verify) and changes loss counts: *Basidiobolus* is placed differently by concatenation and by ASTRAL.
2. **Marker coverage filter.** Drop taxa or markers below a coverage cutoff and report topology stability.
3. **Dereplicate Coemansia**: done with fastANI (62 representatives at >= 99 % ANI; `derep/` in the results folder). Mash sketches from the
   PopSeq project were not checked for reuse.
4. **Fixture test** for the selection SQL (a small parquet fixture) and a check that every selected LOCUSTAG has a
   proteome in `locustag_to_pep.tsv` (the script already fails loudly when one is missing).
