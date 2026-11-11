# Kickxellomycota phylogeny (Analysis E follow-up)

Status: **drafted, not yet run.** Nothing here has been submitted on the HPCC.

## Why

Every cross-clade result in Analysis E and the Kickxellomycota report is tree-free, so genus medians are treated as
independent and the Mucorales baseline cannot be checked. A tree is the gate for ancestral reconstruction of genome
size and intron density, Dollo or Count gain-loss, CAFE, PGLS, and relative-rate tests of selection vs drift.

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

Outputs go to `~/projects/BFD_dashboard/bench/kickx_tree/` (purgeable). `MARKERSET` defaults to `fungi_odb10`, the set the
shared `phyling_db` already contains; BFD's BUSCO is `fungi_odb12`, so decide whether that mismatch matters.

## Still to do before this answers the questions

1. **Second method.** FastTree is a first pass. Run IQ-TREE on the phyling concatenated alignment with a
   site-heterogeneous model (for example LG+C60+F+G with ultrafast bootstrap) and a gene-tree summary method
   (ASTRAL on per-marker trees). Exact IQ-TREE module and the phyling alignment file names must be confirmed on the
   HPCC, so no script is drafted for this step yet. Compare backbones: the order of the zygomycete groups is
   contested (citation to verify) and changes loss counts.
2. **Marker coverage filter.** Drop taxa or markers below a coverage cutoff and report topology stability.
3. **Dereplicate Coemansia** by ANI or Mash distance so near-clonal strains count once (the sketches may already
   exist in the PopSeq project; check before computing new ones).
4. **Fixture test** for the selection SQL (a small parquet fixture) and a check that every selected LOCUSTAG has a
   proteome in `locustag_to_pep.tsv` (the script already fails loudly when one is missing).
