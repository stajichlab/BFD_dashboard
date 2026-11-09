# BFD summary atlas

Static summary of the BFD fungal genome database, grouped by hypothesis. Built from the parquet snapshot in `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables` (species table built 2026-11-06) through the intermediates in `BFD_dashboard/atlas_data/`. Design and open decisions: [PLAN.md](PLAN.md).

**Reports:** [H2 report: are yeast proteins shorter than filamentous orthologs?](H2_report.md)

Genomes: 23,680. Species representatives (one per species, best BUSCO then N50): 7,435. Genus representatives: 1,669.

## How to read the statistics

- The unit of replication is stated on every view. Most cross-clade tests use genus medians of representative genomes, because strains of one species and species of one genus are not independent.
- Quality filters: Q-ann (annotated), Q-rep (one genome per species), Q-busco-rel (BUSCO complete >= the 25th percentile of its class), Q-dup (BUSCO duplicated <= 10 %), Q-n50 (N50 >= 50 kb), Q-model (complete gene models: not partial, start and stop codon present).
- Effect sizes come first (median differences with bootstrap 95 % CI, Cliff's delta, log2 ratios). Kruskal-Wallis is a screening statistic only. Mixed models use random intercepts for CLASS and ORDER within CLASS. PGLS is marked "pending tree" until the genus-level phylogeny finishes.
- Growth form and lifestyle labels are genus-level traits from FungalTraits 1.2 plus an override file for genera it lacks; "dimorphic_yeast" is shown as a separate group and never pooled.

## Hypotheses

| Page | Hypothesis | Views |
|---|---|---|
| [H1_quality](H1_quality.md) | H1. Annotation and assembly quality differ by clade | V01 (ok), V02 (ok), V03 (ok), V04 (ok), V05 (ok), V06 (ok) |
| [H2_protein_length](H2_protein_length.md) | H2. Among orthologous proteins, yeast proteins are shorter | V07 (ok), V08 (ok), V09 (ok), V10 (ok), V11 (ok), V12 (ok), V13 (ok), V14 (ok), V31 (ok) |
| [H3_gene_structure](H3_gene_structure.md) | H3. Yeast genomes have fewer introns and denser, shorter genes | V15 (ok), V16 (ok), V17 (ok), V18 (ok) |
| [H4_secretome_lifestyle](H4_secretome_lifestyle.md) | H4. Secretome and CAZy content track lifestyle | V20 (ok), V21 (ok), V22 (ok), V23 (ok), V24 (ok), V25 (ok) |
| [H5_domains_tf](H5_domains_tf.md) | H5. TF and domain repertoires scale sub-linearly; Zn2Cys6 expansion is Pezizomycotina-specific | V27 (ok), V28 (ok), V29 (ok), V30 (ok) |
| [H6_codon_usage](H6_codon_usage.md) | H6. Codon usage tracks taxonomy more than growth form; GC3 follows genomic GC | V33 (ok), V34 (ok), V35 (ok), V36 (ok), V37 (ok), V38 (ok) |
| [H7_paralogs](H7_paralogs.md) | H7. Recent paralog retention is lower in yeasts; large identical-copy groups mark assembly artefacts | V39 (ok), V40 (ok) |
| H8 (P3) | Telomere repeat monomer is conserved within classes | V41, V42 (P3, not built) |

## Data coverage

See V01 on the H1 page. Tables and figures: `tables/`, `figures/`.
