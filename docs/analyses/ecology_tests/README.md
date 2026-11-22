# Does ecology explain genome properties and protein content, once taxonomy is taken into account?

Status: **design, written before any test was run (2026-11-21).** Results are added below the line "Results".
This replaces the many ordination figures of the Ecology page with a small set of tests. It does not change `bin/summary/summary_ecology.py`.

## Questions

1. **Variance partitioning.** How much of the variation in a genome property or protein-content measure is explained by taxonomy, by ecology, and by both together?
2. **Matched contrasts.** Among closely related genera (same family) that differ in ecology, are genome properties or protein content different, in the same direction across families?
3. **Where is the signal?** Which protein domains and CAZy families differ in those matched contrasts?

## Unit and labels

- Unit: one representative genome per genus (`v_genome_labels.is_genus_rep`, annotated), from the labels file `db/BFD_labels.duckdb`. Genera without a FungalTraits 1.2 `primary_lifestyle` are excluded, never guessed. The suppress list is applied.
- Ecology factor: `primary_lifestyle` (genus-level, FungalTraits 1.2). Lifestyles with fewer than 8 genera are dropped from the fine model. A coarse factor groups them as: saprotroph (all `*_saprotroph`), plant_pathogen, animal_parasite, mycoparasite, mycorrhizal (ectomycorrhizal, arbuscular_mycorrhizal), lichenized, other (dropped).
- Taxonomy factor: ORDER (and, as a finer check, FAMILY). Genera with a missing taxon level are dropped from that model.
- Genera are not independent (shared ancestry beyond order or family). No phylogenetic correction is possible here (no tree). Every p-value below is a test of association within the stated stratum, not of cause.

## Responses (fixed in advance)

- **Genome properties** (log10 except GC and fractions): total length, n genes, genes per Mb, coding fraction, GC %, median intron length, median exons per gene, fraction of genes with an intron, median protein length.
- **Protein content**, per 1,000 proteins of the genome, for genomes with function annotation: CAZy classes GH, GT, PL, CE, AA, CBM; strict secretome; MEROPS proteases; transcription factors. (Count tables: `summary_tables/*_genome.parquet`.)
- **Multivariate profiles** (Bray-Curtis on rates): Pfam domain profile (domains present in at least 5% of genera), CAZy-family profile, TF-family profile.

## Test A: variance partitioning

For each univariate response y, with order as the stratum: R2_adj of taxonomy alone, of ecology alone, of taxonomy + ecology, and the unique part of ecology (ecology after taxonomy, partial R2_adj) and of taxonomy (taxonomy after ecology). Shared part = R2_adj(ecology alone) minus unique(ecology). Significance of unique ecology: F statistic of lifestyle after order, with p from 999 permutations of lifestyle labels **within order** (Frisch-Waugh residualisation). BH over responses. For the multivariate profiles: PERMANOVA of lifestyle with permutations within order (the same stratified test as the existing page), plus the unrestricted PERMANOVA for comparison.
**Reading:** ecology "explains" a response if unique partial R2_adj >= 0.02 and q < 0.05. Taxonomy explaining most and ecology a small unique part is an expected and reportable result.

## Test B: matched contrasts within families

Focal contrasts (fixed): plant_pathogen vs saprotroph; animal_parasite vs saprotroph; mycorrhizal vs saprotroph; lichenized vs saprotroph; mycoparasite vs saprotroph; and, by growth form, yeast vs filamentous_mycelium.
A **clade** is a family with at least 3 genera in each of the two groups. A contrast is evaluable if it has at least 4 clades. For each response: per-clade effect = median(y, group A) - median(y, group B); pooled effect = mean of clade effects weighted by min(nA, nB); p from 999 permutations of group labels within family; q = BH over responses within the contrast.
**Rule:** supported if q < 0.05 **and** the clade effect has the same sign as the pooled effect in at least 75% of clades; mixed if q < 0.05 but fewer than 75% agree; no evidence otherwise. A per-clade table is written so that it can be seen which families carry the effect.

## Test C: domain scan (exploratory, labelled as such)

Same statistic and within-family permutation as Test B, for every Pfam domain and CAZy family present in at least 10% of the genera of the contrast (log1p of rate per 1,000 proteins), BH over all tested features of a contrast. Report the number with q < 0.05 and the top 15 with names. Nothing here is a mechanism; the list generates hypotheses.

## Limits stated in advance

Genus-level ecology label (one label per genus; real species vary); shared ancestry beyond family; annotation quality differs between clades; protein-content tables exist only for function-annotated genomes; lifestyle groups differ in how many families they share, so some contrasts are evaluable in few clades; many tests (BH is applied within each family of tests, not across them).
