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
A **clade** is a family with at least 3 genera in each of the two groups. A contrast is evaluable if it has at least 4 clades. For each response: per-clade effect = median(y, group A) - median(y, group B); pooled effect = mean of clade effects weighted by min(nA, nB); p from 19,999 permutations of group labels within family; q = BH over responses within the contrast.
**Amendment before any run (2026-11-21):** the statistic is the weighted mean over clades of the difference of group means (of log-transformed values), not of medians. The mean difference is linear in the labels, so the within-family permutation of all features is one matrix product; a median over 4,000 features and 20,000 permutations is not feasible. Effects in the per-clade table are reported as mean differences.
**Rule:** supported if q < 0.05 **and** the clade effect has the same sign as the pooled effect in at least 75% of clades; mixed if q < 0.05 but fewer than 75% agree; no evidence otherwise. A per-clade table is written so that it can be seen which families carry the effect.

## Test C: domain scan (exploratory, labelled as such)

Same statistic and within-family permutation as Test B, for every Pfam domain and CAZy family present in at least 10% of the genera of the contrast (log1p of rate per 1,000 proteins), BH over all tested features of a contrast. Report the number with q < 0.05 and the top 15 with names. Nothing here is a mechanism; the list generates hypotheses.

## Limits stated in advance

Genus-level ecology label (one label per genus; real species vary); shared ancestry beyond family; annotation quality differs between clades; protein-content tables exist only for function-annotated genomes; lifestyle groups differ in how many families they share, so some contrasts are evaluable in few clades; many tests (BH is applied within each family of tests, not across them).

---
## Results (2026-11-21; jobs 32361877 and 32361904; 1,489 genus representatives, 190 orders, 533 families)

Order of events: design committed first (`95bd8a5`); one statistic amendment before any run (mean difference instead of median, see Test B). The first run of Test B and C at family level showed that most contrasts could not be evaluated, so an **order-level tier was added after seeing that result (post hoc)**. The first PERMANOVA run failed (all-zero profiles break Bray-Curtis); the fix was to drop genera with an empty profile, and it was rerun. Tables are in `results/`.

### Test A: how much do taxonomy and ecology explain?

Univariate responses (18: nine genome properties, nine protein-content measures; `results/A_variance_partition.tsv`):

| Model | median R2_adj taxonomy | median R2_adj ecology alone | median unique ecology (after taxonomy) | responses where ecology explains (unique >= 0.02 and q < 0.05) |
|---|---|---|---|---|
| order + coarse lifestyle (6 groups) | 0.60 (0.36 to 0.93) | 0.17 | 0.021 (max 0.050) | 9 of 18 |
| order + lifestyle (16 to 17 groups) | 0.61 | 0.27 | 0.022 | 10 of 18 |
| family + coarse lifestyle | 0.76 | 0.16 | 0.006 (max 0.026) | 1 of 18 (median protein length) |

- Taxonomy alone explains 36% to 93% of the variance. Ecology alone explains about 17% (coarse), and about 88% of that (median share) is shared with taxonomy: lifestyle follows ancestry.
- After order, ecology adds a median 2.1 percentage points. Responses that pass the pre-set threshold with order as the stratum: genes per Mb, coding fraction, median protein length, CAZy GH, PL, CE, AA, strict secretome and MEROPS.
- After the finer stratum (family) the unique part of ecology falls to a median 0.6 points, and only median protein length passes. So most of the ecology signal at order level is differences between families inside an order.
- Profiles (PERMANOVA of lifestyle, `results/A_permanova_profiles.tsv`; permutations within order, 999, so p = 0.001 is the floor):

| Profile | coarse lifestyle, R2_adj unrestricted / within order | fine lifestyle, R2_adj unrestricted / within order |
|---|---|---|
| Pfam domains (1,515) | 0.118 / 0.107 | 0.235 / 0.239 |
| CAZy families (179) | 0.145 / 0.135 | 0.273 / 0.276 |
| TF families (26) | 0.157 / 0.162 | 0.231 / 0.235 |

  All p = 0.001. The pseudo-F and R2 are computed on all genera in both tests; the within-order test changes only the permutations. Unlike the unique-part analysis above, this test does not remove the order effect from the statistic, so R2 stays similar.

### Test B: matched contrasts among relatives

Families that contain at least 3 genera of each of two ecologies are few (`results/B_matched_contrasts_family.tsv`):

| Contrast | families (>= 3 genera each side) | genera (A / B) | verdict |
|---|---|---|---|
| plant_pathogen vs saprotroph | 6 | 33 / 37 | no evidence in all 18 responses |
| animal_parasite vs saprotroph | 1 | 3 / 3 | not evaluable |
| mycoparasite vs saprotroph | 1 | 5 / 10 | not evaluable |
| mycorrhizal vs saprotroph | 0 | | not evaluable |
| lichenized vs saprotroph | 0 | | not evaluable |
| yeast vs filamentous | 0 | | not evaluable |

- Closely related genera that differ in ecology are rare: within a family the lifestyle is mostly one value. This is itself a result, and it limits what can be tested at the family level.
- Plant pathogen vs saprotroph at family level: no response reaches q < 0.05. The closest are gene count (pathogens +0.064 log10, same sign in 6 of 6 families, q 0.053) and genome length (+0.059, 6 of 6, q 0.106). Not a finding by the rule.

**Post hoc order-level tier** (`results/B_matched_contrasts_order.tsv`; clades are orders, so relatives are further apart and families inside an order are not controlled):

| Contrast | orders | genera (A / B) | supported | mixed | no evidence |
|---|---|---|---|---|---|
| mycorrhizal vs saprotroph | 5 (Agaricales, Boletales, Cantharellales, Pezizales, Russulales) | 89 / 125 | 16 | 0 | 2 |
| animal_parasite vs saprotroph | 4 | 42 / 89 | 2 (CAZy PL, MEROPS) | 5 | 11 |
| plant_pathogen vs saprotroph | 11 | 155 / 304 | 0 | 7 | 11 |
| mycoparasite vs saprotroph | 2 | 17 / 31 | not evaluable | | |

- Mycorrhizal genera have fewer CAZy genes than saprotrophs in the same order, in 5 of 5 orders for most classes: GH -6.8, AA -3.4, CE -1.2, PL -0.8 per 1,000 genes; strict secretome -22.7; MEROPS -4.0; TFs -3.3; also lower genes per Mb and median exons per gene (q <= 0.04 for all 16). Effect sizes per order are in `results/B_per_clade_effects_order.tsv` (GH: -5.5 to -8.6 in each of the 5 orders).
- Plant pathogens vs saprotrophs: fewer GH (-1.0) and AA (-0.8), more GT (+0.31) per 1,000 genes, with only 55% to 64% of orders in the same direction (mixed by the rule).
- **Caveat that matters:** within an order, ectomycorrhizal genera often sit in different families from the saprotrophs (Boletales: 60 mycorrhizal against 4 saprotrophic genera). The permutation shuffles labels within order, so it does not separate ecology from family. This tier is anti-conservative; the family-level tier, which can control it, could not evaluate these contrasts.

### Test C: domain scan (exploratory; 19,999 permutations; BH over about 1,700 features per contrast)

| Level | Contrast | features with q < 0.05 | of these, same sign in >= 75% of clades |
|---|---|---|---|
| family | plant_pathogen vs saprotroph | 0 of 1,701 | |
| order | mycorrhizal vs saprotroph | 388 of 1,673 | 379 |
| order | animal_parasite vs saprotroph | 145 of 1,704 | 66 |
| order | plant_pathogen vs saprotroph | 634 of 1,707 | 34 |

- The mycorrhizal list is led by loss of plant-cell-wall CAZy families (GH3, GH5, GH16, GH27, GH28, GH31, GH35, GH43, GH51, GH53, GH55), all lower and same sign in 100% (GH16: 80%) of the 5 orders. That is a pattern, not a mechanism.
- Plant-pathogen hits are mostly weak (small effects, 55% to 64% of orders agree); only 34 of 634 agree in at least 75% of orders. Treat that list as noise-prone.
- Pfam names are in `results/C_domain_scan_*.tsv` (CAZy and TF features have no name column).

### What this says

1. Taxonomy explains most of the variance in genome properties and protein content. Ecology adds a small unique part (about 2 points after order, under 1 after family).
2. The one clear within-order signal is ectomycorrhizal vs saprotrophic genera, with large and consistent loss of CAZy, secretome and protease content. It is not confirmed at family level.
3. Plant pathogens, animal parasites and the other contrasts show weak or inconsistent differences among relatives.

### Limits

One genus label per genus; shared ancestry above family is not modeled (no tree); permutations are exchangeable only within the stated clade; order-level tier is post hoc and does not control family; function-annotated genomes only for protein content; the domain scan has no independent replication; no mechanism or function is tested.
