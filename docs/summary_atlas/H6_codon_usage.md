# H6. Codon usage tracks taxonomy more than growth form; GC3 follows genomic GC

Expected: CLASS explains more RSCU variance than growth form.

[Back to index](index.md)

### V33 Codon usage ordination

![](figures/V33_rscu_pca.png)

- Unit = genus (median RSCU over representative genomes; 1,355 genera). RSCU computed with the standard code for all genomes (CTG-clade code-12 genomes keep CTG in the Leu family; see caveat). Filters: Q-ann, Q-rep (one genome per species), Q-busco-rel (BUSCO complete >= class 25th percentile).
- PC1 explains 74.4 %, PC2 5.9 % of RSCU variance.
- PERMANOVA (Euclidean, 1,121 genera in classes with >= 5 genera, 499 permutations):

| model | pseudo_F | R2 | p | n | n_groups | n_perm |
|---|---|---|---|---|---|---|
| CLASS | 55.194 | 0.567 | 0.002 | 1,121 | 27 | 499 |
| growth_form | 91.857 | 0.141 | 0.002 | 1,121 | 3 | 499 |
| growth_form | CLASS (permutations within class) | 91.857 | 0.141 | 0.856 | 1,121 | 3 | 499 |

- Caveat: whole-proteome codon frequencies, not expression-weighted; GC content drives much of PC1 (see V34). Tables: `tables/V33_permanova.tsv`, `tables/V33_genus_pcs.tsv`.


### V34 GC3 versus genomic GC

![](figures/V34_gc3_vs_gc.png)

- n = 5,120 genomes / 5,120 species / 1,164 genera. Filters: Q-ann, Q-rep (one genome per species), Q-busco-rel (BUSCO complete >= class 25th percentile).

| growth_form | n_genomes | slope_gc3_on_gc | intercept | r2 | spearman_rho |
|---|---|---|---|---|---|
| yeast | 1,109 | 1.893 | -35.100 | 0.931 | 0.961 |
| dimorphic_yeast | 207 | 1.688 | -28.564 | 0.804 | 0.944 |
| filamentous_mycelium | 3,804 | 1.616 | -19.821 | 0.688 | 0.809 |

**Residual GC3 (yeast vs filamentous):**
- Unit = genus (median of representative genomes). Yeast genera n = 153, filamentous genera n = 960.
- Median GC3 residual from the filamentous line (%): yeast -4.576, filamentous -1.374; difference -3.201 (bootstrap 95 % CI -3.889 to -2.049); Cliff's delta -0.420; Mann-Whitney p = 6.85e-17.
- Dimorphic yeast genera n = 51, median -5.124 (shown, not pooled).
- Mixed model `gc3_resid ~ gf_yeast + (1|CLASS/ORDER_)` (n = 1,091 genera, 49 classes, converged = True): yeast effect -3.067 (95 % CI -5.112 to -1.022), p = 0.003. PGLS pending tree.
- Caveat: assembly GC includes repeats and masked regions. Table: `tables/V34_gc3_slopes.tsv`.


### V35 Codon-usage distance versus taxonomic rank

![](figures/V35_distance_by_rank.png)

- Unit = genus pair (594,595 pairs among 1,091 yeast or filamentous genera). No p-value: pairs are not independent.

| rank | same_growth_form | median | size |
|---|---|---|---|
| same family | 1 | 0.926 | 4,090 |
| same order | 0 | 1.328 | 6 |
| same order | 1 | 1.215 | 19,580 |
| same class | 0 | 1.065 | 175 |
| same class | 1 | 1.700 | 66,351 |
| same phylum | 0 | 4.479 | 76,249 |
| same phylum | 1 | 2.188 | 115,122 |
| different phylum | 0 | 3.864 | 65,510 |
| different phylum | 1 | 2.434 | 247,512 |

- Reading: if distance rises with rank at both colors, taxonomy structures codon usage; the gap between colors within a rank is the growth-form signal at that rank. Table: `tables/V35_distance_by_rank.tsv`.


### V36 tRNA gene copy number versus codon usage

![](figures/V36_trna_codon_rho_by_class.png)

- n = 5,023 genomes / 5,023 species / 1,237 genera. Filters: Q-ann, Q-rep (one genome per species), Q-busco-rel (BUSCO complete >= class 25th percentile), Q-n50 (N50 >= 50 kb). Per genome: Spearman correlation across 61 sense codons between codon frequency (codon_freq) and the number of tRNA genes whose `codon` field matches.
- Median tRNA genes per representative genome: 158; median codons with >= 1 tRNA gene: 45 of 61. Most frequent `codon` values in gene_trnas: GGC (231,802), GAC (208,962), AAG (200,057), GTT (188,899), GCT (182,782).
- Unit = genus (median of representative genomes). Yeast genera n = 151, filamentous genera n = 884.
- Median Spearman rho: yeast 0.388, filamentous 0.230; difference 0.158 (bootstrap 95 % CI 0.131 to 0.175); Cliff's delta 0.630; Mann-Whitney p = 3.21e-35.
- Dimorphic yeast genera n = 51, median 0.161 (shown, not pooled).
- Mixed model `rho ~ gf_yeast + log_n50 + (1|CLASS/ORDER_)` (n = 1,014 genera, 47 classes, converged = True): yeast effect 0.061 (95 % CI 0.006 to 0.117), p = 0.031. PGLS pending tree.
- Caveat: the `codon` column of gene_trnas was taken as the decoded codon (the sampled rows are consistent with codons, not anticodons; not verified against tRNAscan output); wobble decoding is ignored, so rho is bounded below 1. Table: `tables/V36_genus_medians.tsv`.


### V37 CTG-clade codon usage

![](figures/V37_ctg_by_code.png)

- Pichiomycetes representative genomes: 491 (285 with TRANSL_TABLE 12, 205 with code 1); unit = genus (32 code-12 genera, 26 code-1 genera). Filters: Q-ann, Q-rep (one genome per species), Q-busco-rel (BUSCO complete >= class 25th percentile).
- CTG per 1,000 codons, genus medians: code 12 5.23, code 1 9.20; difference -3.97 (95 % CI -5.20 to -1.17), Cliff's delta -0.483, p = 0.002. All-fungi genus median: 17.70 per 1,000.
- Caveat: the genetic code comes from the samples sheet, not from the sequence; the single code-26 genome (*Pachysolen*) is not shown. Table: `tables/V37_genus_medians.tsv`.


### V38 Amino-acid composition ordination

![](figures/V38_aa_pca.png)

- Unit = genus (1,355 genera). Filters: Q-ann, Q-rep (one genome per species), Q-busco-rel (BUSCO complete >= class 25th percentile). `*` and `X` symbols excluded; frequencies renormalized to the 20 amino acids.
- PC1 75.0 %, PC2 7.2 %. Spearman rho(PC1, genomic GC) = -0.857. Strongest PC1 loadings: A -0.56, R -0.28, P -0.28 ... I +0.33, K +0.37, N +0.40.
**PC1 score, yeast vs filamentous (GC as covariate in the mixed model):**
- Unit = genus (median of representative genomes). Yeast genera n = 152, filamentous genera n = 954.
- Median PC1 score: yeast 0.035, filamentous -0.006; difference 0.041 (bootstrap 95 % CI 0.032 to 0.048); Cliff's delta 0.601; Mann-Whitney p = 9.73e-33.
- Dimorphic yeast genera n = 51, median -0.012 (shown, not pooled).
- Mixed model `pc1_score ~ gf_yeast + genomic_gc + (1|CLASS/ORDER_)` (n = 1,091 genera, 49 classes, converged = True): yeast effect 0.006 (95 % CI 8.96e-04 to 0.012), p = 0.022. PGLS pending tree.
- Tables: `tables/V38_genus_pcs.tsv`, `tables/V38_loadings.tsv`. Points are rasterized in the SVG.

## Verdict

**Verdict: supported.** CLASS explains 57 % of the variance in genus-median RSCU (PERMANOVA R2 = 0.567) against 14 % for growth form, and growth form has no effect once permutations are restricted within class (p = 0.86, V33). GC3 follows genomic GC (Spearman rho 0.81 filamentous, 0.96 yeast, V34); at equal genomic GC, yeast genera have GC3 3.2 points lower than the filamentous regression line. Codon-usage distance rises with taxonomic rank; within a phylum, pairs with different growth forms are about twice as distant as pairs with the same growth form (V35), which reflects the Saccharomycotina/Pezizomycotina split. Amino-acid composition PC1 (75 % of variance) tracks genomic GC (rho = -0.86, V38). tRNA gene counts correlate with codon frequency more strongly in yeasts (median rho 0.39 vs 0.23, V36). Code-12 Pichiomycetes use CTG at 5.2 per 1,000 codons against 9.2 in code-1 Pichiomycetes and 17.7 across fungi (V37). Caveat: whole-proteome codon frequencies, not expression-weighted; code-12 CTG kept in the Leu family for RSCU.
