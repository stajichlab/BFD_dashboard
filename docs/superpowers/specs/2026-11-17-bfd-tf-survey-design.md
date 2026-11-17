# Transcription-factor survey across BFD: design

Date: 2026-11-17. Status: DRAFT, awaiting approval. Repo: `stajichlab/BFD_dashboard`.
Source paper: Lax et al. 2025, Nature Communications 16:10277 (`resources/s41467-025-65177-1.pdf`). DAP-seq of 58 TFs in *Rhizopus microsporus*, linked to 6mA DNA methylation.

## Intent

The paper surveys TF families in 62 fungal genomes (43 early-diverging fungi, EDF) and maps the binding sites of 58 TFs in one species. BFD has more than 23,000 genomes. This work enumerates TFs across all of them, tests the paper's lineage contrasts at that scale, follows three TFs the paper highlights across fungi, and searches for expansions and contractions the paper did not cover. Results go on the Quarto site.

User scope: A (landscape and the paper's contrasts), B (ortholog tracking), C (broader expansion and contraction search). All three.

## What the paper claims that we can test (from its Results, "The fungal TF landscape")

- Zn cluster TFs are expanded in Dikarya, especially Ascomycota, and abundant in Entomophthoromycota.
- HTH/Homeodomain-like TFs are more expanded in EDF than in Ascomycota.
- C2H2/CCHC/CH/C5HC2 zinc fingers are expanded in Monoblepharomycota, Neocallimastigomycota and Chytridiomycota.
- HMG-box TFs are expanded in Glomeromycota.
- HLH and GATA TFs are expanded in Mucoromycota, Mortierellomycota and some Chytridiomycota.
- E2F and Tubby are absent in Dikarya. RFX and STAT are present in Dikarya but absent or rare in EDF.
- WRKY and E2F are almost exclusive to Mucoromycota. Mucoromycota has the highest TF diversity (Shannon index).
- Three TFs with functional data in *R. microsporus*: WC-2C (GATA, white-collar, motif GATC, 1,394 targets, 427 of 2,871 light-regulated genes), Znf1 (C2H2, binds upstream of the zinc transporter gene `zrt1`, Zap1-like), bZip5 (binds the `sexP` promoter).

## What BFD has (checked 2026-11-17)

- `tf_inventory`: 10.9 million protein-domain rows, 23,394 genomes, 21 Pfam-based families (22 Pfam domains, list in `nextflow/assets/fungal_tf_pfam_domains.csv`). This is coarser than the paper's 51 families. It has no E2F, Tubby, RFX, STAT, WRKY, CP2, Copper-fist, TEA or CBF/NF-Y.
- The full `pfam` table holds those domains for the same genomes. Genome counts: RFX_DNA_binding 21,484, CP2 13,649, TEA 21,215, Copper-fist 21,331, CBF/NF-Y 23,334, KilA-N 23,033, Tubby 839, E2F_TDP 604, WRKY 452. No new annotation run is needed.
- `seqdb/proteins` (bgzip + faidx) for all proteomes.
- The paper's anchor species are in BFD: *R. microsporus* (ATCC 11559 is `FBF4BF29`, plus several other strains), *Mucor lusitanicus*, *Phycomyces blakesleeanus*.
- Not available: the paper's supplementary tables (the 88 tested TFs with JGI protein IDs, the 62-species list, the WC-2C target list). The PDF has the main text only. Binding sites and motifs are experimental and cannot be extended from BFD data.

## Goals

- A: one TF family table for all genomes with an extended family set; landscape figures by phylum, class and order; a table that tests each claim above at BFD scale.
- B: for WC-2 (white-collar), Zap1-like C2H2 and bZip5, find orthologs across fungi; report presence, copy number and domain architecture by lineage; give a tree.
- C: a ranked list of lineage-by-family expansions and contractions, with effect sizes and corrected p-values, labelled as hypotheses.
- Outputs follow the site rules: generated files go into the release asset, not git; pages read tables at render time; the suppress list applies.

## Non-goals

- No regulation or binding claims for any genome except as quoted from the paper.
- No new annotation or domain-search run on all proteins. B searches only a few seed families.
- No full-genome orthology inference (OrthoFinder on 23,000 proteomes).

## Anti-goals

- A "supported / not supported" verdict from a threshold nobody agreed on. We report the numbers and the direction of the effect, and the reader judges.
- Counting artefacts as biology: low gene counts from incomplete annotations, and non-independent genomes from the same genus.

## Constraints

- Heavy work runs in SLURM with small requests (the account CPU limit has delayed large jobs by hours; the figure job needed 2 CPUs). Use `$SCRATCH:?` for temp files, results on `/bigdata`.
- Pfam family definitions differ from the paper's InterPro-based ones. The mapping from Pfam to the paper's family names is a hand-curated table with its source, and unmapped families are listed, not hidden.
- Group tests use one value per genus (genus representatives, as in the ordination) to limit pseudo-replication. Phylogeny is not otherwise controlled. This is stated on every page.
- Multiple tests use Benjamini-Hochberg correction. Effect sizes are always shown.
- BUSCO completeness is a covariate or a filter (default: BUSCO complete >= 80%, with a sensitivity run without the filter), because annotation completeness changes TF counts.

## Design

### A. Landscape and contrasts

1. TF table builder (SLURM, DuckDB over `pfam.parquet` plus `tf_inventory.parquet`): per genome and per family: number of proteins with at least one domain of the family; proteome size; TF fraction. Extended family set = the 21 BFD families + E2F_TDP (PF02319), Tub (PF01167), RFX (PF02257), WRKY (PF03106), CP2 (PF04516), CBF/NF-Y (PF00808), Copper-fist (PF00649), TEA (PF01285), KilA-N/APSES (PF04383), STAT (Pfam accessions verified against `pfam_names`). The builder verifies every accession against `pfam_names` and fails on an unknown one.
2. Mapping table `config/tf_family_map.tsv`: Pfam accession, BFD family, paper family name, note. Unmapped paper families are listed in the data-notes page.
3. Diversity: Shannon index of the family distribution per genome, as the paper's reference does (the builder cites the definition it uses).
4. Figures: TF fraction of proteome by lineage; family composition by lineage; presence/absence heatmap by order for the rare families (E2F, Tubby, RFX, STAT, WRKY); TF count against proteome size; Shannon diversity by lineage. All use the global phylum colour map.
5. Claims table: for each claim, the groups, n genera, the median per-genus value in each group, effect size, bootstrap interval over genera, BH-corrected p, and a plain "direction agrees with the paper: yes / no / mixed" column derived from the sign only.

### B. Ortholog tracking

1. Seeds: WC-2 (GATA + PAS), Zap1-like C2H2 (zinc-responsive), bZip5 (yeast Cad1/Yap2-like). Two seed sources, in this order: (a) the paper's *R. microsporus* IDs if the user adds the supplementary tables to `resources/`; (b) literature seeds from UniProt (*N. crassa* WC-2, *S. cerevisiae* Zap1, *S. cerevisiae* Cad1/Yap2, *M. lusitanicus* WC proteins). The seed list and its source are written to a table.
2. Search: the proteins carrying the relevant TF Pfam domain in genus representatives (one genome per genus, best BUSCO). Seed-vs-candidates with MMseqs2 or DIAMOND (whichever module exists on HPCC; the builder checks), keep hits with coverage and identity thresholds stated in the output; verify domain architecture with the existing Pfam hits (for example WC-2 needs GATA plus PAS).
3. Align (MAFFT), trim, build a tree (FastTree; IQ-TREE only for the final small trees). Report: presence and copy number by order, architecture counts, the tree with phylum colours, and the R. microsporus / M. lusitanicus / Phycomyces members marked.
4. Caveat on every page: sequence similarity and tree position are evidence for orthology, not proof; function is not inferred.

### C. Expansions and contractions

1. Unit: genus representatives with BUSCO complete >= 80%. Family count per genus; scaling model log(count + 1) against log(proteome size) across all genera (global fit per family); residual per genus.
2. Lineage test: for each order (and class) with at least 5 genera, compare residuals with the rest of its phylum (Mann-Whitney, BH across all family-by-lineage tests, effect size = rank-biserial correlation). Flag |effect| and corrected p above stated cutoffs, and only if at least 3 genera agree in direction.
3. Presence/absence: genus-level presence of each family by order; flag lineages where a family that is present in at least 80% of genera is absent from at least 80% of the lineage's genera (and the reverse).
4. Outputs: ranked table, a lineage-by-family effect heatmap, per-hit boxplots, all marked "hypothesis". A sensitivity run without the BUSCO filter is shown next to the main run.

### Site integration

New pages after Explore: TF landscape, TF contrasts, TF orthologs, TF expansions. Figures, tables and the index are added to the release asset (stacked on the release-asset work). A small `tf_genomes.json` for the Explore page is optional and deferred. The suppress list applies to all inputs.

## Approaches considered

- **1. Extend the existing summary pipeline and the release-asset flow (chosen).** One new builder script and figure modules; B as a separate small SLURM workflow under `bin/tf/`; pages read tables like the existing ones. Fast, consistent with what is already reviewed.
- 2. Put everything into the BFD Nextflow pipeline (`MERGE_TF` etc.). The right long-term home, but slow to iterate and needs a pipeline release. Deferred.
- 3. Full orthology with OrthoFinder on genus representatives. Too heavy for the question; B's seed-based search answers the three TFs.

## Decisions I made (not yet confirmed)

- Stack the work on branch `quarto-site` as a new branch `tf-survey` after the release-asset work lands. One PR for the site, one for this.
- Genus representatives for tests; BUSCO >= 80% default with a sensitivity run.
- Pfam-based family set with an explicit mapping table; the 21 BFD families are kept, 10 families are added.
- B uses seed search plus trees, not OrthoFinder.
- No verdict column beyond the sign of the effect.

## Left to the builder

- MMseqs2 versus DIAMOND (whichever module exists), identity and coverage cutoffs (stated in the output).
- Exact effect-size and p-value cutoffs for flagging in C (stated on the page).
- Page layout and figure styles (reuse the existing helpers).

## Open items

- The paper's supplementary tables would improve B. If the user adds them to `resources/`, B uses the paper's own protein IDs. The work does not wait for them.
- Paper family names that have no Pfam equivalent in the extended set (for example STAT, if no accession exists in `pfam_names`) are reported as unmapped.
