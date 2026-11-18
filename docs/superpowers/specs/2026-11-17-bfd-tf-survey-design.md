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
- Supplementary material, added by the user on 2026-11-17 in `resources/s41467-025-65177-1_suppl/`:
  - `MOESM3` (xlsx): sheet `TF count` = the paper's TF counts, 51 family rows plus totals, for 62 species (R. microsporus: 605 TFs of 10,700 genes); sheet `Genomes` = the 62 species with their sources (63 JGI MycoCosm, 3 NCBI). This gives the paper's exact family names and a direct way to calibrate BFD's Pfam-based counts.
  - `MOESM4` (xlsx): 57 sheets, one per TF (for example `Znf1`, `GATA Znf10`, `bZip5`), each with the JGI protein IDs of its DAP-seq targets (21,313 rows, 8,290 unique IDs).
  - `MOESM5` (pdf, 5 pages): supplementary figures (bZip target overlap and similar).
  - `MOESM6` (xlsx): 96 protein sequences of the tested TFs, IDs `Rhimi59_2_<JGI id>` (some are DNA-binding-domain fragments with `a`/`b` suffixes).
  - `MOESM7` (xlsx): primers and crRNAs. Not used.
  - `MOESM1` (pdf, 21 pages; added 2026-11-17): Supplementary Information. Supplementary Table 1 lists the 88 TFs tested by DAP-seq with JGI ID, name, peaks, FRIP and family. It is parsed into `config/paper_tf_names.tsv` (88 rows; 58 have FRIP above 5% in at least one replicate, as the paper reports; 40 pass in both). All 88 IDs are in `MOESM6` and in the JGI protein file; the 57 sheet names of `MOESM4` match the table names (`CP2` is `CP2_TF1`). So every tested TF has a name, a full-length sequence and a target list. The three TFs the paper highlights: Znf1 = JGI 1799646 (214 aa), WC-2C = `GATA_Znf10` = JGI 1921537 (369 aa), bZip5 = JGI 1909793 (270 aa). The Supplementary Table 1 replicate labels (Exp 1 and Exp 2) appear swapped relative to the paper's text; it does not matter here.
- JGI download, added by the user in `resources/jgi_atcc_11559/` (JGI MycoCosm `Rhimi59_2`, *R. microsporus* ATCC 11559 v2.0): all protein models (`Rhimi59_2_all_proteins_20220409.aa.fasta.gz`, 258,092 models, headers `jgi|Rhimi59_2|<proteinId>|<name>`) and the gene catalog of the filtered best models (GFF3 and GFF, with `product_name` and `proteinId`). All 92 unique tested-TF IDs from `MOESM6` and all 8,290 unique target IDs from `MOESM4` are in the protein file. So the full-length sequences of the tested TFs and of all targets are available.
- BFD's `FBF4BF29` is a different assembly of the same strain (`GCA_002083735.1_Rhimi_ATCC11559_1`, funannotate, 10,256 proteins; JGI reports about 10,700 genes). Coordinates do not transfer. JGI and BFD proteins are mapped by sequence (reciprocal best hit).
- Binding sites and motifs are experimental. They are not extended to other genomes. What can be extended is the gene sets: the targets mapped into BFD and checked for conservation across fungi.

## Goals

- A: one TF family table for all genomes with an extended family set; landscape figures by phylum, class and order; a calibration of BFD's Pfam-based family counts against the paper's own counts (`MOESM3`) for the species present in both; a table that tests each claim above at BFD scale and, side by side, on the paper's own 62-species table.
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

1. Seeds: WC-2 (GATA + PAS), Zap1-like C2H2 (zinc-responsive), bZip5 (yeast Cad1/Yap2-like). Seed sources, in this order: (a) the paper's own sequences (`MOESM6`, 96 proteins): search them against the BFD *R. microsporus* ATCC 11559 proteome (`FBF4BF29`) to find the matching BFD proteins (near-identical hits), then name them from Supplementary Table 1 (Znf1 = JGI 1799646, WC-2C = `GATA_Znf10` = JGI 1921537, bZip5 = JGI 1909793; the domain check, GATA plus PAS for WC-2C, is a sanity test of the match); (b) literature seeds from UniProt (*N. crassa* WC-2, *S. cerevisiae* Zap1, *S. cerevisiae* Cad1/Yap2). The seed list, source and match result are written to a table.
2. Search: the proteins carrying the relevant TF Pfam domain in genus representatives (one genome per genus, best BUSCO). Seed-vs-candidates with MMseqs2 or DIAMOND (whichever module exists on HPCC; the builder checks), keep hits with coverage and identity thresholds stated in the output; verify domain architecture with the existing Pfam hits (for example WC-2 needs GATA plus PAS).
3. Align (MAFFT), trim, build a tree (FastTree; IQ-TREE only for the final small trees). Report: presence and copy number by order, architecture counts, the tree with phylum colours, and the R. microsporus / M. lusitanicus / Phycomyces members marked.
4. Caveat on every page: sequence similarity and tree position are evidence for orthology, not proof; function is not inferred.
5. Target sets and their conservation: map the JGI protein IDs (the 92 tested TFs and the 8,290 targets) to BFD `FBF4BF29` proteins by reciprocal best hit (MMseqs2; thresholds and the mapping rate are reported). For each TF's target set in BFD, search the targets against the genus-representative proteomes and report, per target, the number of phyla and orders with a hit; compare each TF's targets with all genes of the genome (are targets more conserved or more lineage-specific than the genome average?). Add the Pfam and SwissProt profile of each target set from the BFD tables. These are hypotheses about what the targets are, not about regulation in other species.

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

- No input is missing now. `MOESM1` ties names to IDs. The JGI download resolves the ID mapping: JGI and BFD proteins map by sequence, not by coordinates.
- The JGI protein file holds 258,092 models, filtered and not. The "best" filtered models (about 10,700) come from the GFF3. The mapping uses the filtered set as the reference and reports how many tested TFs and targets fall outside it.
- Paper family names that have no Pfam equivalent in the extended set (for example STAT, if no accession exists in `pfam_names`) are reported as unmapped.
