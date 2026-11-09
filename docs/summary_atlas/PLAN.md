# BFD summary atlas: design plan

Date: 2026-11-08. Author: Claude (design only; no figures built, no SLURM jobs submitted).
Scope: a static Markdown atlas of hypothesis-driven summary views over the BFD tables.
Status of every number in this file: **measured** on the login node with `/usr/bin/python3.12` + DuckDB on 2026-11-08 unless marked "not measured" or "estimate".

Data root: `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables` (parquet snapshot, 134 GB).
DuckDB file with views over the same parquet: `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/db/BFD.duckdb`.

---

## 1. Data coverage (measured)

Genome keys: `species.parquet` has 23,680 rows, 23,680 distinct `LOCUSTAG`, 23,680 distinct `ASMID`, 7,450 distinct `SPECIES`, 1,671 distinct `GENUS`. `LOCUSTAG` == `species_prefix` in the functional tables. `ASMID` joins `asm_stats`, `busco_genome`, `telomere_*`.

| Table | Key | Genomes covered | Rows | Size | Note |
|---|---|---|---|---|---|
| species | LOCUSTAG / ASMID | 23,680 | 23,680 | 1 MB | `SUBPHYLUM` is NULL in every row. `CLASS` NULL in 218 rows, `ORDER` NULL in 300, `PHYLUM` NULL in 3. |
| samples | LOCUSTAG / ASMID | 23,680 | 23,680 | 2 MB | Adds `TRANSL_TABLE`: 1 = 22,525; 12 = 1,154 (all Pichiomycetes); 26 = 1. |
| asm_stats | ASMID | 23,584 | 23,584 | 3 MB | `annotation_status`: annotated 23,582, unannotated 2. N50 quantiles (5/25/50/75 %): 12 kb / 90 kb / 307 kb / 1.07 Mb. |
| busco_genome | ASMID | 23,582 | 23,582 | 1 MB | All `fungi_odb12`, 1,122 markers. Genome-mode BUSCO. |
| telomere_summary | ASMID | 23,582 | 23,582 | 2 MB | 8,369 genomes have >= 1 telomere tract. |
| telomere_tracts | ASMID | 8,369 | 73,242 | 18 MB | |
| aa_freq | species_prefix | 23,452 | 477,535 | 2 MB | 20 amino acids per genome; `*` in 6,254 genomes, `X` in 2,241. 21 prefixes are not in `species`. |
| codon_freq | species_prefix | 23,452 | 1,500,925 | 5 MB | 64 codons per genome; frequencies sum to 1.0. |
| gene_trnas | gene_id prefix | 23,446 | 4,558,316 | 35 MB | Columns: gene_id, amino_acid, codon. |
| tf_inventory | species_prefix | 23,394 | 10,946,927 | 167 MB | One row per TF domain hit; `tf_family` has 21 values. |
| cazy.overview | species_prefix | 23,401 | 9,310,657 | 122 MB | `toolcount` 1: 867,612 rows; 2: 8,443,045 rows. Has `substrate`. |
| cazy.cazymes_hmm | species_prefix | 23,394 | 9,164,068 | 246 MB | |
| signalp.signal_peptide | species_prefix | 23,410 | 18,767,284 | 265 MB | Per-genome SP count quantiles (10/50/90 %): 230 / 783 / 1,470. |
| targetP | species_prefix | 125 | 178,535 | 5 MB | Almost no coverage. Exclude. |
| infernal_rfam | species_prefix | 14 | 3,925 | 1 MB | Exclude. |
| wgd.ks.summary | species_prefix | 7,384 | 7,384 | 1 MB | Interim release (manifest says 5 giant genomes still running). |
| wgd.ks.nodes | species_prefix | 7,384 | 28,982,469 | 570 MB | |
| wgd.ks | species_prefix | 7,384 | 25,831,001 | 2.4 GB | |
| swissprot_subject | sp_id | n/a | 575,503 | 26 MB | Reference proteins. `org_code` gives the source organism. |
| swissprot_subject_summary | sp_id | n/a | 333,632 | 4 MB | 114,004 subjects hit in >= 1,000 genomes; 76,403 in >= 5,000; 58,179 in >= 10,000. |
| swissprot_best | species_prefix | not measured (>= 158 from row-group stats) | 157,614,964 | 4.9 GB | One best hit per protein. Columns `q_len`, `s_len`, `pident_x10`, `qcov_x10`, `scov_x10`. |
| swissprot_hits | species_prefix | not measured | 2,192,440,076 | 37 GB (256 files) | All ranked hits (`rk`). |
| pfam | protein_id only | not measured (>= 2,790) | 368,209,721 | 19 GB | No genome column. genome = `split_part(protein_id,'_',1)`. |
| gene_info | LOCUSTAG | not measured (>= 1,920) | 249,546,610 | 3.8 GB | `gene_type` values seen: protein_coding, tRNA_gene. |
| gene_transcripts | gene_id | not measured | 249,546,610 | 6.5 GB | `is_partial` BOOLEAN; `has_start_codon`/`has_stop_codon` are strings 'TRUE'/'FALSE'/'NULL'. In one 200,000-row sample: 80.6 % TRUE/TRUE and not partial; 2.2 % NULL/NULL. |
| gene_proteins | gene_id / protein_id | not measured | 244,988,294 | 5.2 GB | `length` (aa), `md5checksum`. |
| gene_CDS | transcript_id | not measured | 753,622,745 | 16 GB | |
| gene_exons | transcript_id | not measured | 758,181,061 | 16 GB | |
| gene_introns | transcript_id | not measured | 508,634,451 | 6.5 GB | Has `splice_5`, `splice_3`, `GC_content`. |
| gene_intergenic_distances | species_prefix | not measured (>= 1,820) | 233,406,624 | 2.9 GB | |
| merops | species_prefix | not measured (>= 352) | 43,462,876 | 1.4 GB | BLAST-style rows; several per protein. |
| tmhmm | species_prefix | not measured (>= 370) | 45,707,260 | 1.8 GB | `PredHel` helix count. |
| wolfpsort | species_prefix | not measured (>= 1,985) | 244,790,951 | 1.8 GB | One row per protein, `localization`. |
| predgpi | species_prefix | not measured (>= 1,828) | 227,238,177 | 1.9 GB | GFF-like; `feature` = Chain or GPI-anchor. |
| idp | species_prefix | not measured | 233,329,228 | 4.0 GB | Per disordered region. |
| idp_summary | species_prefix | not measured (>= 1,170) | 144,924,594 | 2.1 GB | **Key form differs**: `protein_id` looks like `F0000039_000002` (gene_id form, no `-T1`). Verify in the build before joining. |

"Genomes covered: not measured" means a full scan is needed. Build job J0 (section 6) computes `count(distinct genome)` for every table and writes `tables/coverage.tsv`; view V01 reads that file.

Physical layout (measured from parquet metadata): all large tables are written in genome order. Row groups hold about 125,000 rows (swissprot_best: about 1,000,000). A filter on a genome prefix range (`gene_id >= 'X_' AND gene_id < 'X`'`) or `species_prefix = 'X'` prunes to a few row groups. Pilot: 16 genomes from `gene_proteins` read in 12-16 s with 4 threads; the same genomes from `swissprot_best` in 0.5-1.0 s. Full per-genome aggregation of a table is a streaming group-by and does not need the whole table in memory.

Protein key forms: `protein_id` = `<prefix>_<6 digits>-T<n>` in pfam, swissprot_best, signalp, tmhmm, wolfpsort, predgpi, cazy, merops, tf_inventory, gene_proteins. `gene_id` = `<prefix>_<6 digits>` in gene_info, gene_transcripts, gene_trnas, gene_intergenic_distances and (sample-verified) idp_summary.

### 1.1 Taxon sampling (measured)

- Phyla: Ascomycota 19,412 genomes (5,176 species); Basidiomycota 3,275 (1,666); Mucoromycota 294; Kickxellomycota 190; Microsporidia 163; Mortierellomycota 100; Chytridiomycota 89; Glomeromycota 60; others < 30 each.
- Largest classes: Sordariomycetes 6,016; Eurotiomycetes 3,220; Pichiomycetes 2,967; Dothideomycetes 2,731; Saccharomycetes 2,647; Agaricomycetes 1,856; Tremellomycetes 559.
- Strains per species: 5,575 species have 1 genome; 859 have 2; 321 have 3. Top: *Saccharomyces cerevisiae* 1,309; *Candidozyma auris* 952; *Fusarium oxysporum* 814; *Pyricularia oryzae* 495; *Macrophomina phaseolina* 466; *Aspergillus fumigatus* 416.
- Over-sampling is therefore extreme. Every cross-clade statistic in this atlas uses a representative set (section 4.3).

### 1.2 BUSCO completeness differs by class (measured; drives the quality filter)

`busco_genome.complete_pct` (fungi_odb12, genome mode), median by class: Sordariomycetes 99.3; Eurotiomycetes 99.1; Dothideomycetes 99.2; Agaricomycetes 98.1; Tremellomycetes 96.4; Microbotryomycetes 96.1; **Pichiomycetes 80.8; Saccharomycetes 75.4; Schizosaccharomycetes 83.6; Dipodascomycetes 89.4; Lipomycetes 93.0; Malasseziomycetes 89.9; Wallemiomycetes 87.3**. Duplicated median: Pucciniomycetes 56.6 (dikaryotic assemblies), Mucoromycetes 12.6, Mortierellomycetes 7.4.

Consequence: a fixed cutoff such as "complete >= 95 %" removes essentially all yeasts (a pilot that asked for Saccharomycetes/Pichiomycetes genomes at >= 95 % returned none). The quality filter must be class-relative (section 4.4).

---

## 2. What exists for orthology, lifestyle and phylogeny (searched 2026-11-08)

### 2.1 Orthology: no table. Proxies and their biases

Searched: `Fungi_BFD_runs` (excluding `genome_annotation/` and `work/`), `Fungi_BFD`, `BFD`, `BFD_Comparative`, `DeltaGain*`. No OrthoFinder, no OrthoMCL, no BUSCO `full_table.tsv`. `Fungi_BFD/nextflow/modules/BFD/BUSCO_PEP/main.nf` runs BUSCO in `$SCRATCH` and emits only `<locustag>.BUSCO_summary.<lineage>.txt`; per-gene tables were not kept. `results/genome_stats/BUSCO_protein/` and `BUSCO_genome/` hold summaries only. `wgd.ks` families are within-genome paralog families, not cross-genome orthogroups.

| Proxy | Family label | Available now | Main biases |
|---|---|---|---|
| P-A: shared SwissProt best hit | `swissprot_best.sp_id` with `qcov_x10 >= 800`, `scov_x10 >= 800`, `pident_x10 >= 400` | Yes (4.9 GB table) | (1) Reference dominance: in the pilot (13 yeast-class genomes) 50-60 % of best hits of Saccharomycotina genomes went to `org_code = 'YEAST'` (*S. cerevisiae*), versus 8-9 % for Tremellomycetes. A yeast and a filamentous genome often hit *different* reference proteins for the same gene family, so `q_len/s_len` is confounded by which reference was hit. Compare `q_len` across genomes **within the same sp_id** so `s_len` cancels; do not compare ratios to different references. (2) The coverage gates remove truncated models and also remove real length variants (long isoforms, fused domains). The test is conservative. (3) Paralogs collapse onto one sp_id; use the per-genome median per sp_id. |
| P-B: model-organism anchor | best hit to each of YEAST, SCHPO, NEUCR, EMENI from `swissprot_hits` (DuckDB view `v_swissprot_model_fungi` already encodes this) | Yes, but 2.2 G rows / 37 GB scan | Removes bias (1): every genome is anchored to the same *S. cerevisiae* protein. Still biased by gates (2). One SLURM job (J5). |
| P-C: Pfam domain architecture | ordered string of `pfam_acc` per protein (domain_i_evalue <= 1e-5, env coordinates), from `pfam.parquet` | Yes (19 GB scan) | A truncated model loses domains and lands in a different architecture class, so truncation drops out of the family instead of making it shorter. Linker/terminal length differences are what remains. Single-domain architectures are heterogeneous. |
| P-D: BUSCO single-copy orthologs | BUSCO protein-mode `full_table.tsv` per genome | **No**; needs a re-run of BUSCO protein mode on 23,400 proteomes with the full table kept (about 23,400 jobs; cost not measured) | Cleanest families (1,122 markers), but biased to conserved housekeeping genes that have conserved length by construction. |

Pilot numbers (8 non-yeast and 13 yeast-class genomes chosen as best-BUSCO per genus; **sampled, not a result**): raw median protein length per genome was 358-466 aa in Sordariomycetes/Eurotiomycetes and 394-421 aa in Saccharomycetes/Pichiomycetes, 377 in *Schizosaccharomyces osmophilus*, 407-451 in Tremellomycetes. Strict-gate best hits (qcov, scov >= 80 %, pident >= 40 %) were 1,400-4,000 per genome (about 30-85 % of best hits). The pilot is too small to support a direction. It shows the gates keep enough proteins per genome for a family-level test.

### 2.2 Lifestyle and growth form: FungalTraits 1.2 is on disk

`/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/refs/lifestyle/FungalTraits_1.2.csv` (2.1 MB, latin-1 encoding; must be converted to UTF-8 before DuckDB reads it). 10,771 rows, 10,765 genera. Columns include `GENUS`, `primary_lifestyle`, `Secondary_lifestyle`, `Growth_form_template`, `Decay_substrate_template`, `Plant_pathogenic_capacity_template`, `Animal_biotrophic_capacity_template`, `Ectomycorrhiza_lineage_template`. One genus has two growth forms; two genera have two primary lifestyles.

Genus join to `species.parquet` (measured): 21,644 of 23,680 genomes, 6,783 of 7,450 species, 1,496 of 1,671 genera matched.

Growth form of matched genomes: filamentous_mycelium 14,512 (4,586 species, 1,011 genera); yeast 5,507 (1,305 species, 135 genera); dimorphic_yeast 811 (285 species, 69 genera); thallus_photosynthetic 535; unicellular-aflagellate (non-yeast; Microsporidia and relatives) 134; zoosporic forms 134; NULL 10.

Primary lifestyle of matched genomes: plant_pathogen 6,979 (1,489 species); nectar/tap_saprotroph 3,861 (798 species; mostly *Saccharomyces*); unspecified_saprotroph 3,439; animal_parasite 1,781; wood_saprotroph 1,134; soil_saprotroph 1,125; litter_saprotroph 852; ectomycorrhizal 663; lichenized 508; mycoparasite 410; others < 200.

Unmatched genera with the most genomes: *Candidozyma* 987 (Pichiomycetes), *Lizanozyma* 112, *Hansenula* 40, *Petasospora* 33, *Maudiozyma* 29, *Monosporozyma* 23, *Linnemannia* 21, *Sungouiella* 21, *Podila* 16, *Entelexis* 13, *Taiozyma* 11, *Dujonia* 11; plus about 300 genomes with NULL genus. Most unmatched names are post-2023 yeast genus names that FungalTraits 1.2 predates. A manual override file of about 30 genera fixes most of the gap (decision D2).

Taxonomy-only fallback for "yeast" (no external data): classes Saccharomycetes, Pichiomycetes, Dipodascomycetes, Lipomycetes, Trigonopsidomycetes, Alloascoideomycetes, Sporopachydermiomycetes (Saccharomycotina; 6,100 genomes in the five largest), Schizosaccharomycetes (68), and the basidiomycete yeast classes Tremellomycetes (559), Microbotryomycetes (318), Cystobasidiomycetes (28), Malasseziomycetes (85). Ustilaginomycetes and Exobasidiomycetes are dimorphic. Pezizomycotina dimorphic yeasts (FungalTraits: Eurotiomycetes 283 genomes, Dothideomycetes 190) are not separable by taxonomy alone.

Views that need labels that do not exist yet: none are blocked. All lifestyle views (H4, H5) use the FungalTraits genus join and report the unmatched fraction. The FungalTraits label is a genus trait, not a genome trait; strain-level lifestyle variation (e.g. *Fusarium oxysporum* f. sp.) is invisible.

### 2.3 Phylogeny: none for this genome set

Found only example scripts: `BFD_dashboard/workflows/run_phyling.sh`, `Fungi_BFD/nextflow/run_phyling.sh` (nf_phyling on a small protein set). No `*.treefile`, `*.nwk`, or phyling result directory exists under `Fungi_BFD_runs/results`, `BFD_Comparative`, or `BFD_dashboard`. The dashboard config has a `treefile` option but no file. PGLS is therefore not possible now. Controls that do not need a tree are in section 4.3. Building a tree is decision D5.

---

## 3. Hypotheses

| ID | Statement | Expected direction (PI or prior) | Views |
|---|---|---|---|
| H1 | Annotation and assembly quality differ by clade and by growth form, and must be modelled before any cross-clade comparison. | Yeast classes have lower fungi_odb12 BUSCO (measured, section 1.2). Fragmented assemblies have more partial models (expected, not measured). | V01-V06 |
| H2 | Among orthologous proteins, yeast proteins are shorter than those of filamentous fungi. | PI hypothesis: yeast shorter. No measurement yet supports or refutes it. | V07-V14 |
| H3 | Yeast genomes have shorter genes, fewer introns and shorter intergenic distances than filamentous genomes; gene density scales with genome size. | Yeast: fewer introns, denser genes (textbook expectation; not measured here). | V15-V19 |
| H4 | Secretome size and CAZy/peptidase content track lifestyle: plant-associated and wood/litter saprotrophs carry more secreted proteins and plant-cell-wall CAZymes than yeasts and animal parasites. | Plant pathogen, wood saprotroph > yeast, animal parasite. Direction from prior literature, unverified here. | V20-V26 |
| H5 | Transcription-factor and domain repertoires scale sub-linearly with gene count, and Zn2Cys6 expansion is specific to Pezizomycotina. | Zn2Cys6 share highest in Pezizomycotina (measured counts show it is the largest TF family overall; clade split not measured). | V27-V32 |
| H6 | Codon usage and amino-acid composition track taxonomy more than growth form or lifestyle; GC3 follows genomic GC; tRNA gene pools match codon usage. | Taxonomy explains more variance than lifestyle (expected; no measurement). | V33-V38 |
| H7 | Recent paralog retention per 1,000 genes is lower in yeasts than in filamentous fungi; large identical-copy groups mark assembly artefacts rather than biology. | Yeast lower (expected; not measured). Framing follows `paralogoscope_pilot/eval/explore/fable_context.md`. | V39-V40 |
| H8 | Telomere repeat monomer is conserved within classes and departs from TTAGGG in Saccharomycotina. | Saccharomycotina have degenerate repeats (prior knowledge; not measured). | V41-V42 |

Total: 8 hypotheses, 42 views. Priority counts: P1 = 22, P2 = 15, P3 = 5.

---

## 4. Shared definitions used by every view

### 4.1 Label table `genome_labels` (intermediate A1)

One row per genome (LOCUSTAG, ASMID) with: taxonomy (PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES), `growth_form` (FungalTraits + override file), `primary_lifestyle`, `yeast_origin` (categorical: Saccharomycotina, Taphrinomycotina, Tremellomycetes, Pucciniomycotina-yeast, Ustilaginomycotina-dimorphic, Pezizomycotina-dimorphic, filamentous, other), `transl_table`, `is_representative`, `quality_tier`, and the quality covariates (busco_complete, busco_duplicated, N50_bp, contig_count, gc_pct, partial_fraction, gene_count).

### 4.2 Yeast / filamentous definition

Default: `growth_form` = 'yeast' versus 'filamentous_mycelium'. 'dimorphic_yeast' is shown as a third group and never pooled. Genomes with no label are excluded from H2-H4 and counted in V01. Sensitivity run: taxonomy-only definition (section 2.2). Both are reported for H2.

### 4.3 Units of replication and phylogenetic controls (no tree)

1. **Representative set (default unit)**: one genome per species: highest `complete_pct` within the species, tie broken by `N50_bp`. Gives at most 7,450 genomes. Removes the strain oversampling listed in 1.1.
2. **Genus-level aggregation**: median of the representative-genome values per genus (1,671 genera; 1,496 with lifestyle). Used for every lifestyle test, because the lifestyle label is a genus trait.
3. **Independent yeast origins as replicates**: Saccharomycotina, Schizosaccharomycetes, Tremellomycetes, Microbotryomycetes + Cystobasidiomycetes, Malasseziomycetes, Pezizomycotina dimorphic. A claim about "yeasts" holds only if the sign is the same in the majority of origins with n >= 5 genera.
4. **Paired sister-clade contrasts**: fixed list, decided with the PI (D3): Schizosaccharomycetes vs Taphrinomycetes + Pneumocystomycetes; Tremellomycetes yeast genera vs Agaricomycetes; Saccharomycotina vs Pezizomycotina (class medians); Malasseziomycetes vs Ustilaginomycetes; Eurotiomycetes dimorphic genera vs other Eurotiomycetes; Microbotryomycetes vs Pucciniomycetes. Test: Wilcoxon signed-rank on the pair differences of genus medians is not valid with 6 pairs; report the sign and effect per pair and a sign test across pairs.
5. **Nested random-effects model**: response per representative genome; random intercepts for CLASS / ORDER / FAMILY / GENUS; fixed effects for growth form and the quality covariates (log N50, busco_complete, partial_fraction). Implemented with `statsmodels` MixedLM (variance components). This is the headline model for H2-H6.
6. **PGLS**: not possible without a tree (section 2.3).

Kruskal-Wallis across classes is used only as a screening statistic and is labelled as such on every page.

### 4.4 Quality filters

| Filter | Definition | Used by |
|---|---|---|
| Q-ann | `annotation_status = 'annotated'` and genome present in `gene_proteins` | all |
| Q-busco-rel | `complete_pct` >= 25th percentile of its CLASS (classes with >= 20 genomes; otherwise >= 25th percentile of its PHYLUM) | all length and count views |
| Q-dup | `duplicated_pct` <= 10 (drops most dikaryotic/unpurged assemblies; Pucciniomycetes largely excluded, stated on the page) | H2, H3, H7 |
| Q-n50 | `N50_bp` >= 50,000 | H2, H3 |
| Q-model | protein-level: `is_partial = false`, `has_start_codon = 'TRUE'`, `has_stop_codon = 'TRUE'` | H2, H3, H4 counts |
| Q-rep | representative set (4.3.1) | all cross-clade statistics |

Every page states which filters it applied and the resulting n (genomes, species, genera).

---

## 5. Views

Format for each view: ID, title, question, tables/columns, SQL sketch (over the intermediates of section 6 where possible), statistic and unit, figure type, filters, caveats, compute, priority.

### H1. Annotation and assembly quality by clade

**V01 Data coverage table** (P1). Question: which genomes are in which table? Tables: all. SQL: `SELECT 'pfam' AS t, count(DISTINCT split_part(protein_id,'_',1)) FROM pfam` per table (job J0). Statistic: counts. Output: Markdown table (this file's section 1 plus the measured large-table counts). Filters: none. Compute: one pass over every table (134 GB). Caveat: counts change when tables are republished; record `built_at` from `_manifest`.

**V02 BUSCO completeness and duplication by class** (P1). Question: how does fungi_odb12 completeness vary by class and growth form? Tables: busco_genome, species, genome_labels. SQL: `SELECT CLASS, growth_form, quantile_cont(complete_pct,[0.1,0.5,0.9]), median(duplicated_pct), count(*) FROM genome_labels GROUP BY 1,2`. Statistic: quantiles per class; no test. Figure: horizontal strip/box plot, classes with >= 20 genomes, colored by growth form; second panel duplicated_pct. Filters: Q-ann. Caveat: fungi_odb12 is a single lineage set; low yeast values partly reflect marker set bias, not annotation quality. Compute: small. Already measured in section 1.2.

**V03 Gene count versus assembly size** (P1). Question: do yeast genomes have fewer genes per Mb? Tables: genome_gene_stats (A2), asm_stats, genome_labels. SQL: `SELECT g.n_protein_coding, a.total_length_bp, l.growth_form, l.CLASS FROM A2 g JOIN asm_stats a USING (ASMID) JOIN genome_labels l USING (LOCUSTAG) WHERE is_representative`. Statistic: log-log OLS slope per growth form; genus medians. Figure: scatter (log-log), points = representative genomes, color = growth form; marginal density. Filters: Q-rep, Q-ann, Q-busco-rel. Caveat: gene count depends on the funannotate model set and on masking (`masked_pct` available as a covariate). Compute: reads A2 (small).

**V04 Partial gene models versus assembly contiguity** (P1). Question: how much of the proteome is truncated, and does it depend on N50 and clade? Tables: A2 (`partial_fraction`, `frac_no_start`, `frac_no_stop`), asm_stats. SQL: `SELECT LOCUSTAG, partial_fraction, N50_bp, contig_count, CLASS FROM A2 JOIN asm_stats JOIN genome_labels`. Statistic: Spearman rho of partial_fraction vs log N50, per class. Figure: scatter + per-class rho table. Filters: Q-ann. Caveat: `has_start_codon` has a 'NULL' string class (2.2 % in a sample); report it separately. Compute: small after A2.

**V05 Representative set and quality tiers** (P1). Question: what does each filter remove? Tables: genome_labels. SQL: `SELECT CLASS, growth_form, count(*) FILTER (WHERE is_representative), count(*) FILTER (WHERE is_representative AND q_busco_rel AND q_dup AND q_n50) FROM genome_labels GROUP BY 1,2`. Statistic: counts. Output: table (per class) + waterfall bar. Filters: shows them. Compute: small.

**V06 Strain oversampling** (P2). Question: which species dominate the raw genome set? Tables: species. SQL: `SELECT SPECIES, count(*) FROM species GROUP BY 1 ORDER BY 2 DESC`. Figure: histogram of genomes per species (log x) + top-20 table. Numbers already measured (1.1). Compute: trivial.

### H2. Yeast proteins are shorter among orthologs

**V07 Raw protein length distribution by growth form and class** (P1). Question: is there a raw length difference before any ortholog control? Tables: A2 (length quantiles per genome), genome_labels. SQL: `SELECT growth_form, CLASS, median(med_protein_len) FROM A2 JOIN genome_labels WHERE is_representative AND q_filters GROUP BY 1,2`. Statistic: genus-median lengths; Mann-Whitney yeast vs filamentous on genus medians; nested mixed model (4.3.5). Figure: ridge/violin of per-genome median protein length by class, colored by growth form. Filters: Q-rep, Q-ann, Q-busco-rel, Q-dup, Q-n50, Q-model. Caveat: raw lengths mix gene content (families present) with length within families; this view only sets the baseline. Compute: small after A2.

**V08 Within-family length contrast, SwissProt proxy (P-A)** (P1). Question: for the same reference family, are yeast proteins shorter? Tables: protein_core (A3), family_length (A4), genome_labels. SQL sketch:
```sql
WITH fam AS (
  SELECT sp_id, genome, median(length) AS len            -- one value per family per genome
  FROM A3 JOIN genome_labels USING (genome)
  WHERE qcov_x10>=800 AND scov_x10>=800 AND pident_x10>=400
    AND is_partial=false AND has_start='TRUE' AND has_stop='TRUE' AND is_representative
  GROUP BY 1,2),
per_fam AS (
  SELECT sp_id,
    median(len) FILTER (WHERE growth_form='yeast')       AS yeast_len,
    median(len) FILTER (WHERE growth_form='filamentous_mycelium') AS fil_len,
    count(DISTINCT GENUS) FILTER (WHERE growth_form='yeast') AS n_yeast_genera,
    count(DISTINCT GENUS) FILTER (WHERE growth_form='filamentous_mycelium') AS n_fil_genera
  FROM fam JOIN genome_labels USING (genome) GROUP BY 1)
SELECT sp_id, log2(yeast_len/fil_len) AS log2_ratio FROM per_fam
WHERE n_yeast_genera>=10 AND n_fil_genera>=10;
```
Statistic: unit = family (sp_id). Distribution of log2 ratio; median with bootstrap CI; sign test (fraction of families with ratio < 1). Repeat with genus-median instead of genome-median inside each family. Figure: histogram of log2 ratios with median line; inset: fraction < 0 by family size bin. Filters: as in SQL plus Q-busco-rel, Q-dup, Q-n50. Caveats: reference-dominance bias (2.1); families where yeast and filamentous genomes hit different references are under-represented after the within-sp_id requirement; the coverage gates make the test conservative. Compute: A3 rows for representatives (about 7,450 genomes; estimate 70-80 M rows), done inside job J3; the figure reads A4 (small).

**V09 Per-genome length offset by yeast origin** (P1). Question: is the effect the same in every independent yeast origin? Tables: A3, A4, genome_labels. SQL: for each (genome, sp_id) compute `log2(len / family_ref_len)` where `family_ref_len` is the median across filamentous representative genera; then `SELECT genome, median(offset) FROM ... GROUP BY genome`. Statistic: unit = genus (median of genome offsets); one panel per yeast origin, each compared with its sister filamentous group; nested mixed model over all. Figure: dot plot, x = genus median offset, rows = yeast origin, with sister groups. Filters: same as V08. Caveat: Pezizomycotina dimorphic genera are few (69 genera, labels at genus level). Compute: small after A4.

**V10 Model-organism-anchored contrast (P-B)** (P2). Question: does the result hold when every protein is anchored to the same *S. cerevisiae* (or *N. crassa*) reference? Tables: swissprot_hits (2.2 G rows), swissprot_subject (`org_code IN ('YEAST','NEUCR','SCHPO','EMENI')`), gene_proteins. SQL: as DuckDB view `v_swissprot_model_fungi` (`QUALIFY row_number() OVER (PARTITION BY protein_id, org_code ORDER BY rk)=1`), then same contrast as V08 with family = YEAST accession, and again with family = NEUCR accession to test reference-side bias. Statistic: as V08; compare the two anchors. Figure: paired histogram. Filters: as V08. Caveat: requires job J5 (37 GB scan); hits are already restricted to the top ranks per protein (`rk`), so anchors may be missing for divergent proteins. Compute: J5, 2.2 G rows.

**V11 Within-architecture length contrast, Pfam proxy (P-C)** (P2). Question: same question with domain architecture as the family label. Tables: A3 (`pfam_arch`, `n_domains`, `length`). SQL: as V08 with `pfam_arch` in place of `sp_id`, restricted to architectures with >= 2 domains present in >= 10 yeast and >= 10 filamentous genera. Statistic: as V08. Figure: histogram + scatter of V08 vs V11 family-level ratios for families that map to both (via sp_id majority architecture). Filters: Q-model and the genome filters. Caveat: architecture classes depend on domain e-value gate; truncation moves proteins between classes (conservative). Compute: A3 (job J2 builds architectures).

**V12 Sensitivity of the length contrast to quality filters** (P1). Question: does the effect change with BUSCO threshold, N50 threshold, exclusion of partial models, or the yeast definition? Tables: A4 built under 6-8 filter settings. SQL: parameterised version of V08. Statistic: median log2 ratio and fraction < 0 for each setting. Output: table + line plot (effect vs BUSCO percentile cutoff). Filters: the grid. Caveat: none beyond V08. Compute: repeats the V08 aggregation on A3 (job J3, all settings in one pass with `GROUPING SETS`).

**V13 Which families drive the difference** (P2). Question: is the shortening uniform, or concentrated in some functional classes or length classes? Tables: A4, swissprot_subject (`go_terms`, `pfam`), A3 (`n_domains`). SQL: join family log2 ratio to reference length bins (`s_len` decile), to domain count bins, and to GO slim terms parsed from `go_terms`. Statistic: median ratio per bin, with family counts. Figure: bar/strip per bin. Filters: as V08. Caveat: GO terms are from the reference, not the query. Compute: small after A4.

**V14 Paired sister-clade contrasts** (P1). Question: does each independent yeast clade differ from its filamentous sister in the same direction? Tables: A4 per-genus offsets (V09). SQL: `SELECT pair_id, median(offset) FILTER (yeast side) - median(offset) FILTER (sister side)`. Statistic: per-pair difference with bootstrap CI over genera; sign test across pairs (n = 6 pairs; low power stated). Output: table + forest plot. Filters: as V08. Caveat: sister pairs are chosen from taxonomy, not from a tree (D3). Compute: small.

### H3. Gene structure

**V15 Introns per gene by class and growth form** (P1). Question: how many introns per protein-coding gene? Tables: A2 (`introns_per_gene`, `frac_intronless`). SQL: `SELECT CLASS, growth_form, median(introns_per_gene), median(frac_intronless) FROM A2 JOIN genome_labels WHERE is_representative GROUP BY 1,2`. Statistic: genus medians; nested mixed model. Figure: strip plot by class. Filters: Q-rep, Q-ann, Q-busco-rel. Caveat: ab initio predictors trained per genome; intron under-prediction in genomes without RNA-seq. Flag `training/` presence if the build can read it (not in tables; D8). Compute: gene_introns (508 M rows) in job J1.

**V16 Gene span versus CDS length** (P2). Question: how much of a gene is intron? Tables: A2 (`median_gene_span`, `median_cds_len`, `intron_fraction`). Figure: scatter of median CDS length vs median gene span per representative genome, color = growth form. Statistic: ratio per genome; genus medians. Filters: as V15 + Q-model. Caveat: UTRs are not in the gene span for most ab initio models. Compute: gene_transcripts in J1.

**V17 Intergenic distance by growth form** (P1). Question: are yeast genes more tightly packed? Tables: A2 (`median_intergenic`, `q10_intergenic`), from gene_intergenic_distances. SQL: `SELECT growth_form, CLASS, median(median_intergenic) ... GROUP BY 1,2`. Statistic: genus medians; mixed model with log genome size as covariate. Figure: box by class. Filters: Q-rep, Q-n50 (short contigs truncate distances). Caveat: distances span contig ends only within contigs; fragmented assemblies bias downward; `gene_intergenic_distances` includes tRNA genes if gene_info does (check in J1). Compute: 233 M rows in J1.

**V18 Gene density versus genome size** (P1). Question: genes per Mb as a function of assembly size and repeat content. Tables: A2, asm_stats (`total_length_bp`, `masked_pct`). Figure: scatter (log genome size vs genes per Mb), color = growth form; second panel vs masked_pct. Statistic: Spearman per growth form; mixed model. Filters: Q-rep, Q-busco-rel. Caveat: masking method differs by wave (not verified). Compute: small after A2.

**V19 Intron length and splice-site classes** (P3). Question: intron length distribution and GT-AG fraction by class. Tables: A2 (`median_intron_len`, `frac_gt_ag`, `frac_gc_ag`), from gene_introns `splice_5`/`splice_3`. Figure: strip by class. Filters: Q-rep. Caveat: non-canonical splice sites may be model errors. Compute: J1.

### H4. Secretome, CAZy and peptidases versus lifestyle

**V20 Secretome size per genome** (P1). Question: how many predicted secreted proteins, by growth form? Tables: A5 (`n_signalp`, `n_secreted_strict` = SignalP prob >= 0.9 AND tmhmm PredHel <= 1 AND no GPI-anchor, `n_proteins`). SQL: `SELECT growth_form, CLASS, median(n_secreted_strict), median(n_secreted_strict*1.0/n_proteins) FROM A5 JOIN genome_labels WHERE is_representative GROUP BY 1,2`. Statistic: genus medians; mixed model with log gene count. Figure: box by class (count and fraction panels). Filters: Q-rep, Q-busco-rel, Q-model for the denominator. Caveat: SignalP cannot see truncated N-termini (partial models reduce counts; report partial_fraction as a covariate). tmhmm covers an unmeasured genome count; where tmhmm is missing the strict count is NULL, not 0. Compute: J4.

**V21 Secretome fraction by primary lifestyle** (P1). Question: do plant pathogens and saprotrophs carry larger secretomes than yeasts and animal parasites? Tables: A5, genome_labels (`primary_lifestyle`). SQL: genus medians of `n_secreted_strict/n_proteins` grouped by lifestyle (lifestyles with >= 10 genera). Statistic: unit = genus; Kruskal-Wallis as screen; pairwise Dunn with BH; mixed model with CLASS random effect to separate lifestyle from clade. Figure: ordered box plot by lifestyle, colored by dominant phylum. Filters: Q-rep, Q-busco-rel. Caveat: lifestyle is a genus trait; `nectar/tap_saprotroph` is dominated by *Saccharomyces* and is confounded with growth form; `plant_pathogen` mixes biotrophs and necrotrophs. Compute: small after A5.

**V22 CAZy class profile by lifestyle** (P1). Question: which CAZy classes (GH, GT, PL, CE, AA, CBM) expand in which lifestyle? Tables: A5 (`cazy_class_counts` map), cazy.overview (`cazyme_fam` parsed to class prefix; count rows with `toolcount >= 2`). SQL: `SELECT species_prefix, regexp_extract(cazyme_fam,'^([A-Z]+)',1) AS cls, count(DISTINCT protein_id) FROM cazy_overview WHERE toolcount>=2 GROUP BY 1,2`. Statistic: genus medians per class; log-ratio versus the filamentous saprotroph median; mixed model per class. Figure: heatmap (lifestyle x CAZy class, cell = median genus count per 1,000 genes). Filters: Q-rep. Caveat: `cazyme_fam` may hold several families per protein ("GH27(12-262)+..."); split on `+` and count a protein once per class. toolcount = 1 rows (868 k) are single-tool calls; the default drops them (D7). Compute: cazy.overview 122 MB (login-node scale; done in J4 for consistency).

**V23 Plant-cell-wall-degrading enzymes versus plant association** (P1). Question: do plant pathogens, endophytes and wood/litter saprotrophs carry more PCW-active families than yeasts and animal parasites? Tables: cazy.overview (`substrate` column: cellulose, xylan, pectin, lignin, ...), A5 (`cazy_substrate_counts`). SQL: `SELECT species_prefix, s, count(DISTINCT protein_id) FROM (SELECT species_prefix, protein_id, unnest(string_split(substrate,'+')) AS s FROM cazy_overview WHERE toolcount>=2) GROUP BY 1,2`. Statistic: genus medians by `Plant_pathogenic_capacity_template` and `Decay_substrate_template` (FungalTraits) and by growth form. Figure: dot plot per substrate group. Filters: Q-rep. Caveat: substrate annotation comes from dbCAN's family-level substrate map, not from the protein. Compute: J4.

**V24 Secreted CAZymes** (P2). Question: what fraction of CAZymes carry a signal peptide, by lifestyle? Tables: cazy.overview joined to signalp and tmhmm on protein_id (A5 `n_cazy_secreted`). Statistic: genus medians. Figure: box by lifestyle. Filters: Q-rep, Q-model. Caveat: as V20. Compute: J4.

**V25 Peptidase (MEROPS) families by lifestyle** (P2). Question: do animal parasites and plant pathogens differ in secreted peptidase family content? Tables: merops (best hit per protein by bitscore; family = first 3 characters of the MEROPS identifier after mapping `merops_id` to a family, which needs the MEROPS id-to-family table; D7), signalp. SQL: `SELECT species_prefix, protein_id, arg_max(merops_id, bitscore) FROM merops WHERE evalue<=1e-10 GROUP BY 1,2`. Statistic: genus medians of family counts per 1,000 genes, by lifestyle. Figure: heatmap (lifestyle x top-15 families). Filters: Q-rep. Caveat: `merops_id` (MER0003338) is a sequence id, not a family; the family mapping file is not in `tables/` (not found; D7). Compute: merops 1.4 GB in J4.

**V26 Subcellular localization profile by class** (P3). Question: does the WoLF PSORT composition differ by growth form? Tables: wolfpsort (A5 `wolfpsort_counts`). Statistic: compositional (CLR) PCA; PERMANOVA at genus level. Figure: stacked bars per class. Filters: Q-rep. Caveat: WoLF PSORT fungal model is coarse. Compute: 245 M rows in J4.

### H5. Domain and transcription-factor repertoires

**V27 TF family counts versus gene count** (P1). Question: does each TF family scale sub-linearly with proteome size? Tables: tf_inventory (count DISTINCT protein_id per `tf_family` per genome; A5 `tf_counts`), A2 `n_proteins`. SQL: `SELECT species_prefix, tf_family, count(DISTINCT protein_id) FROM tf_inventory GROUP BY 1,2`. Statistic: log-log OLS slope per family on representative genomes, with 95 % CI; mixed model with CLASS random intercept. Figure: small multiples (one per family, 21 panels), slope annotated. Filters: Q-rep, Q-busco-rel. Caveat: a protein with two TF domains is counted once per family; family definitions come from the pipeline's Pfam list. Compute: 167 MB.

**V28 Zn2Cys6 share of the TF repertoire by class** (P1). Question: is Zn2Cys6 expansion restricted to Pezizomycotina? Tables: A5 `tf_counts`. Statistic: share per genome; genus medians by class. Figure: box by class, ordered by phylum. Filters: Q-rep. Caveat: `Zn(II)2Cys6 zinc cluster` and `Fungal-specific transcription factor domain` often co-occur on one protein; count proteins, not domains. Compute: small.

**V29 Pfam domain richness versus gene count** (P2). Question: how many distinct Pfam families per genome, and does it saturate? Tables: A7 `pfam_genome_matrix` (genome, pfam_acc, n_proteins). SQL: `SELECT genome, count(*) AS n_families, sum(n_proteins) FROM A7 GROUP BY 1`. Statistic: log-log slope; genus medians by growth form. Figure: scatter. Filters: Q-rep, Q-busco-rel. Caveat: domain_i_evalue gate (1e-5) fixed in A7. Compute: pfam 19 GB in J2.

**V30 Clade-enriched Pfam families, yeast versus filamentous** (P2). Question: which domain families are present in most yeast genera but few filamentous genera, and vice versa? Tables: A7, genome_labels. SQL: presence per genus (present in >= 50 % of representative genomes of the genus), then Fisher exact per family on genus counts, BH-corrected. Statistic: unit = genus. Output: two tables (top 30 each direction) + volcano plot. Filters: Q-rep. Caveat: shared ancestry inflates significance (genera in one order are not independent); report the number of orders in which each family is present. Compute: small after A7.

**V31 Domains per protein and multi-domain fraction** (P2). Question: are yeast proteins shorter because they carry fewer domains? Tables: A3 (`n_domains`, `length`). Statistic: genus medians of multi-domain fraction; within-family (V08) ratio split by domain count bin. Figure: box + binned ratio. Filters: as V08. Caveat: links to V13. Compute: A3.

**V32 Intrinsic disorder fraction by class** (P3). Question: does the disordered fraction differ by growth form? Tables: idp_summary (`IDP_fraction`, keyed by gene_id form), A5 `idp_median`. Statistic: genus medians; mixed model. Figure: box by class. Filters: Q-rep, Q-model. Caveat: key form must be verified (section 1); coverage not measured. Compute: 145 M rows in J4.

### H6. Codon usage and amino-acid composition versus taxonomy

**V33 Codon usage ordination** (P1). Question: do genomes cluster by class or by growth form in RSCU space? Tables: codon_freq (64 codons per genome), samples (`TRANSL_TABLE`), genome_labels. SQL: RSCU = `frequency / avg(frequency) OVER (PARTITION BY species_prefix, amino_acid)` with a codon-to-amino-acid map (code 1; code-12 genomes have CTG = Ser, handled by a second map). Statistic: PCA on genus-median RSCU (1,671 genera); PERMANOVA (`skbio` or `vegan`-style permutation) with CLASS, then growth form | CLASS, then lifestyle | CLASS; report R2 for each. Figure: PCA biplot (PC1-2) colored by class; second colored by growth form. Filters: Q-rep. Caveat: codon frequencies are whole-proteome values, dominated by low-expression genes; this is not the "optimal codon" signal. Compute: 5 MB, login-node scale (done in J4 for provenance).

**V34 GC3 versus genomic GC** (P1). Question: does third-position GC follow genomic GC, and do yeasts deviate? Tables: codon_freq (GC3 = sum of frequencies of codons ending in G or C), asm_stats `gc_pct`. Figure: scatter, color = growth form, with class centroids. Statistic: slope per growth form; residuals by class. Filters: Q-rep. Caveat: assembly GC includes repeats; `masked_pct` as covariate. Compute: small.

**V35 Codon-usage distance versus taxonomic rank** (P1). Question: does codon usage track phylogeny (distance grows with rank) more than growth form? Tables: genus-median RSCU (V33). Statistic: Euclidean distance between genus pairs binned by lowest shared rank (family, order, class, phylum); compare with distance between same-growth-form and different-growth-form pairs within one rank. Figure: box plot of distance by rank, split by same/different growth form. Filters: Q-rep. Caveat: taxonomy rank is not branch length; pairs are not independent (report medians and resampling CI, no p-value). Compute: 1,671 x 1,671 distances (small).

**V36 tRNA gene copy number versus codon usage** (P2). Question: does the tRNA pool match codon frequency within a genome? Tables: gene_trnas (count per `codon` per genome prefix), codon_freq. SQL: `SELECT split_part(gene_id,'_',1) AS g, codon, count(*) FROM gene_trnas GROUP BY 1,2` joined to codon_freq on reverse-complemented anticodon (verify whether `codon` in gene_trnas is the codon or the anticodon; the sample shows 'AGC' for Ser and 'CGT' for Arg, which are codons, but confirm against tRNAscan output in the build). Statistic: Spearman rho per genome; genus medians by class. Figure: box of rho by class + one example scatter. Filters: Q-rep, Q-n50 (tRNA genes on short contigs may be missed). Caveat: tRNA counts in fragmented assemblies are undercounted. Compute: 35 MB.

**V37 CTG-clade codon usage** (P2). Question: is CTG frequency depleted in code-12 genomes relative to code-1 relatives? Tables: codon_freq, samples `TRANSL_TABLE` (12: 1,154 genomes, all Pichiomycetes). Figure: strip of CTG frequency by order within Pichiomycetes, code-1 vs code-12. Statistic: Mann-Whitney on genus medians within Pichiomycetes. Filters: Q-rep. Caveat: table assignment comes from the samples sheet; the *Pachysolen* (code 26) genome is a single point. Compute: small.

**V38 Amino-acid composition ordination** (P2). Question: does proteome composition follow GC and taxonomy? Tables: aa_freq (20 aa; drop `*` and `X`), asm_stats gc_pct. Statistic: PCA on genus medians; correlation of PC1 with GC and with growth form. Figure: PCA colored by GC, inset by class. Filters: Q-rep. Caveat: `*` present in 6,254 genomes suggests internal stops or terminal stop counted; exclude the symbol, state it. Compute: small.

### H7. Paralog retention (7,384 genomes; interim table)

**V39 Young paralog nodes per 1,000 genes by growth form** (P2). Question: do yeasts retain fewer recent duplicates? Tables: wgd.ks.nodes (`kind='tree'`, `dS_mean` in 0.05-1.0, `n_valid>=1`), wgd.ks.summary, A2 `n_proteins`. SQL: `SELECT species_prefix, count(*) FILTER (WHERE kind='tree' AND dS_mean BETWEEN 0.05 AND 1.0) * 1000.0 / n_proteins FROM wgd_ks_nodes JOIN A2 GROUP BY 1`. Statistic: genus medians; mixed model. Figure: box by class. Filters: Q-rep, Q-dup (strict: duplicated_pct <= 5), exclude `hybrid = 1`. Caveat: per `fable_context.md`, dS < 0.05 is heterozygosity/haplotigs and dS > 1.5 is saturated; the window 0.05-1.0 is a choice (D9). Coverage is 7,384 genomes, selection rule for that subset is not documented in the manifest. Compute: 29 M rows (J4).

**V40 Identical-copy groups versus BUSCO duplication** (P2). Question: are large identical-CDS groups assembly artefacts? Tables: wgd.ks.summary (`n_genes_in_identical_groups`, `largest_identical_group`), busco_genome `duplicated_pct`, A2 `n_proteins`. Figure: scatter (genes in identical groups per 1,000 genes vs duplicated_pct), color = phylum; label outliers. Statistic: Spearman; outlier list. Filters: none (artefact detection view). Caveat: TE-derived models also form identical groups (fable_context, Q1 rank 2). Compute: small.

### H8. Telomere repeats

**V41 Telomere monomer by class** (P3). Question: which repeat unit dominates per class? Tables: telomere_summary (`telomere_top_monomer`, `telomere_tracts`), species. SQL: `SELECT CLASS, telomere_top_monomer, count(*) FROM telomere_summary JOIN species USING (ASMID) WHERE telomere_tracts>0 GROUP BY 1,2`. Output: table (top monomer per class with fraction) + bar. Filters: genomes with >= 2 tracts. Caveat: only 8,369 genomes have any tract; the detector's monomer alphabet is unknown here (check `do_telomeres/`). Compute: small.

**V42 Telomere-capped scaffolds versus contiguity** (P3). Question: does the fraction of scaffolds with both-end telomeres track N50? Tables: telomere_summary (`telomere_scaffolds_both_ends`), asm_stats (`contig_count`, `N50_bp`, `t2t_scaffolds`). Figure: scatter, color = phylum. Statistic: Spearman. Filters: none. Caveat: `t2t_scaffolds` in asm_stats and `telomere_scaffolds_both_ends` may come from different detectors; reconcile first. Compute: small.

---

## 6. Build plan

### 6.1 Intermediate tables (built once on SLURM; every figure reads these)

Output directory (decision D1; default proposed): `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/atlas/`. All parquet written with `COMPRESSION ZSTD`, ordered by genome. Each file gets a `_manifest/<name>.json` with row count, source table `built_at` values, filters and the SQL file hash. Large TSV side outputs are `zstd -T0` compressed.

| ID | File | Grain | Sources | Content |
|---|---|---|---|---|
| A1 | genome_labels.parquet | genome | species, samples, asm_stats, busco_genome, FungalTraits 1.2 (UTF-8 copy), `atlas_overrides/growth_form_overrides.tsv` (hand-made), A2 | section 4.1 columns, `is_representative`, quality flags |
| A2 | genome_gene_stats.parquet | genome | gene_info, gene_transcripts, gene_proteins, gene_introns, gene_intergenic_distances | n_genes, n_protein_coding, n_trna, n_transcripts, n_partial, partial_fraction, frac_no_start, frac_no_stop, frac_start_null, protein length quantiles (10/25/50/75/90, mean), median CDS length, median gene span, introns per gene, frac_intronless, median intron length, frac_gt_ag, median and q10 intergenic distance |
| A3 | protein_core.parquet | protein | gene_proteins, gene_transcripts, pfam (A3a), swissprot_best | protein_id, genome, gene_id, length, is_partial, has_start, has_stop, n_domains, pfam_arch, sp_id, pident_x10, qcov_x10, scov_x10, q_len, s_len, org_code |
| A3a | protein_pfam_arch.parquet | protein | pfam (domain_i_evalue <= 1e-5, ordered by env_from) | protein_id, n_domains, pfam_arch (string of pfam_acc without version joined by `;`), has_arch |
| A4 | family_length.parquet | (sp_id or pfam_arch, filter_set) | A3, A1 | per family: n genomes/genera per growth form, median length per growth form, log2 ratio; built for the V12 filter grid; plus `genome_family_offset.parquet` (genome, sp_id, log2 offset) for V09/V14 |
| A5 | genome_function_counts.parquet | genome | signalp, tmhmm, predgpi, wolfpsort, cazy.overview, merops, tf_inventory, idp_summary, codon_freq, aa_freq, gene_trnas, wgd.ks.* | n_signalp (prob >= 0.5 and >= 0.9), n_secreted_strict, n_tm (PredHel >= 2), n_gpi, wolfpsort counts (map), cazy class counts, cazy substrate counts, n_cazy_secreted, merops best-hit family counts (map; needs id-to-family file, D7), tf_counts (map), idp median fraction, RSCU (64 values), GC3, tRNA codon counts (map), young paralog nodes, identical-group stats |
| A7 | pfam_genome_matrix.parquet | (genome, pfam_acc) | pfam | n_proteins, n_domains |
| A9 | protein_model_anchor.parquet | protein | swissprot_hits, swissprot_subject | best YEAST, NEUCR, SCHPO, EMENI hit per protein (sp_id, pident, qcov, scov, s_len) — P2 only |
| coverage.tsv | table | all | count(distinct genome) per table (V01) |

### 6.2 SLURM jobs (sizes are estimates, not measurements)

Each job: `module load` the DuckDB-capable Python (`/usr/bin/python3.12` or a pixi environment; D10), `SET temp_directory='${SCRATCH:?}/duckdb_tmp'`, `SET memory_limit` to 80 % of the allocation, `SET threads` = CPUs, `preserve_insertion_order=false`. Scripts take the project root from an environment variable (never `$(dirname "${BASH_SOURCE[0]}")`). Outputs go to `/bigdata`. Each job verifies its inputs' `_manifest` SHA before running and writes its own manifest after.

| Job | Builds | Input scanned | Request | Estimated wall time | Order |
|---|---|---|---|---|---|
| J0 coverage | coverage.tsv | every table once (134 GB) | 16 CPU, 64 GB, 4 h | 1.5-3 h (I/O bound; estimate) | 1 (can run with J1) |
| J1 gene structure | A2 | gene_info 3.8 GB, gene_transcripts 6.5 GB, gene_proteins 5.2 GB, gene_introns 6.5 GB, gene_intergenic 2.9 GB | 16 CPU, 96 GB, 4 h | 1-2 h (estimate) | 1 |
| J2 pfam | A3a, A7 | pfam 19 GB (368 M rows); `string_agg` ordered per protein | 16 CPU, 128 GB, 4 h | 1-2 h (estimate; string_agg over 368 M rows is the risk; if memory is short, process by `protein_id` prefix hex ranges in 4 chunks inside one job) | 1 |
| J3 protein core + families | A1, A3, A4 | A3a, gene_proteins, gene_transcripts, swissprot_best 4.9 GB; then the V08/V12 aggregations | 16 CPU, 128 GB, 4 h | 1-2 h (estimate) | 2 (after J1, J2) |
| J4 function counts | A5 | signalp 0.3 GB, tmhmm 1.8 GB, predgpi 1.9 GB, wolfpsort 1.8 GB, cazy 0.1 GB, merops 1.4 GB, tf_inventory 0.2 GB, idp_summary 2.1 GB, wgd.ks.nodes 0.6 GB, small tables | 16 CPU, 64 GB, 4 h | 1-1.5 h (estimate) | 1 |
| J5 model anchor (P2) | A9 | swissprot_hits 37 GB, 2.2 G rows | 32 CPU, 192 GB, 8 h | 2-4 h (estimate) | 3, only if D4 approves |
| J6 render | all figures/tables | A1-A7 (small) | 4 CPU, 16 GB, 2 h | < 1 h | 4; can also run on the login node per page |

Five jobs cover P1. This follows the 1-1.5 h sizing rule as far as the table sizes allow; J0 and J5 are single I/O-bound passes and are not split further because job count, not CPU, is the cost.

### 6.3 Build order

1. J0, J1, J2, J4 in parallel (independent).
2. J3 (needs J1, J2).
3. Render P1 views: V01-V05, V07-V09, V12, V14, V15, V17, V18, V20-V23, V27, V28, V33-V35.
4. PI review of the P1 pages, decisions D1-D10.
5. J5 if approved; render P2 views (V06, V10, V11, V13, V16, V24, V25, V29-V31, V36-V40), then P3 (V19, V26, V32, V41, V42).

### 6.4 Code layout (to be created by the build step; nothing created now)

```
BFD_dashboard/
  bin/atlas/
    build_labels.py          # A1 (login-node scale) + override file reader
    sql/J1_gene_stats.sql ... # one SQL file per intermediate, run by run_duckdb_job.py
    run_duckdb_job.py        # generic: reads SQL file, sets SCRATCH temp, writes manifest
    render_views.py          # reads A1-A7, writes figures/ and tables/, one function per view ID
    stats.py                 # mixed model, bootstrap, sign test, PERMANOVA helpers
  workflows/atlas/J0_coverage.sh ... J6_render.sh   # sbatch scripts
  docs/summary_atlas/
    PLAN.md                  # this file
    index.md
    H1_quality.md ... H8_telomeres.md
    figures/<VID>_<slug>.png and .svg
    tables/<VID>_<slug>.tsv  (small tables only; large outputs stay in tables/atlas/)
```

### 6.5 Static Markdown layout

`docs/summary_atlas/index.md`: title, one paragraph on the data snapshot (`built_at` from manifests, genome counts), the coverage table (V01), the hypothesis table (section 3) with links to `H1_quality.md` ... `H8_telomeres.md`, and a "how to read the statistics" box (units of replication, filters, what Kruskal-Wallis means here).

Each `H<n>_<slug>.md`: statement and expected direction; filters applied with the n they leave; one section per view with the figure (`![](figures/V08_family_ratio.png)`), the summary table (inline if <= 20 rows, else link to `tables/`), the statistic with effect size and CI, and the caveats. A final "verdict" section states supported / not supported / not testable with this data, in one sentence each, and lists the next data needed.

Figures: matplotlib, PNG at 200 dpi and SVG, colorblind-safe palette from `lib/figures.py` (already in the repo). Every figure title carries the view ID. Every table has a header comment with filters and n.

---

## 7. Decisions (PI answers received 2026-11-08)

| ID | Decision | Outcome |
|---|---|---|
| D1 | Intermediate location | **Decided by PI**: `/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/` as zstd parquet. Only figures, small tables and Markdown go into git (PI commits). |
| D2 | Yeast definition | **Decided by PI**: FungalTraits 1.2 growth form + override file `atlas_data/labels/growth_form_overrides.tsv` (about 30 genera). `dimorphic_yeast` is a separate third group, never pooled. |
| D3 | Sister-clade pairs for V14 | decided: default (the six pairs in 4.3.4) |
| D4 | SwissProt proxy gates; J5; BUSCO re-run | **Decided by PI**: P1 first. Gates as stated (qcov, scov >= 80 %, pident >= 40 %). J5 (swissprot_hits, 37 GB) runs only after P1 review. No BUSCO re-run. |
| D5 | Phylogeny | **Decided by PI**: build a tree on one representative genome per genus (about 1,671 tips) with nf_phyling + FastTree, for PGLS on genus-level summaries. No 20k-tip tree. Representative rule = D6. |
| D6 | Representative rule | decided: default (highest `complete_pct` within species, tie `N50_bp`; the genus representative is the species representative with the highest `complete_pct`, tie `N50_bp`) |
| D7 | CAZy toolcount; MEROPS family map | decided: default (`toolcount >= 2`; V25 deferred until a MEROPS id-to-family file is found) |
| D8 | RNA-seq training flag | decided: default (not in the first build) |
| D9 | Young-paralog dS window; interim wgd table | decided: default (0.05-1.0; interim table used and marked) |
| D10 | Python environment | decided: DuckDB jobs J0-J4 run on `/usr/bin/python3.12` (duckdb 1.5.5). The render job uses `BFD_dashboard/atlas_env/venv/` (a `--system-site-packages` venv on `/usr/bin/python3.12` with statsmodels 0.15.0 added), because the pixi environment in `atlas_env/pixi.toml` (`[workspace]` table) did not finish installing on the network file system within the build window (`atlas_data/logs/pixi_install.log`). |
| D11 | View trimming | decided: default (keep all; P1 built first) |
| D12 | FungalTraits citation | decided: default (cite FungalTraits 1.2 with a "verify citation" note until checked) |

Operational constraint added 2026-11-08: SLURM maintenance reservation `20261110-Shutdown` starts 2026-11-10 08:00 PST. Every job carries a `--time` that ends before that; intermediates are written atomically (tmp then rename) and skipped when present, so reruns resume; nf_phyling runs with `-resume` and a work dir on `/bigdata`.

## 8. Known limits of this plan

- Genome coverage of the large tables is not measured; J0 measures it.
- No runtime of any full-table job is measured. Estimates in 6.2 come from file sizes and a 16-genome pushdown pilot only.
- The length-hypothesis pilot (21 genomes) does not show a direction and is not evidence either way.
- `SUBPHYLUM` is empty in `species.parquet`; subphylum groupings in this plan are derived from CLASS lists written by hand and must be checked against the taxonomy source.
- `idp_summary` key form and `gene_trnas.codon` semantics (codon vs anticodon) are sample-based observations; the build must assert them.

## 9. Build log (P1, 2026-11-08)

All jobs ran on partition `epyc` (16 CPU) except the render (`short`, 4 CPU) and the phyling driver (`stajichlab`). Login-node work is limited by a 1 GB memory cgroup; every scan and every render runs under SLURM.

| Job | SLURM ID | State | Elapsed | Output |
|---|---|---|---|---|
| labels (A1) | login node | done | < 1 min | `atlas_data/labels/genome_labels.parquet` (23,680 rows; 7,435 species reps; 1,669 genus reps) |
| J0 coverage | 31870349 | COMPLETED | 1 min 51 s | `atlas_data/coverage/*.parquet` (16 tables) |
| J1 gene stats | 31870350 | COMPLETED | 1 min 14 s | `atlas_data/A2/*.parquet` (23,452 genomes) |
| J2 pfam | 31870351 | FAILED (OOM in `string_agg` over 368 M rows at 104 GB) | 2 min 22 s | `A7/pfam_genome_matrix.parquet` only |
| J2b pfam, chunked by protein_id prefix | 31870461 | COMPLETED | 6 min 59 s | `A3/protein_pfam_arch.parquet` (173,294,381 proteins with >= 1 domain), `A5/tmhmm_predhel_check.parquet` |
| J3 protein core | 31870462 | COMPLETED | 2 min 44 s | `A3/protein_core.parquet` (244,988,294 rows), `A4/fam_sp_genome.parquet` (14,384,792), `A4/fam_arch_genome.parquet` (16,032,538), `A4/genome_protein_summary.parquet` |
| J4 function counts | 31870352 | COMPLETED | 21 s | `atlas_data/A5/*.parquet` |
| J6 render (no mixed models) | 31870463, 31870579 | COMPLETED | 25 s, 2 min 52 s | all 22 P1 views |
| J6 render (with mixed models, venv) | 31870629 (V21 failed on a formula-term bug), 31870788 (final, 22/22 views) | COMPLETED | 2 min 47 s | `docs/summary_atlas/*.md`, 46 figure files, 23 tables |
| phyling driver | 31870356 | RUNNING (align task 31870362) | | `atlas_data/phylo/results/` |

The wall-time estimates in section 6.2 were too high by one to two orders of magnitude: DuckDB reads only the needed columns, and the tables are genome-ordered.

Measured in the build: `idp_summary.protein_id` has no `-T` suffix (0 of 144,924,594 rows), so it joins on `gene_id`. 14,577,766 of 18,767,284 SignalP proteins have no TMHMM row; the TMHMM table holds only proteins with >= 1 predicted helix (0 of 45,707,260 rows have PredHel = 0; minimum 1), so a missing row means 0 helices and the secretome counts treat it so.
