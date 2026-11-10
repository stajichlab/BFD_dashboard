# Functional and comparative analyses: plans and findings

Started: 2026-11-10. This page tracks the exploratory analyses run against the full BFD database in MotherDuck, and the plan to turn the useful ones into HPCC jobs. Results are added as each analysis finishes.

## Data access and conventions

- **Source:** MotherDuck database `bfd` (36 tables, all of `Fungi_BFD_runs/tables` except `swissprot_hits`). Loaded by `bin/md_load_bfd.py` / `workflows/md_load_bfd.sh`. Every table has a comment giving its genome key and caveats.
- **Scratch:** intermediate tables go in the MotherDuck database `bfd_scratch`, never in `bfd`. Prefix by analysis:

  | Prefix | Analysis |
  |---|---|
  | `kw_` | Kingdom-wide domain survey (B) |
  | `ao_` | *Antarctolithica obscura* characterization (C) |
  | `rf_` | Functional classifier (D) |
  | none | Kickxellomycota and research agenda (A) |

- **Statistics:** follow the atlas conventions in [summary_atlas/index.md](summary_atlas/index.md):
  - species representatives (7,435) and genus representatives (1,669);
  - genus medians as the replication unit for cross-clade tests;
  - effect sizes first;
  - quality filters Q-busco-rel, Q-dup, Q-n50, Q-model.
- **Known data caveats** (also in the table comments):

  | Item | Caveat |
  |---|---|
  | `tmhmm` | Lists only proteins with at least one TM helix (about 20% of proteins) |
  | `targetP` | 125 genomes; may be replaced by DeepLoc |
  | `infernal_rfam` | 14 genomes |
  | `pfam` | 23,410 genomes; 42 genomes with proteins have no Pfam rows |
  | BUSCO (`fungi_odb12`) | About 75–80% in budding yeasts, from real gene loss |
  | `SUBPHYLUM` | Mostly empty; use CLASS and ORDER |
  | `gene_introns.codon_position`, `codon_frame` | **Bug:** off by one exon. They describe the next CDS junction, and the last intron gets the protein length and frame 0. Recompute from cumulative `gene_CDS` lengths. Confirmed 2026-11-10; the upstream table builder needs fixing |

- **Reporting standard (every analysis):**
  - Include figures (PNG, embedded in the report) and tables.
  - Explain every analytical choice: thresholds, normalization, distance metrics, tests, representative selection, CV design, hyperparameters.
  - Wherever a statistical or ML model or an established method is used, justify it and cite the literature. Never invent a citation; mark uncertain ones "citation to verify".
  - Keep the scripts that produce each figure.
  - Detailed reports and figures go in `docs/analyses/<analysis>/`; this page keeps the summary.
- **End state:** the analyses that prove useful become restartable SLURM jobs in `workflows/` that read the parquet tables on HPCC directly, in the same way as the atlas jobs J0–J10. MotherDuck is used for exploration and for MCP access, not as the production compute.

## A. Review of current ideas and Kickxellomycota deep dive

**Status:** running (Fable agent).

**Questions**
- How do Kickxellomycota reduce genome size and intron content while staying mostly filamentous?
- Is the reduction concentrated in particular orders or lifestyles: Harpellales (gut symbionts), Dimargaritales (mycoparasites), or saprotrophs?
- Is intron reduction a matter of fewer introns or shorter introns?
- Which gene families are lost? Candidates: spliceosome and NMD, RNAi, DNA methylation, transposable-element domains, DNA repair.
- Are these losses convergent with yeast reductions?

The agent also reviews the current dashboard and atlas ideas and proposes a ranked research agenda: gene-family gain, loss and expansion; candidate horizontal transfers; how to define lineage-restricted Pfam domains robustly; and what the dashboard and MCP path need.

**Status:** done 2026-11-10. Full report with 7 figures, tables, methods and references: [analyses/kickxellomycota/](analyses/kickxellomycota/README.md).

**Scope:** 169 Kickxellomycota species representatives (15 genera, 4 orders; Asellariales not sampled). Genera are treated as independent, because no tree is available yet.

**Findings**
1. **Reduction is concentrated in core Kickxellales** (*Coemansia*, *Kickxella*, *Linderina* and relatives). Genus medians: 18.0 Mb, 7,363 genes and 0.65 introns per transcript, against 34.2 Mb, 11,460 genes and 3.83 for Mucorales (Cliff's δ = −1.00 for all three). Harpellales, the gut symbionts, are not reduced (27.6 Mb, 2.06 introns per transcript, AT- and repeat-rich), and neither are Dimargaritales. So the free-living saprotrophs are the reduced group, not the host-associated ones.
2. **The intron deficit is real.** It persists in 30 *Coemansia* assemblies with N50 ≥ 50 kb (0.43 introns per transcript). Within SwissProt-anchored ortholog families shared with Mucorales, *Coemansia* keeps 8.5% of introns (0.43 vs 5.05 per gene; 2,317 families). That is the level of Dipodascomycetes (7.2%) and close to Saccharomycetes (3.5%). CDS length is unchanged.
3. **Intron number fell, not intron length.** Introns per kb of CDS fell 8-fold. The surviving introns are longer, and they are biased toward the 5′ end of genes (39% of sole introns are in the first CDS decile, against 14% in Mucorales), which is the signature of loss by reverse transcription. Zoopagales reached small genomes the other way: they kept intron number and shortened introns.
4. **Gene content: losses unlike the yeasts.** Ku70/80 is absent from all 134 *Coemansia* genomes and from Harpellales; canonical RdRP, histone H1, chitin deacetylase and photolyase are also lost. Unlike budding yeasts, they keep the spliceosome, NMD (Upf1–3, Mago), RNAi (Piwi, Dicer) and DNMT. DNA-transposon domains are depleted. Only 4 of 41 core losses overlap the budding-yeast losses: the intron phenotype converged with yeasts, but gene content mostly did not. Ku loss is a candidate permissive factor for intron loss, but Harpellales (no Ku, 2 introns per transcript) argue against it being sufficient.
5. **Atlas review.**
   - *Missing:* a phylogeny (genus-level p-values are anticonservative), orthogroups, an annotation-gap control for presence/absence calls, and a repeat-class table.
   - *Hidden by pooling:* the pooled "Zoopagomycota s.l." row averages over very different groups.
   - *Over-filtering:* Q-n50 removes 79% of Kickxellomycota.
   - *Bug found:* the `gene_introns` phase bug above.

**Top next studies**
1. Orthogroups plus a species tree, with CAFE, Dollo and Count gain-loss analysis. Load the orthogroups into the release.
2. A robust definition of lineage-restricted domains: model false absences with BUSCO, N50 and partial-gene covariates, and check absences in six frames.
3. Fix and publish shared tables: the `gene_introns` phase bug, `genome_labels`, a `pfam_genome_matrix`, and a taxon synonym table.
4. Ancestral reconstruction of intron positions and Ku/RdRP presence in Kickxellomycota, and a kingdom-wide test of Ku loss against intron density.
5. RNA-seq validation of intron models, and a repeat-class table.

All absences need tblastn confirmation.

## B. Kingdom-wide survey of domain and functional-class change

**Status:** running (Fable agent).

**Approach:** a genome × Pfam count matrix over species representatives (`bfd_scratch.kw_*`), used to find:
- lineage-restricted domains;
- loss modules (domains lost together independently in several lineages);
- family expansions measured against proteome size, and whether each tracks lifestyle or taxonomy;
- shifts in CAZy, MEROPS, TF, secretome and disorder classes.

From these it distills 6–10 patterns worth deeper study.

**Status:** done 2026-11-10. Full report with 9 figures and cited methods: [analyses/kingdom_domains/](analyses/kingdom_domains/README.md).

**Scope:** 7,411 representative genomes, collapsed to 1,173 quality-filtered genera, across 15,052 Pfam domains plus CAZy, MEROPS and TF sets. Nothing is phylogenetically corrected yet (the genus tree is pending), so class and order fixed effects and counts of independent lineages stand in for it. Loss calls in lineages with BUSCO below 80% are provisional.

**Findings**
1. **Taxonomy dominates lifestyle.** After normalizing to proteome size, class explains a median 44% of between-genus variance per Pfam family; lifestyle within order explains 2.4%. CAZy families are the exception (median 5.4%; GH10 19%).
2. **Dikarya-stem loss.** 203 ancestral domains are in ≤5% of Dikarya genera, and 191 of them are retained in ≥50% of Mucoromycota genera. Examples: B12-dependent enzymes, Integrator, the Fanconi core, KICSTOR, SMG1, DNA-PKcs.
3. **Convergent reductive loss.** Fraction of 3,568 core domains absent: Microsporidia 62%, *Hanseniaspora* 28%, *Pneumocystis* 23%. Loss sets overlap 3–11× more than chance. 58 co-absence modules, including complex I, β-oxidation plus glyoxylate cycle, molybdenum cofactor (≥9 lineages), urease, RNAi and sulfate assimilation.
4. **Fatty-acid synthase is absent** in all *Malassezia*, *Pneumocystis*, Glomeromycota and Microsporidia species, against 93–100% elsewhere.
5. **Pezizomycotina expansion** of HET, NACHT, ankyrin, CFEM and rhodopsin-like 7TM families. These follow taxonomy (class explains 49–91% of variance) rather than lifestyle (≤8%).
6. **Ectomycorrhizal contraction** extends beyond plant-cell-wall CAZymes to SDR, GST, aldo-keto reductase and methyltransferase families (20–30% fewer within orders).
7. **Glomeromycota:** median 2,020 kinase proteins per genome, plus Sel1 and TLDc expansions.
8. **Lineage-restricted domains:** 1,041 in total, including large unknown families (Basidiomycota DUF6534; Pezizomycotina DUF8004 and DUF8035) and WRKY in Mucoromycota. *S. cerevisiae*-biased Pfam models produce artefactual "Saccharomycetales-restricted" domains.
9. **Kickxellomycota:** Ku70/80 is in only 4% of 169 species, against 85–100% elsewhere. Needs tblastn confirmation, because BUSCO is about 71% in this group.

**Next (all need the genus tree):**
- Dollo/Count gain-loss reconstruction for the Dikarya-stem and convergent losses;
- phylogenetic logistic regression for FAS loss;
- Pagel correlated-evolution tests for the MoCo, sulfate and urease modules;
- PGLS and sister-pair tests for ECM contraction;
- CAFE for the NLR, 7TM and ankyrin expansions;
- tblastn validation of key absences.

## C. *Antarctolithica obscura* (taxid 3698444)

**Status:** running (Fable agent).

**The genomes:** 8 MAGs from the `FM_MARS` samples. They are unplaced within Ascomycota (no class, order or family), 7.3–10.5 Mb, with GC about 49%.

| LOCUSTAG | Strain | Size (Mb) | Contigs | N50 | BUSCO complete |
|---|---|---|---|---|---|
| FF1FC3B3 | Mars-82.9 | 9.76 | 89 | 181,073 | 93.2% |
| F46C92EF | Mars-17.30 | 9.60 | 124 | 127,517 | 93.0% |
| F7B4B958 | Mars-52.25 | 9.51 | 166 | 91,831 | 92.8% |
| F7E47E53 | Mars-62.17 | 10.47 | 363 | 46,728 | 92.1% |
| FA7C05DE | Mars-55.28 | 9.25 | 620 | 20,473 | 85.7% |
| F1D0D74C | Mars-14.71 | 7.82 | 213 | 63,415 | 77.9% |
| FD4728EA | Mars-9.1 | 7.46 | 644 | 14,431 | 63.5% |
| F3948469 | Mars-19.38 | 7.25 | 735 | 11,512 | 62.5% |

**Approach**
1. **MAG consistency and contamination.** Rely on the four MAGs at about 93% BUSCO. Count a gene as absent only if it is absent from all four.
2. **Per-feature comparison** against Ascomycota species representatives, grouped by class. Feature blocks:
   - genome and gene structure;
   - codon and amino-acid usage;
   - tRNA genes;
   - Pfam, CAZy, MEROPS and TF profiles;
   - secretome, disorder and localization;
   - paralogs and telomeres;
   - SwissProt best-hit taxonomy.
3. **Statistics for each block:**
   - ordination (PCoA) with PERMANOVA;
   - percentiles within each class, with Kruskal–Wallis and ANOVA tests;
   - a random-forest class prediction with feature importance.
4. **Distinctive content:** present in every good MAG but rare in other Ascomycota, or near-universal in Ascomycota but missing from every good MAG.
5. **Phylogenomic placement:** a plan, using `workflows/run_phyling.sh`.

**Status:** done 2026-11-10. Full report with 8 figures and cited methods: [analyses/antarctolithica/](analyses/antarctolithica/README.md). HPCC was in maintenance, so everything ran locally.

**Data quality**
- Only 6 of the 8 MAGs have annotations; Mars-9.1 and Mars-19.38 were skipped as too small.
- They are one species (median amino-acid identity 99.2–99.7%) but distinct genotypes.
- Mars-55.28 carries about 4.4% high-GC contigs from a second fungus.
- The annotation under-calls genes: 33 of 44 "absent" near-universal families are found in the assemblies by blastx. Single-gene absences are therefore weak evidence.

**Findings**
1. **A new deep lineage, not any sampled class.** Phylogenomics (300 BUSCO markers, 159 taxa, FastTree LG+G) places the MAGs as a monophyletic sister to Saccharomycotina + Pezizomycotina, branching after Taphrinomycotina. SH-like support is 1.0, and the placement is stable across 7 runs (marker jackknifes, WAG, long-branch taxa removed; 0.996–1.0). Its root-to-tip length is short, so this is not a long-branch artefact.
2. **The classifier agrees: it's out of distribution.** The highest random-forest class probability is 0.32, against a median of 0.94 for correctly classified genomes; classes left out of training score 0.20–0.49. Dothideomycetes, Eurotiomycetes and the lichen classes each score ≤0.04. The signal is a mosaic across feature blocks: gene structure looks like Taphrinomycetes, gene content like early Saccharomycotina or Taphrinomycotina, and codon, tRNA and TF usage like Pezizomycotina.
3. **The compaction is real, not missing data.** Extrapolated to 100% BUSCO the genome is about 10–10.4 Mb. It has 555 genes per Mb (99.8th percentile of Ascomycota genera); median intergenic distance (178 bp) and median intron length (34.5 bp) are both in the lowest 0.5%. It averages 0.56 introns per transcript, carries about 5,400 genes, and has a minimal tRNA set of 44–46 genes (one per anticodon class).
4. **No Pezizomycotina toolkit.** It lacks HET/NLR, NRPS, reducing PKS, DHN-melanin genes, hydrophobins, AA9, DNA methyltransferase and transposable-element domains, and has few P450s (12.5) and Zn2Cys6 TFs (33).
5. **What it retains looks like adaptation to UV and desiccation:** photolyases (4), UVDE, opsins, the trehalose pathway, expanded PhyH dioxygenases (21 vs median 6) and arylsulfatases (9–10 vs 1).
6. **Side finding:** "Saccharomycotina sp. YJ-2025a/b" group with Basidiomycota, so their phylum label in BFD looks wrong.

**Next:** a full nf_phyling / IQ-TREE run with site-heterogeneous models and AU tests (plan in the report's section 12), better gene annotation of the MAGs, and adding Taphrinomycotina-rich and *Neolecta*-like taxa to the placement analysis.

## D. Functional classifier for unknown genomes (planned)

**Status:** queued (Fable agent). It starts after C finishes, reusing C's feature matrices.

**Goal:** given a genome annotated by the BFD pipeline, predict:
- (a) its taxonomic placement from functional content;
- (b) its ecology: yeast vs. filamentous, saprotroph, plant pathogen, animal parasite, animal symbiont.

**Features**, built once for species representatives and stored in `bfd_scratch.rf_*` and as parquet on HPCC:
- Pfam presence, and counts normalized to proteome size;
- CAZy families and substrates, MEROPS families, TF families;
- secretome fraction, disorder and localization;
- gene structure: introns, gene density, intergenic distance, protein length from `gene_proteins`;
- codon and amino-acid usage (CLR-transformed) and tRNA counts.

**Design**

| Step | Labels | Cross-validation | Output |
|---|---|---|---|
| Taxonomy | PHYLUM, CLASS, ORDER | Grouped by genus, so the same genus is never in both training and test data | Accuracy by rank, plus a novelty score (distance to nearest training genomes, or forest proximity) so new lineages are flagged rather than forced into a class |
| Ecology | FungalTraits genus-level labels plus the atlas override file | **Leave one order or class out**, compared with a taxonomy-only baseline | Shows whether the ecology signal generalizes beyond ancestry |
| Incompleteness | n/a | Randomly drop 10–40% of genes during training | Models that handle MAGs and fragmented assemblies |
| Explanation | n/a | n/a | Permutation importance or SHAP for each prediction |

- **Which ecology traits are testable:**

  | Trait | Origins | Testable? |
  |---|---|---|
  | Yeast vs. filamentous | Many independent origins | Testable |
  | Animal parasite | Several origins | Testable |
  | Saprotroph vs. plant pathogen | n/a | Feasible, but the labels are noisy |
  | Animal symbiont | Few origins | Probably too few to generalize |

- **First unknowns to test:** the 8 *Antarctolithica* MAGs, including the two at about 63% BUSCO.
- **Pilot deliverables:**
  - a script that builds the feature matrix;
  - training and evaluation scripts;
  - a per-trait accuracy table (blocked vs. random cross-validation);
  - a CLI that classifies a new genome.
- **Production path:** an HPCC SLURM job (`workflows/`) that rebuilds the matrix from the parquet tables and retrains, with the model artifact versioned with each data release. Later the classifier could be exposed as an MCP tool or a dashboard view.
- **Later:** add phylogenetic eigenvectors once the genus-level tree from atlas PLAN is finished, to separate ecology from ancestry explicitly.

**Findings:** pending.
