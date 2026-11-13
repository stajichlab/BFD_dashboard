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
  | `kx2_h1_` … `kx2_h4_` | Kickxellomycota reduction hypotheses H1–H4 (E) |
  | `kx_` or none | Kickxellomycota and research agenda (A) |

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

## E. What was lost, gained and reshaped in Kickxellomycota reduction, and how it compares with yeasts

**Status:** first pass done 2026-11-10 (four parallel agents, exploration in MotherDuck). No narrative report yet. The SQL is saved in [analyses/kx2_kickx_reduction/](analyses/kx2_kickx_reduction/README.md) as `sql_h1/` to `sql_h4/` (numbered in run order, with a README per directory giving dependencies and a verification record); see "Reproducibility" below for what was and was not re-run. Scratch tables: `bfd_scratch.kx2_h1_*` to `kx2_h4_*`.

**Scope and conventions:** genus medians over species representatives. Quality filter BUSCO ≥ 85 and N50 ≥ 50 kb for main comparisons, unfiltered shown alongside. No genus has passed the filter in Harpellales or Dimargaritales, so those rows are unfiltered only. Core Kickxellales is 3–5 genera, almost all *Coemansia* (146 of 169 species). Baselines: Mucorales for Kickxellomycota and Zoopagales; Pezizomycotina for yeasts. Nothing is phylogenetically corrected.

**Hypotheses and outcomes**

| | Hypothesis | Outcome |
|---|---|---|
| H1 | Reduction removed or simplified multidomain proteins | **Not supported** |
| H2 | Kickxellales shrank by gene number and introns, yeasts by non-coding compaction | **Not supported** |
| H3 | Reduction pressure differs in kind: streamlining versus drift | **Weak lean**: core Kickxellales yeast-like, Harpellales drift-like |
| H4 | Kickxellales contracted multi-copy families, yeasts lost whole modules | **Mixed** |

**Findings**
1. **H1, architecture.** Mean protein length is unchanged (455 vs 447 aa in Mucorales, δ +0.09), and distinct architectures per 1,000 proteins are about 800–860 in every group. FAS (fused α+β), AROM and CAD are intact, and PKS, multi-NBD ABC and proteins ≥ 2,000 aa are retained. Proteins with ≥ 2 Pfam domains are 28.2% vs 33.2% in Mucorales (δ −0.58; −0.30 on non-partial genes), but Harpellales and Dimargaritales, which are not reduced, show the same direction, and the difference against 19 unreduced relatives is small. Only Microsporidia simplify (−148 aa; ABC and fusion enzymes lost). Within core Kickxellales the multidomain fraction correlates with BUSCO (r = 0.72) and partial-gene fraction (r = −0.70), so part of the deficit is assembly-dependent.
2. **H2, accounting.** About 70% of the size deficit is gene number and about 22% is tighter intergenic space per gene, for core Kickxellales vs Mucorales (15 Mb) and Saccharomycotina vs Pezizomycotina (27.5 Mb) alike. Introns are 9% and 2%; intron number fell about 8-fold but intron bp per gene only 48% (yeasts 71%). Convergent spacers shrank to the yeast value (189 vs 185 bp), but divergent (promoter-side) spacers stayed wider (825 vs 561 bp). **Not measurable here:** UTR is not annotated (exon length equals CDS length everywhere) and gene overlaps are absent from the table, so "intergenic" includes UTRs.
3. **H3, codon signature (genome-wide data only).** At matched GC3, core Kickxellales show more codon bias than expected (ENC-departure residual +0.054; yeasts +0.039, Mucorales +0.008), but tRNA supply matches codon demand less well than in yeasts (0.47 vs 0.61), and GC3 is high (67%, +14 points off the kingdom trend), so GC-biased gene conversion or mutation bias is not excluded. Harpellales are AT-rich (GC3 35%) with an AT-driven amino-acid profile and no excess bias, a drift-like pattern. No per-gene codon data exist in `bfd`, so there is no per-gene ENC or CAI and no ribosomal-protein subset, and no effective-population-size estimates are made.
4. **H4, loss composition.** Against the same Pezizomycotina baseline, Kickxellales lose 41% of gene copies (59% by contraction) and Saccharomycetes 44% (61%); about 30% of single-copy families are lost in each. Observed complete-loss families over the count expected under uniform thinning is 0.62–0.89 for Kickxellales and 0.75–0.92 for yeasts, so only Microsporidia look like targeted whole-family loss. Overlap of lost sets with Saccharomycetes is 0.93 of expected after stratifying the null by how many other lineages lost each family, so gene-content convergence disappears once dispensability is controlled. Overlap with Harpellales is 1.87, but shared ancestry cannot be excluded without a tree. Scaling: Kickxellales are near the proportional curve for CAZy, kinase, TF total and secretome; AKR is over-contracted; serine peptidases, GATA, HMG and Myb TFs are expanded. Yeasts are off the curve in P450 (0.07× the Pezizomycotina fraction), CE (0.08×), AA (0.25×) and SDR (0.37×).
5. **Gains (copy number only, no tree).** [PF00089](https://www.ebi.ac.uk/interpro/entry/pfam/PF00089/) trypsin: 21.5 copies in core Kickxellales (Harpellales 72.5) vs 2 in Mucorales; architecture is the same as Mucorales, so more copies, not new architectures. Also [PF00188](https://www.ebi.ac.uk/interpro/entry/pfam/PF00188/) (10 vs 3), [PF02458](https://www.ebi.ac.uk/interpro/entry/pfam/PF02458/) (17 vs 2–4), [GH46](https://cazy.org/famGH46) (6 vs 2). AA9, laccase ([PF00394](https://www.ebi.ac.uk/interpro/entry/pfam/PF00394/)) and tyrosinase ([PF00264](https://www.ebi.ac.uk/interpro/entry/pfam/PF00264/)) are **not** gained.

**Reading (interpretation)**
- Reduction in both lineages is mostly fewer genes with each gene's content preserved. Intron loss is the visible phenotype but a small share of the size deficit.
- The Mucorales baseline may overstate "reduction": Mucorales carry 2.5× the repeat fraction and 2.4× the intron bp per gene, and Zoopagales already have Kickxellales-sized genomes (about 19.7 Mb). Reduction is relative to an uncertain ancestral state until a tree exists.

**Caveats:** BUSCO of about 71% in Kickxellomycota inflates complete-loss calls, so those are upper bounds. The family decomposition (H1) and Pfam function assignment (H4, GO keyword argmax; 6,310 of 20,084 Pfams, 31%, are in SwissProt subjects, and the first report of 41% unclassified was wrong) are crude. In H1, `arch_per1000` and `bigram_per1000` are per 1,000 Pfam-annotated proteins in the 2,000-protein rarefied sample (annotated count 53–1,867), so groups with low Pfam hit rates such as Microsporidia rest on fewer proteins. Pfam e-value was not varied beyond 1e-5, and the Mucorales matrix includes recent expansions. Codon indices are genome-wide and confounded by GC. The tAI wobble weights and the ENC formula are marked citation to verify.

**Next, in priority order**
1. Species tree and orthogroups (Analysis A, item A1), then ancestral reconstruction of genome size, intron density and gene count against a phylogenetic baseline. This gates H3, H4 and the baseline issue. **Status 2026-11-11:** the atlas genus tree (1,654 genera) is *not* running: its Nov 9 nf_phyling run failed (`PHYLING_ALIGN` exit 137 after 4 h 25 min, then terminated at the maintenance drain) and has not been resubmitted. A cheaper Kickxellomycota-focused pair of trees (65-taxon backbone with outgroups, and 107 *Coemansia* genomes) is drafted in `workflows/phylo/` (taxa in `bin/phylo/kickx_tree_taxa.sql`); the first-pass backbone (63 taxa, `fungi_odb12`, FastTree and ASTRAL) and a *Coemansia* tree (121 taxa, FastTree) ran on the HPCC on 2026-11-11. In the *Coemansia* tree the genus is not monophyletic (*Kickxella* nests inside it, weakly placed) and the strains are redundant (20 of 29 closest pairs differ by under 0.005 substitutions per site), so lineages need ANI dereplication. In both trees Kickxellomycota is monophyletic with Dimargaritales earliest and Kickxellales sister to the unreduced Harpellales, so Ku loss may be a single stem event preceding the intron and size reduction (a hypothesis, untested); *Linderina* and *Basidiobolus* are unstable, and the zygomycete backbone is not resolved. Details and caveats are in `workflows/phylo/README.md`. The atlas genus tree is being built on another machine and must be reconciled with it (tip labels are LOCUSTAGs here), and IQ-TREE with a site-heterogeneous model and a marker-sensitivity run are still to do.
2. tblastn checks of the lost single-copy families in good *Coemansia* assemblies (the check on H4's complete-loss calls).
3. Per-gene codon usage from CDS on the HPCC (ribosomal proteins as the high-expression set) and GC4 tests, to separate selection from GC bias.
4. Trypsin tandem-array and repeat check, and a MEROPS family map.
5. RNA-seq UTR annotation for a few genomes, and repeat coordinates, to open the UTR and masked terms of the accounting.

**Reproducibility** (SQL in `docs/analyses/kx2_kickx_reduction/sql_h1/` to `sql_h4/`, MotherDuck dialect, plain SQL, nothing done in Python or by hand; each README has the dependency graph). Verification was done by re-running the later steps from the existing upstream scratch tables into temporary tables and diffing against the originals:

| | Re-run and matched | Saved but **not** re-run |
|---|---|---|
| H1 | steps 03–07 (genome metrics, arch metrics, genus medians, Cliff's δ; 0 of 300 δ rows differ) | 01, 02 (genome list, 5M-row protein table), 09–12 (family decomposition, fusion, long families), 13–14 |
| H2 | 05 lineage terms (37 rows), 06 decomposition (both orderings), 08 genus table (row count only), 09 δ (119 rows identical) | 01–04 and 07 (genome list, per-gene table, intergenic distances, orientation) |
| H3 | 07–09 and the `resid` step of 10 (0 differing rows); ENC and GC3 recomputed independently (max difference 3e-14) | 04 (tRNA and tAI), 05 (amino acid, disorder), 06 (protein-rate proxy), `resid2` |
| H4 | loss split, lost-set overlap, gains (0 row differences) | 01, 02, 04–06 (medians, class counts, scaling fits, Pfam categories) |

- H2 and H3 start from the Analysis A scratch table `kx_genomes` (1,919 genomes), whose SQL is in [analyses/kickxellomycota/](analyses/kickxellomycota/README.md); it was not rebuilt.
- H2's lineage list is hand-defined over `kx_genomes.grp`, not derived from taxonomy.
- Some reporting queries (H4 leverage and candidate listings; H3 display queries) are comments, not tables. Exploratory design queries that only informed choices are not saved.
- A full end-to-end rebuild from empty scratch has not been done, so "saved" is not yet "proven to rebuild".

**Follow-up: Ku loss and its pathway (2026-11-11).** First pass, report and SQL in [analyses/kx2_ku_pathway/](analyses/kx2_ku_pathway/README.md). On the first-pass backbone tree Ku is lost once, on the stem of Kickxellales + Harpellales; tblastn on 128 genomes (all four Ku positive controls recovered) confirms Ku is absent from all 20 genomes of the 11 stem genera and from all 107 *Coemansia* strains, and shows three Pfam "absences" were annotation misses. By Pfam, Lig4, XLF, Pol4, aprataxin and DNA-PKcs were already lost earlier, which would make Ku loss a late step in an already-dismantled NHEJ rather than its trigger. tblastn supports this for Pol4 (none in any Kickxellomycota genome), and a ligase gene-family tree (203 sequences, all 11 Lig4 references recovered in leave-one-out) shows the "strong" Lig4 tblastn hits were LIG1-type: no annotated Kickxellomycota protein places in the Lig4 clade, so Lig4 was lost on the Kickxellomycota stem before Ku (limits: *Ramicandelaber* unresolved, rescue layer untested, first-pass trees). XLF, aprataxin and DNA-PKcs remain uninformative by tblastn. A RecA-family tree (327 sites; leave-one-out 19 of 31, Dmc1 5 of 5) shows the meiotic recombinase **Dmc1 is absent from the five core Kickxellales genera** (including all 10 *Coemansia* strains in the tree) but present in *Mycoemilia*, Harpellales, Dimargaritales, Entomophthorales and most Mucoromycota, so it is lost later than Ku, on the stem of the same clade that shows the strongest intron and size reduction (a tblastn locus count is consistent but under-detects). Rad55/Rad57 are absent from Mucoromycota too, so the earlier idea of a Kickxellales-specific loss is not supported (details in `analyses/kx2_ku_pathway/family_trees/RECA_RESULT.md`). MRN, Dna2, TERT, Rad51/52, Mus81 and FEN1 are retained. No compensation was found: POLQ, LIG3 and PARP1-like were absent before Ku loss, HR genes are not expanded, and no Ku-specific telomere change is detectable. The tblastn round 2 and 3 queries were chosen after earlier rounds (post hoc), and the species tree and gene-family trees are first-pass.

**Follow-up: intron cross-check with miniprot (2026-11-12).** Report in [analyses/kx2_kickx_reduction/intron_crosscheck/](analyses/kx2_kickx_reduction/intron_crosscheck/README.md). Introns were called from spliced alignments of one fixed set of 1,122 BUSCO `fungi_odb12` proteins on 83 assemblies (no gene finder; three intron-penalty settings, three identity thresholds). **The pre-registered success rule is met:** core Kickxellales are below Mucorales (Cliff's delta -1.00) in all 9 setting × threshold cells and in every identity bin, and a bootstrap over reference proteins agrees. Default setting, identity ≥ 0.5: 0.30 introns per kb of CDS in core Kickxellales, against 2.10 in Mucorales, 1.2 in Harpellales, 0.16 in intron-poor yeast anchors and 1.06 in intron-rich Pezizomycotina anchors; agreement with the funannotate counts is rho = 0.86 across the 83 genomes. The intron deficit is therefore not a gene-model artifact. A strong dependence on the penalty setting was real, but at the most permissive setting 85% of called introns in core Kickxellales are non-canonical (alignment gaps), and counting canonical introns only (post hoc) keeps delta at -1.00. Dimargaritales + *Ramicandelaber* are as intron-poor as the core clade by this measure (0.54 per kb), so intron loss seems to begin earlier than the core clade (descriptive). Limits: 5 versus 8 genera per group, conserved genes only, miniprot aborts on one genome under the most permissive setting.

**Follow-up: meiotic toolkit (2026-11-12).** Report in [analyses/kx2_ku_pathway/meiosis_toolkit/](analyses/kx2_ku_pathway/meiosis_toolkit/README.md). Spo11 (diverged), Hop2, Mnd1 and Msh4 are retained in core Kickxellales; Msh4 is placed in 14 of 14 core genomes by a MutS tree and Msh5 in 7 of 14 (a partial call, not a loss). **Pch2 is absent from Kickxellales and Harpellales**, on the same stem as the Ku loss, and Mlh3 is not seen in any of 21 Kickxellomycota genomes in the trees (patchy in Mucoromycota too, so not specific to the core clade). Rec8 and Hop1 cannot be called (their trees fail the leave-one-out gate). Together with Dmc1 (core clade), the evidence supports a loss of a few meiosis-associated genes at different nodes, not a general collapse of the toolkit in core Kickxellales; the 7-gene completeness score is lower there (median 4 versus 6-7) but tracks assembly quality. Limits: annotated layer only, 14 core genomes in the trees, annotation misses in fragmented assemblies.

**Production path:** add a `workflows/atlas/J12*.sh` job that rebuilds these metrics from the parquet tables on HPCC (reusing `A3/protein_pfam_arch.parquet` from J2b), with a fixture test and the full rebuild above. This is not done yet. The numbers are exploratory until the end-to-end rebuild agrees with them.

## Upstream follow-ups (BFD/Fungi_BFD_runs)

To do when HPCC is back online. These are also in `TODO.md`.

| Issue | Evidence | Action |
|---|---|---|
| Gene annotation under-calls | Analysis C: 33 of 44 near-universal families absent from all 6 annotated *Antarctolithica* proteomes are found in the assemblies by blastx, including ribosomal proteins L22e, L34e and S21e. Mars-9.1 and Mars-19.38 were not annotated (too small) | Find out why gene prediction missed them: training set, filtering, contig edges, short or single-exon genes. Check other small or fragmented genomes for the same problem. Treat single-gene absences as weak evidence until this is fixed. **Added 2026-11-11:** misses are genome-specific (Ku missed in *Tieghemiomyces*, *Caulochytrium*, *Conidiobolus*; *Massospora* has 415 proteins), while in 107 *Coemansia* strains protein-mode BUSCO equals genome-mode (median gap -0.2 points; [analyses/kx2_kickx_reduction/annotation_check/](analyses/kx2_kickx_reduction/annotation_check/README.md)), so no systematic loss is indicated there for conserved genes. Proposed remedy: a validated, source-tagged miniprot gap-filling step after predict for genomes with a large BUSCO gap, a tiny proteome or a MAG, extending the existing `do_annotation_busco_gap` diagnostic (details in `TODO.md`); not to be adopted until the rescue passes its validation |
| `gene_introns.codon_position`, `codon_frame` | Off by one exon; confirmed 2026-11-10 | Fix the table builder and rebuild `gene_introns` |
| "Saccharomycotina sp. YJ-2025a/b" (FF90B97F, F6EB82EF) | Group with Basidiomycota in the analysis C tree | Check and correct the phylum. Excluded from analyses C and D |

