# Kingdom-wide survey of functional-domain gain, loss, expansion and contraction in BFD

Agent B report, 2026-11-10. Exploratory and data-driven; it builds on the summary atlas (`docs/summary_atlas`) and does not repeat its hypothesis tests. Kickxellomycota is only noted in passing (another agent covers it).

## 1. Main findings

1. **Taxonomy dominates, lifestyle is a second-order signal.** After normalising copy number to proteome size, class membership explains a median 44 % of the between-genus variance of a domain family (2,606 Pfam families); lifestyle within order explains a median 2.4 %. CAZy families are the exception (median 5.4 %, up to 19 % for GH10).
2. **The Dikarya ancestor lost a block of ancestral, animal-like machinery.** 203 Pfam domains are present in at least half the genera of two or more early-diverging groups and in 5 % or fewer of 1,083 Dikarya genera; 191 of them are still in Mucoromycota s.l. They include cobalamin-dependent enzymes, Integrator, the Fanconi anaemia core, KICSTOR and biopterin metabolism (Fig 9).
3. **Reduced lineages lose the same genes.** Of 3,568 "core" domains, Microsporidia lack 62 %, *Hanseniaspora* 28 %, *Pneumocystis* 23 %, Neocallimastigomycetes 17 %, *Malassezia* 13 %, *Schizosaccharomyces* 12 %, *Wallemia* 10 %. Loss sets of independent lineages overlap about 3 to 11 times more than chance (log2 1.4 to 3.5), and co-absence clustering recovers pathway-level modules: complex I, beta-oxidation plus glyoxylate cycle, molybdenum cofactor, urease, glycogen metabolism, RNAi, type I fatty-acid synthase, sulfate assimilation (Fig 2).
4. **Fatty-acid synthase is absent from four independent host-dependent lineages**: 0 of 21 *Malassezia* species, 0 of 7 *Pneumocystis*, 0 of 31 Glomeromycota, 0 of 93 Microsporidia (and 27 % of 11 Zoopagomycota), against 93 to 100 % elsewhere.
5. **Pezizomycotina are defined by a non-self-recognition and surface-sensing expansion**: HET, NACHT, NPHP3-type NTPase, ankyrin, the "Fung_rhodopsin" 7TM family and CFEM. These are zero or near zero in every yeast lineage and in most Basidiomycota. Their copy number tracks class (eta squared 0.49 to 0.91), not lifestyle (partial R2 0.02 to 0.05).
6. **Ectomycorrhizal genera contract far more than plant-cell-wall enzymes.** Within orders they also carry fewer short-chain dehydrogenases, GSTs, aldo-keto reductases and methyltransferases and a smaller secretome, while kinases are unchanged and a transposon integrase domain is higher.
7. **Glomeromycota combine extreme expansions with metabolic loss**: median 2,020 protein-kinase proteins per genome (125 elsewhere), 734 Sel1, 622 TLDc; FAS and catalase absent.
8. **Many "lineage-restricted" domains are annotation history, not biology.** Saccharomycetales has 210 restricted domains from 13 genera because Pfam families were built on *S. cerevisiae* proteins. Restricted domains in zoosporic fungi are retained flagellar genes, not gains.

Caveat that applies to everything below: no test here is phylogenetically corrected (the genus tree is not finished). Clade-balanced weights, class and order fixed effects, and counts of independent lineages are used as stand-ins, and every pattern in section 8 names the tree-based test that should replace them.

## 2. Data and compute

- **Source.** Pfam hits (`bfd.pfam`, 368 M rows), CAZy (`cazy_overview`), MEROPS, `tf_inventory`, SignalP/TMHMM/PredGPI secretome counts, `idp_summary`. I reused the atlas intermediates that already aggregate these per genome (`atlas_data/A7/pfam_genome_matrix.parquet`, `A5/*`, `A4/genome_protein_summary.parquet`, `A11/pfam_trait_counts.parquet`, `labels/genome_labels*.parquet`), as instructed.
- **Compute.** The HPCC was in a full maintenance reservation (`20261110-Shutdown`) during this session, so no SLURM job could start. I copied the small atlas parquet files (about 370 MB) to the local scratchpad and ran everything there with DuckDB, numpy and scipy in a scratch virtualenv. Nothing heavy ran on the login node; my one pending job was cancelled.
- **MotherDuck scratch tables** (database `bfd_scratch`, all prefixed `kw_`; `bfd` was not modified):

| table | rows | content |
|---|---|---|
| `kw_rep_genomes` | 7,435 | one representative genome per species (annotated; best BUSCO complete, then N50, then ASMID); 1,669 genera. Same count as the atlas `is_representative` set |
| `kw_pfam_matrix` | 32,265,977 | genome x Pfam: distinct proteins and domain count, gate `domain_i_evalue <= 1e-5`, version suffix stripped; 7,392 genomes x 15,052 Pfams |
| `kw_genus_pfam` | 7,479,899 | genus x Pfam: share of representative genomes with the domain, mean copies (no quality filter); 1,654 genera |

These mirror the local matrices so the analyses can be rerun in SQL. The quality flags and FungalTraits labels exist only in the atlas parquet, so the statistics below were computed locally, not in MotherDuck.

SQL used for the matrix (one scan of `pfam`):

```sql
CREATE TABLE bfd_scratch.kw_rep_genomes AS
WITH g AS (SELECT s.*, b.complete_pct, a.N50_bp,
   row_number() OVER (PARTITION BY s.SPECIES ORDER BY b.complete_pct DESC NULLS LAST, a.N50_bp DESC NULLS LAST, s.ASMID) rk
  FROM bfd.species s LEFT JOIN bfd.busco_genome b USING (ASMID) LEFT JOIN bfd.asm_stats a USING (ASMID)
  WHERE a.annotation_status = 'annotated')
SELECT * FROM g WHERE rk = 1;

CREATE TABLE bfd_scratch.kw_pfam_matrix AS
SELECT split_part(protein_id,'_',1) AS genome, split_part(pfam_acc,'.',1) AS pfam_acc, any_value(pfam_id) AS pfam_id,
       count(DISTINCT protein_id) AS n_proteins, count(*) AS n_domains
FROM bfd.pfam
WHERE domain_i_evalue <= 1e-5 AND split_part(protein_id,'_',1) IN (SELECT LOCUSTAG FROM bfd_scratch.kw_rep_genomes)
GROUP BY 1, 2;
```

- **Scripts** (kept in `kw/` next to this report): `step1_local.py` (matrices), `lib.py` (loaders, clade groups), `step2_restricted.py`, `step3_loss.py`, `step4_expand.py`, `step5_classes.py`, `step6_checks.py`, `make_tables.py`, `make_figures.py`. Intermediate tables are in `kw/work/` (TSV and parquet). Figures are in `fable_B_figs/`.

## 3. Analytical choices and why

**Unit of replication: the genus.** Strains of a species and species of a genus are not independent, and BFD sampling is very uneven (1,614 Sordariomycetes species against 7 *Pneumocystis*). Following the atlas, one genome represents each species (best BUSCO, then N50), and species are collapsed to genus: a domain is *present in a genus* when at least 50 % of its representative genomes carry it, and the genus copy number is the median across those genomes (zeros included). The 50 % rule tolerates single missed gene models and single contaminated assemblies. Genera are still not independent (Felsenstein 1985, Am Nat 125:1), so I add clade weights or fixed effects everywhere and treat p-values as screening statistics only.

**Quality sets.** "q genera" (n = 1,173, from 4,743 genomes) pass the atlas filters: BUSCO complete at or above the 25th percentile of the genome's class, BUSCO duplication at most 10 %, N50 at least 50 kb. BUSCO is judged relative to class because low BUSCO is real biology in Microsporidia and budding yeasts (Simão et al. 2015, Bioinformatics 31:3210; Manni et al. 2021, Mol Biol Evol 38:4647). "All genera" (n = 1,654, 7,411 genomes) drops the three filters. Presence inside a clade is measured on q genera; absence outside a clade is checked on all genera, so a domain is not called restricted merely because low-quality genomes were removed. Loss analyses of single-genus lineages (*Malassezia*, *Pneumocystis*, *Schizosaccharomyces*, *Wallemia*, *Hanseniaspora*) use all genera and are confirmed at species level (Table 5).

**Domain gate.** Pfam hits with domain i-Evalue at most 1e-5, counting distinct proteins per family per genome (the atlas A7 rule). This is stricter than Pfam gathering thresholds for some families and looser for others; it was kept for consistency with the atlas (Pfam: Mistry et al. 2021, Nucleic Acids Res 49:D412).

**Clade groups.** PHYLUM, CLASS, ORDER as in `species`. SUBPHYLUM is empty, so I assigned subphylum-level groups from class names (Saccharomycotina = the seven budding-yeast classes in the current taxonomy; Taphrinomycotina; Pezizomycotina; Agaricomycotina; Ustilaginomycotina; Pucciniomycotina) and three pooled early-diverging groups (Mucoromycota s.l. = Mucoromycota, Mortierellomycota, Glomeromycota, Calcarisporiellomycota; Zoopagomycota s.l. = Zoopagomycota, Kickxellomycota, Entomophthoromycota, Basidiobolomycota; zoosporic = Chytridiomycota, Blastocladiomycota, Neocallimastigomycota, Monoblepharomycota, Olpidiomycota, Sanchytriomycota). The pooled groups are convenience groups, not clades.

**Lineage restriction (section 4).** For every clade with at least 5 q genera I computed the fraction of clade genera with the domain (`in`) and the number and fraction of genera outside with it (`out`, on all genera). A domain is *restricted* when in >= 0.8 and out <= 0.02; *strictly novel* when out = 0. Each domain is assigned to the clade maximising Youden's J = in - out (Youden 1950, Cancer 3:32), ties to the larger clade, so a domain is not reported at several nested ranks. This is a phylogenetic-profile screen (Pellegrini et al. 1999, PNAS 96:4285). It cannot tell gain from retention; that needs an outgroup or a gain/loss reconstruction.

**Ancestral and core sets (section 5).** Without a tree or non-fungal outgroups, ancestry is approximated by breadth: a domain is *ancestral* if present in at least 50 % of q genera in two or more of the three early-diverging groups, or in one of them plus both Ascomycota and Basidiomycota (4,475 domains). It is *core* if it is also in at least 70 % of Pezizomycotina and 70 % of Agaricomycetes genera (3,568 domains). Core domains are present on both sides of every deep split, so their absence from a nested lineage is a derived loss under Dollo-type reasoning (a complex character is unlikely to be regained; Farris 1977, Syst Zool 26:77). A lineage has *lost* a domain when at most 10 % of its genera carry it.

**Co-loss overlap.** For two lineages, observed shared losses are compared with the expectation under independence, n_A x n_B / 3,568, reported as log2(observed/expected). This is the hypergeometric mean; with these counts every reported enrichment has p far below 1e-10, so I report effect sizes only.

**Co-absence modules.** Among core domains absent from 3 to 70 % of genera (class-weighted), I computed the phi coefficient (Pearson correlation of binary absence vectors) across q genera, with each genus weighted by 1/(number of genera in its class) so that every class contributes equally. Microsporidia were excluded because they lack 62 % of all core domains and would create one giant module. Domains were clustered by average linkage (UPGMA; Sokal and Michener 1958, Univ Kansas Sci Bull 38:1409) on 1 - phi and cut at 0.35 (within-module phi about 0.65 or higher). Correlated presence/absence predicts functional linkage (Pellegrini et al. 1999), but raw correlation across species overstates the evidence because one ancient loss is inherited by many descendants (Barker and Pagel 2005, PLoS Comput Biol 1:e3; Maddison and FitzJohn 2015, Syst Biol 64:127). Class weighting reduces but does not remove this. For each module I therefore list the number of separate lineages in which it is absent; that count, not the correlation, is the evidence for convergence.

**Copy-number normalisation (section 6).** Functional categories scale with genome size as power laws with category-specific exponents (van Nimwegen 2003, Trends Genet 19:479; Molina and van Nimwegen 2009, Trends Genet 25:243), so "copies per 1,000 proteins" is biased whenever the exponent is not 1. For each family present in at least 25 % of genera with mean genus count of at least 1 (2,831 families of all types) I fitted log10(count + 1) = a + b log10(proteome size) across q genera by weighted least squares with class-balanced weights (as above), and used the residual. A residual of +0.3 means twice the copy number expected for that proteome size. The +1 pseudocount keeps zeros; it compresses families with low counts, which is why the ranking requires a median of at least 2 copies. Per-1,000 medians are also given (Table 9) because they are easier to read.

**Taxonomy versus lifestyle.** eta squared (between-group over total sum of squares) of residuals on CLASS (classes with fewer than 5 genera pooled). Lifestyle is the FungalTraits primary lifestyle of the genus (Põlme et al. 2020, Fungal Diversity 105:1); lifestyles with fewer than 10 q genera are pooled as "other"; 1,054 of 1,173 genera have a label. Its effect is the partial R2 after removing class or order fixed effects, obtained by residualising both the response and the lifestyle dummies on the taxon dummies (Frisch-Waugh-Lovell; Frisch and Waugh 1933, Econometrica 1:387; Lovell 1963, J Am Stat Assoc 58:993). Order fixed effects (orders with fewer than 3 genera pooled within class) are the conservative version: only genera that differ in lifestyle within the same order contribute. F-test p-values are corrected by Benjamini-Hochberg (Benjamini and Hochberg 1995, J R Stat Soc B 57:289). Fixed effects for taxonomic rank are a crude substitute for a phylogenetic covariance; they assume independence within orders.

**Functional classes.** CAZy family and substrate from dbCAN (Zhang et al. 2018, Nucleic Acids Res 46:W95; CAZy: Drula et al. 2022, Nucleic Acids Res 50:D571), MEROPS family of the best hit (Rawlings et al. 2018, Nucleic Acids Res 46:D624), TF family from `tf_inventory`, strict secretome (SignalP >= 0.9, at most one TM helix, no GPI anchor), GPI, TM and disorder from the atlas A5 tables. They pass through the same allometry and variance partition as Pfam families.

## 4. Lineage-restricted domains

{n_strict} clade-domain pairs meet in >= 0.8 and out <= 0.02; after assigning each domain to its best clade, 1,041 distinct domains are restricted to some clade (6.9 % of 15,052).

![Fig 1](fable_B_figs/fig1_restricted_heatmap.png)

**Table 1. Restricted domains per clade** (best-clade assignment; `n_zero_outside` = present in >= 80 % of clade genera and in no genus outside, counted per clade without the best-clade rule).

{restr_counts}

**Table 2. Top restricted domains per clade** (lowest outside count first).

{restr_top}

Observations.

- **Dikarya (17), Ascomycota (71), Pezizomycotina (108), Saccharomycotina (55), Basidiomycota (6), Agaricomycetes (4).** The Ascomycota set is mostly named regulators and complex subunits known from yeast and *Aspergillus* (Rox3/Med19, PalH/Rim21, eisosome, DASH Spc34, Fig1). The Pezizomycotina set is 47 strictly novel domains, many still DUFs (DUF8004 and DUF8035 are in 100 % of 666 genera and in none of 766 other genera), plus KpeA, AcrB and Ace1 regulators. Basidiomycota have few restricted domains but two are large families: DUF6534 (97.9 % of 291 genera; median 50.5 proteins per genome in Agaricomycetes, up to 211 in *Gerronema*) and its partner DUF6533, and the Clampless-1 domains (94 %).
- **Microsporidia (24; 14 strictly novel)**: polar-tube proteins PTP2 and PTP3, spore-wall proteins EnP1, SWP7, MICSWaP. These are the infection apparatus; a clean positive control for the method.
- **Zoosporic fungi (38 pooled, 89 Chytridiomycota)**: intraflagellar transport, axonemal dynein, centriole and cilium proteins. This is retention of the ancestral flagellum, lost in the non-flagellated fungi, not a gain.
- **Mucoromycota (21), Mortierellomycota (46), Glomeromycota (58)**: a mix of retained ancestral genes (UFL1 UFMylation ligase, Castor1 arginine sensor, PEPCK-GTP) and real lineage expansions: WRKY DNA-binding domain in 100 % of Mucoromycota genera with a median of 5 copies in Mucoromycetes, RNF213-type AAA domains at 15 to 19 copies and GVIN/URGCP GTPases in all 6 Glomeromycota genera.
- **Small yeast clades dominate the count**: Saccharomycetales 210 from 13 genera, Saccharomycetes 14, Serinales 12. Saccharomycotina-wide strictly novel domains include glucose-signalling regulators (STD1/MTH1, Gsf2), Vhr1 and POC1.

Interpretation and caveats.

- Pfam is built from characterised proteins, so clades containing *S. cerevisiae*, *S. pombe*, *Aspergillus*, *Neurospora* and *Ustilago* have many narrow families and under-studied clades have few. The Saccharomycetales count is this bias. Counts in Table 1 are not comparable between clades as a measure of innovation.
- Fast-evolving proteins fall below the HMM threshold in distant relatives, which looks like restriction (homology-detection failure; Weisman, Murray and Eddy 2020, PLoS Biol 18:e3000862). The Saccharomycotina "restricted" regulators may be of this kind.
- Clades with 5 or 6 q genera (Spizellomycetales, Glomeromycota, Kickxellomycota, Orbiliomycetes, Geoglossomycetes) meet in >= 0.8 by chance more easily; their lists have 6 to 30 outside genera and should be read as enrichment, not restriction.
- Bacterial-type domains restricted to one small clade (YscW in Mortierellomycota, HSDR_N_2 in Orbiliomycetes) may be endosymbiont or contaminant sequence; Mortierellomycota are known to carry *Burkholderia*-related endobacteria. Check scaffold-level taxonomy before interpreting.

## 5. Repeated loss and loss modules

### 5.1 Ancestral domains missing from Dikarya

{n_dl} ancestral domains are present in 5 % or fewer of 1,083 Dikarya q genera; 191 of these are in at least half of the 33 Mucoromycota s.l. genera. Only one is confined to zoosporic fungi, because flagellar domains fail the "two early-diverging groups" rule; this set is therefore not the flagellum.

![Fig 9](fable_B_figs/fig9_dikarya_lost.png)

Recognisable sets (Fig 9): cobalamin-dependent metabolism (methionine synthase B12 domains, B12-binding, cobalamin adenosyltransferase, MMADHC, propionyl-CoA carboxylase BT domain: 67 to 100 % of early-diverging genera, 0 to 0.1 % of Dikarya); Integrator subunits INTS2/3/4/7; Fanconi anaemia core FANCD2/FANCI/FANCL; KICSTOR (kaptin, ITFG2, C12orf66) and LAMTOR5; SMG1 and Smg8/9; DNA-PKcs; biopterin-dependent aromatic amino-acid hydroxylase and PTPS; presenilin; BK potassium channel; Tubby; tripeptidyl peptidase II; type III glutamine synthetase.

Interpretation: one or a few losses on the Dikarya stem removed several animal-type regulatory and repair complexes. B12 loss in Dikarya is published (Orłowska et al. 2021, Genome Biol Evol; citation to verify). I did not find a published kingdom-wide statement for Integrator, Fanconi and KICSTOR as a group; I have not searched the literature systematically. Alternative: these domains are in Dikarya but diverged below threshold. This is testable with profile-profile or structure searches on a few Dikarya proteomes.

### 5.2 Core-domain loss per lineage

**Table 3. Core domains (of 3,568) absent from each lineage** (present in <= 10 % of its genera; all-genera set).

{loss_lineage}

![Fig 3](fable_B_figs/fig3_coabsence_and_loss_vs_size.png)

- Loss scales with proteome reduction but not only with it (Fig 3b). Neocallimastigomycetes have a large proteome (median 14,667) and lack 17 % of core domains (BUSCO 64 %, so part of this is incompleteness; the mitochondrial losses are expected for anaerobes with hydrogenosomes). Glomeromycetes (27,234 proteins) lack 3.3 %; Entomophthoromycetes (18,241) lack 6.1 %.
- *Hanseniaspora* lacks 1,016 core domains, four times the Saccharomycetes median (270), consistent with the published loss of cell-cycle and DNA-repair genes in this genus (Steenwyk et al. 2019, PLoS Biol 17:e3000255). Its median BUSCO is 49 %, which is itself a consequence of the losses; I cannot separate real loss from annotation gaps with this table alone.
- The large filamentous classes lack 0 to 2 core domains each, by construction of the core set.

**Overlap between independent lineages** (log2 observed/expected shared losses; full matrices in `kw/work/step3.out`): Saccharomycetes vs *Schizosaccharomyces* 2.2 (147 shared; 32 expected); *Malassezia* vs *Wallemia* 2.0 (185 shared; 46 expected); *Malassezia* vs Saccharomycetes 2.2; *Pneumocystis* vs *Malassezia* 1.4 (291 shared; 108 expected); Pucciniomycetes vs Ustilaginomycetes 4.0; Kickxellomycota vs Zoopagomycetes 3.4 and vs Entomophthoromycetes 3.3. Lecanoromycetes, a control with 2 losses, shows none. Every pair of reduced lineages overlaps more than expected, so there is a common pool of dispensable domains. 1,072 core domains are never lost in any of the 37 lineages; 850 are lost in exactly one; 114 are lost in 10 or more. The most often lost: MIF (23 lineages), Thaumatin, alkaline ceramidase, MIOX, cyanate lyase, NACHT, tyrosinase, PARP, animal haem peroxidase (12 to 18 lineages each).

### 5.3 Loss modules

1,697 core domains are variably absent; average-linkage clustering of class-weighted co-absence gives 58 modules of 4 or more domains. Modules were named by me from their Pfam content.

![Fig 2](fable_B_figs/fig2_loss_modules.png)

**Table 4. Selected co-absence modules.**

{modules}

**Table 5. Species-level check of marker domains** (fraction of representative species with the domain; all species representatives, no quality filter). FAS = Fas_alpha_ACP PF18325; Piwi PF02171; complex I 49 kDa PF00346; MoaC PF01967; urease alpha PF00449; glycogen synthase PF05693; glucan synthase PF02364; isocitrate lyase PF00463; Ku PF02735; catalase PF00199; ATP-sulfurylase PF01747.

{species_loss}

Observations.

- **Respiratory complex I** (17 co-absent domains): absent from Microsporidia, Neocallimastigomycota, *Schizosaccharomyces* (0 of 7), Saccharomycetaceae (0 of 158) and *Hanseniaspora* (0 of 22), present in 96 to 99 % of Dipodascomycetes and Pichiomycetes. At least four independent losses.
- **Beta-oxidation plus glyoxylate cycle** (17 domains): absent from Microsporidia, Neocallimastigomycota, *Pneumocystis*, *Schizosaccharomyces* (isocitrate lyase 14 %) and *Hanseniaspora* (0 %).
- **Type I FAS** (8 domains): see finding 4. Lipid dependence is documented for *Malassezia* (Xu et al. 2007, PNAS 104:18730), arbuscular mycorrhizal fungi (Wewer et al. 2014, Plant J 79:398; Jiang et al. 2017, Science 356:1172) and *Pneumocystis* (Ma et al. 2016, Nat Commun 7:10740). The BFD matrix recovers all of them in one screen and adds a partial signal in Zoopagomycota (3 of 11 species retain FAS).
- **Molybdenum cofactor** (7 domains) is the most often lost pathway: absent or rare in 9 lineages at module level, and at species level in Lipomycetes (0 %), Pichiomycetes (19 %), Saccharomycetaceae (0 %), Tremellomycetes (25 %; *Cryptococcus* 7 %), *Puccinia* (0 of 11), *Malassezia*, *Wallemia*, *Schizosaccharomyces*, *Pneumocystis*, Neocallimastigomycota. Mo-dependent enzymes (nitrate reductase, sulfite oxidase type) form the adjacent module.
- **Urease** is absent from all Saccharomycotina classes (0 to 1 %), consistent with replacement by urea amidolyase, and also from Erysiphaceae (6 %), *Pneumocystis*, *Batrachochytrium*, Blastocladiomycota and Neocallimastigomycota.
- **Sulfate assimilation** (ATP-sulfurylase): absent from Erysiphaceae (0 of 17), *Puccinia* (0 of 11), Zoopagomycota (0 of 11), *Pneumocystis*, *Batrachochytrium*, Blastocladiomycota (11 %), Neocallimastigomycota and Microsporidia. The powdery-mildew and rust losses are published (Spanu et al. 2010, Science 330:1543; Duplessis et al. 2011, PNAS 108:9166); the zoosporic and Zoopagomycota cases extend the pattern.
- **RNAi** (Dicer plus Argonaute domains): Dicer dimerisation domain is absent from Saccharomycetes, Pichiomycetes and Dipodascomycetes (0 %), *Pneumocystis*, *Malassezia* (5 %), *Wallemia* and *Hanseniaspora*; Piwi is retained in 54 % of Saccharomycetaceae and 63 % of *Candida*, matching the non-canonical Dicer of budding yeasts (Drinnenberg et al. 2009, Science 326:544).
- **Glycogen metabolism** (9 domains): absent from Microsporidia, *Malassezia* (5 %), *Wallemia* and *Schizosaccharomyces* (0 %).
- **Catalase** is absent from Glomeromycota (0 of 31), *Pneumocystis*, Neocallimastigomycota, and rare in Chytridiomycota (35 %) and Blastocladiomycota (11 %).
- **Ku70/Ku80** is in 4 % of 169 Kickxellomycota species and in 85 to 100 % of almost every other group (noted for the Kickxellomycota agent; BUSCO there is 71 %, so confirm with a tblastn search).
- Two modules are not losses in reduced lineages but late gains: the DASH kinetochore complex and the SUR7/KAR9/ATG11 set are absent from Chytridiomycetes, Blastocladiomycota and Neocallimastigomycota and present in Mucoromycota and Dikarya. They passed the ancestral filter through Mucoromycota s.l. and Zoopagomycota s.l. They date innovations to after the zoosporic grade.

Caveats. Absence in lineages with BUSCO below about 80 % (Neocallimastigomycota 73 %, Zoopagomycota 73 %, Kickxellomycota 71 %, Blastocladiomycota 80 %) can be missing assembly or gene models. A whole pathway of 7 to 17 domains missing together while housekeeping modules are present is hard to explain by random incompleteness, but single-domain calls in these groups are not reliable. Mitochondrial-genome-encoded subunits are not in the protein sets and do not affect the nuclear complex I domains used here.

## 6. Expansions and contractions relative to proteome size

![Fig 4](fable_B_figs/fig4_allometry_scatter.png)

**Table 6. Variance partition by feature type** (medians over families).

{var_summary}

Slopes are mostly well below 1 (median 0.24 for Pfam families): most families do not grow with the proteome. Families with slopes above 1 are the ones that make large proteomes large: ankyrin (1.3 to 1.4), F-box-like (1.41), GDSL lipase (1.33), carboxylesterase (1.24), GMC oxidoreductase (1.23), NACHT (1.18), protein kinase (1.16), CAZy AA (1.46) and CE (1.60). These slopes come from a class-balanced fit and differ from the genome-weighted slopes in atlas view V27 and page S1; the two answer different questions (the kingdom-wide trend across classes here, the trend among sampled genomes there).

![Fig 5](fable_B_figs/fig5_class_residual_heatmap.png)

**Table 7. Families whose residual copy number differs most between classes** (Pfam, median >= 2 copies; SD of class-mean residuals over the 24 classes with >= 5 q genera; residuals in log10 units).

{exp_class}

Observations.

- **NLR and allorecognition parts** (HET, NACHT, NPHP3_N, NB-ARC, TPR_10, Ank) lead the list. Per 1,000 proteins: HET 6.0 to 7.3 in Sordariomycetes, Leotiomycetes and Dothideomycetes, 0.55 in Eurotiomycetes, 0.10 in Agaricomycetes and Pezizomycetes, 0 in all yeasts. Pezizomycetes carry the largest ankyrin load of any class: median 352 Ank_2 proteins (42 elsewhere), *Choiromyces* 1,236 and *Tuber* 1,036. This agrees with atlas page C1 and adds the Pezizomycetes ankyrin outlier. Fungal NLR diversity is reviewed in Dyrka et al. 2014 (Genome Biol Evol 6:3137).
- **"Fung_rhodopsin" (PF20684)**: eta squared 0.91, the most class-structured large family. 4.6 to 11.5 per 1,000 proteins in Pezizomycotina (Lecanoromycetes median 118 proteins, *Parmelia* 257), 0 to 0.13 in every other group. With CFEM this is the PTH11-type receptor repertoire (name-to-family mapping: citation to verify against the Pfam entry).
- **Secondary metabolism**: PKS domains (ketoacyl-synt, PS-DH, acyl transferase) are highest in Lecanoromycetes (median 26 ketoacyl-synt proteins against 7 elsewhere; *Parmelia* 79), then Sordariomycetes, Leotiomycetes and Eurotiomycetes. UstYa (RiPP-associated) is highest in Leotiomycetes. P450 is 7.8 to 11.7 per 1,000 in Pezizomycotina classes and Agaricomycetes and 1.3 to 1.6 in yeasts and Tremellomycetes.
- **Zn2Cys6 and Fungal_trans**: high in Sordariomycetes, Leotiomycetes, Eurotiomycetes and Saccharomycotina (10 per 1,000) and in Tremellomycetes for Fungal_trans (10.2), low in Agaricomycetes (2.2). So the "Pezizomycotina-specific" expansion of atlas H5 is shared with budding yeasts and with Tremellomycetes when normalised.
- **Glomeromycetes are the lowest class for most metabolic families and the highest for a few**: TLDc (+1.55), Sel1 (+1.23), TE-derived domains (+1.05), protein kinase (+0.41). Raw medians: 2,020 Pkinase proteins, 734 Sel1, 622 TLDc. Kinase and Sel1 expansion in *Rhizophagus* is published (Tisserant et al. 2013, PNAS 110:20117); here it holds for all 6 q genera.

**Table 8. Families with the largest lifestyle effect within orders** (partial R2 after order fixed effects; effects are mean within-order residuals in log10 units).

{exp_life}

![Fig 6](fable_B_figs/fig6_taxonomy_vs_lifestyle.png)

![Fig 7](fable_B_figs/fig7_lifestyle_effects.png)

- 1,163 of 2,831 families have a lifestyle effect at BH q < 0.01 within orders; 496 with partial R2 >= 0.05. The largest are plant-cell-wall enzymes: GH10 (0.19), CBM1 (0.18), GH43, AA9 (0.17), GH7, pectinesterase, cellobiose dehydrogenase.
- **In nearly every top family the lowest lifestyle is ectomycorrhizal** and the highest is litter saprotroph. The ECM deficit is not confined to CAZymes: aldo-keto reductase, GST, alcohol dehydrogenase, SDR, methyltransferase and MFS-type efflux families are 0.10 to 0.15 log10 units lower (20 to 30 % fewer copies than order-mates). In Agaricomycetes absolute medians, ECM (44 genera) versus wood saprotroph (100 genera): cellulose-active CAZymes 5 vs 25, AA9 3 vs 15, GH7 0 vs 3.5, SDR (PF00106) 56.5 vs 98.5, GST_N_2 15 vs 28, P450 112.5 vs 149, strict secretome 2.8 % vs 5.1 % of proteins; protein kinases 161 vs 177 (no change); proteome 10,838 vs 12,227; fraction of proteins with any Pfam domain 0.595 vs 0.630. Integrase_H2C2 is the one family that is higher in ECM. The PCWDE loss and transposon gain are established (Kohler et al. 2015, Nat Genet 47:410; Miyauchi et al. 2020, Nat Commun 11:5125); the breadth of the oxidoreductase and detoxification contraction is the part worth testing.
- **Mycoparasites** are higher for chitin-active CAZymes and GH55 and lower for pectin enzymes (Fig 7). **Animal parasites** are lower for plant-cell-wall families.
- **Lichenised** genera show almost no within-order effect. Lichenisation is fixed at order or class level (Lecanoromycetes), so the model has no contrast. The Lecanoromycetes PKS and 7TM expansions (Table 7) may be lichen biology or Lecanoromycetes history; this design cannot say which.
- NLR families have small lifestyle effects (HET 0.02, NACHT 0.04).

## 7. Functional-class shifts

![Fig 8](fable_B_figs/fig8_functional_class_shifts.png)

**Table 9. Copies per 1,000 proteins, median of genus medians (q genera), selected families and classes.** Group sizes are in Table 10.

{selected}

**Table 10. Proteome-level traits by group** (median of genus medians; fractions of the proteome).

{scalar}

Which classes shift most (eta squared by class; partial R2 of lifestyle within order):

| trait | eta2 class | pR2 lifestyle within order | lowest lifestyle (within-order effect) |
|---|---|---|---|
| fraction of proteins with Pfam | 0.66 | 0.108 | ectomycorrhizal (-0.047) |
| log proteome size | 0.59 | 0.018 | dung saprotroph |
| median protein length | 0.54 | 0.060 | ectomycorrhizal (-17 aa) |
| strict secretome fraction | 0.47 | 0.087 | ectomycorrhizal (-0.0096) |
| median disorder fraction | 0.45 | 0.033 | wood saprotroph (ECM highest, +0.012) |
| TM (>= 2 helices) fraction | 0.40 | 0.096 | ectomycorrhizal (-0.011) |
| GPI-anchored fraction | 0.39 | 0.072 | ectomycorrhizal (-0.0017) |

- **CAZy.** Substrate classes for cellulose, xylan and pectin vary most between classes (SD of class means 0.39 to 0.45 log10) and have the largest lifestyle effects (0.14 to 0.17). Per 1,000 proteins, cellulose-active enzymes are 1.8 to 3.7 in Sordariomycetes, Dothideomycetes, Leotiomycetes, Pezizomycetes, Orbiliomycetes and Agaricomycetes, 0.3 in Eurotiomycetes and Lecanoromycetes, and 0 in Saccharomycotina, Tremellomycetes and Ustilaginomycotina. Lignin-active (AA) enzymes peak in Agaricomycetes (2.2). GT is the only class enriched in yeasts (12.3 per 1,000 in Saccharomycotina against 5.0 in Agaricomycetes), because the cell-wall and glycosylation GT set is fixed while the proteome shrinks. PL is absent (median 0) from Lecanoromycetes, Saccharomycotina and Taphrinomycotina.
- **Peptidases.** Sedolisin S53 varies most (Leotiomycetes +0.74; 0 in Saccharomycotina, Tremellomycetes, Pucciniomycotina and early-diverging groups). Subtilisin S08A is highest in Orbiliomycetes (2.0 per 1,000; nematode-trapping fungi), Zoopagomycota s.l. and Sordariomycetes. Lifestyle effects on peptidase families are small (median partial R2 0.028; largest S33 0.12 and S09X 0.11, both lower in ECM).
- **Transcription factors.** Zn2Cys6 and its partner domain separate Dikarya classes; the generic families (C2H2, bHLH, GATA, homeodomain, HMG, bZIP, HSF) are all highest in Mucoromycetes after normalisation (+0.4 to +0.5 log10) and per 1,000 proteins in Mucoromycota s.l. (for example GATA 2.1 against 0.5 in Sordariomycetes). Mucoromycetes have a known whole-genome-duplication history, which may explain retained TF paralogs (the atlas `wgd_ks_summary` could test this).
- **Secretome.** Highest in Orbiliomycetes (9.6 %), Zoopagomycota s.l. (8.1 %), Sordariomycetes and Dothideomycetes (7.0 %); lowest in Saccharomycotina (2.5 %), Taphrinomycotina (2.6 %) and Microsporidia (2.8 %). This agrees with atlas H4.
- **Disorder.** Highest in Ustilaginomycotina (0.48), Pucciniomycotina (0.46), Pezizomycetes (0.45) and Tremellomycetes (0.44); lowest in Microsporidia (0.32) and Saccharomycotina (0.33). The basidiomycete yeasts are therefore not low-disorder, consistent with atlas P2.
- **Fraction of proteins with any Pfam domain** runs from 0.86 (Saccharomycotina) to 0.59 to 0.61 (Pezizomycetes, Agaricomycetes) and 0.44 (Microsporidia). Every per-proteome normalisation in this report inherits this: large proteomes contain more domain-less models, real or spurious.

## 8. Patterns worth deeper investigation

Ranked by my judgment of novelty times tractability. Each lists the observation, alternative explanations, and a concrete design. The genus tree (`atlas_data/phylo`, FastTree on fungi_odb10 markers) is required for most designs. R packages already on the system: `phylolm`, `ape`, `phytools`, `corHMM`.

**P1. A Dikarya-stem loss of animal-type regulatory and repair complexes.**
- Observation: 203 ancestral domains in <= 5 % of 1,083 Dikarya genera; 191 retained in >= 50 % of 33 Mucoromycota s.l. genera. Includes Integrator, Fanconi core, KICSTOR, B12 metabolism, SMG1, DNA-PKcs (Fig 9).
- Alternatives: (a) divergence below HMM threshold in Dikarya; (b) the losses are spread over several nodes (Zoopagomycota s.l. already lack some, for example SMG1 and LETM1_C at 0 %); (c) the "ancestral" call is wrong and these were acquired in Mucoromycota; unlikely for domains shared with animals.
- Design: add non-fungal opisthokont outgroups (animal, choanoflagellate, nucleariid proteomes) to the matrix. Reconstruct gains and losses on the genus tree with Dollo parsimony and with the phylogenetic birth-death model in Count (Csűrös 2010, Bioinformatics 26:1910) and report the number of losses per branch. Test whether the Dikarya stem carries more losses than expected from its length. For 20 focal complexes run HHsearch or Foldseek on three Dikarya proteomes to rule out (a). Tables: `pfam`, `swissprot_best` for ortholog confirmation.

**P2. A common pool of dispensable genes: convergent loss in every reduced lineage.**
- Observation: loss sets of independent reduced lineages overlap at log2 1.4 to 3.5 over expectation; 114 core domains lost in >= 10 of 37 lineages, 1,072 never lost.
- Alternatives: (a) shared annotation failure for short or fast-evolving proteins; (b) low-BUSCO lineages inflate overlap; (c) nested lineages counted as independent (the four Saccharomycotina classes share ancestry).
- Design: on the genus tree, map each core domain's losses by Dollo parsimony and count independent loss events per domain. Model loss count per domain with a negative-binomial regression on protein properties (length, disorder, copy number, essentiality of the *S. cerevisiae* or *S. pombe* ortholog from `swissprot_best`, pathway membership). Validate absences with tblastn against assemblies for a sample of 200 domain-lineage pairs. Gene loss as an adaptive or neutral force is reviewed in Albalat and Cañestro 2016 (Nat Rev Genet 17:379).

**P3. Lipid auxotrophy by FAS loss in at least four host-dependent lineages.**
- Observation: Table 5. 0 of 21, 7, 31 and 93 species in *Malassezia*, *Pneumocystis*, Glomeromycota, Microsporidia; 3 of 11 in Zoopagomycota.
- Alternatives: the 8-domain FAS module is two giant genes, so one unassembled locus removes all domains at once (relevant for Zoopagomycota, BUSCO 73 %). Not an issue for the four main lineages, where every species agrees.
- Design: phylogenetic logistic regression of FAS presence on obligate host association (Ives and Garland 2010, Syst Biol 59:9; `phylolm::phyloglm`), with the caveat that four origins give little power (Maddison and FitzJohn 2015). More informative: ask what is gained or retained alongside the loss (lipases, fatty-acid transporters, CD36-like or FAT domains) by comparing each lineage with its sister. Tables: `pfam`, `merops` (lipase S09/S33), `cazy_overview`.

**P4. Loss of inorganic sulfur and nitrogen assimilation and of the molybdenum cofactor.**
- Observation: MoCo biosynthesis absent in 9 or more lineages; ATP-sulfurylase absent in Erysiphaceae, *Puccinia*, Zoopagomycota, *Pneumocystis*, *Batrachochytrium*, Neocallimastigomycota, Microsporidia; urease absent in Saccharomycotina and several parasites.
- Alternatives: MoCo loss follows loss of its client enzymes (nitrate reductase) in any lineage with reduced nitrogen sources, not parasitism as such (it is also lost in free-living Lipomycetes, Saccharomycetaceae and Tremellomycetes). The order of loss cannot be read from presence/absence.
- Design: test correlated evolution of MoCo synthesis and each Mo-enzyme with Pagel's discrete model (Pagel 1994, Proc R Soc B 255:37; `phytools::fitPagel` or `corHMM`), including the dependent-order models (does enzyme loss precede cofactor loss?). Cross with atlas hypothesis T2 (inorganic-nutrient transporters in obligate biotrophs).

**P5. Ectomycorrhizal contraction extends to oxidoreductase and detoxification families.**
- Observation: within orders, ECM genera have 20 to 30 % fewer SDR, GST, aldo-keto reductase, ADH and methyltransferase proteins, on top of the known PCWDE loss; kinases unchanged; 51 ECM q genera.
- Alternatives: (a) ECM genomes are transposon-rich and harder to annotate, lowering the Pfam-annotated fraction (0.595 vs 0.630) and biasing every family; kinases being unchanged argues against a uniform bias; (b) few independent ECM origins within orders (Agaricales, Boletales, Russulales, Pezizales dominate); (c) FungalTraits labels are genus-level and some genera are mixed.
- Design: PGLS of log copy number on ECM status with log proteome size as covariate (Grafen 1989, Phil Trans R Soc B 326:119; Freckleton, Harvey and Pagel 2002, Am Nat 160:712; `phylolm`, Ho and Ané 2014, Syst Biol 63:397), per family, BH-corrected. Then sister-pair contrasts for each independent ECM origin using `ecm_lineage` in `genome_labels_ext`, with a sign test across origins (as atlas V14 did for yeasts). Check whether the contracted oxidoreductases are secreted or cytosolic (`signalp_signal_peptide`).

**P6. The Pezizomycotina recognition and sensing complex: NLR, HET, ankyrin, PTH11-type 7TM, CFEM.**
- Observation: eta squared by class 0.49 to 0.91; lifestyle partial R2 <= 0.08; Pezizomycetes ankyrin median 352 (elsewhere 42); Lecanoromycetes 7TM median 118.
- Alternatives: (a) repeat domains inflate protein counts through fragmented models in repeat-rich genomes (*Tuber* is 58 % repeats; Martin et al. 2010, Nature 464:1033); (b) TE-associated ankyrin proteins; (c) the families expand together only because all scale steeply with proteome size (slopes 1.1 to 1.4).
- Design: count domain architectures instead of domains using `A3/protein_pfam_arch.parquet` (NACHT-Ank, HET-NACHT-WD40, 7TM-CFEM) and exclude proteins co-occurring with TE domains. Model family size with a birth-death process along the tree (CAFE 5; Mendes et al. 2020, Bioinformatics 36:5516) to locate the branches with rate shifts. Test co-expansion of 7TM and NLR families by PGLS on residuals.

**P7. Glomeromycota: signalling expansion with metabolic loss.**
- Observation: 2,020 Pkinase, 734 Sel1, 622 TLDc proteins (genus medians, 6 q genera); RNF213-type AAA 15 to 19 copies and GVIN/URGCP GTPases in all genera and almost nowhere else; FAS and catalase absent.
- Alternatives: transposon-linked gene duplication and allelic copies in dikaryotic or heterokaryotic assemblies inflate counts (BUSCO duplication should be checked per genome; `wgd_ks_summary` identical-copy statistics); only 6 q genera.
- Design: use `wgd_ks_nodes` to date the kinase and Sel1 paralogs (dS distribution) and compare with TE-domain load. Architecture analysis of the kinases (TKL with or without receptor domains). Compare with the other obligate plant biotrophs (Pucciniomycetes, Erysiphaceae) to separate "obligate biotroph" from "Glomeromycota".

**P8. Large Basidiomycota-specific families of unknown function.**
- Observation: DUF6534 in 97.9 % of 291 Basidiomycota genera, median 50.5 copies in Agaricomycetes (IQR 30 to 74; *Gerronema* 211, *Mycena* 199, *Armillaria* 183), 2 genera outside; DUF6533 similar; Clampless-1 domains in 94 %.
- Alternatives: a fast-evolving family whose Ascomycota homologs fall below threshold; counts inflated by fragmented models.
- Design: architecture, TM and signal-peptide profile of DUF6533/6534 proteins (`tmhmm`, `signalp_signal_peptide`), genomic clustering (`gene_info` coordinates), copy number against fruit-body and lifestyle traits by PGLS within Agaricomycetes (222 q genera gives real power). Structure-based homology search for remote relatives.

**P9. Timing of fungal cell-biology innovations after the zoosporic grade.**
- Observation: modules such as the DASH kinetochore complex (8 domains) and the SUR7/KAR9/ATG11 set (19 domains) are absent from Chytridiomycetes (35 genera), Blastocladiomycota and Neocallimastigomycota but present in Mucoromycota and Dikarya; the Ascomycota and Dikarya restricted sets (71 and 17 domains) add later layers.
- Alternatives: homology-detection failure with distance from *S. cerevisiae*, the species these families were defined in (Weisman et al. 2020). This is the main risk for every "innovation" call here.
- Design: stratigraphic assignment of each domain to its earliest node by Dollo parsimony on the tree, with a per-family detection-failure simulation (absence expected from evolutionary rate and protein length, as in Weisman et al. 2020) to flag families whose apparent origin is explained by divergence alone.

**P10. Lichen secondary metabolism versus Lecanoromycetes history.**
- Observation: Lecanoromycetes carry the most PKS domains per proteome (median 26 ketoacyl-synt proteins, 74 q genera) and the most PTH11-type 7TM proteins; within-order lifestyle contrast is nil because lichenisation is fixed at high rank.
- Alternatives: phylogenetic inertia; metagenome-derived lichen assemblies (BUSCO median 89.7 %, only 74 of 164 genera pass the quality filters) can include photobiont or other fungal contigs.
- Design: use the independent lichenised lineages outside Lecanoromycetes (Arthoniomycetes, Lichinomycetes, Candelariomycetes, Coniocybomycetes, lichenised Eurotiomycetes and Dothideomycetes; `primary_photobiont` in `genome_labels_ext`) and the secondarily non-lichenised Lecanoromycetes as natural contrasts in a phylogenetic logistic regression or PGLS. Screen contigs for taxonomic consistency first.

## 9. Limits of this survey

- No phylogenetic correction; class weights and rank fixed effects are stand-ins (Felsenstein 1985; Uyeda, Zenil-Ferguson and Pennell 2018, Syst Biol 67:1091).
- Presence/absence comes from predicted proteins. Unannotated genes, fragmented models and contamination all enter. Loss calls in lineages with BUSCO below 80 % are provisional.
- Pfam coverage is biased toward model organisms; counts of restricted domains measure curation as much as innovation.
- Ancestral and core sets are defined by breadth, without outgroups.
- FungalTraits lifestyle is one label per genus, and 119 of 1,173 q genera have none.
- Thresholds (50 % for genus presence, 0.8/0.02 for restriction, 10 % for loss, 0.35 linkage cut) were set once and not tuned; I did not run a sensitivity grid. The specific counts would move with them; the named modules and the top-ranked families are well separated from the thresholds (for example FAS at 0 % against 93 to 100 %).
- Literature statements marked "citation to verify" were written from memory and not checked online. Other citations are ones I am confident of in author, year and journal; volume and page numbers should still be checked before publication.
