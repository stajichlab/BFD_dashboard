# Hypotheses: micro-to-macro shifts in fungal organization

Date: 2026-11-08. Author: Claude (design). Build status: G1, G3, P2, S1, C1 and C3 are built (PI approval 2026-11-08; pages `HS_*.md`, code `bin/atlas/hs_labels.py`, `bin/atlas/sql/J9_hs.sql`, `bin/atlas/render_hs.py`, jobs `workflows/atlas/J9_hs.sh`, `J10_hs_render.sh`; intermediates `atlas_data/A11/`, labels `atlas_data/labels/genome_labels_ext.parquet`). A "Results" link sits under each built hypothesis. The other hypotheses are design only.
Scope: new, testable hypotheses that connect shifts in organism size and organization to genome and proteome traits in BFD. They build on atlas pages H1-H7 and do not repeat views that exist.

Every number in this file is **measured** on the login node (`/usr/bin/python3.12`, DuckDB, memory limit 500 MB) from `atlas_data/labels/genome_labels.parquet`, the small A2/A4/A5 per-genome tables, and a filtered scan of `atlas_data/A7/pfam_genome_matrix.parquet`, unless marked **not measured** or **pilot**. "Pilot" means a quick screen used to check that a signal is worth a view. Pilot values are not results.

[Atlas index](index.md) · [PLAN.md](PLAN.md)

**Decision (PI, 2026-11-08): Microsporidia can be skipped from some analyses.** Default rule applied to every hypothesis below:

- **Exclude Microsporidia** from kingdom-wide scaling, density and distribution views (S1, S3, S4, C1, C2, C3, P2, P3, G2, G3, G4, T1, T3), from the early-diverging outgroup sets used for polarity and ancestral-state reconstruction (their long branches and extreme reduction distort both), and from the H2/H2b re-framing (Theme F).
- **Keep Microsporidia only where the hypothesis is about reduction or intracellular parasitism**, as one replicate among several: T2 (transporter loss in obligate biotrophs), S2 (signaling-sensor loss in parasites), P1 (protein shortening in reduced genomes) and G1 (genome size under obligate parasitism). Use only genomes < 10 Mb (removes the 34 metagenome-derived bins, section 1.4). Report each result with and without Microsporidia.
- Microsporidia rows stay in the descriptive tables of section 1, for reference.

---

## 1. Contrast groups and how many genera BFD has for each

Counts below use `is_representative AND q_ann AND q_busco_rel` (one genome per species, class-relative BUSCO filter) unless stated. Raw counts (no filter) are given where the filter removes a group.

### 1.1 Groups from existing labels (growth form, taxonomy)

| Group | Definition in BFD | Genera (filtered) | Genomes (filtered) | Note |
|---|---|---|---|---|
| Yeast, Saccharomycotina | `yeast_origin = 'Saccharomycotina'` | 144 raw | 6,019 raw | 118 genera used in H2 |
| Yeast, Taphrinomycotina | `yeast_origin = 'Taphrinomycotina'` | 3 raw | 80 raw | *Schizosaccharomyces* (5 filtered genomes), *Taphrina* (9), *Pneumocystis* (6) |
| Yeast, Tremellomycetes | `yeast_origin = 'Tremellomycetes'` | 15 | | 9 dimorphic genera also in the class |
| Yeast, Pucciniomycotina | `yeast_origin = 'Pucciniomycotina-yeast'` | 14 raw | 281 raw | Microbotryomycetes 5, Cystobasidiomycetes 2 filtered |
| Yeast, Ustilaginomycotina | `yeast_origin = 'Ustilaginomycotina-yeast'` | 6 raw | 101 raw | *Malassezia* 19 genomes filtered |
| Dimorphic, Pezizomycotina | `yeast_origin = 'Pezizomycotina-dimorphic'` | 21 raw | 474 raw | Eurotiomycetes: 15 dimorphic vs 32 filamentous genera filtered |
| Dimorphic, Ustilaginomycotina | `yeast_origin = 'Ustilaginomycotina-dimorphic'` | 23 raw | 199 raw | Ustilaginomycetes 11, Exobasidiomycetes 9 |
| Filamentous, not Agaricomycetes | `growth_form = 'filamentous_mycelium'`, CLASS <> Agaricomycetes | 584 | 2,638 | Mostly Pezizomycotina |
| Agaricomycetes | CLASS | 291 | 845 | |
| Lichenized (thallus) | `growth_form = 'thallus_photosynthetic'` | 156 | 327 | Lecanoromycetes 122, Arthoniomycetes 7, Lichinomycetes 4, Candelariomycetes 3, Dothideomycetes 3, Eurotiomycetes 2, Coniocybomycetes 2 |
| Zoosporic (flagellated) | `growth_form LIKE 'zoosporic%'` | 43 | 63 | Chytridiomycetes 29, Neocallimastigomycetes 5, Blastocladiomycetes 4, Aphelidea 2, others 1 each |
| Mucoromycota s.l. (coenocytic) | PHYLUM in Mucoromycota, Mortierellomycota, Glomeromycota | 59 | 169 | Mucoromycota 35, Mortierellomycota 12, Glomeromycota 11 |
| Zoopagomycota s.l. | PHYLUM in Kickxellomycota, Zoopagomycota, Entomophthoromycota, Basidiobolomycota | 28 | 160 | *Coemansia* alone is 112 of the 160 genomes |
| Microsporidia | PHYLUM | 23 | 73 | See the quality warning in 1.4 |

### 1.2 Groups from FungalTraits 1.2 columns not yet used by the atlas

FungalTraits 1.2 (`/bigdata/stajichlab/shared/projects/BFD/DeltaGain_Fungi/refs/lifestyle/FungalTraits_1.2.csv`, latin-1, 10,771 rows) has 24 columns. The atlas uses `Growth_form_template`, `primary_lifestyle`, `Secondary_lifestyle`, `Plant_pathogenic_capacity_template`, `Decay_substrate_template`, `Animal_biotrophic_capacity_template`. The columns below are on disk and not yet in `genome_labels`:

| Column | Values (whole table) | Use here |
|---|---|---|
| `Fruitbody_type_template` | perithecium 2,114; apothecium 1,561; none 1,249; agaricoid 533; corticioid 467; polyporoid 238; zoosporangium 191; rust 184; other 148; gasteroid-hypogeous 127; gasteroid 123; cleistothecium 108; smut 107; clavarioid 78; thyrothecium 51; cyphelloid 46; tremelloid 39; mazaedium 32; chasmothecium 28; hysterothecium 18; phalloid 11; zygosporangium 9; clathroid 7; blank 3,294 | Fruiting-body complexity (C2, S1, S2, P3) |
| `Hymenium_type_template` | closed 2,770; smooth 2,045; none 1,255; gills 425; poroid 382; other 348; hydnoid 98; gel 28; blank 3,420 | Fruiting-body complexity (C2) |
| `Ectomycorrhiza_exploration_type_template` | long-distance 126; short-distance_delicate 37; short-distance_coarse 34; mat 33; unknown 32; contact 31; medium-distance_smooth 21; medium-distance_fringe 18; blank 10,439 | Rhizomorph proxy (T3) |
| `Ectomycorrhiza_lineage_template` | 70 values | Independent ECM origins (T3) |
| `primary_photobiont` | chlorococcoid 617; trentepohlioid 230; cyanobacterial 111 | Lichen C/N source (T3) |
| `Decay_type_template` | white_rot 552; soft_rot 147; brown_rot 98; mold 94; keratinolytic 40; chitinolytic 11 | Covariate for T1 |
| `Endophytic_interaction_capability_template` | foliar_endophyte 307; root_endophyte 77; root-associated 63; others < 25 | Not used here |
| `Aquatic_habitat_template` | non-aquatic 7,180; partly_aquatic 410; freshwater 227; marine 191 | Not used here |
| `Comment_on_lifestyle_template` | 157 values; includes `obligate_plant_pathogen` and `obligate_animal_parasite` | Obligate biotrophy (T2, G1) |

Measured genus counts in BFD after the filter:

- Agaricomycetes by fruit body: agaricoid 144 genera (9 orders), polyporoid 45 (6), corticioid 49 (10), gasteroid 9, clavarioid 10, gasteroid-hypogeous 7 (4 orders), tremelloid 2, unlabelled 19.
- Agaricomycetes by hymenium: gills 100 genera (8 orders), poroid 94 (6), smooth 46 (13), closed 16 (4), hydnoid 10 (4).
- Ascomycota by fruit body: perithecium: Sordariomycetes 117, Dothideomycetes 94; apothecium: Lecanoromycetes 122, Leotiomycetes 39, Pezizomycetes 11, Arthoniomycetes 7, Orbiliomycetes 5; cleistothecium: Eurotiomycetes 10, Sordariomycetes 8; none: Sordariomycetes 41, Pichiomycetes 36, Saccharomycetes 19, Dipodascomycetes 13; gasteroid-hypogeous: Pezizomycetes 5.
- Ectomycorrhizal genera inside orders that also hold non-ECM genera: Boletales 45 vs 7, Agaricales 9 vs 95, Pezizales 7 vs 9, Russulales 5 vs 11, Cantharellales 3 vs 5, Gomphales 3 vs 2, Thelephorales 3 vs 0. Exploration type: long-distance 45, contact 6, short-distance delicate 6, medium smooth 6, mat 6, short coarse 4, medium fringe 4.
- Lichenized by photobiont: chlorococcoid 124 genera, trentepohlioid 13, cyanobacterial 9.
- Obligate biotrophs (raw genus counts, filter in brackets): Pucciniales 8 genera, 131 genomes (5 genera, 18 genomes); Erysiphaceae (powdery mildews; BFD places them in ORDER_ = Helotiales) 7 genera, 45 genomes (2 genera, 3 genomes: the Leotiomycetes class threshold removes them); Glomeromycota 15 genera, 60 genomes (11, 22); Taphrinales 3 genera (2); Pneumocystidales 1 genus, 11 genomes (6 genomes); Entomophthorales 6 genera (6); Neocallimastigales 5 (5); Microsporidia 37 genera (23).

### 1.3 Labels that do not exist on disk

- Sclerotium formation: no column in FungalTraits 1.2. Needs a hand-made list.
- Rhizomorph formation: no column. The ECM exploration type (long-distance = rhizomorphic, to check) is a partial proxy and covers only ECM genera.
- Septate vs coenocytic hyphae: no column. Taxonomy proxy only (Mucoromycota s.l. and Zoopagomycota s.l. = coenocytic; Dikarya = septate).
- Flagellum: `Growth_form_template` zoosporic values only. The number of independent flagellum losses is a literature question (to check).
- Mating system (bipolar / tetrapolar) for S4: no column. Needs a hand-made list.
- Intracellular vs extracellular (haustorial) obligate biotrophy for T2 and G1: no column. Needs a hand-made list of about 15 genera.
- Hypogeous fruiting: only `gasteroid-hypogeous`; sequestrate but epigeous forms are not separated.

### 1.4 Quality warnings specific to these groups

- **Microsporidia**: 163 genomes; 129 are < 10 Mb. The rest (up to 56 Mb, up to 16,855 proteins, BUSCO fungi_odb12 about 10 %) are metagenome-derived bins (`gmOTU` names) and likely hold host or co-assembled sequence. `q_busco_rel` for Microsporidia is at phylum level with a 25th percentile near 10 %, so it removes almost nothing. Every Microsporidia view must add `total_length_bp < 10e6` (32 genera remain, raw) and should state that fungi_odb12 is not a valid completeness measure for them (a `microsporidia_odb10` set may exist in `funannotate_db`; to check).
- **Powdery mildews**: `q_busco_rel` removes 42 of 45 genomes because the Leotiomycetes threshold is high. Use a family-level relaxed filter (BUSCO >= 90 %) for T2 and G1 and say so.
- **Pucciniales and Glomeromycota**: BUSCO duplication is high (dikaryotic / multinucleate assemblies). `q_dup` removes most of them. Views on these groups must not use `q_dup`, and must treat copy-number traits with care.
- **Zoopagomycota s.l.**: *Coemansia* is 112 of 160 filtered genomes. Genus medians are mandatory.
- **Lichens**: BUSCO median 93.2 in the filtered set (metagenome-derived assemblies). Report `busco_complete` as a covariate.

### 1.5 Transporter annotation: what exists

- **TCDB is not in the production tables.** `Fungi_BFD_runs/tables/` has no `tcdb.parquet`, and `BFD_function/results/function/` has no `tcdb/` directory.
- A one-genome test exists: `Fungi_BFD_runs/BFD_function_tcdb_tf_cgc_test/` (run 2026-11-06, *Aaosphaeria arxii* CBS 175.79, LOCUSTAG FF5840CF). It produced `tables/tcdb.parquet` with 11,251 rows for 2,821 proteins (columns: species_prefix, protein_id, tc_number, percent_identity, aln_length, mismatches, gap_openings, q_start, q_end, s_start, s_end, evalue, bitscore; up to 5 hits per protein). The module is `Fungi_BFD/nextflow/modules/BFD/TCDB` (blastp against `lib/tcdb/tcdb.fasta`, 24,281 sequences, `-max_target_seqs 5`, e-value from `params.tcdb_evalue`). `MERGE_TCDB` and `SETUP_TCDB_DB` exist. The `cazy_cgc.parquet` from the same test has 0 rows.
- There is no TC-family description or substrate table on disk; `tc_number` (for example `2.A.1.1.x`) gives class, subclass and family by position, which is enough for family-level counts. Descriptions need a download from tcdb.org (not done).
- A full TCDB run is about 23,400 blastp jobs. Runtime per genome is **not measured**.
- Until TCDB runs, transporter traits use Pfam families from A7 (`pfam_genome_matrix`: genome, pfam_acc, pfam_id, n_proteins, n_domains; 15,610 distinct families) and SwissProt best hits (A3 `protein_core.sp_id` joined to `swissprot_subject.go_terms`, which holds GO ids such as `GO:0055085` transmembrane transport). Family-level Pfam proxies are listed per hypothesis and marked "to check" for specificity.

### 1.6 Phylogeny status

`atlas_data/phylo/` has no tree file. The nf_phyling run log ends at 17:39 on 2026-11-08 without a result; `pep_fa/` holds 1,217 staged proteomes (20:36). PGLS is not possible yet. Every hypothesis below lists a no-tree control.

### 1.7 Pilot screen used to pick families (A7, representative genomes, Q-ann + Q-busco-rel; median copies per 1,000 proteins per genome; genome level, not genus level)

| Pfam | Name | Yeast (n=1,109) | Dimorphic (207) | Filamentous non-Agarico. (2,638) | Agaricomycetes (845) | Lichen (327) | Microsporidia (73) | Zoosporic (63) | Mucoromycota s.l. (169) | Zoopagomycota s.l. (160) |
|---|---|---|---|---|---|---|---|---|---|---|
| PF07690 | MFS_1 | 19.9 | 24.9 | 28.8 | 12.4 | 17.6 | 1.2 | 5.1 | 10.5 | 10.3 |
| PF00083 | Sugar_tr | 10.0 | 12.4 | 14.5 | 5.5 | 8.9 | 0.5 | 2.2 | 4.9 | 5.1 |
| PF00324 | AA_permease | 4.8 | 2.7 | 3.8 | 1.5 | 3.3 | 0.0 | 1.0 | 1.8 | 1.1 |
| PF00153 | Mito_carr | 5.9 | 4.4 | 3.3 | 2.9 | 3.5 | 0.0 | 4.0 | 4.8 | 4.8 |
| PF00005 | ABC_tran | 5.0 | 4.6 | 4.4 | 3.5 | 3.4 | 4.2 | 7.2 | 4.0 | 7.8 |
| PF00909 | Ammonium_transp | 0.54 | 0.32 | 0.31 | 0.19 | 0.31 | 0.00 | 0.25 | 0.26 | 0.30 |
| PF00069 | Pkinase | 17.0 | 13.3 | 10.4 | 13.8 | 13.0 | 11.1 | 17.6 | 19.8 | 13.5 |
| PF00071 | Ras | 5.3 | 4.1 | 2.9 | 2.7 | 3.0 | 4.5 | 4.2 | 5.2 | 4.6 |
| PF00512 | HisKA | 0.56 | 1.07 | 1.00 | 0.28 | 0.81 | 0.00 | 0.43 | 0.90 | 0.16 |
| PF00989 | PAS (fraction of genomes with >= 1) | 0.59 | 1.00 | 0.98 | 0.78 | 1.00 | 0.00 | 0.95 | 1.00 | 0.19 |
| PF00360 | PHY phytochrome (fraction present) | 0.14 | 0.91 | 0.91 | 0.62 | 0.98 | 0.00 | 0.11 | 0.05 | 0.08 |
| PF01036 | Bac_rhodopsin (fraction present) | 0.56 | 1.00 | 0.88 | 0.42 | 1.00 | 0.00 | 0.54 | 0.01 | 0.00 |
| PF00503 | G-alpha | 0.36 | 0.42 | 0.31 | 1.21 | 0.49 | 0.00 | 0.82 | 0.90 | 0.30 |
| PF02076 | STE3 pheromone receptor | 0.18 | 0.11 | 0.08 | 0.44 | 0.10 | 0.00 | 0.00 | 0.00 | 0.00 |
| PF00211 | Guanylate_cyc | 0.18 | 0.14 | 0.09 | 0.08 | 0.10 | 0.00 | 1.13 | 0.17 | 0.15 |
| PF00027 | cNMP_binding | 0.50 | 0.41 | 0.34 | 0.23 | 0.40 | 0.40 | 1.62 | 0.50 | 0.44 |
| PF00520 | Ion_trans | 0.32 | 0.29 | 0.27 | 0.28 | 0.26 | 0.00 | 1.23 | 0.81 | 0.48 |
| PF03028 | Dynein_heavy | 0.17 | 0.13 | 0.08 | 0.08 | 0.10 | 0.00 | 0.98 | 0.08 | 0.15 |
| PF06985 | HET | 0.00 | 0.00 | 4.69 | 0.23 | 3.71 | 0.00 | 0.00 | 0.00 | 0.00 |
| PF05729 | NACHT | 0.00 | 0.00 | 1.47 | 1.91 | 1.83 | 0.00 | 0.07 | 0.08 | 0.00 |
| PF01185 | Hydrophobin (fraction present / per 1,000) | 0.01 / 0 | 0.37 / 0 | 0.52 / 0.06 | 0.94 / 1.29 | 0.12 / 0 | 0 / 0 | 0 / 0 | 0.03 / 0 | 0.01 / 0 |
| PF07249 | Cerato-platanin (fraction present) | 0.00 | 0.05 | 0.93 | 0.92 | 0.04 | 0.00 | 0.17 | 0.02 | 0.18 |
| PF05730 | CFEM | 0.36 | 0.46 | 1.08 | 0.30 | 0.82 | 0.00 | 0.00 | 0.00 | 0.00 |

Pilot genus-median proteome traits (same filter; group median of genus medians):

| Group | Genera | Genome Mb | Proteins | Median protein aa | Fraction < 100 aa | Introns per transcript | Intergenic bp | Fraction with Pfam | Multi-domain fraction | TM proteins / proteome | Strict secretome / proteome | GPI / proteome | Median IDP fraction | TF / proteome |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Filamentous non-Agarico. | 584 | 37.9 | 11,177 | 405 | 0.037 | 1.79 | 911 | 0.73 | 0.34 | 0.196 | 0.065 | 0.0151 | 0.401 | 0.0234 |
| Agaricomycetes | 291 | 44.4 | 12,562 | 380 | 0.040 | 4.91 | 687 | 0.60 | 0.27 | 0.172 | 0.047 | 0.0117 | 0.414 | 0.0123 |
| Lichen | 156 | 30.4 | 9,562 | 402 | 0.032 | 2.01 | 811 | 0.69 | 0.32 | 0.193 | 0.057 | 0.0124 | 0.418 | 0.0167 |
| Yeast (all origins) | 153 | 13.1 | 5,727 | 413 | 0.034 | 0.24 | 408 | 0.85 | 0.40 | 0.192 | 0.026 | 0.0158 | 0.335 | 0.0243 |
| Mucoromycota s.l. | 59 | 37.7 | 12,231 | 364 | 0.065 | 3.36 | 677 | 0.67 | 0.34 | 0.178 | 0.035 | 0.0103 | 0.384 | 0.0252 |
| Dimorphic | 51 | 20.9 | 7,473 | 436 | 0.029 | 1.56 | 591 | 0.75 | 0.36 | 0.192 | 0.042 | 0.0137 | 0.442 | 0.0221 |
| Zoosporic | 43 | 31.5 | 10,299 | 398 | 0.038 | 3.57 | 506 | 0.62 | 0.32 | 0.210 | 0.047 | 0.0184 | 0.378 | 0.0130 |
| Zoopagomycota s.l. | 28 | 25.7 | 8,692 | 371 | 0.063 | 1.90 | 781 | 0.69 | 0.33 | 0.181 | 0.075 | 0.0166 | 0.396 | 0.0181 |
| Microsporidia | 24 | 8.5 | 3,185 | 251 | 0.125 | 0.21 | 488 | 0.44 | 0.18 | 0.151 | 0.028 | 0.0057 | 0.316 | 0.0092 |

Obligate biotrophs with large genomes (genus medians, same filter): Pucciniomycetes 5 genera, 166 Mb, 18,356 proteins, median protein 302 aa; Glomeromycetes 8 genera, 156 Mb, 28,571 proteins, median protein 227 aa. These short medians are a warning (see P1).

---

## 1.8 Evolutionary framing: shared elements

Every hypothesis below carries an "Evolutionary framing" block with six parts: polarity, placement, mode, methods, order of events, caveats. This section holds what is common to all of them.

**Outgroups available in BFD for polarizing changes** (genera after `is_representative AND q_ann AND q_busco_rel`; raw in brackets where different):

| Shift to polarize | Derived lineage | Sister / outgroups in BFD | Genera |
|---|---|---|---|
| Saccharomycotina yeast | Saccharomycotina (144 raw) | Pezizomycotina (about 1,000 raw); Taphrinomycotina: *Schizosaccharomyces*, *Taphrina*, *Pneumocystis* | 3 Taphrinomycotina genera (5, 9, 6 filtered genomes); Basidiomycota 432 raw genera |
| Basidiomycete yeasts | Tremellomycetes 15, Pucciniomycotina 7, Ustilaginomycotina 6 + 20 dimorphic | Agaricomycetes 291; Pucciniomycetes 5; Ustilaginomycetes 11 | as listed |
| Dikarya-level shifts (septation, NLRs, hydrophobins, loss of flagellum) | Dikarya | Mucoromycota 35, Mortierellomycota 12, Glomeromycota 11; Kickxellomycota 14, Zoopagomycota 7, Entomophthoromycota 6, Basidiobolomycota 1; Chytridiomycota 32, Blastocladiomycota 5, Neocallimastigomycota 5, Monoblepharomycota 2; Aphelidiomycota 2, Rozellomycota 1; Microsporidia 23 | 156 non-Dikarya genera in total |
| Complex fruit bodies | Agaricomycetes fruit-body types | Corticioid Agaricomycetes 49; Tremellomycetes; Dacrymycetes (count not measured); Pucciniomycotina | as listed |
| Obligate biotrophy | each lineage in T2 | nearest free-living sister in BFD (listed in T2) | as listed |

Non-fungal outgroups (animals, Nucleariida, Oomycota) are not in BFD. A Dikarya-wide polarity therefore rests on the early-diverging fungal lineages, and the root position among them is taken from the literature (to check; Rozellomycota + Microsporidia + Aphelidiomycota as the earliest branch).

**The genus tree.** `atlas_data/phylo/genus_representatives.tsv` lists 1,654 genus representatives (one genome per genus; fungi_odb10 markers; phyling alignment; FastTree). Not finished (section 1.6). Its branch lengths will be amino-acid substitutions per site, not time. Consequences: (i) PGLS, phylogenetic signal (lambda, K), phylogenetic logistic regression, ancestral state reconstruction and stochastic character mapping run on this tree without change, because they use relative branch lengths. (ii) Rate-shift and OU models (bayou, OUwie, l1ou) estimate rates per unit of molecular change, not per million years; rate comparisons between clades then mix evolutionary rate with substitution-rate differences. Yeasts and Microsporidia have fast molecular rates (to check), so their branches are long and rate estimates on them are biased. A dated tree (penalized likelihood with `ape::chronos` or treePL, with fossil or secondary calibrations) is needed before any statement about rate per time. (iii) FastTree topology at 1,654 tips with fungi_odb10 markers is adequate for class- and order-level placement; deep nodes among early-diverging lineages and the position of long-branch taxa (Microsporidia, Saccharomycotina) are uncertain, so transitions placed on those nodes must be checked on an alternative topology (constrain to the published consensus, to check, and rerun). (iv) Gene-family birth-death models (CAFE) need an ultrametric tree; they cannot run on the raw FastTree branch lengths.

**Software on this system** (checked 2026-11-08): R 4.5.2 (`/opt/linux/rocky/8.x/x86_64/pkgs/R/4.5.2/bin/Rscript`) with `ape`, `phytools`, `geiger`, `nlme`, `phylolm`, `bayou`, `corHMM`, `OUwie`, `treeio` installed; `caper`, `l1ou`, `PhylogeneticEM`, `motmot`, `RPANDA` not installed. Python `/usr/bin/python3.12` has Biopython, pandas, scipy; `statsmodels` is in `atlas_env/venv`; `ete3`, `dendropy`, `scikit-bio` are missing. CAFE5 was not looked for. Method to package map:

| Method | Package | Use |
|---|---|---|
| Phylogenetic signal (Pagel's lambda, Blomberg's K) | `phytools::phylosig`, `geiger::fitContinuous` | Every continuous trait, before any contrast |
| PGLS | `nlme::gls` with `ape::corPagel`/`corBrownian`; `phylolm::phylolm` (fast at 1,654 tips) | Genus-median contrasts, growth form as a factor |
| Phylogenetic logistic regression | `phylolm::phyloglm` | Presence/absence traits (photoreceptors, transporter families, hydrophobins) |
| Ancestral state, continuous | `phytools::fastAnc`, `anc.ML`; `geiger` | Order-of-events tests along a path |
| Ancestral state, discrete; stochastic character mapping | `phytools::make.simmap`, `corHMM::corHMM` (hidden rates) | Growth-form transitions; count independent origins |
| Correlated evolution of two binary traits (Pagel 1994) | `phytools::fitPagel`, `corHMM` | Loss order (for example flagellum loss vs signaling loss) |
| OU and rate-shift models | `OUwie` (regimes from simmap), `bayou` (shift detection); `l1ou` and `PhylogeneticEM` would need installation | Shifts in optimum at growth-form transitions |
| Gene-family birth-death (CAFE-style) | CAFE5 (not on disk; to install) or `badirate`; needs an ultrametric tree | Family expansion/contraction per branch |
| Dating | `ape::chronos` (penalized likelihood), treePL (not checked) | Needed for any per-time rate statement |

**Lineage sampling caveats that apply to all blocks.** BFD is 82 % Ascomycota by genome; non-Dikarya are 156 filtered genera. Taphrinomycotina is 3 genera. Transition counts from stochastic mapping depend on the sampled genera; a growth-form change that appears "once" may be several with denser sampling. Annotation quality differs by clade (H1): BUSCO median 75-84 % in Saccharomycotina, about 10 % (invalid marker set) in Microsporidia, 84 % in zoosporic lineages, 99 % in large filamentous classes. Any ancestral-state estimate of a count trait inherits these biases on the derived branches. Long branches (Microsporidia, Saccharomycotina, Pneumocystis) attract spurious gene losses from Pfam and SwissProt detection failure; presence traits on those branches need a second detector before a loss is called.

---

## 2. Hypotheses

Format: statement and expected direction; rationale; traits and how to compute them; contrast and unit; confounders and controls; views; compute; priority; novelty; evolutionary framing. Pfam accessions are named from the A7 `pfam_id` column and need a specificity check before use. Literature claims are marked "to check" unless stated otherwise.

### Theme A: Transport

#### T1. Uptake replaces depolymerization: the transporter-to-secretome ratio rises in every independent yeast origin, but absolute transporter density does not

**Statement.** In each yeast origin, the ratio (plasma-membrane transporter proteins) / (secreted CAZymes + secreted peptidases) is higher than in the filamentous sister group. The absolute density of MFS transporters per 1,000 proteins is *not* higher in yeasts. Expected: ratio up in >= 4 of 5 origins; MFS density down or equal.

**Rationale.** A single cell in a liquid or sugar-rich niche takes up small molecules and does not need to depolymerize polymers at a distance. Hyphae both secrete enzymes and import products. The pilot shows yeasts have fewer MFS and sugar transporters per 1,000 proteins (19.9 vs 28.8; 10.0 vs 14.5) but more amino-acid permeases (4.8 vs 3.8) and a secretome fraction less than half of the filamentous value (0.026 vs 0.065, H4). So the ratio, not the transporter count, is the candidate trait.

**Traits.** Numerator: proteins with >= 1 domain from a transporter set (PF07690 MFS_1, PF00083 Sugar_tr, PF00324 AA_permease, PF13520 AA_permease_2, PF01490 Aa_trans, PF03169 OPT, PF00909 Ammonium_transp, PF00916 Sulfate_transp, PF00860 Xan_ur_permease, PF02133 Transp_cyt_pur, PF00005 ABC_tran with PF00664 ABC_membrane, PF00230 MIP, PF01554 MatE, PF00999 Na_H_Exchanger; to check) AND tmhmm PredHel >= 4 AND no signal-peptide-only secreted call. Denominator: A5 `n_cazy_secreted` + `n_merops_secreted`. Also report each term per 1,000 proteins. When TCDB is available: TC class 2.A (porters) and 3.A.1 (ABC) counts replace the Pfam set.

**Contrast and unit.** Genus medians. Five yeast origins each against its sister (the six pairs of PLAN 4.3.4), plus Pezizomycotina dimorphic vs filamentous Eurotiomycetes (15 vs 32 genera). Sign test across origins.

**Confounders.** Proteome size (ratios cancel most of it; also fit log-log allometry of each term and use residuals); secretome undercount from partial models (Q-model; `partial_fraction` covariate); Pfam family size scaling (MFS_1 is a large family; report per-1,000 and residual from allometry); lifestyle (nectar/tap saprotroph is almost all *Saccharomyces*; stratify by `primary_lifestyle` within filamentous group); Agaricomycetes have low MFS density too (12.4), so the filamentous reference must be the taxonomic sister, not all filamentous fungi.

**Views.** (1) Forest plot: log2 ratio difference per origin pair with bootstrap CI. (2) Scatter: transporter count vs secreted-enzyme count (log-log), points = genera, colour = growth form, with allometric lines. (3) Table: per-family density by origin (MFS, Sugar_tr, AA permeases, ABC). (4) Same plot restricted to `primary_lifestyle` in (unspecified_saprotroph, soil_saprotroph) to remove the nectar confound.

**Compute.** A7 filtered scan (2-3 s on the login node), A5 (small), tmhmm join needs A3-scale protein rows (`tmhmm.parquet` 1.8 GB, 45.7 M rows; SLURM, minutes by the build log).

**Priority.** P1. **Novelty.** The secretome side is known (H4). The ratio framing with replicated origins is less explored (to check).

**Evolutionary framing.** (1) Polarity: the ancestral Dikarya state is a hyphal, secreting organism (outgroups: Mucoromycota 35 genera, Zoopagomycota s.l. 28, zoosporic lineages 43; all have secretome fractions of 0.035-0.075 in the pilot). The derived change is in each yeast lineage. For Saccharomycotina the outgroups are Taphrinomycotina (3 genera; *Taphrina* is dimorphic and *Pneumocystis* is a parasite, so *Schizosaccharomyces* is the only free-living yeast outgroup) and Basidiomycota. (2) Placement: 5 independent yeast origins plus 2 dimorphic transitions; stochastic mapping of growth form on the genus tree gives the branches and the number of transitions. (3) Mode: expected to be loss of secreted enzyme families (contraction) on the yeast stem branch, with transporter families kept (no expansion); a burst at the transition is expected, not gradual decline; no link to WGD except that the Saccharomycetaceae WGD post-dates the yeast transition (to check; wgd.ks tables can show the WGD peak but not the order). (4) Methods: PGLS of the log ratio on growth form (regime from simmap); `OUwie` with two regimes (yeast, hyphal) to test for a shift in optimum rather than a rate change; CAFE-style birth-death on secreted CAZy family counts after dating. (5) Order of events: reconstruct ancestral secretome fraction and transporter density along the path Pezizomycotina root -> Saccharomycotina stem -> crown; if the secretome drop maps to the stem and the transporter ratio change is spread along crown branches, uptake tuning came after secretion loss. (6) Caveats: Taphrinomycotina sampling is 3 genera; secretome counts in Saccharomycotina are reduced by lower BUSCO completeness (section 1.4); the long Saccharomycotina branch inflates apparent family loss.

#### T2. Obligate biotrophs and intracellular parasites lose the same inorganic-nutrient transporter families, in independent origins

**Statement.** Transporter families for inorganic nitrogen, sulfate and nucleobase uptake (ammonium transporters, sulfate permeases, NRT2 nitrate transporters, Xan_ur permeases, Transp_cyt_pur) are absent or reduced to one copy in >= 4 independent obligate biotroph lineages, while sugar and amino-acid transporters are retained. Intracellular parasites (Microsporidia, *Pneumocystis*) show the deeper loss and gain bacterial-type nucleotide transporters (PF03219 TLC; to check).

**Rationale.** Host cells provide reduced nitrogen and organic sulfur. The pilot shows Microsporidia have 0 Ammonium_transp, 0 Xan_ur_permease, 0 OPT, 0 TrkH, 0 Transp_cyt_pur and 0.01 Mito_carr presence. Rusts and powdery mildews are reported to have lost nitrate and sulfate assimilation (to check); the transporter side across many lineages is not in the atlas.

**Traits.** Presence and copy number per genome of: PF00909 Ammonium_transp, PF00916 Sulfate_transp, PF00860 Xan_ur_permease, PF02133 Transp_cyt_pur, PF03169 OPT, PF02386 TrkH, PF03219 TLC (nucleotide transporter), PF00083 Sugar_tr, PF00324/PF13520 amino-acid permeases, PF00153 Mito_carr. NRT2 has no specific Pfam (MFS_1): use SwissProt best hits to NRT2 family entries (`swissprot_subject.protein_name LIKE '%nitrate transporter%'`; to check). Assimilation enzymes as a paired trait: nitrate reductase (PF00174 Oxidored_molyb; to check), sulfite reductase (PF01077 NIR_SIR; to check).

**Contrast and unit.** Genus presence (present in >= 50 % of representative genomes of the genus). Lineages: Pucciniales (5 genera), Erysiphaceae (7 genera raw; relaxed filter), Glomeromycota (11), Microsporidia (23; size < 10 Mb), *Pneumocystis* (1), Entomophthorales (6), Neocallimastigomycetes (5; anaerobic, not biotroph, as a control for "nutrient-rich gut"). Each lineage vs its nearest free-living sister in BFD (Pucciniales vs Microbotryomycetes + Agaricomycotina; Erysiphaceae vs other Helotiales (67 genera); Glomeromycota vs Mortierellomycota + Mucoromycota; Microsporidia vs Rozellomycota (2) and chytrids; Entomophthorales vs Basidiobolus + Kickxellomycota). Report the number of lineages with loss.

**Confounders.** Assembly completeness (rusts and AM fungi have high duplication and fragmentation; use presence, not counts; relaxed BUSCO; note Microsporidia bins); Pfam sensitivity on divergent Microsporidia proteins (report `frac_with_pfam` per genome and use SwissProt hits as a second detector); genome size (large rust genomes inflate counts of everything; use presence); metagenome bins with host contamination give false presence (size filter).

**Views.** (1) Presence heatmap: lineage x transporter family (rows = genera, grouped by lineage; cell = copies). (2) Dot plot: loss count per lineage with sister-group values. (3) Table: families lost in >= 4 lineages vs families lost in 1. (4) Scatter: transporter family richness vs assimilation enzyme presence.

**Compute.** A7 filtered scan (small). SwissProt hits: A3 `protein_core` filtered to about 60 genera (SLURM, minutes).

**Priority.** P1. **Novelty.** Loss of nitrate/sulfate assimilation is known for some lineages (to check). A kingdom-wide replicated test with transporters is less explored.

**Evolutionary framing.** (1) Polarity: presence of ammonium, sulfate, nitrate and nucleobase transporters is ancestral; the pilot shows them in 95-100 % of free-living genomes in every phylum with >= 5 genera. The loss is on the parasite branch. Outgroups per lineage: Pucciniales vs Microbotryomycetes (5 genera) and Agaricomycotina; Erysiphaceae vs other Helotiales (67); Glomeromycota vs Mortierellomycota (12) and Mucoromycota (35); Microsporidia vs Rozellomycota (1) and Chytridiomycota (32); Entomophthorales vs Basidiobolomycota (1) and Kickxellomycota (14). (2) Placement: >= 5 independent losses expected; stochastic mapping of a binary "obligate biotroph" character (hand list) gives the branches; losses that map to the same stem branch as the lifestyle change support the link. (3) Mode: gene loss (family contraction to zero), expected as a burst on the stem branch; for Microsporidia the loss is part of genome reduction (G1); no WGD link expected. (4) Methods: Pagel 1994 correlated evolution (`phytools::fitPagel`) between "obligate biotroph" and "transporter family present", one family at a time, with the dependent model where transporter loss depends on lifestyle; phylogenetic logistic regression (`phyloglm`) of presence on lifestyle; `corHMM` with hidden rates for the lifestyle character. (5) Order of events: the Pagel transition-rate matrix gives which state changes first; if the rate of transporter loss is high only in the biotroph state and the reverse transition (lifestyle change after loss) is near zero, lifestyle change came first. Compare with assimilation enzymes: enzyme loss before or after transporter loss along the same path (ancestral presence per node). (6) Caveats: Microsporidia and *Pneumocystis* are long branches with weak Pfam detection, so "absence" needs a SwissProt check; Erysiphaceae are 2 filtered genera; rust assemblies are fragmented; a single origin inherited by a clade (Pucciniales 5 genera from one origin) counts once.

#### T3. Symbiosis changes which transporters expand: ECM fungi expand ammonium and amino-acid/peptide transporters as they lose PCW CAZymes, and lichen mycobionts with cyanobacterial photobionts carry fewer nitrogen transporters than those with green algae

**Statement.** (a) Within each order that holds both ECM and non-ECM genera (Boletales 45 vs 7, Agaricales 9 vs 95, Pezizales 7 vs 9, Russulales 5 vs 11, Cantharellales 3 vs 5, Gomphales 3 vs 2), ECM genera have higher density of ammonium transporters, amino-acid permeases and oligopeptide transporters per 1,000 proteins, and lower PCW CAZyme density (the second part is H4 V23). (b) Within Lecanoromycetes, genera with cyanobacterial photobionts (9) have fewer ammonium and nitrate transporters than genera with chlorococcoid photobionts (124), and both have more polyol/sugar transporters than non-lichen Pezizomycotina. (c) ECM genera with long-distance exploration type (rhizomorphs; 45 genera, mostly Boletales) have more hydrophobin and cell-wall structural protein genes than contact or short-distance types (6 + 6 + 4 genera).

**Rationale.** ECM fungi take organic N from soil and give it to the plant; cyanolichens receive fixed N from the photobiont; green-algal lichens receive polyols. Rhizomorphs are differentiated multicellular organs with a hydrophobic outer layer (to check).

**Traits.** Pfam families as in T1 and T2; PF00909 Ammonium_transp; PF03169 OPT; PF00324 AA_permease; PF00083 Sugar_tr; hydrophobins PF01185, PF06766, PF29465, PF29785; CFEM PF05730. CAZy PCW from A5 `cazy_substrate_counts`.

**Contrast and unit.** Genus medians within order (a) or within class (b), (c) within ECM genera. For (a) a per-order sign test over >= 5 orders is the replicate. For (b) and (c) the counts are small; report effect sizes with bootstrap CI, no p-value claim.

**Confounders.** Boletales dominate ECM (45 of 70 genera); order random intercept required; genome size and TE content (ECM genomes are large; use per-1,000 and allometric residuals); lichen assemblies are metagenomic (BUSCO covariate; photobiont contamination can add algal transporters: require a fungal SwissProt best hit for counted proteins); exploration type is mostly long-distance, and the other classes are 4-6 genera each (state low power).

**Views.** (1) Paired dot plot per order: ECM vs non-ECM genus medians, one panel per transporter family. (2) Lichen: box of transporter density by photobiont type. (3) ECM exploration type vs hydrophobin density, strip plot.

**Compute.** A7 (small), A5 (small), A3 join for the fungal-hit filter (SLURM, minutes).

**Priority.** P2 for (a); P3 for (b) and (c) because of counts. **Novelty.** ECM CAZyme loss is well known. Transporter gain across replicated ECM origins and the photobiont contrast are less explored (to check).

**Evolutionary framing.** (1) Polarity: saprotrophy with a full PCW CAZyme set is ancestral in Agaricomycetes and Pezizales (outgroups: non-ECM genera in the same orders, listed above; corticioid saprotrophs 49 genera). ECM is derived in each order. For lichens, the ancestral Lecanoromycetes state is uncertain (lichenization may be ancestral to the class, with losses; to check), and the within-class contrast is 122 vs 3 genera, so polarity needs Arthoniomycetes (7) and Lichinomycetes (4) and non-lichen Pezizomycotina outgroups. (2) Placement: ECM origins in BFD by order: Boletales, Agaricales, Russulales, Cantharellales, Gomphales, Thelephorales, Pezizales, Eurotiales (*Elaphomyces*, 2 genera), Endogonales (1 genus): >= 8 origins if each order holds one (Agaricales holds several; to check with `Ectomycorrhiza_lineage_template`, 25 lineages seen in BFD). (3) Mode: CAZyme loss (contraction) is expected on each ECM stem; transporter change may be expansion (gene duplication, young paralogs from wgd.ks nodes in the transporter families) or retention; for cyanolichens, loss of N transporters. (4) Methods: PGLS with the ECM lineage as a regime; `OUwie` multi-regime on transporter density; `fitPagel` between ECM state and transporter copy class; CAFE-style for ammonium transporter family size after dating. (5) Order of events: ancestral reconstruction of PCW CAZyme count and ammonium transporter count along each ECM stem; test whether CAZyme loss precedes transporter gain (CAZyme drop on the stem, transporter gain on crown branches) or the two change together. (6) Caveats: Boletales ECM is largely one origin (45 genera inherit it); rhizomorph exploration types other than long-distance have 4-6 genera; lichen assemblies are metagenomic and may carry algal transporter genes.

### Theme B: Signaling

#### S1 (falsification test). Multicellular complexity does not require more signaling genes: kinase, small-GTPase and GPCR counts scale sub-linearly with proteome size, and residuals are not higher in fruiting-body-forming lineages

**Results:** [HS_S1_signaling_scaling.md](HS_S1_signaling_scaling.md) (built 2026-11-08).

**Statement.** The popular idea: complex multicellular fungi (fruiting bodies) need a larger signaling repertoire. Prediction that would support it: positive residuals of kinase, Ras-family, RhoGEF/GAP and GPCR counts (after proteome-size allometry) in agaricoid/polyporoid Agaricomycetes and in perithecial/apothecial Pezizomycotina relative to genera with no fruit body in the same class. Prediction that falsifies it: no residual difference, or negative. The pilot points to the second: kinases per 1,000 proteins are 17.0 in yeasts, 10.4 in non-Agaricomycetes filamentous fungi, 13.8 in Agaricomycetes, 19.8 in Mucoromycota; Ras-family 5.3 (yeast) vs 2.9 (filamentous). Absolute counts (pilot, from medians): about 97 kinases per yeast proteome, 116 per filamentous Pezizomycotina proteome, 173 per Agaricomycetes proteome, 242 per Mucoromycota s.l. proteome.

**Rationale.** The kinome is near-constant in Dikarya and expanded in Mucoromycota and zoosporic lineages (to check against published kinome surveys). If complexity needed more signaling genes, the yeast proteome would not keep the full set.

**Traits.** Proteins with PF00069 Pkinase or PF07714 PK_Tyr_Ser-Thr (count a protein once); PF00071 Ras and PF00025 Arf; PF00620 RhoGAP, PF00621 RhoGEF; GPCR set (PF00001 7tm_1, PF02116 STE2, PF02076 STE3, PF01036 Bac_rhodopsin, PF10192 and other fungal GPCR classes; to check); two-component (PF00512 HisKA, PF00072 Response_reg); PF00503 G-alpha, PF00615 RGS; calcium (PF00036, PF13499 EF-hand, PF00168 C2); PF00018 SH3; PF00169 PH. From A7. Kinase subfamilies (MAPK, PKA, CDK) would need a classifier; not in tables.

**Contrast and unit.** Log-log OLS of count vs `n_proteins` on genus medians, all fungi, slope with CI; then residuals by (i) Agaricomycetes fruit-body type (agaricoid 144, polyporoid 45, corticioid 49, gasteroid 9, clavarioid 10 genera) and hymenium type (gills 100, poroid 94, smooth 46, closed 16, hydnoid 10); (ii) Ascomycota fruit-body type within class (Sordariomycetes perithecium 117 vs none 41; Eurotiomycetes cleistothecium 10 vs unlabelled; Leotiomycetes apothecium 39 vs none 5); (iii) yeast origins vs sisters.

**Confounders.** Proteome size (the allometry is the control); Pfam family size and detection (Pkinase is well detected; GPCR families are small and clade-specific); gene-model fragmentation splits kinases (Q-model; `partial_fraction` covariate); phylogenetic non-independence (class/order random intercepts now; PGLS later); FungalTraits "none" fruit body is partly "unknown".

**Views.** (1) Allometry panel per signaling class (count vs proteome size, log-log; genus points coloured by group; slope and CI). (2) Residual strip plot by fruit-body type (Agaricomycetes) and by hymenium type. (3) Residual strip by yeast origin with sisters. (4) Table: absolute kinase count per group with CI.

**Compute.** A7 filtered scan (small). All on the login node after a one-time extract.

**Priority.** P1. **Novelty.** Sub-linear kinase scaling is reported in some studies (to check); a kingdom-wide test against fruit-body type is less explored.

**Evolutionary framing.** (1) Polarity: the pilot shows the largest kinase densities in the early-diverging lineages (Mucoromycota s.l. 19.8 per 1,000, zoosporic 17.6) and the lowest in filamentous Pezizomycotina (10.4). The polarized question is therefore not "did complex forms gain kinases" but "did Dikarya hyphal lineages lose or dilute them, and did yeasts keep the ancestral count". Outgroups: Chytridiomycota 32, Blastocladiomycota 5, Mucoromycota 35, Mortierellomycota 12, Kickxellomycota 14 genera. (2) Placement: ancestral kinase count at the Dikarya node, the Agaricomycetes node and each yeast stem; a single Dikarya-wide shift inherited by all classes versus separate shifts at each fruit-body origin. Fruit-body origins in Agaricomycetes are several (to check; `corHMM` on fruit-body type with hidden rates). (3) Mode: expected mode is proteome growth by non-signaling families (secreted enzymes, orphans) that dilutes a near-constant kinome, so the count shows no burst; kinome expansion in Mucoromycota may be linked to WGD or large-scale duplication in Mucorales (to check with wgd.ks: young paralog nodes inside kinase families per genome). (4) Methods: PGLS of log count on log proteome size with fruit-body type as a factor; `phylosig` to show whether kinase count carries phylogenetic signal beyond proteome size; `bayou` on the allometric residual to find shift branches; CAFE-style on kinase family size after dating. (5) Order of events: along the path root -> Dikarya -> Agaricomycetes -> agaricoid genera, reconstruct proteome size and kinase count; if proteome size rises before any kinase change, dilution is confirmed. (6) Caveats: kinase families are well detected, but subfamily turnover (loss of one MAPK, gain of another) is invisible to counts; the deep-node placement of Mucoromycota vs Zoopagomycota is uncertain in a FastTree genus tree; gene-model fragmentation in low-BUSCO clades splits kinases.

#### S2. Environmental sensing modules (photoreceptors, phytochrome, PAS/GAF, two-component) are lost in yeasts, in intracellular parasites and in hypogeous fruiters, in independent origins

**Statement.** Presence of phytochrome (PF00360 PHY + PF01590 GAF), opsin (PF01036 Bac_rhodopsin), cryptochrome/photolyase (PF00875 DNA_photolyase, PF03441 FAD_binding_7), WC-1-type PAS+GATA proteins (PF00989 PAS with PF00320 GATA on one protein; architecture from A3 `pfam_arch`), and two-component histidine kinases (PF00512 HisKA + PF00072 Response_reg) drops in each independent yeast origin, in intracellular parasites, and in hypogeous fruiters (gasteroid-hypogeous Agaricomycetes 7 genera, Pezizales truffles 5 genera, *Elaphomyces* in Eurotiales) relative to their epigeous sisters. Expected: phytochrome and opsin lost first; PAS kinases reduced; two-component reduced by copy.

**Rationale.** Surface-living mycelia use light and humidity to time sporulation. Pilot: PHY present in 14 % of yeast genomes vs 91 % of filamentous; Bac_rhodopsin 56 % vs 88 %; PAS 59 % vs 98 %; HisKA 0.56 vs 1.00 per 1,000; Microsporidia 0 for all. Agaricomycetes: PHY 62 %, opsin 42 %. The hypogeous contrast is a new convergence test with at least 3 origins in BFD.

**Traits.** Presence and copies from A7; WC-1 architecture from A3 `pfam_arch` (PAS;PAS;GATA pattern; to check). Two-component: count of HisKA proteins; subfamilies (osmosensor, light, phytochrome-linked) need a classifier (not available; use architecture with PHY or PAS).

**Contrast and unit.** Genus presence. Replicates: 5 yeast origins vs sisters; parasites (Microsporidia, *Pneumocystis*, Pucciniales, Erysiphaceae, Glomeromycota; see T2); hypogeous vs epigeous within Boletales, Agaricales, Pezizales (counts: 7 Agaricomycetes genera in 4 orders; 5 Pezizomycetes genera). Sign test over origins.

**Confounders.** Detection (small families; use presence); contamination in lichen or metagenome bins; FungalTraits labels fruit body at genus level (sequestrate species inside epigeous genera are invisible); N50 (presence of a 1-copy gene depends on assembly completeness; use `busco_complete` covariate and require >= 90 % for presence/absence claims).

**Views.** (1) Presence matrix (families x genera) ordered by taxonomy with group colour bar. (2) Loss-count per origin with sister values (dot plot). (3) Hypogeous vs epigeous paired plot per order. (4) Copies of HisKA per 1,000 by group.

**Compute.** A7 filtered scan (small); A3 architecture filter for WC-1 (SLURM, minutes).

**Priority.** P1. **Novelty.** Photoreceptor loss in yeasts is known in parts (to check); the hypogeous and parasite convergence test is less explored.

**Evolutionary framing.** (1) Polarity: phytochrome, opsin and PAS/GATA photoreceptors are present in the pilot in 91-100 % of filamentous Pezizomycotina and Mucoromycota genomes (PHY 91 % / 5 %: Mucoromycota lack phytochrome but keep PAS and opsin partly), so presence is ancestral for Dikarya; for the fungal root, zoosporic lineages have PAS (95 %) and opsin (54 %) but PHY in only 11 %, so phytochrome may be a Dikarya gain (to check) and its absence in yeasts is a loss within Dikarya. Outgroups for each yeast origin as in T1; for hypogeous fruiters the epigeous genera of the same order. (2) Placement: losses on each yeast stem (5), on each intracellular-parasite stem (>= 3), on each hypogeous branch (>= 3 orders); a loss shared by a whole clade (all Saccharomycotina) counts once. (3) Mode: gene loss, expected in a fixed order (phytochrome first, then opsin, then WC-1/PAS, then two-component copies); a burst at the transition would show as all losses on one stem branch; gradual loss would spread them across crown branches. (4) Methods: `fitPagel` for growth form (yeast / hyphal) vs each photoreceptor presence; `phyloglm` with `busco_complete` as a covariate; `make.simmap` on each presence character to count independent losses and place them. (5) Order of events: ancestral presence at successive nodes along the path to each derived clade; the first node where phytochrome is absent vs the first where opsin is absent gives the loss order; compare the order across the yeast origins. (6) Caveats: single-copy genes are missed on fragmented assemblies (require BUSCO >= 90 % for loss calls); hypogeous labels are genus-level and the groups are 7 + 5 genera; the hypogeous test cannot separate loss of light sensing from loss of other above-ground traits.

#### S3. Zoosporic lineages keep a cyclic-nucleotide and ion-channel signaling toolkit that is lost with the flagellum

**Statement.** Guanylate cyclase, cNMP-binding, phosphodiesterase (PDEase_I), voltage-gated ion channel (Ion_trans), and axonemal/IFT genes form a co-lost module. Genera with zoospores (43 filtered genera in Chytridiomycota, Blastocladiomycota, Neocallimastigomycota, Monoblepharomycota, Aphelidiomycota) have 5-10x the density of these families; non-flagellated lineages have lost IFT/BBS/radial-spoke genes and show a parallel drop in cyclic-nucleotide families. Expected: strong positive association of the two gene sets across genera; the signaling set is retained partly in some non-flagellated lineages (Mucoromycota Ion_trans 0.81 per 1,000 in the pilot) and lost completely in Microsporidia.

**Rationale.** Zoospore chemotaxis and phototaxis use cGMP and calcium signaling in Blastocladiella and Allomyces (to check). Pilot: Guanylate_cyc 1.13 per 1,000 in zoosporic vs 0.09 in filamentous; cNMP_binding 1.62 vs 0.34; PDEase_I 0.82 vs 0.08; Ion_trans 1.23 vs 0.27; Dynein_heavy 0.98 vs 0.08. IFT57 (PF10498), IFT46 (PF12317), MIP-T3/IFT54 (PF10243) and Radial_spoke_3 (PF06098) are present in 89-94 % of zoosporic genomes and 0 % elsewhere, so the flagellar set is a clean marker.

**Traits.** Flagellar set: PF10498, PF12317, PF10243, PF06098, PF04712 Radial_spoke, PF14779 BBS1, PF14783 BBS2, PF03028 Dynein_heavy count, PF01221 Dynein_light. Signaling set: PF00211 Guanylate_cyc, PF00027 cNMP_binding, PF00233 PDEase_I, PF00520 Ion_trans, PF08016 PKD_channel, PF00036/PF13499 EF-hand. All from A7.

**Contrast and unit.** Genus medians. Zoosporic (43 genera; 3 phyla as replicates: Chytridiomycota 31, Blastocladiomycota 5, Neocallimastigomycota 5) vs the nearest non-flagellated groups (Zoopagomycota s.l. 28 genera; Mucoromycota s.l. 59; Microsporidia 23 as the extreme loss). Within zoosporic: Neocallimastigomycetes (anaerobic gut, 5 genera) vs aerobic chytrids. The number of independent flagellum losses in the fungal tree is a literature question (to check; PGLS needed to count transitions).

**Confounders.** Low BUSCO in zoosporic genomes (median 83.8 in the filtered set; Neocallimastigomycetes 49.7): use per-1,000 values and `busco_complete` covariate; chytrid gene models are AT-rich and may be fragmented (Q-model); Pfam detection on divergent proteins (report `frac_with_pfam`); small genus counts in Blastocladiomycota (5).

**Views.** (1) Heatmap: flagellar set and signaling set presence per genus, ordered by phylum. (2) Scatter: flagellar-set score vs signaling-set density (genus points). (3) Density per 1,000 by group (strip) for each signaling family.

**Compute.** A7 filtered scan (small).

**Priority.** P2. **Novelty.** Flagellar gene loss is known. Linking it to a specific signaling module across all BFD zoosporic genera is less explored (to check).

**Evolutionary framing.** (1) Polarity: the flagellum and its signaling module are ancestral (the opisthokont root has a flagellum; in BFD the earliest-branching sampled lineages, Chytridiomycota 32, Blastocladiomycota 5, Aphelidiomycota 2, Rozellomycota 1 genus, have it). The loss is the derived state in the terrestrial lineages. (2) Placement: the number of losses is a literature question (to check; between 1 and about 6 depending on the topology among Zoopagomycota, Mucoromycota, Olpidiomycota and Dikarya). BFD can test it: stochastic mapping of the flagellar-gene score on the genus tree, with the deep topology constrained to the published consensus as an alternative. (3) Mode: block loss of >= 10 families on one branch (burst); the signaling module may be lost on the same branch or decay gradually afterward (Mucoromycota keep Ion_trans at 0.81 per 1,000); no WGD link. (4) Methods: `fitPagel` between flagellar-set presence and each signaling family presence; `make.simmap` to count flagellum losses; PGLS of signaling density on flagellar score with `busco_complete` covariate. (5) Order of events: with the Pagel dependent model, estimate whether the rate of signaling loss rises only after flagellar loss (flagellum first) or whether the two rates are symmetric; ancestral density of Guanylate_cyc at nodes before and after each flagellar loss. (6) Caveats: zoosporic BUSCO is 84 % median (Neocallimastigomycetes 50 %); AT-rich chytrid genomes give fragmented models; the deep-node order is the weak point of a FastTree genus tree, and a change in it changes the loss count.

#### S4. Pheromone-receptor and G-alpha copy number expands in Agaricomycetes and tracks mating-system complexity, not fruit-body complexity

**Statement.** STE3-type pheromone receptors (PF02076) and G-alpha (PF00503) copies per genome are higher in Agaricomycetes (pilot: 0.44 and 1.21 per 1,000) than in any Ascomycota group (0.08-0.18 and 0.31-0.49). Within Agaricomycetes, copy number is higher in tetrapolar than bipolar genera and does not differ by fruit-body type after the mating-system split. STE2 (PF02116) is absent from Basidiomycota (pilot: present in 1 % of Agaricomycetes genomes), consistent with the known Ascomycota restriction (to check).

**Rationale.** Tetrapolar Agaricomycetes carry multi-allelic pheromone-receptor loci with several receptor genes; this is mating biology, not multicellular signaling. The test separates two explanations of the Agaricomycetes expansion.

**Traits.** PF02076 STE3 count; PF08015 Pheromone (precursor; present in 360 genomes); PF00503 G-alpha; PF00615 RGS; PF05722 Ustilago_mating; PF04769 MATalpha_HMGbox; PF17043 MAT1-1-2. From A7. Mating-system label: not on disk; needs a hand-made list from the literature (about 100 Agaricomycetes genera; to build).

**Contrast and unit.** Genus medians. Agaricomycetes (291 genera) vs Pezizomycotina (584) for the class-level claim; within Agaricomycetes, mating system (label needed) and fruit-body type (FungalTraits) in a two-factor model with order random intercept.

**Confounders.** Dikaryon assemblies carry both alleles (BUSCO duplication; use haploid-filtered genomes: `busco_duplicated <= 10` where possible, and report the duplication covariate); STE3 Pfam also hits non-pheromone GPCRs in some clades (to check); without the mating-system label only the class-level claim is testable.

**Views.** (1) Strip plot of STE3 and G-alpha per genome by class. (2) Within Agaricomycetes: STE3 count by mating system (when the label exists) and by fruit-body type. (3) Scatter STE3 vs busco_duplicated to show the allele effect.

**Compute.** A7 filtered scan (small).

**Priority.** P2 (P3 until the mating-system list exists). **Novelty.** Multi-allelic receptors are well known in model Agaricomycetes; the kingdom-wide copy-number scan and the fruit-body separation are less explored.

**Evolutionary framing.** (1) Polarity: one or two STE3-type receptors per genome is the state in Pezizomycotina (0.08 per 1,000), Saccharomycotina (0.18) and Mucoromycota (0, which needs a check: Mucoromycota use trisporic acid, not peptide pheromones, so the family may be absent ancestrally or lost). The Agaricomycetes value (0.44) is derived expansion on the Agaricomycotina or Agaricomycetes stem (Tremellomycetes and Pucciniomycotina values give the placement; not yet extracted). STE2 absence in Basidiomycota is a loss on the Basidiomycota stem if Mucoromycota have it (pilot: 0 %, so STE2 may be an Ascomycota gain; to check). (2) Placement: expansion placed at one node (inherited by the class) vs repeated expansions in tetrapolar lineages; bipolar reversals (to check: several bipolar Agaricomycetes derived from tetrapolar ancestors) are the within-class test. (3) Mode: tandem duplication (expansion) at the mating-type locus; wgd.ks young paralog nodes within the STE3 family per genome show whether copies are recent; no WGD link. (4) Methods: `bayou` on STE3 count to find the shift branch; `OUwie` with mating-system regimes once the label exists; `make.simmap` on mating system. (5) Order of events: ancestral STE3 count at the Agaricomycotina node vs the Agaricomycetes node vs fruit-body origins; if the expansion predates fruit-body origins, it is not a fruit-body trait. (6) Caveats: dikaryotic assemblies double counts (duplication covariate); the mating-system label is missing; STE3 Pfam specificity.

### Theme C: Communication, adhesion, cell wall

#### C1. NLR-type allorecognition repertoires (HET, NACHT, NB-ARC) track septate heterokaryotic hyphae and collapse in every yeast origin, in coenocytic lineages and in dikaryon-forming Agaricomycetes

**Results:** [HS_C1_nlr_allorecognition.md](HS_C1_nlr_allorecognition.md) (built 2026-11-08).

**Statement.** The density of HET-domain, NACHT and NB-ARC proteins per 1,000 proteins is: filamentous Pezizomycotina > lichenized Pezizomycotina > Agaricomycetes > Mucoromycota s.l. ≈ yeasts ≈ Microsporidia ≈ 0. Within Pezizomycotina, dimorphic Eurotiomycetes (15 genera) have fewer HET proteins than filamentous Eurotiomycetes (32). Expected: an order-of-magnitude drop at each shift to a non-networking growth form.

**Rationale.** Heterokaryon incompatibility genes act at hyphal fusion between genetically distinct septate hyphae. Yeasts do not fuse somatically; coenocytic Mucoromycota fuse rarely (to check); Agaricomycetes fuse but regulate the dikaryon through mating type. Pilot: HET 4.69 per 1,000 in filamentous Pezizomycotina, 3.71 in lichens, 0.23 in Agaricomycetes, 0 in yeasts, Mucoromycota and Microsporidia (present in 33 % of Mucoromycota genomes); NACHT 1.47 / 1.83 / 1.91 / 0.08 / 0 / 0.

**Traits.** PF06985 HET, PF05729 NACHT, PF17100 NACHT_N, PF00931 NB-ARC, PF17107 SesA, PF11558 HET-s, PF17108 HET-S, plus the repeat partners PF12796/PF00023 Ank, PF00515/PF13414 TPR, PF00400 WD40 on the same protein (A3 `pfam_arch`): count NLR architectures (NACHT or NB-ARC with a C-terminal repeat) separately from bare HET proteins.

**Contrast and unit.** Genus medians. Replicates: 5 yeast origins vs sisters; Pezizomycotina dimorphic vs filamentous within Eurotiomycetes; Mucoromycota s.l. vs Dikarya; Agaricomycetes by fruit-body type (does complexity add NLRs? pilot says NACHT is as dense in Agaricomycetes as in Pezizomycotina, HET is not).

**Confounders.** Proteome size (per-1,000 and allometry); HET/NACHT families are large and fast-evolving (Pfam detection varies by clade; report presence per genome and e-value sensitivity); strain-level variation is high in *Neurospora* and *Podospora* (genus medians); phylogeny (class random intercept).

**Views.** (1) Strip plot: HET, NACHT, NB-ARC and NLR-architecture density by group. (2) Paired origin plot (yeast vs sister; dimorphic vs filamentous). (3) Architecture table: top NLR domain combinations per phylum.

**Compute.** A7 filtered scan (small); A3 architecture scan (3.7 GB, SLURM, minutes).

**Priority.** P1. **Novelty.** NLR diversity in filamentous Ascomycota is described (to check: Dyrka et al. 2014). Replicated loss across yeast origins and the dimorphic within-class contrast are less explored.

**Evolutionary framing.** (1) Polarity: the pilot shows HET at 0 per 1,000 in Mucoromycota s.l. (present in 33 % of genomes at low copy), Zoopagomycota s.l. (1 %), zoosporic (0, present in 33 %) and Microsporidia (0). So HET-domain expansion is a derived gain on the Pezizomycotina stem, not a loss elsewhere; NACHT is present at low copy in the early-diverging lineages (0.07-0.08 per 1,000, present in about 60 % of genomes) and expanded in both Pezizomycotina and Agaricomycetes (1.5-1.9), which suggests either one expansion on the Dikarya stem with loss in yeasts, or two expansions. Outgroups: Mucoromycota 35, Kickxellomycota 14, Chytridiomycota 32 genera. (2) Placement: HET gain on the Pezizomycotina stem (single origin inherited by the subphylum); NACHT ancestral reconstruction at the Dikarya, Ascomycota and Basidiomycota nodes decides one vs two expansions; losses on 5 yeast stems and on the Onygenales dimorphic branches. (3) Mode: family expansion by duplication and domain shuffling (NLR architectures), likely gradual along Pezizomycotina branches (counts vary 10-fold among genera, to be measured); losses in yeasts as bursts; wgd.ks young paralog nodes inside HET/NACHT families give the recent duplication rate per genome. (4) Methods: `phylosig` (expect high lambda); PGLS on growth form; `bayou` for shift branches on log density; CAFE-style on family size after dating; `corHMM` on growth form to place yeast transitions. (5) Order of events: along Eurotiomycetes branches leading to dimorphic Onygenales, reconstruct HET count and growth form; if HET contraction maps to the branch where dimorphism arises, the two changed together; if it precedes, NLR loss is not caused by the yeast phase. (6) Caveats: HET and NACHT are fast-evolving and Pfam detection varies by clade; counts in *Neurospora* and *Podospora* vary by strain; the Pezizomycotina stem is one branch, so the gain is a single event and gives no replication.

#### C2. Surface-protein repertoires swap at each organization shift: hydrophobins and cerato-platanins mark aerial multicellular structures; lineage-specific adhesins mark unicellular forms; neither tracks fruit-body complexity once presence is set

**Statement.** (a) Hydrophobin class I (PF01185) and class II (PF06766) density per 1,000 proteins is near zero in every yeast origin (pilot: 1 % of yeast genomes have any), present in about half of filamentous Pezizomycotina genomes (0.06 per 1,000), and highest in Agaricomycetes (94 % present, 1.29 per 1,000). (b) Within Agaricomycetes, hydrophobin count is higher in agaricoid and polyporoid genera than in corticioid genera, and higher for gills/poroid than for smooth hymenia; the same sign is expected in Pezizomycotina for apothecium/perithecium vs none. (c) In yeasts, GPI-anchored proteins with no Pfam domain, or with Flo11/Flocculin/PA14/Candida_ALS_N domains, replace hydrophobins; the GPI-anchored fraction of the proteome itself does not differ (pilot: 0.0158 yeast vs 0.0151 filamentous). (d) Falsification part: after presence is set, hydrophobin copy number does not scale with hymenium complexity (gills vs poroid vs smooth).

**Rationale.** Hydrophobins coat aerial hyphae and fruit bodies; yeasts have none and use cell-wall adhesins for flocculation and host attachment. Pilot hydrophobin presence in lichens is only 12 %, which is lower than expected from the literature on lichen hydrophobins (to check); this is itself a view.

**Traits.** A7: PF01185, PF06766, PF29465 Hydrophobin_like, PF29785 Hydrophobin_D, PF07249 Cerato-platanin, PF05730 CFEM, PF12296 HsbA, PF10182 Flo11, PF00624 Flocculin, PF13928 Flocculin_t3, PF07691 PA14, PF11766 Candida_ALS_N, PF10528 GLEYA, PF29914 EGT2_adhesin. GPI-anchored proteins without Pfam: A3 `protein_core` (has_arch = false) joined to `predgpi` (feature = GPI-anchor) and signalp (A5 `secreted_proteins.parquet` holds the strict secreted set; GPI flag is in A5 `secretome.n_gpi` only at genome level; a protein-level join needs `predgpi.parquet`, 1.9 GB). Ser/Thr-rich composition needs sequences, which are not in the tables (not available).

**Contrast and unit.** Genus medians. (a) 5 yeast origins vs sisters. (b) Agaricomycetes fruit-body type (agaricoid 144, polyporoid 45, corticioid 49 genera) and hymenium (gills 100, poroid 94, smooth 46) with order random intercept; Ascomycota by class. (c) yeast origins. Lichens (156 genera) as a separate panel.

**Confounders.** Hydrophobins are short (about 100 aa) and cysteine-rich; gene prediction misses them without RNA-seq (D8 flag not available; use `frac_short_lt100` and `busco_complete` as covariates; report presence rather than count where prediction is doubtful); Pfam detection of divergent class I hydrophobins (to check: HMM sensitivity differs for Agaricomycetes); FungalTraits fruit-body type is a genus label; proteome size.

**Views.** (1) Strip plot: hydrophobin density by group and by yeast origin. (2) Agaricomycetes: hydrophobin count by fruit-body type and by hymenium type (box with genus points). (3) Stacked bar: surface-protein family composition (hydrophobin, cerato-platanin, CFEM, Flo/PA14/ALS, GPI-no-domain) per group. (4) Lichen hydrophobin presence vs busco_complete (artefact check).

**Compute.** A7 filtered scan (small); predgpi + A3 join (SLURM, minutes).

**Priority.** P1. **Novelty.** Hydrophobin absence in yeasts is known. The fruit-body/hymenium contrast at 291 genera and the GPI-no-domain replacement test are less explored.

**Evolutionary framing.** (1) Polarity: hydrophobins are absent from zoosporic lineages (0 %), Zoopagomycota s.l. (1 %) and Mucoromycota s.l. (3 %) in the pilot, so they are a Dikarya gain (to check: reported as Dikarya-specific), and their absence in yeasts is a loss within Dikarya. Class I (PF01185) dominates Agaricomycetes, class II (PF06766) Pezizomycotina (pilot: 71 % present vs 4 %); which class is ancestral in Dikarya is a polarity question that the genus tree can address by ancestral presence at the Dikarya node (both classes present in some Pezizomycotina genomes). Adhesin families (Flo11, Flocculin, ALS) are Saccharomycotina gains. (2) Placement: hydrophobin gain on the Dikarya stem; losses on each yeast stem (5 independent); adhesin gains on the Saccharomycotina stem or crown; fruit-body-related expansion in Agaricomycetes placed by `bayou`. (3) Mode: gain = family birth (no detectable homolog outside Dikarya); loss = contraction to zero at the yeast transition (burst); expansion in Agaricomycetes by tandem duplication (wgd.ks young nodes within the family). (4) Methods: `make.simmap` on hydrophobin presence; `fitPagel` with growth form; PGLS of log count on fruit-body type with proteome size; CAFE-style after dating. (5) Order of events: along the Saccharomycotina path, reconstruct hydrophobin count and adhesin count; the hypothesis predicts hydrophobin loss on the stem and adhesin family birth later on crown branches (a swap, not a simultaneous replacement). For Agaricomycetes, hydrophobin expansion before or after the first fruit-body origin. (6) Caveats: hydrophobins are short and often unannotated without RNA-seq, so an apparent loss on a poorly annotated branch (Saccharomycotina BUSCO 75-84 %) needs a sequence-level check; the 12 % presence in lichens may be an annotation artefact; family birth is defined by Pfam detection limits.

#### C3. Septal-pore, Woronin-body and cell-fusion gene sets are present only in lineages with septate, fusing hyphae; Agaricomycetes replace the Pezizomycotina set rather than lack it

**Results:** [HS_C3_woronin_fusion.md](HS_C3_woronin_fusion.md) (built 2026-11-08; gene set changed to SwissProt entries from several fungi because SwissProt lacks the *N. crassa* SPA and HAM proteins; the septal-pore module is not testable).

**Statement.** A curated *Neurospora crassa* gene set for Woronin bodies (HEX-1, LAH-1/LAH-2, WSC), septal-pore proteins (SPA-1..SPA-10 to check), and cell-fusion/chemotropism (SO/HAM-1, HAM-5..HAM-11, MAK-1/MAK-2 MAPK modules, PRM-1) has anchors in >= 90 % of filamentous Pezizomycotina genera, is lost as a block in Saccharomycotina and in Pezizomycotina dimorphic genera partly, and is absent from Mucoromycota s.l.; Agaricomycetes lack HEX-1 but keep the fusion/MAPK modules.

**Rationale.** Woronin bodies and septal plugging are Pezizomycotina innovations; Agaricomycetes have dolipore septa with different proteins (to check). Cell fusion machinery is shared more widely. The test asks which sub-module tracks septation and which tracks fusion.

**Traits.** A9 `protein_model_anchor.parquet` gives the best NEUCR hit per protein (235.8 M rows, protein x reference). Build the gene list as NEUCR SwissProt accessions (to check each entry; HEX1_NEUCR exists; SPA and HAM proteins may be in TrEMBL only, not SwissProt: then the test needs a custom BLAST, which is not in tables). Presence = strict anchor (qcov, scov >= 80 %, pident >= 40 %; relaxed to pident >= 25 % for divergent clades).

**Contrast and unit.** Genus presence. Groups as in C1.

**Confounders.** Anchor coverage drops with distance from *N. crassa* (A9 shows 620 strict NEUCR anchors in a median filamentous genome vs 282 in yeasts); report relaxed-gate presence and the clade's overall anchor rate as a covariate; absence in Basidiomycota may be divergence, not loss. SwissProt coverage of the gene set is the main risk.

**Views.** (1) Presence matrix (genes x groups, cell = fraction of genera with an anchor). (2) Module-level score by group.

**Compute.** A9 filtered scan on about 40 accessions (3.2 GB file; SLURM, minutes).

**Priority.** P2. **Novelty.** Woronin body phylogenetic distribution is known (to check); the sub-module decomposition across all BFD groups is less explored.

**Evolutionary framing.** (1) Polarity: Woronin-body proteins are a Pezizomycotina gain (HEX-1 absent outside; to check); septal-pore proteins follow septation, which is a Dikarya gain relative to coenocytic Mucoromycota (35 genera) and Zoopagomycota s.l. (28); cell-fusion MAPK modules are ancestral and shared with yeasts (the pheromone MAPK pathway). Each sub-module therefore has a different polarity, and the test is to confirm these from ancestral presence on the genus tree. (2) Placement: Woronin set gained on the Pezizomycotina stem; septal set on the Dikarya stem, with a Basidiomycota replacement (dolipore) that will show as absence of the Pezizomycotina proteins on the Basidiomycota stem; fusion set lost on yeast stems partly. (3) Mode: gene birth for Woronin proteins; loss on yeast stems; no expansion expected. (4) Methods: `make.simmap` per gene; `fitPagel` for growth form vs module presence; phylogenetic logistic regression with the clade anchor rate as covariate. (5) Order of events: along the Saccharomycotina path, which sub-module is lost first (septal pore proteins before fusion proteins?); along the Onygenales dimorphic path the same question with a within-class control. (6) Caveats: anchors to *N. crassa* fade with distance (A9 strict anchors: 620 median in filamentous vs 282 in yeasts), so "absence" in Basidiomycota and early-diverging lineages may be divergence; SwissProt may lack several SPA and HAM proteins.

### Theme D: Protein architecture

#### P1. Genome reduction shortens proteins within orthologs by removing disordered and linker segments; in large-genome obligate biotrophs the short raw median is a gene-model artefact

**Statement.** (a) Within ortholog families (H2 proxies P-B anchored to SCHPO or YEAST, and P-C same architecture), Microsporidia proteins are shorter than those of free-living relatives by more than the Saccharomycotina effect (which is about 1.6 %), and the loss is concentrated in disordered regions and inter-domain linkers (domain count unchanged). (b) In Pucciniales and Glomeromycota, the low raw median protein length (302 and 227 aa against about 400 aa) is not seen within ortholog families; it comes from many short, domain-less, often TE-derived gene models.

**Rationale.** Pilot: Microsporidia median protein 251 aa, 12.5 % of proteins < 100 aa, multi-domain fraction 0.18, median IDP fraction 0.316. Published Microsporidia proteins are shorter than homologs (to check). Rust and AM fungal genomes are large and TE-rich with fragmented gene models; the within-family test separates biology from annotation.

**Traits.** H2 machinery: A4 `fam_anchor_genome.parquet` (family = reference accession, per-genome median length), A3 `protein_core` (length, n_domains, has_arch), `idp_summary` (per-protein disorder fraction; joins on gene_id form), `swissprot_best` alignment coordinates (q_start, q_end) to locate the missing segment (N-terminal, internal, C-terminal; the H2 report lists this as next step 3).

**Contrast and unit.** Family (log2 ratio group / sister) and genus (offset). Groups: Microsporidia (23 genera; size < 10 Mb), Pucciniales (5), Glomeromycota (11), Erysiphaceae (relaxed filter), *Pneumocystis* (1) each vs the H2 filamentous reference. Add the per-family disorder ratio and domain-count ratio.

**Confounders.** Microsporidia divergence gives few strict anchors (relax pident to >= 25 %, report family counts); gene-model truncation (Q-model; the filter grid of V12); the ortholog proxies mix paralogs (as in H2); IDP predictor bias with amino-acid composition and GC (control by within-family comparison and GC covariate); for rusts, dikaryotic duplicates (do not use Q-dup; use per-genome median per family).

**Views.** (1) Forest plot: within-family length ratio per reduced lineage with CI, beside the raw median ratio. (2) Per-family scatter: domain-count ratio vs length ratio (Microsporidia). (3) Location of the deleted segment: histogram of alignment start/end offsets. (4) Rust/AM: fraction of proteome that is domain-less and < 150 aa, vs masked_pct.

**Compute.** A4 anchor families (small after filter); A3 + idp_summary join for about 60 genera (SLURM, minutes).

**Priority.** P1. **Novelty.** Microsporidia protein shortening is reported (to check). The disorder/linker localization and the rust/AM artefact test are less explored, and the artefact part protects every other view on these groups.

**Evolutionary framing.** (1) Polarity: the ancestral protein length per family is the filamentous and early-diverging median (zoosporic 398 aa, Mucoromycota 364 aa, filamentous Pezizomycotina 405 aa raw medians); shortening is derived on the Microsporidia stem. Outgroups: Rozellomycota (1 genus, *Rozella*), Aphelidiomycota (2), Chytridiomycota (32). For rusts and AM fungi the sisters are Microbotryomycetes (5) and Mortierellomycota (12). (2) Placement: one shift on the Microsporidia stem inherited by all 23 genera (single origin) vs continued shortening along crown branches (gradual); the rust/AM case predicts no within-family shift at all. (3) Mode: loss of protein regions (disorder, linkers), not domain loss (7 of 301 Saccharomycotina shrinkers lose a domain in H2b; the Microsporidia prediction is the same); coupled with genome reduction (G1); no WGD. (4) Methods: PGLS of the genus offset on genome size and lifestyle; `bayou` on the offset to find whether one shift (stem) or several (crown) explain the Microsporidia pattern; `phylosig` on the offset. (5) Order of events: along the path chytrids -> Rozellomycota -> Microsporidia stem -> crown, reconstruct genome size, gene density and the protein-length offset; the hypothesis predicts genome compaction (gene density up) at or before the offset change, as H2b found inside Saccharomycotina (compaction predicts offset within the clade). (6) Caveats: Microsporidia is the longest branch in the kingdom, so strict anchors are few and the offset rests on conserved families only; the 23-genus crown gives little within-clade variation in genome size (most are 2-6 Mb); rust and AM assemblies are dikaryotic or multinucleate and fragmented.

#### P2. Intrinsic disorder is lower in every independent yeast origin and in reduced parasites than in filamentous sisters, after GC control

**Results:** [HS_P2_disorder.md](HS_P2_disorder.md) (built 2026-11-08; Microsporidia excluded per PI decision, so the "reduced parasites" part is not tested).

**Statement.** The median per-protein disorder fraction (idp_summary) is lower in each yeast origin than in its filamentous sister (pilot pooled yeasts 0.335 vs 0.401), lowest in Microsporidia (0.316), and highest in dimorphic Pezizomycotina (0.442) and lichens (0.418). The effect survives a within-ortholog-family comparison and a GC covariate. Expected: sign consistent in >= 4 of 5 yeast origins.

**Rationale.** Disordered regions carry regulatory and interaction sites; a larger, more differentiated organism may need more of them (to check: disorder correlates with organismal complexity in eukaryotes). But disorder predictors also respond to amino-acid composition, which tracks GC (H6: aa composition PC1 vs GC rho = -0.86), so GC control is mandatory.

**Traits.** A5 `idp_genome.parquet` (median_idp_fraction, frac_idp_ge50; 23,412 genomes); per-protein `idp_summary` joined to A3 families for the within-family test; `gc_pct` from labels; GC3 from codon_freq.

**Contrast and unit.** Genus medians; 5 yeast origins vs sisters; Eurotiomycetes dimorphic vs filamentous; reduced parasites vs sisters; mixed model with class/order random intercepts and `gc_pct`, `busco_complete`, `frac_with_pfam` as covariates. PGLS when the tree exists.

**Confounders.** GC and amino-acid composition (covariate; and a within-family test); predictor (`idp` tool not named in the manifest; check which predictor; AIUPred appears in `BFD_function/results/function/aiupred`); gene-model truncation (truncated proteins lose termini, which are often disordered: Q-model); proteome composition (secreted small proteins are often disordered: compute also on the intracellular subset only).

**Views.** (1) Strip plot: median IDP fraction by group and by yeast origin with sisters. (2) Scatter: IDP fraction vs gc_pct (genus points, colour = growth form) with residuals. (3) Within-family IDP ratio histogram (as V08). (4) Dimorphic Eurotiomycetes vs filamentous Eurotiomycetes paired plot.

**Compute.** A5 (small) for (1)-(2); idp_summary 2.1 GB join (SLURM, minutes) for (3). V32 in PLAN is a P3 version of (1) only; this hypothesis replaces it.

**Priority.** P1. **Novelty.** Disorder vs complexity is argued across kingdoms; the replicated within-fungi test with GC control is less explored.

**Evolutionary framing.** (1) Polarity: early-diverging lineages have median IDP fractions of 0.378 (zoosporic) and 0.384-0.396 (Mucoromycota s.l., Zoopagomycota s.l.), close to filamentous Pezizomycotina (0.401) and Agaricomycetes (0.414); yeasts (0.335) and Microsporidia (0.316) are below, dimorphic Pezizomycotina (0.442) above. So the ancestral state is about 0.38-0.40, and both yeast origins and reduced parasites are derived decreases; the dimorphic increase is a derived gain that needs a GC check first. (2) Placement: decreases on each yeast stem (5 origins) and on the Microsporidia stem; `bayou` finds whether the decrease is one shift per stem or spread along crown branches. (3) Mode: loss of disordered sequence inside proteins (H2b: 188 of 301 Saccharomycotina shrinker families lose >= 10 disordered residues), so the mode is region deletion, not gene loss; possibly coupled to compaction (G1, P1). (4) Methods: `phylosig` on IDP fraction and on GC (both expected to carry strong signal); PGLS with GC as covariate and growth form as factor; `OUwie` two-regime; ancestral reconstruction of IDP fraction and GC on the same tree to test whether they changed on the same branches. (5) Order of events: along each yeast path, reconstruct GC, genome size and IDP fraction; if IDP drops on a branch where GC does not change, the disorder loss is not a composition artefact; if compaction precedes the IDP drop, disorder loss follows genome reduction. (6) Caveats: the disorder predictor responds to composition; Saccharomycotina and Microsporidia are long branches with lower annotation completeness; truncated models lose termini (Q-model applied).

#### P3 (falsification test). Complex multicellularity does not raise domain complexity: domain-architecture richness and domain promiscuity scale with proteome size only, and Agaricomycetes with fruiting bodies sit on or below the allometric line

**Statement.** The popular idea: organismal complexity goes with more multi-domain proteins and more domain combinations ("domain promiscuity"). Prediction supporting it: positive residuals of (i) distinct Pfam architectures per genome, (ii) distinct domain-pair combinations, (iii) multi-domain fraction, after proteome-size allometry, in fruit-body-forming Agaricomycetes and Pezizomycotina relative to genera with no fruit body and to yeasts. Prediction that falsifies it: residuals near zero or negative. Pilot: multi-domain fraction is 0.27 in Agaricomycetes, 0.34 in other filamentous fungi, 0.40 in yeasts; fraction of proteins with any Pfam domain 0.60 / 0.73 / 0.85 (also in H2 report, section 6).

**Rationale.** Agaricomycetes proteomes carry many lineage-specific proteins without domains (small secreted, effector-like, TE-derived), which dilutes domain-based complexity measures. The test must therefore use both the whole proteome and the Pfam-annotated subset.

**Traits.** A3 `protein_pfam_arch` (pfam_arch string, n_domains): per genome, count distinct architectures, count distinct unordered domain pairs, mean domains per annotated protein, and the promiscuity index (number of distinct partner domains per domain, averaged over domains present in >= 2 architectures). A4 `genome_protein_summary` (frac_multidomain, median_domains_if_any) already holds the first two summaries.

**Contrast and unit.** Log-log allometry on genus medians; residuals by Agaricomycetes fruit-body type and hymenium type; by Ascomycota fruit-body type within class; yeast origins. Random intercepts for class/order.

**Confounders.** Proteome size (allometry); Pfam coverage differences by clade (compute on annotated proteins and report `frac_with_pfam`); gene-model fragmentation reduces domains per protein (Q-model); TE-derived models (masked_pct covariate; exclude architectures made only of TE domains such as PF00078 RVT_1, PF03184 DDE_1); phylogeny.

**Views.** (1) Allometry panels: architectures, domain pairs, promiscuity vs proteome size. (2) Residual strips by fruit-body type and hymenium type. (3) Residual strip by yeast origin.

**Compute.** A3 `protein_pfam_arch` (1.4 GB; one SLURM pass, minutes).

**Priority.** P2. **Novelty.** Promiscuity-complexity links are published for animals and plants (to check); a fungal test against fruit-body type is less explored.

**Evolutionary framing.** (1) Polarity: multi-domain fraction is 0.32-0.34 in zoosporic, Mucoromycota s.l., Zoopagomycota s.l. and filamentous Pezizomycotina, so that is the ancestral value; Agaricomycetes (0.27) and Microsporidia (0.18) are derived decreases, yeasts (0.40) a derived increase. The Agaricomycetes decrease must be polarized against Tremellomycetes and Pucciniomycotina (values not yet extracted) to place it on the Agaricomycetes stem or deeper. (2) Placement: `bayou` on the allometric residual of architecture richness; the falsification prediction is that no shift maps to fruit-body origins; a shift on the Agaricomycetes stem would be a class-wide dilution event. (3) Mode: dilution by birth of domain-less proteins (orphan gain, TE-derived models) rather than loss of multi-domain proteins; test by reconstructing the absolute count of multi-domain proteins (expected constant) and the count of domain-less proteins (expected to rise) separately. (4) Methods: PGLS of log richness on log proteome size and fruit-body type; `phylosig`; `bayou` on residuals; CAFE-style on architecture classes after dating. (5) Order of events: along the path to agaricoid genera, reconstruct domain-less protein count and architecture count; if the domain-less count rises before fruit-body origins, complexity did not drive it. (6) Caveats: Pfam coverage differs by clade (Agaricomycetes 0.60 vs Pezizomycotina 0.73 with any domain); TE-derived models inflate domain-less counts on repeat-rich branches; fruit-body origins are genus-level labels.

### Theme E: Genome architecture

#### G1. Obligate parasitism moves genome size in two directions set by the host-cell boundary: intracellular parasites shrink, haustorial extracellular biotrophs expand through repeats, and both lose gene content

**Results:** [HS_G1_genome_size_parasitism.md](HS_G1_genome_size_parasitism.md) (built 2026-11-08).

**Statement.** Lineages that live inside host cells (Microsporidia 23 genera; Rozellomycota 2; *Pneumocystis* 1; *Olpidium* 1; endobiotic chytrids such as *Synchytrium* 2; to check which are truly intracellular) have genomes < 15 Mb, high gene density, few introns and short intergenic distances. Lineages that feed through haustoria from outside the host protoplast (Pucciniales 5, Erysiphaceae 7 raw, Glomeromycota 11, Entomophthorales 6; to check) have genomes > 50 Mb with high masked fraction, low gene density, and intron counts like their sisters. Both groups have fewer Pfam families than sisters. Expected: a bimodal genome-size distribution among obligate biotrophs with no intermediate.

**Rationale.** Intracellular life imposes a small cell and nutrient import, which favors deletion; extracellular biotrophs lose selection on repeat removal (to check: the "two-speed" and repeat-expansion literature for rusts and mildews). BFD gives several independent origins for each side.

**Traits.** `total_length_bp`, `masked_pct`, genes per Mb (A2 n_protein_coding / size), introns per transcript, median intergenic (A2), Pfam family richness (A7 count per genome, V29), BUSCO completeness as covariate. All small tables.

**Contrast and unit.** Genus medians; each lineage vs its free-living sister (as in T2); two-group comparison (intracellular vs haustorial) with lineage as the replicate (about 5 vs 4).

**Confounders.** Metagenome-derived Microsporidia bins (size filter); assembly fragmentation shortens intergenic distances (Q-n50); masking method may differ by annotation wave (not verified; PLAN V18 caveat); rust assemblies are dikaryotic (size can be doubled: use `busco_duplicated` to flag and halve where > 50 %); Erysiphaceae fail the class BUSCO filter (relaxed filter).

**Views.** (1) Scatter: genome size vs genes per Mb, log-log, all genera grey, obligate lineages coloured by intracellular/haustorial. (2) Dot plot per lineage: size ratio to sister, with masked_pct. (3) Intron and intergenic panels per lineage. (4) Pfam family richness residual per lineage.

**Compute.** Small tables only (login node).

**Priority.** P1. **Novelty.** Both directions are known for single lineages (to check). The paired, replicated contrast across the kingdom in one figure is less explored, and it is cheap.

**Evolutionary framing.** (1) Polarity: the ancestral fungal genome is about 25-40 Mb with 8,000-12,000 genes (BFD medians: zoosporic 31.5 Mb, Mucoromycota s.l. 37.7, Zoopagomycota s.l. 25.7, filamentous Pezizomycotina 37.9). Microsporidia (8.5 Mb median, 2-6 Mb for most) are a derived reduction; rusts (166 Mb) and AM fungi (156 Mb) are derived expansions. Outgroups: for Microsporidia, Rozellomycota (1 genus, 9.3 Mb raw median) and Aphelidiomycota (2, 15.7 Mb) are already small, so part of the reduction may predate Microsporidia (to check; ancestral reconstruction at the Rozellomycota + Microsporidia node); for rusts, Microbotryomycetes (5 genera); for AM fungi, Mortierellomycota (12) and Mucoromycota (35). (2) Placement: each lineage is one origin; the test is across lineages (about 5 intracellular, about 4 haustorial), so the unit is the stem branch of each lineage; `bayou` on log genome size locates the shift branches without a lifestyle label. (3) Mode: reduction by gene loss and intergenic deletion (burst on the stem, then stasis); expansion by TE accumulation (masked_pct up) with gene count roughly constant; gene content loss in both (Pfam richness down); no WGD, although rust and AM genomes have high BUSCO duplication for assembly reasons. wgd.ks `n_genes_in_identical_groups` separates haplotig duplication from true duplication. (4) Methods: `bayou` and `OUwie` on log genome size, genes per Mb, masked_pct; PGLS with the hand-made `obligate_biotroph_mode` factor; `fitPagel` for lifestyle vs size class (binary: < 15 Mb, > 50 Mb). A dated tree is needed to say whether expansion was fast. (5) Order of events: along each path, reconstruct genome size, gene count, intron density and Pfam richness; the hypothesis predicts gene-content loss on the stem in both directions, with size change after (intracellular: size falls after gene loss; haustorial: size rises after gene loss). (6) Caveats: Microsporidia bins (size filter); rust assemblies may be diploid-sized; Erysiphaceae fail the class filter; a 1-genus outgroup (Rozellomycota) gives a weak polarization; long branches in Microsporidia bias rate estimates.

#### G2. Yeast proteome shrinkage is family loss, not copy loss: the yeast proteome has fewer Pfam families and fewer domain-less proteins, but the same copies per retained family

**Statement.** Decompose the proteome-size difference between each yeast origin and its sister into three parts: (i) Pfam families absent, (ii) copies per shared family, (iii) proteins without any Pfam domain. Expected: parts (i) and (iii) explain most of the deficit; part (ii) is near zero or positive. The same decomposition for Microsporidia is expected to show copy loss as well (families reduced to one copy).

**Rationale.** V29 shows yeast genera have 4,062 vs 4,755 families, and H2 shows yeasts have more multi-domain and fewer domain-less proteins (85 % with a domain vs 70 %). V39 shows no difference in young paralog retention. A formal decomposition tells which process built the small proteome.

**Traits.** A7 per genome: n_families, sum(n_proteins) over families, per-family copies; A4 `frac_with_pfam`; A2 `n_proteins`. Shared-family copy ratio = median over shared families of log2(copies_yeast / copies_sister) with family size bins (families with 1, 2-5, > 5 copies in the sister).

**Contrast and unit.** Genus medians; 5 yeast origins vs sisters; dimorphic Eurotiomycetes vs filamentous; Microsporidia vs chytrids; Agaricomycetes vs Pezizomycotina (same decomposition answers where the Agaricomycetes proteome is larger: pilot says domain-less proteins).

**Confounders.** Pfam detection (families with weak HMMs appear "lost" in divergent clades; use the A7 e-value gate and report sensitivity at 1e-3); gene-model fragmentation inflates copies (Q-model, Q-dup); TE-derived families (exclude TE domains); proteome size (the decomposition is in absolute proteins and in fractions).

**Views.** (1) Stacked bar per origin: proteome deficit split into the three parts. (2) Histogram of per-family copy log-ratio by family size bin. (3) Table: top 30 families lost in >= 4 yeast origins (extends V30, which pools yeasts).

**Compute.** A7 full scan grouped per genome (355 MB, 1.6 M rows for 100 families in 2.6 s; the full file needs a SLURM pass of minutes; the per-genome aggregates are small).

**Priority.** P1. **Novelty.** Less explored as an explicit decomposition; cheap, and it reuses A7.

**Evolutionary framing.** (1) Polarity: the ancestral Dikarya proteome holds about 10,000-12,000 proteins and about 4,700 Pfam families (V29 filamentous value; early-diverging lineages not yet measured for family richness). The yeast proteome (5,700 proteins, 4,062 families) is a derived reduction; Agaricomycetes (12,562) a small derived increase. Outgroups as in T1. (2) Placement: 5 yeast stems; family loss events can be placed per family with `make.simmap` on presence, which gives for each lost family the branch of loss and whether the loss is shared by the whole origin (stem) or repeated in crown lineages. (3) Mode: family death (contraction to zero) and orphan loss on stems; copy number per retained family constant (no contraction) and, in Saccharomycetaceae, a WGD (wgd.ks peak; to check) that did not restore proteome size, which argues for strong loss after duplication. The Saccharomycetaceae WGD is the one place in this file where the wgd.ks tables answer a direct question: the fraction of WGD-age paralog nodes retained vs the family-loss count on the same branch. (4) Methods: CAFE-style birth-death on Pfam family counts (needs a dated tree) with yeast stems as separate rate classes; `make.simmap` on presence per family; PGLS of family richness on proteome size. (5) Order of events: along the Saccharomycotina path, ancestral family richness and ancestral copies per family; the hypothesis predicts richness falls on the stem while copies per family stay flat until the WGD. (6) Caveats: family "loss" on a long branch may be detection failure (second detector needed for Saccharomycotina-specific divergence); gene-model fragmentation raises copies; Pfam families born in Saccharomycotina (V30 lists POC1, STD1, MPM1, Gsf2) are gains that the decomposition must count separately.

#### G3 (falsification test). Unicellular growth does not by itself drive intron loss: basidiomycete yeast origins keep filamentous-like intron densities

**Results:** [HS_G3_intron_loss_by_origin.md](HS_G3_intron_loss_by_origin.md) (built 2026-11-08).

**Statement.** The popular reading of H3: yeasts have few introns because they are unicellular. Prediction that supports it: lower introns per transcript than the sister in each yeast origin. Prediction that falsifies it: Tremellomycetes (15 genera), Pucciniomycotina yeasts (7), Malassezia (1), Ustilaginomycotina dimorphic (20) and Pezizomycotina dimorphic (15 in Eurotiomycetes) keep intron densities like their sisters, so the loss is a Saccharomycotina (and Schizosaccharomyces) trait. Also: coenocytic Mucoromycota (3.36 per transcript, pilot) are intron-rich, and Microsporidia (0.21) are intron-poor, which separates genome reduction from unicellularity.

**Rationale.** V15 reports yeast vs filamentous by class and pooled (0.24 vs 2.03), and the mixed model with class random intercepts gave a yeast effect of 0.001 (CI -0.39 to 0.39), which already hints that the pooled effect is clade-bound. The by-origin forest plot is not built.

**Traits.** A2 `introns_per_transcript`, `frac_intronless`, `median_intron_len`. Small.

**Contrast and unit.** Genus medians; per-origin difference to sister with bootstrap CI; sign test across origins; Microsporidia and Mucoromycota as outgroups on the same axis.

**Confounders.** Ab initio intron under-prediction without RNA-seq (D8; not available; `busco_complete` and `frac_gt_ag` as proxies); Tremellomycetes assemblies are good (median BUSCO 96.4), Saccharomycotina lower (75-84); report both.

**Views.** (1) Forest plot: introns per transcript, yeast origin minus sister, with CI. (2) Strip: introns per transcript for all groups in section 1.7 (adds Microsporidia, zoosporic, coenocytic to V15). (3) Intron length by group.

**Compute.** Small tables (login node).

**Priority.** P2 (cheap; extends V15). **Novelty.** The Saccharomycotina intron loss is well known. The replicated-origin framing with basidiomycete yeasts is a clean falsification test and is cheap.

**Evolutionary framing.** (1) Polarity: the ancestral fungal gene is intron-rich (zoosporic 3.57 per transcript, Mucoromycota s.l. 3.36, Agaricomycetes 4.91); Pezizomycotina (1.79) and Zoopagomycota s.l. (1.90) are already reduced; Saccharomycotina (0.24 pooled yeast value) and Microsporidia (0.21) are extreme reductions. Polarity for Saccharomycotina: Taphrinomycotina (*Schizosaccharomyces* about 1 per gene, to check with A2) and Pezizomycotina. (2) Placement: intron loss is expected on the Saccharomycotina stem (one event inherited by 144 genera) and the Microsporidia stem; the falsification test predicts no shift on the Tremellomycetes, Pucciniomycotina or Ustilaginomycotina yeast stems. (3) Mode: intron loss in Saccharomycotina is gradual within the clade according to the literature (to check); `bayou` separates a stem burst from crown decline; no WGD link (the Saccharomycetaceae WGD post-dates most loss, to check). (4) Methods: `phylosig` (expect lambda near 1); PGLS on growth form; `bayou` for shift branches; `fitPagel` between growth form and intron class (binary: < 0.5 per transcript). (5) Order of events: along the Saccharomycotina path, reconstruct growth form (simmap), genome size and intron density; if the intron drop maps to the same stem as the growth-form change but basidiomycete yeast stems show no drop, the loss is a Saccharomycotina event, not a yeast event. (6) Caveats: intron prediction depends on ab initio training without RNA-seq (D8); Taphrinomycotina 3 genera; the Saccharomycotina long branch.

#### G4. Codon-usage selection is stronger in fast-growing unicellular forms and weaker in reduced parasites, in independent origins

**Statement.** The departure of codon usage from the GC-expected neutral curve (effective number of codons vs GC3, Wright's curve; to check) and the match between tRNA gene copies and codon frequencies (V36: yeast rho 0.39 vs filamentous 0.23) are larger in each yeast origin than in its sister, and smallest in Microsporidia and other obligate parasites (small effective population size), independent of genome GC.

**Rationale.** Translational selection scales with growth rate and effective population size. Yeast origins and parasites give replicated contrasts in both directions.

**Traits.** `codon_freq` (whole-proteome, 64 codons per genome) gives a genome-level ENC and GC3; A5 `trna_codon_counts`. A per-gene ENC (which separates highly expressed genes) needs `gene_CDS` (16 GB) and is a SLURM job (not measured); the genome-level version is on the login node.

**Contrast and unit.** Genus medians; yeast origins vs sisters; parasite lineages vs sisters; GC3 as covariate (V34 shows yeasts sit 3.2 points below the GC3 regression line).

**Confounders.** Whole-proteome codon frequencies are dominated by low-expression genes (per-gene ENC needed for a strong test); GC (the ENC-GC3 curve is the control); tRNA undercount on fragmented assemblies (Q-n50).

**Views.** (1) ENC vs GC3 with Wright's curve, genus points by group. (2) Forest plot: ENC departure per yeast origin and per parasite lineage.

**Compute.** Small tables for the genome-level version; `gene_CDS` pass for the per-gene version (SLURM; not estimated).

**Priority.** P3. **Novelty.** Codon selection in *Saccharomyces* and *Candida* is well known; replicated origins and the parasite side are less explored (to check).

**Evolutionary framing.** (1) Polarity: codon usage tracks class (H6: CLASS explains 57 % of RSCU variance), so the ancestral state is clade-specific and must be reconstructed per clade; the trait to polarize is the ENC departure from the GC3 curve, with early-diverging lineages (Mucoromycota 35, zoosporic 43 genera) and Taphrinomycotina (3) as outgroups for the Saccharomycotina change. (2) Placement: shifts on yeast stems (selection up) and parasite stems (selection down); H6 V37 already places the CTG reassignment inside Pichiomycetes (code 12, 1,154 genomes). (3) Mode: gradual (selection acts per gene over long periods), so a stem burst is not expected; the test is a rate or optimum shift under an OU model, not a step. (4) Methods: PGLS of ENC departure on growth form with GC3 as covariate; `OUwie` regimes from the growth-form simmap; `phylosig`. (5) Order of events: along the Saccharomycotina path, reconstruct GC3, ENC departure and tRNA-codon correlation; if tRNA pool matching changes after the ENC departure, selection on codons preceded tRNA gene adjustment. (6) Caveats: whole-proteome frequencies; tRNA counts depend on assembly contiguity; the Microsporidia branch is long and GC-poor.

### Theme F: Existing results, re-framed as evolutionary test cases

The results in H2, H2b, H3, H4 and H6 are trait-vs-growth-form contrasts. Each becomes a polarity question once the genus tree exists. The tests below use the same genus-median tables that those pages already wrote (`tables/V09_*`, `V15_*`, `V20_*`, `V33_*`, `V50_*`, `V55_*`, `V56_*`), so each is one PGLS or ancestral-state run.

**H2 / H2b: ortholog shortening and compaction in Saccharomycotina.** Result: Saccharomycotina proteins are 1.6 % shorter within families (118 genera); basidiomycete yeasts are not shorter (+0.2 to +1.5 %); inside Saccharomycotina a more compact genome predicts shorter orthologs (rho -0.43 with genes per Mb); 301 consistent shrinker families lose mostly internal or N-terminal disordered sequence (188 of 301 lose >= 10 disordered residues; 7 lose a domain); 54 consistent gainers carry N-terminal extensions with no domain gain. Polarity question: did Saccharomycotina proteins shrink, or did Pezizomycotina proteins grow? The filamentous reference in V08/V09 is Pezizomycotina plus other filamentous genera, which assumes the Pezizomycotina length is ancestral. Test: reconstruct the per-genus offset at the Ascomycota node using Taphrinomycotina (3 genera; *Schizosaccharomyces* offset -0.043 in H2 report section 5, which is more negative than Saccharomycotina) and Basidiomycota (432 raw genera) and the early-diverging lineages as outgroups. If the ancestral Ascomycota offset is near the Saccharomycotina value, Pezizomycotina proteins lengthened (gain of N-terminal and linker sequence) and the "yeast shortening" is a Pezizomycotina extension. The Taphrinomycotina value points that way and must be checked per family. Placement: `bayou` on the genus offset locates the shift branch (Saccharomycotina stem vs Pezizomycotina stem vs both). Mode and order: along the Saccharomycotina path, reconstruct genes per Mb, genome size and the offset (H2b V50 shows the within-clade coupling); the order test is whether compaction (genes per Mb) shifts on an earlier branch than the offset. For the shrinker families, map the per-family loss on the tree with `fastAnc` to find whether loss is stem-wide (one deletion) or repeated (several deletions in crown lineages; `n_origins_shorter` in V55 is the current proxy and is 1 for most families). For the gainers, the N-terminal extensions are candidate gains on the Saccharomycotina stem (or losses in Pezizomycotina); the Taphrinomycotina state decides. Caveats: Taphrinomycotina is 3 genera; Saccharomycotina is a long branch (fewer strict anchors); the compaction signal inside the clade is small (0.049 log2 per log10 bp).

**H3: intron loss.** Result: yeast genera 0.24 introns per transcript vs 2.03 filamentous; the class/order mixed model gives a yeast effect near zero, so the pooled effect is clade-bound. Polarity: ancestral fungi are intron-rich (early-diverging lineages 3.4-3.6 per transcript; Agaricomycetes 4.9), Pezizomycotina are intermediate (1.8), Saccharomycotina extreme. Test: ancestral intron density at the Ascomycota, Pezizomycotina and Saccharomycotina nodes; `bayou` shift branches; G3 above adds the basidiomycete yeast origins as the replicated falsification. Order: intron loss vs genome compaction vs growth-form change on the Saccharomycotina stem (H2b V50 found no within-Saccharomycotina association between intron density and the length offset, rho -0.01, so introns and protein length changed on different branches or by different processes).

**H4: secretome loss in yeasts.** Result: strict secretome 2.6 % in yeasts vs 5.9 % in filamentous; the mixed model keeps a yeast effect below zero (fraction -0.013, CI -0.023 to -0.004), so the deficit is not only class membership. Polarity: hyphal secretion is ancestral (early-diverging lineages 3.5-7.5 %). Test: `fitPagel` between growth form and secretome class across the 5 yeast origins; `OUwie` with yeast and hyphal regimes on secretome fraction; ancestral secretome fraction at each yeast stem. Order: T1 above (secretome loss before or after transporter change); and secretome loss vs intron loss vs compaction on the Saccharomycotina stem from the same ancestral-state runs. A single answer for the Saccharomycotina stem is the most useful product: which of compaction, intron loss, secretome loss and protein shortening map to the stem, and which to crown branches.

**H6: codon usage and the CTG reassignment.** Result: CLASS explains 57 % of RSCU variance, growth form 14 % and nothing within class; code-12 Pichiomycetes use CTG at 5.2 per 1,000 codons vs 9.2 in code-1 Pichiomycetes. Polarity: codon usage is a clade trait inherited along the tree, so the test is phylogenetic signal (`phylosig` on the RSCU PC scores; expect lambda near 1), and the growth-form effect becomes a PGLS term expected to be near zero. The CTG reassignment is one event on the stem of the code-12 clade (1,154 genomes, all Pichiomycetes); the test is ancestral CTG frequency along the path into the code-12 clade: the "codon capture" model predicts CTG frequency fell to near zero before the reassignment (ancestral nodes below the clade already CTG-poor), the "ambiguous intermediate" model predicts no prior drop (to check both). V37 gives the tips; the tree gives the nodes. Caveat: the genus tree holds one genome per genus, so within-genus code variation (if any) is invisible; `TRANSL_TABLE` comes from the samples sheet.

---

## 3. Summary table

| ID | Theme | Statement (short) | Replicates (genera, filtered) | Priority | Falsification test |
|---|---|---|---|---|---|
| T1 | Transport | Transporter-to-secretome ratio rises in each yeast origin; MFS density does not | 5 yeast origins (144, 3, 15, 14, 6 raw) + dimorphic Eurotiomycetes 15 vs 32 | P1 | |
| T2 | Transport | Same inorganic-N/S transporter families lost in independent obligate biotrophs | Pucciniales 5, Erysiphaceae 2 (7 raw), Glomeromycota 11, Microsporidia 23, Pneumocystis 1, Entomophthorales 6 | P1 | |
| T3 | Transport | Symbiosis changes which transporters expand (ECM within orders; lichen photobiont; rhizomorph proxy) | ECM in 6 orders (45/7, 9/95, 7/9, 5/11, 3/5, 3/2); photobiont 124 vs 13 vs 9 | P2 / P3 | |
| S1 | Signaling | Kinase, GTPase, GPCR counts scale sub-linearly; no residual for fruit-body forms | Fruit body 144/45/49 genera; hymenium 100/94/46; 5 yeast origins | P1 | yes |
| S2 | Signaling | Light and two-component sensing lost in yeasts, intracellular parasites, hypogeous fruiters | Yeast origins; parasites; hypogeous 7 + 5 genera in >= 4 orders | P1 | |
| S3 | Signaling | Cyclic-nucleotide and ion-channel module co-lost with the flagellum | Zoosporic 43 genera in 3 phyla vs 28 + 59 non-flagellated | P2 | |
| S4 | Communication | STE3 and G-alpha expansion tracks mating system, not fruit body | Agaricomycetes 291 vs Pezizomycotina 584; mating-system label missing | P2 / P3 | partly |
| C1 | Communication | NLR allorecognition repertoire tracks septate heterokaryotic hyphae | Yeast origins; dimorphic 15 vs 32; Mucoromycota s.l. 59; Agaricomycetes 291 | P1 | |
| C2 | Cell wall | Hydrophobins vs adhesins swap at each shift; no scaling with hymenium complexity | Yeast origins; fruit body 144/45/49; hymenium 100/94/46; lichens 156 | P1 | partly |
| C3 | Communication | Woronin/septal/fusion gene sets present only with septate fusing hyphae | Groups as C1; needs NEUCR accession list | P2 | |
| P1 | Architecture | Reduced genomes shorten proteins in linkers; rust/AM short medians are artefacts | Microsporidia 23; Pucciniales 5; Glomeromycota 11; Erysiphaceae; Pneumocystis | P1 | |
| P2 | Architecture | Disorder lower in each yeast origin and reduced parasites after GC control | 5 yeast origins; dimorphic 15 vs 32; parasites | P1 | |
| P3 | Architecture | Domain-architecture richness and promiscuity do not rise with fruit-body complexity | Fruit body 144/45/49; Ascomycota by class; yeast origins | P2 | yes |
| G1 | Genome | Two directions of genome size under obligate parasitism | Intracellular about 5 lineages; haustorial about 4 | P1 | |
| G2 | Genome | Yeast proteome shrinkage is family loss and orphan loss, not copy loss | 5 yeast origins; Microsporidia; Agaricomycetes | P1 | |
| G3 | Genome | Basidiomycete yeast origins keep intron density | Tremellomycetes 15, Pucciniomycotina 7, Ustilaginomycotina 20 + 6, Pezizomycotina dimorphic 15 | P2 | yes |
| G4 | Genome | Codon selection stronger in yeasts, weaker in parasites | Yeast origins; parasite lineages | P3 | |

Falsification tests: S1, P3, G3 (and parts of S4, C2).

## 4. Shared build notes

- New label columns to add to `genome_labels` (A1 rebuild, login-node scale): `fruitbody_type`, `hymenium_type`, `ecm_exploration_type`, `ecm_lineage`, `primary_photobiont`, `decay_type`, `lifestyle_comment` from FungalTraits; `obligate_biotroph_mode` (intracellular / haustorial / none; hand list); `flagellated` (from growth form); `coenocytic` (from phylum); `mating_system` (hand list; S4 only). Also a `size_ok_microsporidia` flag (< 10 Mb) and a `relaxed_busco_90` flag for Erysiphaceae, Pucciniales and Glomeromycota.
- One new intermediate, A10 `pfam_trait_sets.parquet` (genome x trait set x n_proteins), built from A7 with a curated Pfam set file (`atlas_data/labels/pfam_trait_sets.tsv`: set name, pfam_acc, note, "checked" flag). T1, T2, T3, S1, S2, S3, S4, C1, C2 read it. Build time: one A7 pass (seconds to minutes).
- One A3-scale pass (SLURM, minutes by the build log) serves C1 (NLR architectures), C2 (GPI without domain), P1 and P2 (within-family length and disorder), P3 (architecture richness). Chunk it like J2b.
- PGLS: when `atlas_data/phylo/results/` holds a tree, every genus-median contrast above gets a PGLS line; until then the no-tree controls (independent origins, within-order pairs, class/order random intercepts) are the stated result.
- TCDB: decide whether to run `run_tcdb` on all genomes (about 23,400 blastp jobs; per-genome time not measured). T1-T3 can start on Pfam proxies and be re-run on TC families later.

## 5. Recommended build order

1. **G1 + G3 + P2 (genome-level part)**: small tables only; one login-node session; they also validate the new label columns and the Microsporidia size filter that every other view needs.
2. **S1 + C1 + C2 (A7 part)**: one A10 build from A7, then three pages; S1 is a falsification test with large genus counts (144/45/49 and 100/94/46), C1 and C2 have the largest pilot effects.
3. **P1 + G2**: reuse the H2 family machinery (A4 anchors) and A7; P1 protects every later claim about rusts, AM fungi and Microsporidia from gene-model artefacts.
4. Then T1, T2, S2 (need the tmhmm / A3 joins and the hand-made obligate-biotroph list), then P3, S3, C3, T3, S4, G4.
