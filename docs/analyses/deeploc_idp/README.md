# DeepLoc and intrinsic-disorder (IDR) patterns across the BFD proteomes

Status: **design, controls and decision rules, written before any test was run** (this file is committed first; results are added after).
Data: the local `BFD.duckdb` (`deeploc`, `idp_summary`, `idp`, `signalp`, `tmhmm`, `predgpi`, `wolfpsort`, `gene_proteins`, `pfam`,
`tf_inventory`, `v_secretome_proteins`, `v_protein_annotation`). Predictions are computational throughout.

## What the data audit showed (job 32112426, read-only, kept in `audit/`)

- Keys `(species_prefix, protein_id)` join cleanly to `gene_proteins` (244,812,918 proteins, 23,412 genomes) for every table; no mismatches.
- `deeploc` covers every protein (one model, no duplicate keys). `idp_summary` has a row for 59.8% of proteins. `signalp` and `tmhmm` hold
  positive calls only. `targetp` covers 125 genomes and is not used.
- The share of proteins with an `idp_summary` row rises steeply with length (32.7% under 100 aa, 89.2% at 1,000 aa or more).
  **Every IDR comparison below is therefore made within length bins** (under 100, 100 to 299, 300 to 999, 1,000 and over).
- Summed segment lengths differ from `IDP_residues` for most proteins in a sample (cause not investigated), so **residue counts are not
  used; only `IDP_fraction`** is.

## Unit and measures

- Unit: **one representative genome per genus** (`is_genus_rep`, annotated, 1,654 with proteins): the genomes are not independent but
  this removes within-genus pseudo-replication. Class-level summaries use classes with at least 8 genus representatives (20 classes).
  A phylogenetic correction is not possible (no tree here), so class is used as a stratum and contrasts are made **within genome**.
- IDR measure: `IDR_all0` = sum of `IDP_fraction` over all proteins of a category / number of proteins of the category (a protein
  with no row counts as 0); alongside it the share with a row and the mean among proteins with a row. Primary = `IDR_all0`.
- Matched contrasts: for a pair of categories in one genome, the difference of `IDR_all0` is computed in each length bin and averaged
  with weights equal to the smaller category count in the bin; the genome's contrast is that weighted mean.
- Localization: primary DeepLoc label = class with the highest probability among cytoplasm, nucleus, extracellular, cell membrane,
  mitochondrion, plastid, endoplasmic reticulum, lysosome/vacuole, Golgi, peroxisome.
- Statistics: effect sizes first (median of genome contrasts per class, bootstrap 95% interval over genomes); each hypothesis has one
  primary test (sign test of the class medians); Benjamini-Hochberg over the primary tests of H1 to H7.

## Controls (pre-set pass criteria; a failed control is reported and caveats every result that depends on it)

| Control | Pass criterion |
|---|---|
| C1 transcription factors (`tf_inventory`) are DeepLoc nuclear | at least 90% of TF proteins in at least 90% of classes |
| C2 mitochondrial carriers ([PF00153](https://www.ebi.ac.uk/interpro/entry/pfam/PF00153/)) are DeepLoc mitochondrion | at least 70% of carrier proteins |
| C3 SignalP-positive proteins are DeepLoc extracellular or membrane | at least 60% of strict-secretome proteins |
| C4 known small-secreted effector-like families (CFEM [PF05730](https://www.ebi.ac.uk/interpro/entry/pfam/PF05730/), LysM [PF01476](https://www.ebi.ac.uk/interpro/entry/pfam/PF01476/), CAP [PF00188](https://www.ebi.ac.uk/interpro/entry/pfam/PF00188/), Chitin_bind_1 [PF00187](https://www.ebi.ac.uk/interpro/entry/pfam/PF00187/)) are recovered by the secretion filter | at least 70% of their proteins pass `secreted_strict` |
| C5 predictor sanity: nuclear proteins have higher `IDR_all0` than mitochondrial proteins at matched length | positive in at least 90% of classes |

## Hypotheses and kill rules (primary test per hypothesis; all matched within genome and length bin)

| # | Hypothesis | Prediction | Supported if | Killed if |
|---|---|---|---|---|
| H1 | TFs are more disordered than other nuclear proteins | matched contrast TF minus other nuclear > 0 | positive in >= 90% of classes and pooled median >= 0.05 | positive in < 60% of classes |
| H2 | Unannotated proteins (no Pfam, no CAZy, no MEROPS hit) are more disordered than annotated ones | contrast unannotated minus annotated > 0 | positive in >= 90% of classes | positive in < 60% of classes |
| H3 | Strict-secretome proteins are less disordered than non-secreted | contrast secreted minus other < 0 | negative in >= 90% of classes | negative in < 60% of classes |
| H4 | Genome-wide `IDR_all0` tracks genomic GC (composition) | class-demeaned Spearman rho between genome `IDR_all0` and GC > 0 | rho >= 0.30 | abs(rho) < 0.15 |
| H5 | The nuclear-predicted share rises with proteome size | class-demeaned Spearman rho (nuclear share, log proteins) > 0 | rho >= 0.20 | rho <= 0.05 |
| H6 | Genome-reduced Microsporidia are less disordered than other genomes at matched length | genome contrast Microsporidia vs the rest, matched by length bin | median contrast < 0 and bootstrap interval excludes 0 | median contrast >= 0 |
| H7 | Effector-like candidates are more frequent in plant-pathogen than saprotroph genera of the same class | candidates per 1,000 proteins, plant pathogen minus saprotroph, within class | positive in >= 75% of evaluable classes (at least 3 evaluable) | positive in < 50% of evaluable classes, or fewer than 3 evaluable classes (then "not evaluable") |

Notes on H4 to H7: H4 asks whether IDR predictions mostly reflect amino-acid composition, which would make H1 to H3 partly a composition
effect; if H4 is supported, results are also reported after adjusting for GC. H6 and H7 depend on annotation quality (gene-model
errors in reduced and fragmented genomes can change IDR and secretome calls); this is stated, not corrected.

## Effector-like candidates (transparent rule, not validated effectors)

A protein is an **effector-like candidate** if all of the following hold: strict secretome (SignalP signal peptide, no transmembrane
helix, not GPI-anchored; `v_secretome_proteins.secreted_strict`), length at most 300 aa, no Pfam hit, no CAZy hit, no MEROPS hit.
This describes small secreted proteins of unknown function; many are not effectors, and real effectors can fail the rule (for
example, those with a Pfam hit or longer than 300 aa; C4 measures the second failure for known families). Lifestyle labels for H7
are FungalTraits 1.2 `primary_lifestyle` at genus level (`plant_pathogen`, and saprotrophs = any `*_saprotroph`), taken verbatim from
`atlas_data/labels/fungaltraits_1.2.utf8.csv`; genera without a label are excluded, never guessed.

## Not done

No phylogenetic correction (no tree); no validation against experimental effector sets; no statement about function. Any
literature link between disorder and function is a citation to verify and none is made here.

## Files (planned)

`01_build_cube.py` (HPCC job: genome x length bin x localization x flags table), `02_analyze.py` (controls, H1 to H7, candidates),
`test_deeploc_idp.py` (fixture test), `audit/` (audit script and log), `results/`.
