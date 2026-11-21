# DeepLoc and intrinsic-disorder (IDR) patterns across the BFD proteomes

Status: **done 2026-11-16 (see Results at the end).** The design, controls and decision rules were written, and committed, before any test was run.
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

## Files

`01_build_cube.py` (HPCC job 32132454, 2 min: genome x length bin x localization x flags table), `02_analyze.py` (controls, H1 to H7,
candidate counts), `test_deeploc_idp.py` (4 fixture tests pass on the HPCC), `audit/` (audit script and log), `results/`.
The cube tables and the list of effector-like candidate proteins (245,547 proteins, `ssp_candidates.parquet`) stay in
`bench/deeploc_idp/cube/` on the HPCC.

---
## Results (2026-11-16; 1,669 genus-representative genomes, 18,424,213 proteins)

Order of events: the design and rules above were committed first (`4575c84`), the scripts next (`dbacaf5`), then the cube was built
and the analysis run once. No rule or threshold was changed after the results were seen (one code bug, a column named `loc` clashing
with the pandas `.loc` attribute, was fixed before the first run).

### Controls

| Control | Observed | Criterion | Result |
|---|---|---|---|
| C1 TFs DeepLoc nuclear (share of classes with median >= 0.90) | 1.00 | 0.90 | pass |
| C2 mitochondrial carriers DeepLoc mitochondrion | 0.964 | 0.70 | pass |
| C3 strict secretome DeepLoc extracellular or membrane | 0.909 | 0.60 | pass |
| C4 effector-like families recovered by the secretion filter | **0.431** | 0.70 | **FAIL** |
| C5 nuclear above mitochondrial IDR_all0 (share of classes) | 1.00 | 0.90 | pass |

C4 failed: only 43% of proteins carrying CFEM, LysM, CAP or Chitin_bind_1 domains pass the strict secretion filter. The reason was not
tested (these domains also occur in non-secreted and multidomain proteins, and the signal-peptide call is a single predictor), so it
is not known whether the filter is too strict or the control too broad. Per the pre-registered rule it is a **caveat on H7 and on the
candidate list**: real effectors may be missed.

### Hypotheses (`results/hypotheses.tsv`; q = Benjamini-Hochberg over the seven primary tests)

| # | Result (primary statistic) | Effect size | Verdict by the pre-set rule |
|---|---|---|---|
| H1 TFs more disordered than other nuclear proteins | positive in 19 of 20 classes (p 2e-05, q 3.5e-05) | pooled median +0.124 IDR_all0 | supported |
| H2 unannotated more disordered than annotated | positive in 20 of 20 classes (p 9.5e-07) | pooled median +0.162 | supported |
| H3 strict secretome less disordered than the rest | negative in 18 of 20 classes (p 2e-04) | pooled median -0.137 | supported |
| H4 genome IDR_all0 tracks GC (class-demeaned) | rho 0.188 (n = 1,491 genomes) | between the thresholds 0.15 and 0.30 | **mixed** |
| H5 nuclear share rises with proteome size | rho **-0.228** (n = 1,491), opposite sign | nuclear share falls with proteome size | **killed** |
| H6 Microsporidia less disordered at matched length | median genome contrast -0.192, bootstrap interval -0.201 to -0.178 (30 genomes) | about -0.19 IDR_all0 | supported |
| H7 effector-like candidates more frequent in plant-pathogen genera | positive in 3 of 4 evaluable classes (p 0.31, q 0.36) | per 1,000 proteins: Dothideomycetes +1.77, Agaricomycetes +0.41, Sordariomycetes +0.33, Leotiomycetes -1.76 | meets the rule (>= 75% of at least 3 classes), **but weak** |

### How to read these

- **H1 to H3** hold in nearly every class with effect sizes of 0.12 to 0.16 in IDR fraction, after matching for protein length and
  within each genome. Because H4 is only mixed (rho 0.19), a share of the IDR signal follows composition; a GC-adjusted version of
  H1 to H3 was not pre-registered and was not run.
- **H5 is opposite to the prediction**, not merely null; why was not investigated.
- **H6** is a large effect, but Microsporidia genomes are small and annotation-sensitive, and the 30 genera are not independent;
  it should be read as "Microsporidia proteomes are strongly depleted in predicted disorder at matched length".
- **H7 meets the pre-registered rule only narrowly.** There are 4 evaluable classes, the sign test is not significant (p 0.31), the
  per-class differences are small against rates of 11 to 15 candidates per 1,000 proteins (one class goes the other way), and C4
  failed. Treat it as no evidence, not as support.
- The class is the unit of the sign tests; classes are not independent (shared ancestry), and no phylogenetic correction was possible.
- The 245,547 effector-like candidates (about 13 per 1,000 proteins) are small secreted proteins of unknown function selected by the
  rule above, not validated effectors.

### Not done

A two-sided test for H5; ~~the cause of the C4 failure~~ and ~~GC-adjusted H1 to H3~~ (both done post hoc on 2026-11-21, see the next section); any phylogenetic correction; any link to function
(literature claims about disorder and function are citations to verify and none is made).

---
## Post hoc checks (2026-11-21; NOT pre-registered; they do not change the verdicts above)

Scripts `03_c4_diagnosis.py` (job 32358756, 1.5 min) and `04_gc_adjusted.py` (job 32358729, 43 s); tables in `results/posthoc/`;
tests `test_gc_adjusted.py`. Both were chosen after the results were known.

### Why C4 failed (29,062 proteins with CFEM, LysM, CAP or Chitin_bind_1 domains at e <= 1e-5)

Each protein gets one reason, in the order no SignalP row, SignalP < 0.9, more than 1 TMHMM helix, GPI anchor, pass.

| Family | proteins | pass strict | no SignalP row | SignalP < 0.9 | > 1 TM helix | GPI anchor |
|---|---|---|---|---|---|---|
| CFEM | 12,294 | 19.2% | 10.7% | 0.9% | 18.0% | **51.1%** |
| CAP | 5,970 | 66.3% | 22.7% | 3.3% | 1.1% | 6.7% |
| Chitin_bind_1 | 5,517 | 67.3% | 25.5% | 3.2% | 1.3% | 2.8% |
| LysM | 5,281 | 47.9% | **47.2%** | 2.0% | 0.6% | 2.3% |
| all | 29,062 | 43.2% | 22.6% | 2.0% | 8.2% | 23.9% |

- CFEM is 42% of the control proteins, and 69% of them fail because of a predicted GPI anchor or more than one TM helix. The strict filter
  removes GPI-anchored proteins by design, so for CFEM the control and the filter disagree by construction.
- Without CFEM the pooled pass rate is 60.8%, still under 70%. The main loss there is "no SignalP row" (22.7% to 47.2%).
- If GPI anchors and TM helices are allowed (SignalP >= 0.9 only), the pooled rate is 75.4%; with the GPI and TM-count rules of
  `secreted_noTM` it is 67.2%.
- DeepLoc calls 92.9% of these proteins extracellular or cell membrane (CAP 93.0%, CFEM 99.5%, Chitin_bind_1 99.2%, LysM 70.6%),
  so the localization predictor and the signal-peptide filter disagree for a large share.
- Not tested: whether the proteins without a SignalP row are real secreted proteins missed by SignalP, or intracellular proteins that
  carry the domain (LysM and CAP occur in both). So the filter is shown to be too strict for GPI-anchored families, but it is not shown to miss
  real effectors.
- Consequence for H7 and the candidate list: the rule excludes GPI-anchored and multi-TM proteins and every protein without a
  SignalP call. The candidate list is therefore a conservative subset. H7 was not re-run with a looser rule.

### GC-adjusted H1 to H3 (classes with at least 8 genomes)

The contrasts are already within genome, so a genome-level GC difference cannot create them directly. The check asks whether the size of
the contrast depends on genome GC. Slopes use class-centred GC.

| | as pre-registered: classes in direction, pooled median | GC-adjusted: classes in direction, pooled median | slope per GC point | Spearman rho (contrast, GC) | low-GC half / high-GC half median |
|---|---|---|---|---|---|
| H1 | 19 of 20, +0.118 | 19 of 20, +0.118 | -0.0003 | -0.03 | +0.120 / +0.116 |
| H2 | 20 of 20, +0.163 | 20 of 20, +0.163 | +0.0017 | +0.17 | +0.159 / +0.166 |
| H3 | 18 of 20, -0.140 | 18 of 20, -0.140 | -0.0012 | -0.09 | -0.139 / -0.142 |

- The verdicts do not change; the adjusted and unadjusted medians agree to three decimals. The pooled medians differ slightly from the
  table above (H1 +0.124) because this check keeps only genomes in classes with at least 8 valid contrasts.
- The H2 contrast rises a little with GC (rho 0.17); the others do not depend on it.
- This adjusts for genome GC only. It does not adjust for the amino-acid composition of each protein, because the cube holds no
  per-protein composition. H4 (rho 0.19) therefore still leaves open that disorder predictions partly follow composition.
