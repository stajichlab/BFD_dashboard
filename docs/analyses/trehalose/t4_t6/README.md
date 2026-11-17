# T4 and T6: entomopathogen trehalose co-expansion, and neutral-trehalase retention in reduced genomes

Status: **done 2026-11-16: T4 mixed (TPS supported, ATH in 67% of the focal genera); T6 mixed with no support for preferential loss** (see Results at the end).
The design and decision rules below were written, and committed, before any result.
Parent hypotheses: [`../README.md`](../README.md) section 7 (T4, T6). Data: the local `BFD.duckdb`; unit = one representative genome per genus
(`is_genus_rep`, annotated). Domain gate for Pfam: `domain_i_evalue <= 1e-5`, version stripped; copies = proteins carrying the domain.

## T4. Entomopathogens co-expand trehalose synthesis and hydrolysis

Hypothesis (from the parent README): entomopathogens take up host trehalose (insect haemolymph is rich in trehalose, **citation to verify**), so
TPS ([PF00982](https://www.ebi.ac.uk/interpro/entry/pfam/PF00982/)), acid trehalase ATH/GH65 ([PF03632](https://www.ebi.ac.uk/interpro/entry/pfam/PF03632/)) and trehalose
transport should co-expand in them.

- **Focal set F:** Sordariomycetes genus representatives whose FungalTraits 1.2 `primary_lifestyle` is `animal_parasite` (verbatim label; there is no
  entomopathogen label, so this set also contains nematode-trapping and other animal parasites, a stated limit).
  **Reference R:** all other Sordariomycetes genus representatives that have a FungalTraits label (unlabelled genera are excluded, not guessed).
- **Measure:** copies per 1,000 proteins per genome. Effect size first: median difference with bootstrap 95% interval (genera resampled), Cliff's delta, Mann-Whitney p.
- **T4a (TPS):** supported if the median difference is positive, its interval excludes 0 and delta >= 0.33; killed if delta <= 0.
- **T4b (ATH), the kill rule of the parent README:** killed if fewer than 50% of F genera carry at least one ATH copy; supported if at least
  75% do **and** ATH per 1,000 is higher in F than R with delta >= 0.33; otherwise mixed.
- **T4c (trehalose transport):** **not tested**. No Pfam family or curated set specific to trehalose transport exists in the tables (MFS_1 is
  far too generic); registered here so its absence is explicit.
- **Positive control:** subtilisin-like proteases ([PF00082](https://www.ebi.ac.uk/interpro/entry/pfam/PF00082/)) and GH18 chitinases ([PF00704](https://www.ebi.ac.uk/interpro/entry/pfam/PF00704/)) per 1,000 proteins
  higher in F than R (delta > 0 and interval excluding 0 for both). If the control fails, F is not capturing the expected animal-parasite
  signal and T4 is read with that caveat (that these families expand in entomopathogens is **citation to verify**).
- **T4 overall:** supported only if T4a and T4b are supported; killed if T4b is killed; otherwise mixed.

## T6. Neutral trehalase is preferentially dropped in reduced genomes

Hypothesis: neutral trehalase (NTH, [PF01204](https://www.ebi.ac.uk/interpro/entry/pfam/PF01204/)) loses copies in reduced genomes more than most families do.

- **Families:** all Pfam families present in at least 90% of genus representatives (the near-universal set), plus the trehalose and glycogen modules of the parent README.
- **Contrasts** (a reduced group against a reference group; each side needs at least 5 genera, otherwise the contrast is "not evaluable"):
  R1 order Kickxellales vs Mucorales; R2 Microsporidia vs Chytridiomycota; R3 Saccharomycotina vs Pezizomycotina (subphylum derived from class,
  `bin/summary/taxonomy.py`); R4 Malasseziomycetes vs Ustilaginomycetes + Exobasidiomycetes.
  **Controls without reduction:** K1 Eurotiomycetes vs Sordariomycetes; K2 Dothideomycetes vs Leotiomycetes. A control is valid only if the
  ratio of median proteome sizes is at least 0.7; a reduction contrast is expected to be below 0.8 (reported).
- **Retention ratio** of family f: (mean copies per genome in the reduced group + 0.01) / (mean copies in the reference group + 0.01).
  **Percentile of NTH** = share of families whose retention ratio is below NTH's; low means NTH is dropped more than most families. TPS, GSY, GPH,
  and the families whose Pfam name starts with `Ribosomal` (a housekeeping set) are reported the same way, as description.
- **Supported** if NTH is at or below the 10th percentile in at least 3 of the 4 reduction contrasts **and** above the 25th percentile in both
  controls. **Killed** if NTH is above the 25th percentile in at least 3 of the 4 reduction contrasts (not preferentially dropped) **or** at or below the 10th
  percentile in both controls (the same pattern without reduction, the kill rule of the parent README). Otherwise mixed.
- The earlier allometry result (NTH at the 85th percentile of the null; not significant) is a different statistic and is not re-used here.

## Limits stated in advance

Genus-level unit and shared ancestry (contrasts are between clades, so each is one comparison); Pfam presence is annotation-dependent
(reduced and fragmented genomes lose calls to gene-model gaps); `animal_parasite` is a coarse genus-level label; retention ratios are
over means that are sensitive to a few high-copy genomes; no function or mechanism is tested.

## Files

`01_build_tables.py` (HPCC job 32135806), `02_analyze.py`, `test_t4t6.py` (5 fixture tests pass on the HPCC), `results/`. The family copy
table (4,414,468 rows, 1,654 genus representatives, 2,780 near-universal families plus the module families) stays in
`bench/trehalose_T4T6/tables/` on the HPCC.

---
## Results (2026-11-16)

Order of events: design and rules committed first (`b22fce1`), scripts next (`e4e4c20`), then one build and one analysis run. No rule or
threshold was changed after the results were seen (the only edit after the first test run was a fixture-test bound).

### T4: mixed by the rule (TPS supported, ATH not in 75% of the focal genera)

Focal set: 27 Sordariomycetes genera labelled `animal_parasite` by FungalTraits; reference: 209 other labelled Sordariomycetes genera.
Per 1,000 proteins (`results/t4_contrasts.tsv`):

| Family | Median focal | Median reference | Difference (95% interval) | Cliff's delta |
|---|---|---|---|---|
| TPS | 0.375 | 0.269 | +0.105 (0.079 to 0.170) | 0.63 |
| ATH (GH65) | 0.091 | 0.000 | +0.091 (0 to 0.098) | 0.44 |
| NTH | 0.289 | 0.239 | +0.050 (0.025 to 0.074) | 0.38 |
| TPP | 0.206 | 0.179 | +0.027 (0.006 to 0.058) | 0.35 |
| GSY | 0.096 | 0.090 | +0.006 (0.001 to 0.015) | 0.22 |
| GPH | 0.096 | 0.091 | +0.006 (0.000 to 0.018) | 0.24 |
| Subtilisin S8 proteases (control) | 2.42 | 1.90 | +0.52 (0.10 to 1.21) | 0.43 |
| GH18 chitinases (control) | 1.84 | 1.40 | +0.44 (0.26 to 0.71) | 0.51 |

- **T4a (TPS): supported.** Delta 0.63 and the interval excludes 0.
- **T4b (ATH): mixed.** 18 of the 27 focal genera (67%) carry at least one ATH copy: more than the 50% below which the parent README's rule
  kills the hypothesis, fewer than the 75% needed for "supported". So the parent README's observation that ATH is absent in
  *Beauveria* is true of that genus but not of the focal set as a whole; ATH is patchy inside it.
- **T4c (trehalose transport): not tested** (stated in advance).
- **Positive control passes**: subtilisin proteases and GH18 chitinases are both higher in the focal set, so the label does pick up an
  animal-parasite signal (that these families expand in entomopathogens is a citation to verify).
- **Pattern:** the focal set is higher in the trehalose synthesis and hydrolysis families (TPS, ATH, NTH, TPP) and only marginally in the glycogen
  arm (GSY, GPH, delta about 0.2), so the signal is specific to the trehalose side. This is an association with a coarse label, not evidence for host
  trehalose uptake.

### T6: mixed by the rule, with no support for preferential loss

| Contrast | Genera (reduced / reference) | Proteome ratio | NTH retention ratio | NTH percentile among 2,780 families |
|---|---|---|---|---|
| R1 Kickxellales vs Mucorales | 7 / 36 | 0.64 | 0.64 | 0.37 |
| R2 Microsporidia vs Chytridiomycota | 30 / 36 | 0.26 | 0.78 | **0.90** |
| R3 Saccharomycotina vs Pezizomycotina | 145 / 870 | 0.52 | 0.77 | 0.20 |
| R4 Malasseziomycetes vs two Ustilaginomycotina classes | 1 / 31 | not evaluable (Malassezia is one genus) | | |
| K1 control: Eurotiomycetes vs Sordariomycetes | 71 / 279 | 0.86 | 0.48 | **0.004** |
| K2 control: Dothideomycetes vs Leotiomycetes | 207 / 89 | 0.92 | 1.89 | **1.00** |

- No evaluable reduction contrast has NTH at or below the 10th percentile (0 of 3); in R1 NTH falls exactly in proportion to the
  proteome (0.64 against 0.64), and in R2 and R3 it falls **less** than the proteome (0.78 against 0.26; 0.77 against 0.52). That is not
  preferential loss.
- The two controls without reduction differ in NTH by about two-fold in opposite directions (ratio 0.48 and 1.89, percentiles 0.004
  and 1.00): neutral-trehalase copy number varies strongly between ordinary classes, which swamps any reduction effect. The rule was not
  triggered as "killed" (only 2 of 3 reduction contrasts are above the 25th percentile, and the controls are not both at or below the 10th),
  so the pre-registered verdict is **mixed**; the evidence leans against the hypothesis.
- R4 failed by design: the genus-representative unit leaves one genus for Malasseziomycetes, which the 5-genera minimum excluded.
- TPS and the glycogen families are reported alongside in `results/t6_retention.tsv` for description (TPS percentiles 0.17, 0.77, 0.88 in R1 to R3).

### Limits

One comparison per contrast (clades are not independent replicates); `animal_parasite` is a coarse genus-level label that includes nematode
parasites; Pfam presence depends on annotation (a reduced or fragmented genome can lose calls); retention ratios use means, which a few
high-copy genomes move; no mechanism or function is tested.

### Not done

Trehalose transport (T4c); a retention test that matches orthologs rather than domain families; a control with more than two non-reduced
contrasts; a replacement for R4 (a rule change after the result, so not made here).
