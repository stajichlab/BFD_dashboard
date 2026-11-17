# T4 and T6: entomopathogen trehalose co-expansion, and neutral-trehalase retention in reduced genomes

Status: **design and decision rules, written before any result** (committed first; results are added after).
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

## Files (planned)

`01_build_tables.py` (HPCC job), `02_analyze.py`, `test_t4t6.py`, `results/`.
