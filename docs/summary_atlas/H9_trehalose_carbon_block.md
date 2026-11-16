# H9. Trehalose synthesis is near-universal and moves with the carbon-storage block (pilot)

Status: **pilot.** The numbers below are from the first pass of analysis G (MotherDuck, computed on a laptop, 2026-11-10).
The reproducible rebuild is the HPCC job `workflows/atlas/J13_trehalose.sh` (inputs `docs/analyses/trehalose/scripts/th0_inputs.sql`,
statistics `th1`-`th5`, fixture tests `test_stats_util.py`); until its outputs replace the numbers here (see the old-versus-new table in
[the report](../analyses/trehalose/README.md), section 0) treat them as exploratory. No figures are built on this page; views V57-V61 are planned.

[Back to index](index.md) · full report: [analyses/trehalose/README.md](../analyses/trehalose/README.md) · related: [H4](H4_secretome_lifestyle.md), [H5](H5_domains_tf.md), [H1](H1_quality.md)

## Hypotheses

| ID | Hypothesis | Supports | Falsifies |
|---|---|---|---|
| T-a | The trehalose pathway is near-universal, so "has trehalose" is not informative; copy number and branch balance are | TPS, TPP and neutral trehalase present in > 95% of quality genomes | clades with a clean, replicated absence |
| T-b | Trehalose synthesis varies together with glycogen, UDP-glucose and glycerol supply (a carbon-partitioning block), not with desiccation-protectant or spore modules | class-demeaned correlations of 0.65-0.86 inside the block, about 0 outside | no block structure |
| T-c | Trehalose gene density follows proteome size like an ordinary family | percentile in the kingdom-wide allometry null near 50 (control families agree) | percentile in the tails |
| T-d | Apparent TPS/TPP absences (Lecanoromycetes, a *Scedosporium* clade) are annotation gaps, not losses | absences cluster in poorly scanned, fragmented genomes; a class-stratified co-absence null removes the signal | the signal survives the class-stratified null and tblastn finds no gene |

## Views (planned, not built)

| View | Content |
|---|---|
| V57 | module prevalence and median copies per 1,000 proteins by class |
| V58 | class-demeaned correlation matrix of module densities |
| V59 | allometry null: percentile of each module's copies-versus-proteome-size slope among 4,954 families |
| V60 | TPS/TPP absence map with the `qcov` scan-completeness control |
| V61 | trehalose-versus-glycogen balance contrasts for curated stress-habitat genera |

## First-pass findings (observed; MotherDuck run)

1. Among 5,107 quality-kept species representatives, TPS is present in 98.2%, TPP in 98.7% and neutral trehalase in 99.6%; GH65 acid trehalase in 57.3%
   (patchy by class: absent in all Mucoromycetes, Tremellomycetes, Ustilaginomycetes, Chytridiomycetes and Glomeromycetes genera). The GH65 family also contains
   non-trehalase enzymes (citation to verify), so it is a family call.
2. Class-demeaned Spearman correlations between module densities across 1,290 genera: TPS-TPP 0.86, GPD-GSY 0.75, GPH-GSY 0.74; TPS against hydrophobin, DIT1,
   HSP12 and sHSP is within +-0.1. Genus medians are treated as independent, which they are not (no phylogenetic correction).
3. Against the kingdom-wide allometry null, TPS sits at the 53rd percentile (the PGI control agrees); glycerol synthesis is at the 13th and neutral trehalase at the 91st percentile.
4. TPS is absent in 46 of 122 Lecanoromycetes genera and 6 *Scedosporium*-clade genera. This is a **candidate**, not a loss call: it was tested with a null that ignored the clade
   and the absent genomes are gap-prone (see the report, section 0, for the q-value and null corrections made on this branch); the tblastn check (T1) has not been run.
5. Desiccation lifestyle is weakly predicted (within-class AUC 0.68-0.69 for 11 curated xerophile genera, small sets, genus-level labels).
6. *Coemansia* keeps the pathway (TPS 0.44 per 1,000 proteins against 0.67 in Mucorales), so the Kickxellomycota reduction analyses can use it as a control region.

## Method references

Cliff's delta: Cliff (1993), *Psychological Bulletin* 114:494-509. Stratified expected counts: Mantel and Haenszel (1959), *J Natl Cancer Inst* 22:719-748.
Benjamini-Hochberg: Benjamini and Hochberg (1995), *J R Stat Soc B* 57:289-300. AUC as Mann-Whitney: Hanley and McNeil (1982), *Radiology* 143:29-36.
(Journal details written from memory; citation to verify.) Non-independence of species (Barker and Pagel 2005; Maddison and FitzJohn 2015) as cited by analysis B: citation to verify.

## Cross-links

- [Analysis C](../FUNCTIONAL_ANALYSES.md): finding 5 there reads the trehalose pathway as adaptation to UV and desiccation; this page supports "retained, not expanded".
- [Analysis B](../FUNCTIONAL_ANALYSES.md): the TPS+TPP pair is already a two-domain module in its `loss_modules` table, below its reporting cut.
- [Analysis E](../FUNCTIONAL_ANALYSES.md): 12 of the 80 low-scan-coverage genomes are Kickxellomycetes, so E's complete-loss calls need a per-genome scan-coverage (`frac_pfam`) check.
