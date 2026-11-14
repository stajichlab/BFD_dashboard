# kx2_h4: Kickxellales vs yeast reduction (hypothesis H4)

SQL that produced the H4 analysis (loss composition, scaling residuals, loss function and overlap, gains).
All steps are MotherDuck SQL run with the query_rw tool in database `bfd_scratch`; reads of `bfd` are read-only.
Nothing was computed in Python or by hand: the OLS/SMA fits (`regr_slope`, `stddev(y)/stddev(x)*sign(corr)`) and the
stratified overlap null are plain SQL. Reporting/diagnostic queries are kept as commented blocks inside each file.

## Run order and table dependencies

| File | Creates | Depends on |
|---|---|---|
| 01_meta_nprot_gen.sql | kx2_h4_nprot, kx2_h4_meta, kx2_h4_gen | bfd.gene_proteins, kw_rep_genomes |
| 02_gp_all_gm.sql | kx2_h4_gp_all, kx2_h4_gm | 01, kw_pfam_matrix |
| 03_loss_split.sql | kx2_h4_loss | 02 |
| 04_class_counts.sql | kx2_h4_cazy_fam, kx2_h4_tf, kx2_h4_sp, kx2_h4_cls | 01, kw_pfam_matrix, bfd.cazy_overview/tf_inventory/signalp_signal_peptide |
| 05_scaling.sql | kx2_h4_gcls, _fit, _resid, _resid_grp, _scaling_summary | 01, 04 |
| 06_pfcat.sql | kx2_h4_pfcat | bfd.swissprot_subject |
| 07_lostset_overlap.sql | kx2_h4_lostset, kx2_h4_overlap | 02, 06 (enrichment query only) |
| 08_gains_trypsin.sql | kx2_h4_gains, kx2_h4_rep_gen, kx2_h4_trypsin_arch | 01, 02, 04, bfd.pfam |

03, 05-07 can run in any order after their dependencies. Each file's statements are run one per call (or the whole file
in one call; the multi-statement form worked).

## Key thresholds (also in each file header)
- Core baseline family: prevalence >= 0.5 and median >= 1 (Mucorales or Pezizomycotina); lostset uses families core in both.
- Lost: group median 0 and prevalence < 0.3. BUSCO>=85 variant = genus median complete_pct >= 85 (flt 'b85').
- Gain: median above the 95th percentile of both baselines and >= 2 copies above the higher baseline median.
- Scaling label: |z|>=1 in both OLS and SMA residuals (over-contracted/expanded), else mixed or on-curve.

## Verification (2026-11-11)
The loss split, lostset/overlap and gains steps were re-run into `kx2_h4_chk_*` tables from the saved SQL and compared
to the originals and to the reported values; all four chk tables were then dropped (none remain).

| Check | Reported | Re-run | Result |
|---|---|---|---|
| Kickxellales copies lost vs Pezizomycotina | 41% | 41.4% | reproduced |
| Contraction share of loss vs Pezizomycotina | 59% | 58.6% | reproduced |
| Kickxellales vs Mucorales (lost / contraction) | 46% / 81% | 45.7% / 81.2% | reproduced |
| BUSCO>=85 Kickxellales vs Mucorales | 42% / 80% | 42.0% / 80.0% | reproduced |
| Overlap obs/exp, Saccharomycetes | 0.93 | 0.93 (189 obs) | reproduced |
| Overlap, Pneumocystis / Harpellales | 1.03 / 1.87 | 1.03 / 1.87 | reproduced |
| Trypsin copies Kickxellales (all / b85), Harpellales | 21.5 / 22, 72.5 | 21.5 / 22, 72.5 | reproduced |
| CAP, BAHD, GH46-Pfam Kickxellales | 10, 17, 4 | 10, 17, 4 | reproduced |
| Row-level diff vs original loss / overlap / gains tables | n/a | 0 rows differ both directions | identical |

## Not reconstructed exactly / caveats
- The checks used the existing `kx2_h4_gm` and `kx2_h4_gp_all` as inputs (steps 01, 02 and 04-06 were not re-run end to end
  because each is a large aggregation); their SQL is copied from what was run, but a full re-run was not performed.
- 01-08 are the final working versions. Failed first attempts (an ambiguous-column error in an enrichment query, a
  misplaced `label` alias, `list_quantile` instead of `list_aggregate(...,'quantile_cont',0.95)`, `x` not in GROUP BY in the
  first residual query, and a non-read-only statement sent to the read-only tool) were corrected; only the corrected form is saved.
- `kx2_h4_fit` keeps a leftover column `sd0` (stddev of y, unused).
- `kx2_h4_pfcat` ties in the category argmax are broken arbitrarily by row_number; 6,310 of 20,084 Pfams in SwissProt subjects (31%) are 'unclassified' (the earlier report figure of 41% was wrong); Pfams absent from SwissProt are also treated as unclassified downstream.
- Group sizes were not saved as a table; they are reproduced by the commented query in 01.
- The `kx2_h4_scaling_summary` leverage query and the CAZy/targeted-candidate listings are stored as comments, not tables.
