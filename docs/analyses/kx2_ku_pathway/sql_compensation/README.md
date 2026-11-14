# Q2: compensatory DNA-repair / telomere changes with Ku loss (Kickxellales + Harpellales)

SQL for the compensation analysis. Run in MotherDuck with default database `bfd_scratch`; `bfd` is read-only. All intermediate tables are `bfd_scratch.kx2_kc_*`.

## Run order and dependencies

```
01_genomes -> 02_pfhits_arch -> 04_calls -> 05_gene_counts -> 06_wide -> 08_genus_group -> 09_replication
              03_swissprot_class -^                              |-> 07_telomere_genome -> 10_telomere_groups
                                                                 |-> 11_artefact_checks (also reads 04)
                                                                 12_verify (reads 01, 06, 07; creates and drops chk tables)
```

Manual steps: (1) the replication lineage assignment in 09 is manual (Kickxellales + Harpellales counted as one lineage; other Ku-negative genera one lineage each; no tree for the full set);
(2) 06 needs the feature map in the same file (`kx2_kc_featmap`). No Python.

## Quality filters
- `q_loose`: BUSCO complete >= 60, >= 3000 genes, BUSCO duplicated <= 10 (same filter for all clades; early-diverging fungi are 60-85% BUSCO). 21,483 genomes.
- `q_strict`: BUSCO >= 85, dup <= 10 (sensitivity, 11).
- Unit of comparison: genus median (`kx2_kc_genus_long`); groups KxHp (Kickxellales+Harpellales, 7 genera), DimRam (Dimargaritales+Ramicandelaberales, 3), ZoopEnto (12), Mucoro (Mucoromycota+Mortierellomycota+Glomeromycota, 54), EarlyOther (43), Dikarya (1430).
- Tieghemiomyces (55% BUSCO duplicated), Kickxella, Linderina, Spiromyces, Furculomyces fail q_loose and are absent from genus tables.

## Genes, identification rules, confidence
Pfam names verified against `swissprot_subject.pfam` of reference proteins; accession = split_part(pfam_acc,'.',1). Domain hits i-evalue <= 1e-5. Best-hit name classes need best SwissProt hit e-value <= 1e-10 AND the Pfam domain.

| Gene / family | Rule | Pfam verified | Confidence |
|---|---|---|---|
| Ku70/80 | Ku domain ([PF02735](https://www.ebi.ac.uk/interpro/entry/pfam/PF02735/)) or Ku_N [PF03731](https://www.ebi.ac.uk/interpro/entry/pfam/PF03731/) + Ku_C [PF03730](https://www.ebi.ac.uk/interpro/entry/pfam/PF03730/) | yes | high presence; absence is weak evidence (annotation) |
| POLQ | DNA_pol_A [PF00476](https://www.ebi.ac.uk/interpro/entry/pfam/PF00476/) + DEAD [PF00270](https://www.ebi.ac.uk/interpro/entry/pfam/PF00270/) or Helicase_C [PF00271](https://www.ebi.ac.uk/interpro/entry/pfam/PF00271/) or POLQ_helical [PF21099](https://www.ebi.ac.uk/interpro/entry/pfam/PF21099/) in one protein; best-hit name 'DNA polymerase theta' NOT used (hits are Mus308/HELQ helicases: 5 of 6763 carry DNA_pol_A) | yes (human POLQ O75417) | high for presence (27 proteins kingdom-wide); absence moderate |
| PARP1-like | PARP [PF00644](https://www.ebi.ac.uk/interpro/entry/pfam/PF00644/) + WGR [PF05406](https://www.ebi.ac.uk/interpro/entry/pfam/PF05406/) / zf-PARP [PF00645](https://www.ebi.ac.uk/interpro/entry/pfam/PF00645/) / PARP_reg [PF02877](https://www.ebi.ac.uk/interpro/entry/pfam/PF02877/) | yes (P09874) | moderate (family-level) |
| XRCC1 | XRCC1_N [PF01834](https://www.ebi.ac.uk/interpro/entry/pfam/PF01834/) | yes | 1 protein kingdom-wide: not informative |
| LIG1 / LIG3 / LIG4 | DNA_ligase_A_M [PF01068](https://www.ebi.ac.uk/interpro/entry/pfam/PF01068/) + best-hit name 'DNA ligase 1/3/4'; family count `ligfam` is name-independent | yes | LIG1 high; LIG3/LIG4 moderate (name subtyping) |
| FEN1 / EXO1 / GEN1 | XPG_I [PF00867](https://www.ebi.ac.uk/interpro/entry/pfam/PF00867/) + name; `xpgfam` = all XPG_I proteins | yes | FEN1 high; EXO1/GEN1 moderate |
| Rad51, Dmc1, paralogs | Rad51 [PF08423](https://www.ebi.ac.uk/interpro/entry/pfam/PF08423/); subtype from best-hit name (Rad51, Dmc1, Rad55 (rhp55), Rad57 (rhp57), XRCC2/3, Rad51B/C/D); `rad51fam` = all Rad51-domain proteins | yes | rad51fam high; subtypes moderate; Rad55 low (divergent, ATPase [PF06745](https://www.ebi.ac.uk/interpro/entry/pfam/PF06745/) route) |
| Rad52 / Rad59 | Rad52_Rad22 [PF04098](https://www.ebi.ac.uk/interpro/entry/pfam/PF04098/); `rad52fam` count is primary; Rad52 vs Rad59 name split unreliable | yes | family high; split low |
| Rad54 / Rdh54 | Rad54_N [PF08658](https://www.ebi.ac.uk/interpro/entry/pfam/PF08658/) (`rad54N`) and, for presence, SNF2-rel_dom [PF00176](https://www.ebi.ac.uk/interpro/entry/pfam/PF00176/) + best hit 'rhp54'/Rad54/Rdh54 (`rad54name`); Rad54 vs Rdh54 not separated | yes | presence high via rad54name; N-domain detection moderate |
| BRCA2 | BRCA2 [PF00634](https://www.ebi.ac.uk/interpro/entry/pfam/PF00634/) or BRCA-2_OB1 [PF09103](https://www.ebi.ac.uk/interpro/entry/pfam/PF09103/) | yes | low (fungal BRCA2 divergent) |
| Mus81, Mms4/Eme1, XPF/Rad1, Rad10/ERCC1 | ERCC4 [PF02732](https://www.ebi.ac.uk/interpro/entry/pfam/PF02732/) + name; Rad10 [PF03834](https://www.ebi.ac.uk/interpro/entry/pfam/PF03834/); Mms4 name class low sensitivity | yes | Mus81, XPF moderate; Mms4 low |
| TERT | Telomerase_RBD [PF12009](https://www.ebi.ac.uk/interpro/entry/pfam/PF12009/) or RVT_1 [PF00078](https://www.ebi.ac.uk/interpro/entry/pfam/PF00078/) + name 'Telomerase reverse transcriptase' | yes | high |
| Telomerase RNA | not available: `infernal_rfam` has 14 genomes (no Kickxellomycota) and only the Ascomycota model (Telomerase_Asco, RF02462) | n/a | not assessable |

## Verification record (2026-11-11)
- 12_verify.sql: independent recount of ku, rad51fam, rad54N, ligfam, parp_any, tertdom, rad52fam, polq from bfd.pfam for 313 genomes: 0 mismatches against `kx2_kc_wide`.
  Telomere detection recomputed from bfd.telomere_summary for N50 >= 500 kb: Harpellales 3/3, Kickxellales 0/2, Mucoromycota+Mortierellomycota+Glomeromycota 12/63; identical to `kx2_kc_tel_genome`. chk tables dropped.
- Pipeline was re-executed once end to end for 03 -> 06 -> 08 after the rhp54/rhp55/rhp57/rhm52 name fix (the first Rad54 name regex missed the S. pombe style names); the numbers in the report are from the re-executed tables.
- NOT independently re-computed: name-based classes (03/04), genus/Cliff tables (08), replication (09), tables 10/11. 07 was built before the last 06 re-run; the columns it uses (q_loose, Ku) are unchanged by that re-run.
- Not verifiable here: any gene absence by tblastn (no sequences accessible); all absences rest on BFD annotation.
- asm_stats.telomere_fwd/rev were not used (definition undocumented, looser than telomere_summary).
