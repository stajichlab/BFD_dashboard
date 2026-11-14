# Ku-interacting machinery in Kickxellomycota: presence/absence and co-loss (Q1)

All computation ran in MotherDuck (`bfd_scratch`, prefix `kx2_ku_`); `bfd` was read only. Tree: first-pass FastTree backbone
(63 genus representatives; `bfd_scratch.kx2_tree_taxa`, tree='backbone'). Losses are placed by **Dollo parsimony on that
first-pass topology, not a likelihood reconstruction**. Compensation (alternative end joining, HR copy number, telomere
sequence) is Q2 and is not covered here.

## Run order and dependencies

| step | file | needs | makes |
|---|---|---|---|
| 01 | `01_genomes.sql` | kx2_tree_taxa, bfd.species/busco/asm_stats/gene_proteins | kx2_ku_genomes (169 = 63 backbone + 106 further Coemansia strains; the backbone Coemansia is the 107th) |
| 02 | `02_pfam_hits.sql` | 01 | kx2_ku_accs, kx2_ku_pfam (panel Pfam hits, i-evalue <= 1e-3) |
| 03 | `03_swissprot_best.sql` | 01 | kx2_ku_sp (SwissProt best hit per protein) |
| 04 | `04_rules.sql` | - | kx2_ku_rules (46 rules) |
| 04b | `04b_rules_controls_patch.sql` | 04 | updates control rules (see caveat below) |
| 05 | `05_calls.sql` | 02, 03, 04b | kx2_ku_calls (protein-level) |
| 06 | `06_matrix.sql` | 01, 04b, 05 | kx2_ku_matrix (gene x genome presence) |
| 07 | `07_pfam_coverage.sql` | 01 | kx2_ku_pfam_all, kx2_ku_panel (2,017 families), kx2_ku_cov |
| 08 | `08_prevalence_controls.sql` | 06, 07 | kx2_ku_prev, kx2_ku_ctrl_cov, kx2_ku_ctrl_grp |
| 09 | `09_keydomain_sensitivity.sql` | 01, 02 | kx2_ku_key, kx2_ku_keymatrix, kx2_ku_key_prev |
| 10 | `10_tree_clades.py` then run its output `10_tree_clades.sql` | the Newick | kx2_ku_tree_edges (manual/Python step; the run pasted the same VALUES on one line) |
| 11 | `11_dollo.sql` | 06, 10 | kx2_ku_node_label, kx2_ku_clade_tips, kx2_ku_dollo, kx2_ku_dollo_sum |
| 11b | `11b_dollo_keydomain.sql` | 09, 11 | kx2_ku_dollo_key |
| 12 | `12_pikk_relaxed.sql` | 02, 03 | kx2_ku_pikk, kx2_ku_pikk_grp |
| 13 | `13_coloss.sql` | 06, 09, 11 | kx2_ku_coloss |
| 14 | `14_tblastn_worklist.sql` | 11, 13 | kx2_ku_queries, kx2_ku_worklist |
| 15 | `15_verify.sql` | 05, 06, 13 | comparison queries only (chk tables dropped) |

Graph: 01 -> (02, 03, 07); 04 -> 04b; (02, 03, 04b) -> 05 -> 06 -> (08 with 07, 13, 11); 02 -> 09 -> (11b, 13); 10 -> 11 -> (11b, 13, 14); 13 -> 14.
Manual steps: step 10 (Python tree text parsing, `python3 -I`); the group labels and node labels in 11 are hand-coded from the tree.

## Gene list and identification rules (table `kx2_ku_rules`)

Pfam accessions below were read from the `pfam` field of the named SwissProt reference in `bfd.swissprot_subject` and found in
`bfd.pfam`; none are from memory. Name rules use the SwissProt best hit (pident >= 25, subject coverage >= 40%, -log10 e >= 10).
Pfam rules use domain i-evalue <= 1e-5 (<= 1e-3 where marked "1e-3"). Reference proteins: Ku (S. cerevisiae YKU70 P32807, YKU80 Q04437;
A. nidulans ku70 Q5AVC7, ku80 Q5B4H8; N. crassa mus-51 Q7SA95, mus-52 Q7RX73).

| gene | rule (domains / name) | conf |
|---|---|---|
| Ku (70/80 not separated) | [PF02735](https://www.ebi.ac.uk/interpro/entry/pfam/PF02735/) Ku core and [PF03731](https://www.ebi.ac.uk/interpro/entry/pfam/PF03731/) Ku_N in one protein | high |
| DNA-PKcs | PIKK [PF00454](https://www.ebi.ac.uk/interpro/entry/pfam/PF00454/) plus any of [PF08163](https://www.ebi.ac.uk/interpro/entry/pfam/PF08163/), [PF20500](https://www.ebi.ac.uk/interpro/entry/pfam/PF20500/), [PF20502](https://www.ebi.ac.uk/interpro/entry/pfam/PF20502/), [PF19704](https://www.ebi.ac.uk/interpro/entry/pfam/PF19704/) (human PRKDC P78527) | high |
| Lig4 | [PF04675](https://www.ebi.ac.uk/interpro/entry/pfam/PF04675/) + [PF01068](https://www.ebi.ac.uk/interpro/entry/pfam/PF01068/) and best hit named DNA ligase 4/lig4/dnl4/mus-53 (separates from Lig1) | high |
| XRCC4 / Lif1 | [PF06632](https://www.ebi.ac.uk/interpro/entry/pfam/PF06632/) (1e-3) or best-hit name; "like"/XLF excluded | low |
| XLF / Nej1 | [PF09302](https://www.ebi.ac.uk/interpro/entry/pfam/PF09302/) (1e-3) or best-hit name (xlf, Nej1, Cernunnos) | low |
| PAXX | [PF15384](https://www.ebi.ac.uk/interpro/entry/pfam/PF15384/) (1e-3) | low |
| SNM1/Pso2/Artemis family | [PF07522](https://www.ebi.ac.uk/interpro/entry/pfam/PF07522/) + [PF00753](https://www.ebi.ac.uk/interpro/entry/pfam/PF00753/); Artemis proper by name only | medium / low |
| Pol4 (X-family) | [PF14792](https://www.ebi.ac.uk/interpro/entry/pfam/PF14792/) + [PF14791](https://www.ebi.ac.uk/interpro/entry/pfam/PF14791/) | medium |
| PNKP/Tpp1 | [PF08645](https://www.ebi.ac.uk/interpro/entry/pfam/PF08645/) | medium |
| TDP1 / TDP2 | [PF06087](https://www.ebi.ac.uk/interpro/entry/pfam/PF06087/); TDP2 = [PF03372](https://www.ebi.ac.uk/interpro/entry/pfam/PF03372/) + name (generic domain) | medium / low |
| Aprataxin (Hnt3) | [PF01230](https://www.ebi.ac.uk/interpro/entry/pfam/PF01230/) HIT + [PF16278](https://www.ebi.ac.uk/interpro/entry/pfam/PF16278/) | medium |
| Mre11 | [PF00149](https://www.ebi.ac.uk/interpro/entry/pfam/PF00149/) Metallophos + [PF04152](https://www.ebi.ac.uk/interpro/entry/pfam/PF04152/) DNA-binding | high |
| Rad50 | [PF04423](https://www.ebi.ac.uk/interpro/entry/pfam/PF04423/) zinc hook + [PF13476](https://www.ebi.ac.uk/interpro/entry/pfam/PF13476/), or best hit Rad50 (1e-3) | high |
| Nbs1/Xrs2 | name only (FHA/BRCT generic) | low |
| Sae2/CtIP | [PF08573](https://www.ebi.ac.uk/interpro/entry/pfam/PF08573/) (1e-3) | medium |
| Exo1 | [PF00867](https://www.ebi.ac.uk/interpro/entry/pfam/PF00867/) + [PF00752](https://www.ebi.ac.uk/interpro/entry/pfam/PF00752/) and best hit Exo1 (FEN1/XPG/GEN1 excluded) | medium |
| Dna2 | [PF08696](https://www.ebi.ac.uk/interpro/entry/pfam/PF08696/) + [PF13086](https://www.ebi.ac.uk/interpro/entry/pfam/PF13086/) | high |
| ATM/Tel1 | PIKK + [PF11640](https://www.ebi.ac.uk/interpro/entry/pfam/PF11640/) TAN | high |
| ATR/Mec1 | PIKK + [PF08064](https://www.ebi.ac.uk/interpro/entry/pfam/PF08064/) UME | high |
| ATRIP/Ddc2 | [PF09798](https://www.ebi.ac.uk/interpro/entry/pfam/PF09798/) (named HEAT_LCD1 in BFD) or name (1e-3) | low |
| 53BP1/Rad9/Crb2 | [PF08605](https://www.ebi.ac.uk/interpro/entry/pfam/PF08605/) or name | low |
| Rad53/Chk2 | [PF00498](https://www.ebi.ac.uk/interpro/entry/pfam/PF00498/) FHA + [PF00069](https://www.ebi.ac.uk/interpro/entry/pfam/PF00069/) kinase and best hit Rad53/Chk2/Cds1 | medium |
| TERT | [PF12009](https://www.ebi.ac.uk/interpro/entry/pfam/PF12009/) | high |
| Est1 (SMG6 as comparator) | [PF10374](https://www.ebi.ac.uk/interpro/entry/pfam/PF10374/) + [PF10373](https://www.ebi.ac.uk/interpro/entry/pfam/PF10373/), without PIN [PF13638](https://www.ebi.ac.uk/interpro/entry/pfam/PF13638/) | medium |
| TPP1/Tpz1/Est3 | [PF10341](https://www.ebi.ac.uk/interpro/entry/pfam/PF10341/) (one family, not separable) (1e-3) | low |
| POT1/Cdc13-OB, Pot1 | [PF02765](https://www.ebi.ac.uk/interpro/entry/pfam/PF02765/) (1e-3); Pot1 = + name | low |
| Rap1 | [PF08914](https://www.ebi.ac.uk/interpro/entry/pfam/PF08914/) (+ PF09197, PF11626, PF21411) or name | low |
| Rif1 | [PF12231](https://www.ebi.ac.uk/interpro/entry/pfam/PF12231/) | medium |
| Taz1/TRF | [PF08558](https://www.ebi.ac.uk/interpro/entry/pfam/PF08558/) (1e-3) | low |
| Stn1, Ten1, CTC1 | [PF10451](https://www.ebi.ac.uk/interpro/entry/pfam/PF10451/)/PF09170; [PF15490](https://www.ebi.ac.uk/interpro/entry/pfam/PF15490/)/PF12658; PF15489/PF15491 (1e-3) | low |
| Cdc13 | name only (Cdc13 Pfam models gave no hits) | low |
| controls (strict) | Rad51 family ([PF08423](https://www.ebi.ac.uk/interpro/entry/pfam/PF08423/) + PF14520); Pol delta (PF24055 + PF00136 + PF03104); Pol epsilon ([PF08490](https://www.ebi.ac.uk/interpro/entry/pfam/PF08490/) + PF22634); Mcm2 ([PF12619](https://www.ebi.ac.uk/interpro/entry/pfam/PF12619/) + PF00493) | high |
| weak controls | Pol alpha (PF12254 + PF00136), Rad54 (PF08658 + PF00176), Rad52 (PF04098) | not used for false-absence |

Pfam hits for PF08599, PF09038, PF12658, PF12659, PF16853, PF18233, PF18691, PF22395, PF22520, PF28409 and PF25385 (and others in `kx2_ku_accs`)
were absent from the 169 genomes at 1e-3 (checked in `kx2_ku_accs` x `kx2_ku_pfam`), so those families contribute nothing.
Honest limits: Ku70 vs Ku80 name split (`Ku70`, `Ku80` rows) is unreliable (best-hit name vs subunit) and is not used; only combined `Ku` is.
Within-gene Lig4/Exo1/Rad53 depend on the best-hit name and therefore on SwissProt coverage.

**Control revision (tuning on the data, stated plainly):** the first Rad51 control used a name rule and gave 67/107 in Coemansia although
the family domain is in 110/111 Kickxellales-core genomes; Pol delta used PF14260 (missed in about 10%). They were switched to family-level
domain rules in `04b`; Pol alpha, Rad54 (Rad54_N absent in 0/4 Harpellales, 32/111 Kickxellales core, a model-sensitivity or real change, not a
clean control) and Rad52 were demoted. 05/06 were re-run after 04b; every number below is from the revised rules.

## Main results

Presence by group (strict rule; genera for backbone groups; Coemansia = 107 strains). Selected genes; full table in `kx2_ku_prev`.

| gene | root (15) | Mucoro. (21) | Zoopag. (8) | Entomo.+Basid. (4) | Dimarg. (3) | Ramic. (1) | Harpell. (4) | Mycoem./Spiro. (2) | Kickx. core (5) | Coemansia (107) |
|---|---|---|---|---|---|---|---|---|---|---|
| Ku | 14 | 21 | 7 | 2 | 2 | 1 | 0 | 0 | 0 | 0 |
| Ku, relaxed key domain | 15 | 21 | 8 | 3 | 3 | 1 | 0 | 0 | 0 | 0 |
| DNA-PKcs | 8 | 17 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Lig4 | 14 | 19 | 4 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| XLF/Nej1 (low) | 10 | 20 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Pol4 | 12 | 20 | 2 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| Aprataxin | 10 | 10 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| Mre11 | 14 | 20 | 7 | 4 | 3 | 1 | 4 | 2 | 3 | 91 |
| Rad50 | 13 | 21 | 7 | 4 | 3 | 1 | 3 | 2 | 4 | 94 |
| Dna2 | 12 | 17 | 6 | 4 | 2 | 1 | 4 | 2 | 5 | 98 |
| Sae2/CtIP | 13 | 19 | 4 | 2 | 3 | 1 | 3 | 0 | 0 | 1 |
| Exo1 | 10 | 7 | 5 | 2 | 3 | 0 | 3 | 1 | 5 | 58 |
| ATM/Tel1 (TAN rule) | 6 | 1 | 3 | 1 | 0 | 0 | 0 | 1 | 3 | 77 |
| ATR/Mec1 (UME rule) | 11 | 8 | 6 | 3 | 3 | 1 | 0 | 1 | 2 | 76 |
| TERT | 13 | 17 | 5 | 3 | 3 | 1 | 3 | 1 | 5 | 102 |
| Est1 | 12 | 18 | 4 | 2 | 3 | 1 | 3 | 0 | 5 | 97 |
| Rif1 | 7 | 5 | 1 | 1 | 3 | 0 | 0 | 0 | 4 | 102 |
| Taz1/TRF (low) | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Stn1 (low) | 7 | 10 | 2 | 1 | 1 | 1 | 0 | 0 | 2 | 63 |
| Ten1 (low) | 11 | 17 | 4 | 2 | 3 | 1 | 3 | 2 | 5 | 92 |

Relaxed PIKK census (12, best-hit names, no coverage filter): ATM-like present in 4/4 Harpellales and 101/107 Coemansia, ATR-like in
4/4 Harpellales and 105/107 Coemansia, so the TAN/UME-rule absences in Harpellales are truncated or diverged proteins (best-hit subject
coverage 11-27%), not clear losses. DNA-PKcs by best-hit name PRKDC occurs in 17/21 Mucoromycota and 6/15 root and in no Kickxellomycota.
ATM-like proteins are rare in Mucoromycota by both rules (1/21), which is why the co-loss class for ATM/Tel1 reads "rare everywhere" (misleading).

Dollo placement (11, strict; 11b relaxed key-domain): **Ku has one robust loss, on the stem of Kickxellales + Harpellales (node 39, 11 genera)**;
Dimargaritales and Ramicandelaber (sister lineages, 4 genera) keep Ku. Strict rule gives 5 losses (2 internal: stem node 39 and the
Neoconidiobolus + Capillidium pair; 3 single-genus: Piptocephalis, Tieghemiomyces, Caulochytrium); relaxed key-domain gives 2 (node 39 plus
Neoconidiobolus alone), so the singletons are probably annotation misses (Tieghemiomyces is the 55%-duplicated caveat taxon).
Lig4, XLF/Nej1, Pol4 and aprataxin are each lost one branch **earlier**, on the stem of Entomophthorales + Kickxellomycota (node 25, 18 genera;
Basidiobolus keeps Lig4), and DNA-PKcs earlier still (node 7: Basidiobolus + Zoopagales + Entomophthorales + Kickxellomycota, 27 genera;
the Mucoromycota keep it). Sae2/CtIP is lost on Kickxellales s.l. (node 47, 7 genera) but kept in Harpellales (3/4) and in the Ku-retaining lineages.
Rif1 is retained in Kickxellales core and Dimargaritales and lost independently several times (Harpellales 0/4, Mycoemilia/Spiromyces, Ramicandelaber, Entomophthorales).
Mre11, Rad50, Dna2, TERT, Est1, Ten1, PNKP, TDP1 are retained everywhere.

Co-loss classes (`kx2_ku_coloss`): lost with/after Ku = Ku only; lost earlier than Ku = DNA-PKcs, Lig4, XLF/Nej1, Pol4, aprataxin; retained = MRN (Mre11, Rad50),
Dna2, TERT, Est1, Ten1, PNKP, TDP1, ATRIP-like; partial/unclear = Sae2, Rif1, Exo1, ATR/Mec1 (Harpellales), Rad53/Chk2, Stn1; uninformative (rare even in outgroups) =
XRCC4/Lif1, PAXX, TDP2, Pot1, Cdc13, CTC1, Rap1, Taz1/TRF, TPP1/Tpz1/Est3, 53BP1.

False absence (08, 09): with the strict multi-domain rules, the four strict controls are missed in about 12.5% (Kickxellomycota backbone), 14.3% (Coemansia strains) and 14.1%
(other backbone) of gene-genome calls; with single key domains the rates are 3.6%, 2.3% and 6.3%. Correlation of control detection with the
Pfam-panel coverage is r = 0.34 and with BUSCO complete r = 0.39 (n = 169). Pfam panel coverage (2,017 families found at i-evalue <= 1e-5 in >= 95% of the 48
non-Kickxellomycota backbone genomes; data-derived, not the curated ~1,400 list; see 07) is 94.1% in Coemansia strains and 96-98% in other Kickxellomycota groups, vs 98-99% elsewhere;
Spiromyces is the low outlier (0.809; BUSCO 44%). A single absent protein therefore carries a false-absence chance of about 2-14% depending on rule stringency; the Ku result
(absent in 107/107 strains and 11/11 backbone genera while present in 4/4 sister-lineage genera by the same rule) is far outside that, but systematic divergence beyond HMM
detection would not show up as a random miss, hence the tblastn worklist.

tblastn worklist (14, not run): 211 gene x genome rows (152 priority 1 across 7 genes and 27 genomes; 59 priority 2), covering Ku, DNA-PKcs, Lig4, XLF/Nej1, Pol4, aprataxin, Sae2, ATR/Mec1, plus
positive controls (Ku in Ramicandelaber, Dimargaris, Rhizopus; Rad51 in the 11 stem genera). Query accessions are fungal SwissProt entries verified in `bfd.swissprot_subject`
(Ku: P32807, Q04437, Q7SA95, Q7RX73, Q5AVC7, Q5B4H8; Lig4: Q5BH83, Q7SB49, Q08387; XLF: Q9P7J2, Q06148; Pol4: Q09693, P25615; aprataxin: Q08702, O74859; Sae2: P46946, O74986; DNA-PKcs has only human P78527 in
SwissProt, so a Mucoromycota BFD protein such as F574EFC1_006650-T1 (Actinomortierella) is a better query).

## What the survey can and cannot show

Observed: Ku (PF02735 + PF03731) is absent from all 107 Coemansia strains and from every Kickxellales and Harpellales backbone genus, and present in Ramicandelaber, Dimargaritales (by key domain, 3/3),
Mucoromycota and almost all outgroups. Lig4, XLF/Nej1, Pol4, aprataxin and DNA-PKcs are already missing in the Ku-retaining Kickxellomycota lineages.
Interpretation (hypotheses, not shown): the stem loss of Ku fell in a lineage that had already lost the rest of the classical NHEJ complex (only a Lig1-type ligase remains), so Ku
was perhaps already functionally orphaned; the loss of Ku in Kickxellales + Harpellales would then be a last step rather than the start of NHEJ decay, and could explain tolerance of its loss. Retention of MRN, Dna2, Mre11-Rad50,
TERT and Est1 argues that the resection and telomerase arms were not lost with Ku. Any link of Ku loss to intron or genome-size reduction is outside this analysis.
Cannot show: that any absence is a true gene loss (annotation under-calls genes; every absence needs tblastn); that losses are single events (first-pass tree, Linderina and Basidiobolus unstable; Entomophthorales
sister to Kickxellomycota in this tree, which changes the node-25 placement if the true position differs); function of Ku-lacking cells; or anything about low-confidence families (XRCC4, XLF, CST, shelterin, Rap1), where
absence mostly reflects fast evolution and no usable domain model. Literature statements on pathway function are not made; any such comparison is "citation to verify".

## Verification record

- `15_verify.sql`: the 05 (calls), 06 (matrix) and 13 (co-loss) statements were re-run into `kx2_ku_chk_calls`, `kx2_ku_chk_matrix`, `kx2_ku_chk_coloss` and compared with `EXCEPT`
  in both directions against the reported tables: 0 differences in each (7,517 call rows both sides; 7,774 matrix rows; 46 co-loss rows). The chk tables were then dropped.
- Limits of that check: it re-uses the same inputs (kx2_ku_pfam, kx2_ku_sp, rules, kx2_ku_keymatrix, clade tips), so it shows determinism and that the saved SQL equals what was run, not an independent re-derivation.
  Not re-run: 01-04b, 07-12 and 14 (the reported tables came from the same statements run once); the edge table in 10 was inserted as a single-line VALUES identical in content to `10_tree_clades.sql`;
  the first run of 05/06 preceded the `04b` control patch and was replaced by the re-run, so the reported tables correspond to the post-patch rules.
