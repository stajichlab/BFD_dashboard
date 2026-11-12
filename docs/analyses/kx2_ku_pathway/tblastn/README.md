# tblastn check of the Ku-pathway absences

Purpose: BFD gene annotation under-calls genes, so "Ku absent" from Pfam is weak evidence. This searches the assemblies directly with
fungal reference proteins. It covers the 8 genes on the Q1 worklist (Ku, DNA-PKcs, Lig4, Pol4, XLF/Nej1, aprataxin, ATR/Mec1,
Sae2/CtIP) plus a Rad51 control, in the 31 worklist genomes (24 Kickxellomycota and sister genomes with absences to check, including
9 *Coemansia* strains, plus positive controls: Ku in *Rhizopus*, *Dimargaris* and *Ramicandelaber*; Rad51 in the stem genera). The
worklist is `../sql_coloss/14_tblastn_worklist.sql`.

| File | Role |
|---|---|
| `genomes.tsv`, `queries.tsv` | inputs (genomes from the worklist; 23 query proteins, 22 from UniProt plus one *Actinomortierella* DNA-PKcs from BFD) |
| `fetch_queries.sh` | builds `queries.faa` (checks each UniProt length against `bfd.swissprot_subject`) |
| `run_tblastn.sh` | SLURM array, one task per genome, `tblastn` against `input_clean_genomes/<ASMID>.fa.gz` (unmasked; masked file as fallback), BLAST+ 2.17.0. The `input/dna/*.scaffolds.fa` links are broken (they point at deleted `genome_annotation/.../predict_results/` copies), which stopped the first submission |
| `summarize_tblastn.py` | per genome and gene: strong / weak / none |

## Settings and thresholds (fixed before the results were seen)

- Search: e-value 1e-3 (permissive), up to 50 targets per query, BLAST defaults for the rest (including low-complexity masking).
- **strong** = best e-value <= 1e-10 and merged HSPs cover >= 30% of the query on one scaffold. **weak** = e-value <= 1e-5 but not
  strong. **none** = otherwise. Coverage merges all HSPs of a query on its best scaffold, so a gene split by introns counts once.
- The best query per gene is the one with the highest coverage x bit score, so one close reference is enough to count.

## How to read it

- A **strong** hit in a genome Pfam called absent means the gene is present in the assembly and was missed by annotation.
- **None** in a contiguous, complete assembly, with the positive controls recovered, supports a true loss. None in a fragmented
  assembly is still weak, because the region could be missing from the assembly. A diverged gene could also fall below the
  thresholds, which is why the **weak** class and the controls matter.
- The positive controls must be strong (Ku in *Rhizopus*, *Dimargaris*, *Ramicandelaber*; Rad51 in the stem genera). If a control
  is not strong, the threshold is too strict for these genomes and the negative calls cannot be trusted.

## Round 1 result (2026-11-11) and why there is a round 2

Saved as `results/round1_calls.tsv` and `results/round1_by_group.tsv`. Positive controls: Rad51 strong in 30 of 31 genomes (the
search works for a conserved protein); Ku strong in *Rhizopus*, *Ramicandelaber* and *Dispira* but **none in *Dimargaris*** (Pfam-positive),
so 3 of 4 Ku positives were recovered. The small or fast-evolving proteins (Pol4, XLF/Nej1, aprataxin, Sae2) fail the control check:
Pfam finds Pol4 and XLF in nearly all Mucoromycota, but yeast and *Aspergillus* queries find neither in *Rhizopus*, so their
negatives are uninformative. DNA-PKcs is strong only in *Rhizopus* (BFD query); elsewhere the weak hits cover about 15% of the
query, consistent with the shared PIKK kinase domain. Ku is none in all 14 core Kickxellales, 4 Harpellales and 2 *Mycoemilia*/
*Spiromyces* genomes. Lig4 has strong hits in several lineages where Pfam called it absent (*Ramicandelaber*, one Dimargaritales
genus, one Harpellales genus, *Mycoemilia*/*Spiromyces*, *Basidiobolus*) but in none of the 14 core Kickxellales; a strong Lig4 hit may
still be the LIG1-type ligase (shared catalytic domain), so this needs a reciprocal-best-hit check.

**Round 2 is a post-hoc change**, made after seeing that the controls did not all pass: the queries are BFD proteins that Pfam
called present in nearby lineages (`queries_round2.tsv`: Ku from *Dimargaris*, *Dispira*, *Ramicandelaber*, *Rhizopus*,
*Basidiobolus*, *Conidiobolus*, *Syncephalis*; Lig4, Pol4, XLF/Nej1, aprataxin, DNA-PKcs and ATR from *Basidiobolus*, *Rhizopus*,
*Syncephalis* and the Kickxellomycota genera). Same genomes, same thresholds, same calls. A query protein from the genome being
searched is dropped (leave-one-out), so no control is recovered by its own protein; a gene whose only query was from that genome
is reported `not_tested`, not `none`. Sae2 and the Rad51 control are not repeated (the Sae2 proteins have no name evidence).

## Round 2 result (2026-11-11)

Saved as `results/round2_calls.tsv` and `results/round2_by_group.tsv`. Assembly check: the scaffold counts of all 31
`input_clean_genomes` files match `bfd.asm_stats.contig_count` (0 differences; no masked-file fallbacks).

**Ku.** With closer-relative queries all four Ku positives are recovered strongly (*Dimargaris*, *Dispira*, *Ramicandelaber*,
*Rhizopus*). tblastn also finds Ku in three genomes where Pfam called it absent (*Tieghemiomyces*, *Caulochytrium*,
*Conidiobolus*, all strong), so those were annotation misses and not losses. Ku remains **none in all 14 core Kickxellales, all 4
Harpellales and both *Mycoemilia*/*Spiromyces* genomes** (20 of 20 genomes in the 11 genera under the stem node), which supports a
real loss on the stem of Kickxellales + Harpellales. *Piptocephalis* and *Neoconidiobolus* are none and *Capillidium* is weak, so
independent Ku loss (or strong divergence) in those lineages is not excluded.

| Gene | What round 2 shows | Reading |
|---|---|---|
| Ku | controls 4 of 4 strong; none in all 20 stem-clade genomes | supported absence in the stem clade |
| ATR/Mec1 | strong in all Kickxellomycota genomes, including all 4 Harpellales | present everywhere; the weak Harpellales hits in round 1 were query divergence |
| Pol4 | strong in *Rhizopus* and Zoopagales (control passes); none in all 24 Kickxellomycota and all 4 Entomophthorales genomes | supported absence, but from two queries only, one an unusually large (1,604 aa) gene model |
| Lig4 | strong or weak hits in nearly every lineage, strong in 2 of 4 Harpellales, **0 strong and 13 weak in 14 core Kickxellales** | **unresolved**: a single query and the ligase catalytic domain is shared with LIG1; needs a reciprocal-best-hit check |
| XLF/Nej1 | none everywhere, including *Rhizopus* | uninformative (control fails; small, fast-evolving) |
| Aprataxin | strong in the root genus, one Dimargaritales and one Entomophthorales genus; none or weak elsewhere | uninformative |
| DNA-PKcs | weak (PIKK-like) hits everywhere; *Rhizopus* untested (its own protein was the only query) | uninformative |

Caveats: round 2 queries were chosen after round 1 (post hoc, declared above); Pol4 and Lig4 each rest on one or two proteins; a
"none" in a fragmented assembly can still mean a missing region; Rad51 was not repeated (round 1: 30 of 31 genomes); weak and strong
hits were not checked by reciprocal best hit against a proteome.

## Round 3 design (2026-11-11; written before the round 3 results were seen)

Round 3 is **post hoc**, like round 2: it adds queries after rounds 1 and 2 showed which genes lacked sensitivity. Same thresholds,
same calls, same leave-one-out rule, same BLAST settings. Files: `queries_round3.tsv` (94 queries = the 44 of rounds 1 and 2 plus 50
new), `genomes_round3.tsv` (128 genomes = the 31 plus the remaining 97 annotated *Coemansia* strains with BUSCO complete >= 60 and
>= 3,000 proteins; column 6 is the expected scaffold count from `bfd.asm_stats`, and a task fails if its assembly disagrees),
`coemansia_strain_table.py`, `test_tblastn_fixture.py`. Run with `ROUND=3` and `--array=1-128%12` (outputs `hits_round3/`).

New queries: Ku (8 BFD proteins from Glomeromycota, Mortierellomycota, Mucorales, Endogonales, Umbelopsidales, chytrids,
Blastocladiomycota and *Syncephalis*); Pol4 (8 BFD proteins, 450-700 aa, from Mucoromycota, Glomeromycota, Zoopagales, chytrids,
Blastocladiomycota, Monoblepharomycota); Lig4 (4 UniProt: *S. pombe*, *Coprinopsis*, *Cryptococcus*, *Candida*, plus 9 BFD proteins,
850-1,200 aa); DNA-PKcs (5 BFD proteins of 3,800-4,500 aa from Mucorales, Umbelopsidales, Glomeromycota and chytrids; the 1,999 aa
*Rhizopus* model of round 2 is kept but is a partial model); Sae2/CtIP (10 BFD proteins that Q1 called present, name evidence weak);
and the Rad51 family (Dmc1, Rad55, Rad57 from *S. cerevisiae* and *S. pombe*, keys `RecA_DMC1`, `RecA_RAD55`, `RecA_RAD57`). **These
RecA-family keys are hit evidence only:** every RecA-like protein hits every other, so a hit under a key does not say which paralog
it is. That needs the gene-family trees. Every UniProt length was checked against `bfd.swissprot_annot`, and every BFD protein must
fall in its length range (`fetch_queries.sh` fails otherwise).

**Positive controls (from the Q1 Pfam calls, fixed before the results):** Ku in *Basidiobolus*, *Conidiobolus*, *Dimargaris*,
*Dispira*, *Ramicandelaber*, *Rhizopus*; Pol4 in *Basidiobolus*, *Rhizopus*; ATR in 13 genomes; Sae2 in 10 (weak name evidence); Lig4 in
*Basidiobolus* only; aprataxin in *Basidiobolus* only; DNA-PKcs and XLF/Nej1 in *Rhizopus* only. **A gene's negatives count as
informative only if at least 75% of its control genomes are strong**; with one control, one strong hit is "one control only" (weak
evidence) and a miss makes the gene uninformative. A control is judged with its own protein left out.

## Round 3 results (2026-11-11)

All 128 genomes ran (every assembly matched its expected scaffold count). Files in `results/`: `round3_calls.tsv`, `round3_by_group.tsv`,
`round3_coemansia_strains.tsv` (per strain, with BUSCO, N50 and the `none_in_fragmented` flag), `round3_lig4_cterminal.tsv`,
`round3_query_lengths.txt`. Evaluated with `round3_evaluate.py`; `test_tblastn_fixture.py` passes.

**Controls (judged as declared above).**

| Gene | Controls strong | Verdict |
|---|---|---|
| Ku | 6 of 6 | pass |
| Pol4 | 2 of 2 | pass |
| ATR/Mec1 | 12 of 12 | pass |
| Lig4, aprataxin, DNA-PKcs | 1 of 1 each | one control only |
| XLF/Nej1 | 0 of 1 | **uninformative** |
| Sae2/CtIP | 1 of 10 (6 weak, 3 none) | **uninformative** |
| Rad51 (completeness proxy) | 126 of 128 strong (106 of 107 *Coemansia*; failures: *Dipsacomyces*, N50 15 kb, and strain RSA 2675, N50 17 kb) | assemblies are largely complete for a core gene |

**Ku: no hit at all in 0 of 107 *Coemansia* strains** (not even an HSP at e-value 1e-3 with any of the 22 Ku queries); controls pass
6 of 6. 75 of the 107 assemblies are fragmented (N50 < 50 kb or BUSCO < 70), 32 are not, and Ku is absent in both sets, including the
two with N50 >= 500 kb (RSA 1278, RSA 71). The 21 other genomes are unchanged from round 2 (Ku strong in *Dimargaris*, *Dispira*,
*Tieghemiomyces*, *Ramicandelaber*, *Rhizopus*, *Caulochytrium*, *Basidiobolus*, *Conidiobolus*; weak in *Capillidium*; none in all
Harpellales, core Kickxellales, *Mycoemilia*/*Spiromyces*, *Piptocephalis*, *Neoconidiobolus*).

**Pol4: none in all 107 *Coemansia* strains and all 14 other Kickxellomycota genomes**, with controls passing (2 of 2). Pol4 hits
exist in *Rhizopus*, *Basidiobolus* (now strong with a *Phlyctochytrium* query) and *Piptocephalis*; *Caulochytrium* is none.

**Lig4: the "strong" calls are not evidence of Lig4.** Round 3 moved 16 original-genome calls from weak to strong, 13 of them
core Kickxellales (10 *Coemansia* strains, *Kickxella*, *Linderina*, *Martensiomyces*) and 14 of the 16 through the same *Zoophagus*
query, and 102 of 107 *Coemansia* strains are strong. These hits are one block
that covers query residues 1-417 of 964 (the N-terminal ligase core that Lig4 shares with LIG1); the 30% coverage rule is met by that
core alone, so **the threshold is not specific for Lig4** (a design flaw I did not foresee). The more specific observation (`lig4_cterminal_evidence.py`, evidence
only) is whether Lig4 query hits also reach the C-terminal third, where Lig4 has its BRCT region: yes in 13 of the 17 non-core
genomes (Dimargaritales 3 of 3, *Ramicandelaber* 1 of 1, Harpellales 2 of 4, *Mycoemilia*/*Spiromyces* 1 of 2, Entomophthorales 3 of 4,
*Rhizopus*, *Caulochytrium*, *Piptocephalis*), but **in 0 of 111 core Kickxellales genomes** (maximum 0.20, median 0.12). BRCT domains
occur in many proteins, so this is supportive, not specific; the identity needs the gene-family trees.

**DNA-PKcs:** *Rhizopus* is strong (full-length *Cunninghamella* query, e = 1e-167, 89% coverage). All Kickxellomycota genomes are weak
(median coverage 0.16, the shared PIKK kinase region), except one *Coemansia* strain (RSA 485, 23.8% BUSCO-duplicated) that is strong with the
partial *Rhizopus* model: its hits are the query's C-terminal kinase region and the ATR queries hit the same scaffold, so it is a PIKK
paralog, not a DNA-PKcs call. **Aprataxin:** weak in 31 and none in 76 strains, one control only. **Sae2, XLF/Nej1: uninformative.**

**RecA-family keys (hit evidence only, not paralog calls):** `RecA_DMC1` and `RecA_RAD57` are strong in all 107 strains and `RecA_RAD55`
weak in 8 (none in 99); a Dmc1 or Rad57 query also matches Rad51, so these do not show that Dmc1, Rad55 or Rad57 exist.

**Changes from round 2 in the original 31 genomes:** 20 calls, all upgrades (16 weak to strong, 2 not_tested to strong, 1 none to strong,
1 none to weak); none downgraded. All Lig4 upgrades are explained by the shared ligase core (above); the others are DNA-PKcs in *Rhizopus*
(full-length query, previously untested) and Pol4 in *Basidiobolus*.

**Still uninformative:** XLF/Nej1, Sae2/CtIP, aprataxin and DNA-PKcs outside Mucoromycota (one control only, or weak PIKK-like hits); Lig4
identity (needs trees); Dmc1/Rad55/Rad57 identity (needs trees); anything about fragmented assemblies (75 of 107 strains) beyond core
genes such as Rad51 being present.

## Run (HPCC)

    BENCH=~/projects/BFD_dashboard/bench/kx2_ku_tblastn REPO=<checkout> bash fetch_queries.sh      # login node, light
    cd $BENCH && BENCH=$BENCH REPO=<checkout> sbatch $REPO/docs/analyses/kx2_ku_pathway/tblastn/run_tblastn.sh
    python3 summarize_tblastn.py $BENCH genomes.tsv results/tblastn

Hit tables live in `$BENCH/hits/` (purgeable); the summaries are copied into `results/`.
