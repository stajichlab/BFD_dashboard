# Is the Dmc1 loss part of a loss of the meiotic toolkit in core Kickxellales?

Date: 2026-11-12. Status: **first pass, annotation level plus tblastn plus first-pass trees** (MotherDuck SQL on `bfd`, HPCC jobs; scratch tables
`bfd_scratch.kx2_mt_*`). Dmc1 is absent from the five core Kickxellales genera (*Coemansia*, *Kickxella*, *Linderina*, *Dipsacomyces*,
*Martensiomyces*) but present in *Mycoemilia*, Harpellales, Dimargaritales, Entomophthorales and most Mucoromycota
([../family_trees/RECA_RESULT.md](../family_trees/RECA_RESULT.md)). Question: do other meiosis-specific genes go with it?

Genomes: the same 169 as the Ku analysis (63 backbone genera + 106 more *Coemansia* strains; groups as before: core Kickxellales 111,
*Mycoemilia*/*Spiromyces* 2, *Ramicandelaber* 1, Harpellales 4, Dimargaritales 3, Zoopagales 8, Entomophthorales + *Basidiobolus* 4,
Mucoromycota 21, root taxa 15).

## Directory

| Path | What |
|---|---|
| `sql/01_rules.sql` | gene and control rules, with the revision note |
| `sql/02_calls_matrix.sql` | protein-level calls and the gene x genome matrix |
| `sql/03_group_summary.sql` | control miss rates, presence by group, independent recount from `bfd.pfam`, Spo11 domain split |
| `sql/04_score.sql` | per-genome completeness score and the assembly-quality effect |
| `tblastn/` | queries, genomes, scripts, positive-control evaluator, fixtures; `results/` (small tables) |
| `trees/` | references, input builder, tree and classification scripts, fixture test |

## Design, and what was fixed when (disclosure)

- **Fixed before any group-level result was viewed:** the gene and control rules, the control rule "a gene's absence in a group is informative only
  if that group's single-domain control is missed in <= 10% of genomes", the tblastn thresholds (strong: e-value <= 1e-10 and >= 30% query
  coverage; weak: <= 1e-5; same as the Ku analysis) and the pass rule (a gene is informative only if tblastn recovers it in >= 75% of the
  genomes where annotation already calls it present in the reference groups), the tree trimming mode, and the definition of the score.
- **Revised after looking at control row counts only (no meiotic-gene group result viewed)**, documented in `sql/01_rules.sql`: the Mlh1 control
  used the wrong Pfam domain; the Msh2 name control had too few name hits; the Rad50 zinc-hook model failed everywhere; the Hop1 name rule collided
  with an unrelated protein also named HOP1.
- **Not blind:** the tblastn query set (Spo11, Msh5, Pch2, Mer3, Hop2, Mnd1) was chosen after the annotation-level group table had been seen, and
  the Spo11 domain split, the SwissProt-name evidence and the BUSCO-bin analysis were added after seeing that table. They are labelled as such below.
- Never invented: every Pfam accession and name, and every reference accession and length, was read from `bfd.swissprot_subject` / `bfd.pfam`
  (references fetched from UniProt by accession, lengths checked). No literature is cited; statements about gene function are marked citation to verify.

## Genes, identification and confidence (`sql/01_rules.sql`)

| Gene | Evidence | Confidence | Note |
|---|---|---|---|
| Spo11 | TP6A_N [PF04406](https://www.ebi.ac.uk/interpro/entry/pfam/PF04406/) + Toprim [PF21180](https://www.ebi.ac.uk/interpro/entry/pfam/PF21180/) (single domain, or name, used for the score) | high for the two-domain rule, which under-detects (see below) | in fungi the only TopoVI-A-like protein |
| Hop2 | WHD_TBPIP [PF07106](https://www.ebi.ac.uk/interpro/entry/pfam/PF07106/) | medium | |
| Mnd1 | [PF03962](https://www.ebi.ac.uk/interpro/entry/pfam/PF03962/) | medium | |
| Msh5 | MutS5_N [PF27035](https://www.ebi.ac.uk/interpro/entry/pfam/PF27035/) or best-hit name | medium | the Pfam model fails in Kickxellales |
| Msh4 | best-hit name only; tree | medium | no specific domain |
| Mer3/Hfm1 | DEAD + Helicase_C + Sec63 with name; or name | medium | Brr2 and other Sec63 helicases cross-hit |
| Pch2/TRIP13 | lid [PF23242](https://www.ebi.ac.uk/interpro/entry/pfam/PF23242/) / N [PF23563](https://www.ebi.ac.uk/interpro/entry/pfam/PF23563/) or name | medium | domains are specific to TRIP13 |
| Mlh3, Rec8, Hop1 | tree (paralogs of Mlh1/Pms1, Rad21, Mad2/Rev7) | see trees | annotation alone cannot separate them |
| Mei4, Mei5, Mer2, Rec114, Red1, Spo16, Mek1, Ndt80 | Pfam | low | **uninformative**: not detectable by Pfam outside yeasts even in the outgroups (Mei4, Mer2 and Red1 in 0 of 21 Mucoromycota; Spo16 in 3; Rec114 in 3; Mei5 in 9) |
| Zip1, Zip2, Zip3 | none | not assessed | yeast proteins have no Pfam; the names collide with zinc transporters |
| Controls | Mre11, Mcm2, Pol delta, Pol epsilon, Rad51 family, MutL, MutS, kleisin, Mlh1, MutS_I+II | high | single- and multi-domain styles |

## Controls and false-absence

Single-domain controls are missed in **0-5.6% of genomes in every group (3.7% in core Kickxellales)**, under the 10% threshold, so single-domain absence calls
are informative. Multi-domain controls are missed in about 12% of core Kickxellales genomes (Mre11 16%, and the Rad51 two-domain rule 38% because HHH_5 PF14520 is not detected
in this lineage), so multi-domain rules carry more false-absence risk, as in the Ku analysis. The score uses single specific domains or names to limit this.

## Annotation-level result by group (genomes with the gene / genomes; `sql/03_group_summary.sql`)

| Gene | core Kx (111) | Myco/Spiro (2) | Rami (1) | Harp (4) | Dimarg (3) | Zoop (8) | Ento (4) | Mucor (21) | Root (15) |
|---|---|---|---|---|---|---|---|---|---|
| Hop2 | 98 | 1 | 0 | 4 | 2 | 6 | 4 | 20 | 15 |
| Mnd1 | 104 | 2 | 1 | 4 | 3 | 6 | 4 | 20 | 14 |
| Swi5/Sae3 (not meiosis-specific) | 94 | 2 | 1 | 4 | 3 | 8 | 1 | 16 | 12 |
| Msh5, two-domain Pfam rule | 0 | 0 | 0 | 4 | 1 | 5 | 2 | 16 | 12 |
| Msh5, best-hit name | 62 | 1 | 0 | 4 | 1 | 6 | 4 | 15 | 12 |
| Msh4, best-hit name | 94 | 1 | 1 | 4 | 1 | 7 | 4 | 18 | 12 |
| Spo11, two-domain Pfam rule | 8 | 1 | 1 | 4 | 1 | 6 | 2 | 18 | 11 |
| Spo11, either domain | 54 | 1 | 1 | 4 | 2 | 7 | 2 | 20 | 15 |
| Mer3/Hfm1, best-hit name | 27 | 2 | 1 | 4 | 2 | 5 | 2 | 14 | 10 |
| Pch2, TRIP13 domains (lid or N) | 0 | 0 | 1 | 0 | 2 | 8 | 1 | 7 | 12 |

Independent recount of Spo11, Msh5, Hop2 and Mnd1 straight from `bfd.pfam` matches these counts in every group (verification, `sql/03_group_summary.sql` part c).

## What the evidence says, gene by gene (observed)

1. **Hop2 and Mnd1 are retained** in core Kickxellales (98/111 and 104/111, against 20/21 and 20/21 in Mucoromycota; the shortfall is in the poorest assemblies). Their tblastn controls pass (Hop2 10/12, Mnd1 13/14).
2. **Spo11 is present but diverged; the first reading of a loss was an artifact.** The two-domain rule finds it in only 8/111 core genomes, because the N-terminal TP6A_N model is detected in 9 of them,
   whereas the Toprim domain is detected in 53/111 (51 of 107 *Coemansia* strains, *Martensiomyces* and *Dipsacomyces*). Toprim-only proteins carry Spo11-family best hits (15 proteins named spo-11, SPO11-1/2 or rec12 at 26-38%
   identity) or no SwissProt hit (35 proteins), the signature of a diverged protein. **By assembly quality, Spo11 is detected in 19 of 21 *Coemansia* strains with BUSCO >= 90, in 17/33 at 80-90, 10/21 at 70-80 and 6/32 below 70**, so the rest is annotation in fragmented assemblies. tblastn fails the control rule for Spo11
   (6 of 10 positives recovered, under 75%), so it is uninformative here.
3. **Msh5 is probably retained; the Pfam rule is wrong in this lineage.** The MutS5_N rule finds 0/111, but the best-hit name finds Msh5 in 62/111 (Harpellales 4/4, Mucoromycota 15/21), Msh4 in 94/111, and **the MutS family has the usual size in core Kickxellales (median 6 proteins per genome; 6 in Mucoromycota)**,
   which a loss of two paralogs would reduce. Detection stays near 50% regardless of assembly quality (11/21 at BUSCO >= 90), so this is not an assembly effect; the tree classifies it (below).
4. **Pch2 is absent from Kickxellales and Harpellales by two annotation methods:** the two TRIP13-specific domains are found in 0/111 core, 0/4 Harpellales and 0/2 *Mycoemilia*/*Spiromyces*, and no best hit is named Pch2/TRIP13 (0/121), against Zoopagales 8/8, root taxa 12/15, Dimargaritales 2/3, *Ramicandelaber* 1/1. It is also 0/21 in the best *Coemansia* assemblies. This fits a loss on the stem of Kickxellales + Harpellales, the same stem as the Ku loss;
   Mucoromycota show 7/21, so Pch2 is not universally retained either. tblastn cannot test it (the AAA ATPase family cross-hits).
5. **Mer3/Hfm1 is reduced in annotation** (27/111 core vs 62-100% elsewhere), but it is a 1,100-1,400 aa helicase in fragmented assemblies; in the best *Coemansia* assemblies it is detected in 15/21 (71%), similar to the outgroups, so it is probably retained.
6. **Uninformative genes** (Mei4, Mei5, Mer2, Rec114, Red1, Spo16, Mek1, Zip1-3): not detected even in the outgroups, so their absence from Kickxellales carries no information.

## tblastn (128 genomes; `tblastn/`, `results/meiosis_*.tsv`)

Positive-control gate: Msh2, Mlh1, Rad51 recovered in 14/14 reference genomes, Rad21 14/14; Hop2 10/12; Mnd1 13/14; **Spo11 6/10 (fails, uninformative)**. Msh5, Pch2 and Mer3 formally pass (8/8, 5/5, 6/6, small n) but are
**strong in all 111 core Kickxellales genomes, which is not specific**: their queries match the ubiquitous Msh2/6, other AAA ATPases and Brr2-like helicases, so the strong hits identify the family, not the gene (rule fixed in advance: a hit alone does not identify a paralog). Msh4, Hop1, Rec8 and Mlh3
are hit evidence only. The tree and the domain-specific evidence above are therefore the basis for those genes, not tblastn.

## Toolkit completeness score (`sql/04_score.sql`; defined before group values were viewed)

Seven informative genes (Hop2, Mnd1, Msh4, Msh5, Spo11, Pch2, Mer3), one point for each detected by any specific evidence. Dmc1 is not in the score.

| Group | n | median score | genomes with score >= 6 |
|---|---|---|---|
| core Kickxellales | 111 | 4 | 16 |
| *Mycoemilia*/*Spiromyces* | 2 | 4 | 0 |
| *Ramicandelaber* | 1 | 5 | 0 |
| Harpellales | 4 | 6 | 4 |
| Dimargaritales | 3 | 5 | 1 |
| Zoopagales | 8 | 6.5 | 5 |
| Entomophthorales + *Basidiobolus* | 4 | 5.5 | 2 |
| Mucoromycota | 21 | 6 | 15 |
| root taxa | 15 | 7 | 11 |

The score is lower in core Kickxellales, but **inside *Coemansia* it depends on assembly quality** (Pearson 0.54 with BUSCO, 0.52 with ln N50; median 5 at BUSCO >= 90 vs 3 below 70), and the best bin loses only Pch2 completely
(Pch2 0/21; Hop2 21/21, Mnd1 21/21, Msh4 20/21, Spo11 19/21, Mer3 15/21, Msh5 11/21).

## Gene-family trees (`trees/`, `trees/results/`; IQ-TREE, annotated layer only, 14 core Kickxellales genomes)

Four trees of the annotated proteins with the relevant domains plus UniProt references: MutS (401 leaves, 955 sites), MutL (190 leaves, 830 sites), kleisin (146 leaves, 570 sites)
and HORMA (144 leaves, 221 sites). Classification uses the same rule as the RecA tree (nearest reference-containing ancestor holds one clade, UFBoot >= 90 and SH-aLRT >= 80).
**Gate (fixed in advance):** a clade's absence calls are used only if every reference of that clade is re-assigned correctly in leave-one-out (`trees/results/classification_summary.txt`).

| Tree | leave-one-out (references correct) | clades that pass the gate | clades that fail |
|---|---|---|---|
| MutS | 24/27 | Msh3 6/6, Msh6 5/5, **Msh4 4/4, Msh5 4/4** | Msh2 3/5, Msh1 2/3 |
| MutL | 10/15 | Pms1 5/5, **Mlh3 3/3** | Mlh1 2/5, Mlh2 0/2 |
| Kleisin | 2/8 | none | Rad21 0/4, Rec8 2/4 (**tree unusable; no Rec8 call**) |
| HORMA | 11/13 | Mad2 6/6, Rev7 5/5 | **Hop1 0/2 (no Hop1 call)** |

Genomes with the clade present / genomes in the tree (`*_by_group.tsv`; "not seen" = no member assigned or named for that clade):

| Clade | core Kx (14) | Harp (4) | Dimarg (3) | Mucor (21) | Root (15) |
|---|---|---|---|---|---|
| Msh4 | 14 | 4 | 2 | 18 | 13 |
| Msh5 | 7 | 4 | 2 | 16 | 12 |
| Mlh3 | **0** | **0** | **0** | 6 | 3 |
| Pms1 (control) | 14 | 4 | 3 | 20 | 15 |
| Mad2 (control) | 14 | 4 | 3 | 20 | 12 |

1. **Msh4 is retained in all 14 core genomes** (the tree confirms the name-based call). **Msh5 is placed in 7 of 14 core genomes**, against 4/4 Harpellales and 16/21 Mucoromycota: reduced but not absent, and not
   a clean loss (a single loss on the core stem would give 0 of 14). The unplaced genomes are unclassified by an annotation-only layer, and the name-based recount and the BUSCO-bin result above (about half detected at every assembly quality) are consistent with a poorly annotated, diverged Msh5, not an absent one.
   The difference between 7/14 and 16/21 was not tested statistically. **Verdict: Msh5 not lost; call remains partial.**
2. **Mlh3 is not seen in any Kickxellomycota genome in the tree (0/21)**, including Harpellales and Dimargaritales, and is present in only 6 of 21 Mucoromycota and 3 of 15 root taxa. Its gate passes, so the absence is informative, but it is
   **not specific to core Kickxellales**: it is lost, or too diverged or rare to annotate, in the whole Kickxellomycota sample and is patchy even in the outgroups. Mlh1 fails the gate (2/5), so it is not called.
3. **Kleisin (Rec8) and HORMA (Hop1) cannot be used**: the trees do not recover their references, and the "uncertain" rows (for example Rec8 14/14 uncertain in core Kickxellales) are not evidence of presence. Mad2 is present everywhere (control).
   Incidentally, Rev7 (the HORMA subunit of translesion polymerase zeta, not a meiotic protein) passes its gate and is placed in 0 of 21 Kickxellomycota genomes but 7/21 Mucoromycota; not pursued.

## Verdict: is the meiotic toolkit eroded in core Kickxellales, or only Dmc1 lost?

**Mostly only Dmc1, plus two earlier or broader losses; no broad erosion specific to the core clade is shown (observed, with limits).**

| Gene | Result | Node (descriptive, not a reconstruction) |
|---|---|---|
| Dmc1 | absent from the five core genera (RecA tree, 5/5 in leave-one-out) | stem of core Kickxellales |
| Pch2 | absent from all Kickxellales and Harpellales (two domain-based methods; 0/111 core, 0/4 Harpellales) | stem of Kickxellales + Harpellales, the Ku stem |
| Mlh3 | not seen in any Kickxellomycota genome in the tree (0/21), patchy outside | before the core, or lost with low annotation sensitivity |
| Spo11, Hop2, Mnd1, Msh4 | retained (Spo11 diverged: 19/21 in the best assemblies; Msh4 14/14 by tree) | not lost |
| Mer3 | probably retained (15/21 in the best assemblies) | not lost |
| Msh5 | placed in 7/14; partial | not called |
| Rec8, Hop1 | not callable (trees fail the gate) | unknown |
| Mei4, Mei5, Mer2, Rec114, Red1, Spo16, Mek1, Zip1-3 | uninformative (undetectable even in the outgroups) | unknown |

The completeness score (median 4 of 7 in core Kickxellales against 6-7 elsewhere) points to erosion, but within *Coemansia* it tracks assembly quality (Pearson 0.54 with BUSCO), and in the best assemblies only Pch2 is entirely missing.
So what is lost is Dmc1 (core clade), Pch2 (earlier, with Ku) and probably Mlh3; the genes that act with or around Dmc1 and Spo11 (Hop2, Mnd1, Msh4, Spo11) are still there. **Hypothesis, not tested:** the retention of Hop2/Mnd1 (Dmc1 partners in
other fungi, citation to verify) without Dmc1 could mean that they serve Rad51-dependent repair in this lineage, or that a Dmc1 is present but not annotated (the tblastn locus count gives one Rad51-like locus in 100/111 core genomes, which supports loss
but under-detects). Whether core Kickxellales reproduce sexually by meiosis is not addressed by BFD data.

## Limits

- Annotated layer only (the rescue layer was not usable at tree level, as for the RecA tree); every "not seen" is subject to annotation misses in fragmented assemblies. The trees contain 14 core genomes, not all 111, so tree-level calls are not per-genome counts for the whole clade.
- Msh2, Msh1, Mlh1, Mlh2, Rad21, Rec8 and Hop1 fail leave-one-out, so those genes have no tree call (of these only Rec8 and Hop1 are meiosis-specific).
- tblastn is uninformative for Spo11, Msh5, Pch2 and Mer3 (control failure or family cross-hits, above); the domain and tree evidence carries these calls.
- Single ML tree per family; no site-heterogeneous model checks.
- The gene list and rules are disclosed above; the score and gate were defined before group values were viewed. No literature is cited; functional statements are citations to verify.
