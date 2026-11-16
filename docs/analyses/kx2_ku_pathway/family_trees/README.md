# Gene-family trees to classify Lig4, Dmc1/Rad55/Rad57 and Pol4 candidates

Date: 2026-11-11. Status: **first pass** (maximum-likelihood trees on the HPCC; scripts and small results here, big files in
`~/projects/BFD_dashboard/bench/kx2_ku_family/`, which is purgeable).

## Why

A tblastn hit cannot tell Lig4 from LIG1 or LIG3 (they share the ligase catalytic domain), or Dmc1 from Rad51 and Rad55/Rad57 (the same
RecA/Rad51 domain), and Pol4 sits in a family with Pol lambda, mu, beta and TdT. The tblastn rounds in `../tblastn/` found "strong" Lig4 hits in
lineages where Pfam found none, so the question is whether those hits are Lig4. A gene-family tree with verified reference sequences answers it
by placement, and extends to the Dmc1/Rad55/Rad57 calls (Q2, name-based) and Pol4.

## Design

| Step | File | What it does |
|---|---|---|
| 1 | `references.tsv`, `01_fetch_references.sh` | 22 ligase, 31 RecA-family and 15 Pol X reference proteins. Accessions were looked up in `bfd.swissprot_subject` (gene, organism, Pfam) and each fetched UniProt length is checked against `s_len`; a mismatch stops the run. Clades: LIG1, LIG3, LIG4, LIG6; Rad51, Dmc1, Rad55, Rad57, Rad51B/C/D, XRCC2, XRCC3, RadA (archaeal outgroup); Pol4, Pol lambda, Pol mu, TdT, Pol beta. |
| 2 | `02_candidates.py/.sh` | Annotated candidates from the BFD Pfam table: proteins with the family domain at domain e-value <= 1e-5 in the 63 backbone genomes and the 9 extra *Coemansia* strains (72 genomes). Domains are Pfam names checked in the data: ligase PF01068 DNA_ligase_A_M; RecA/Rad51 PF08423 Rad51; Pol X PF14792 DNA_pol_B_palm or PF14791 DNA_pol_B_thumb. |
| 3 | `03_extract_proteins.py/.sh` | Sequences of those proteins from the BFD proteomes (read-only). |
| 4 | `04_miniprot.sh`, `04b_translate_models.py`, `04t_test_intron_penalty.sh` | **Rescue of loci annotation missed**: miniprot (0.2-r116) maps the references plus every annotated family protein (closer relatives, found necessary because the distant references alone gave no Pol X model at all) to the 31 worklist assemblies (`input_clean_genomes/<ASMID>.fa.gz`). Settings `-G 5k -J 16`; the best model per locus is translated by a tested, phase-aware translator. |
| 5 | `05_build_inputs.py/.sh` | Two inputs per family, as **layers**: `annotated` (references + annotated BFD candidates; does not depend on the rescue) and `rescue` (the same plus rescued models, kind=`rescued`). A rescued model overlapping an annotated candidate on the same scaffold and strand by >= 50% of the shorter locus is a duplicate and dropped. Genomes whose rescue did not finish are listed in `rescue_failed_genomes.txt`. |
| 5v | `05v_validate_models.py/.sh` | Validation of the rescue against annotated genes (see Results). |
| 6 | `06_trees.sh` (`LAYER=`), `06b_retrim_tree.sh` | mafft 7.505 `--auto`, trimAl 1.5.1 `-automated1`, IQ-TREE 3.0.1 ModelFinder (`-mset LG,WAG,JTT,Q.pfam`), 1000 ultrafast bootstraps and 1000 SH-aLRT. `06b` is a one-off re-trim (`-gt 0.4 -resoverlap 0.5 -seqoverlap 50`) for trees whose automated trim left too few columns. |
| 7 | `07_classify.py`, `07a_ligase_architecture.py`, `07_classify.sh` (`<layer>`) | Placement, leave-one-out, per-genome and per-group tables; ligase domain architecture for annotated proteins. |
| tests | `tests/test_fam_lib.py`, `fam_lib.py` | 14 fixture tests (spliced models on both strands, phase-aware translation and frameshifts, internal stops, intron span, strand normalisation, locus de-duplication, tree parsing, assigned / ambiguous / low-support placement, leave-one-out). Run `python3 -I tests/test_fam_lib.py`. |

### Decisions and thresholds (fixed before the trees existed)

- **Minimum length** for annotated and rescued members: ligase 300 aa, RecA/Rad51 200 aa, Pol X 200 aa (about the core domain). Models with an
  internal stop codon, an intron over 5 kb, or below the minimum are recorded in the tables but not put in the tree.
- **Placement rule** (`fam_lib.classify`): a member is *assigned* to a reference clade when the nearest ancestor that contains references holds only
  that clade's references and has UFBoot >= 90 and SH-aLRT >= 80. Single-clade but below the thresholds is *low_support*; mixed references is
  *ambiguous*. Long branches are flagged (terminal branch above the 95th percentile). The thresholds are conventional cut-offs for these two
  supports, not tuned.
- **Presence per genome**: in the 31 genomes searched with miniprot, *present* = an assigned member (annotated or rescued), *uncertain* = only
  low_support or ambiguous members, *absent* = nothing placed in that clade. The other 41 genomes were not searched for missed loci, so their
  call is "not seen in annotation", which is weaker.
- **Standard genetic code** for translating models. Whether any of these genomes needs an alternative code (BFD's `samples` table has a translation-table field) was **not** checked; see limits.
- **Validation criterion for the rescue** (fixed before the results): at least 80% of the annotated loci that miniprot recovers must translate cleanly, and the median identity of those models to the annotated protein must be >= 0.95.

### Problems found and fixed on the way (the rescue took several rounds; each is why a setting is what it is)

1. **Chimeric models.** miniprot's default maximum intron size (200 kb) stitched distant exons together (the first *Coemansia spiralis* ligase models
   spanned 92 and 101 kb; 20 of 30 ligase loci had internal stops). Introns here are short (median about 96 bp in *Coemansia*), so the cap is 5 kb.
2. **De-duplication matched nothing.** `gene_info` stores the strand as `1` / `-1`, miniprot writes `+` / `-`, so no rescued model was ever recognised as
   an annotated gene (0 duplicates, 22 "rescued" ligase models that were all annotated genes). Fixed with `fam_lib.norm_strand` and a test.
3. **No Pol X models from distant references.** Yeast and animal Pol X references gave no miniprot model in any genome, even *Rhizopus*, where annotation has
   Pol X proteins. The annotated family proteins are now also used as queries (a re-found annotated gene is removed as a duplicate).
4. **Short introns read as coding.** On the annotated *Rhizopus* Pol X gene (507 aa, one intron of about 41 bp) the default intron penalty gives a 1,611 nt model,
   `-J 16`, 13, 10 and 5 give exactly the annotated 1,521 nt. `-J 10` crashes miniprot 0.2-r116 on 13 of 31 genomes (an internal assertion), `-J 13`
   and up ran on the crashing test genome, so `-J 16` is used.
5. **Frameshifts.** The translator ignored the GFF3 phase after the first block, so a frameshift gave a wrong reading frame; it now follows each block's phase,
   counts frame breaks and excludes such models. The earlier tests used a wrong phase and passed only because of this (corrected).
6. **Alignment trimming.** trimAl `-automated1` left the RecA/Rad51 alignment with 41 columns and the rescue-layer ligase and Pol X alignments with 13
   and 65 (see Results); a gentler re-trim (`06b`) was written and queued for RecA.

## Results (2026-11-11, first pass)

### Rescue validation (`results/validation_summary.txt`, `results/validation_models.tsv`)

Final rescue run: miniprot `-G 5k -J 16` on all 31 genomes, **0 genomes failed**. Against annotated genes in those genomes (reciprocal overlap >= 90%):
ligase 35 of 35 annotated loci recovered, all translated cleanly, 32 identical to the annotated protein, median identity 1.000; RecA/Rad51 130 of 130
recovered, 121 clean (93%), median identity 1.000; Pol X only 2 annotated loci in these genomes (both recovered and identical), so Pol X is barely validated.
The criterion was passed for all three, i.e. the rescue **translates real genes correctly**; it does not say that rescued models in other loci are real genes.

### Trees (`results/annotated/*.treefile`, versions in `*_versions.txt`; ModelFinder BIC, IQ-TREE 3.0.1 with 1000 ultrafast bootstraps and 1000 SH-aLRT)

| Tree (layer) | Sequences | Alignment columns after trimming | Model | Leave-one-out (references re-assigned to own clade) |
|---|---|---|---|---|
| ligase, annotated | 203 | 330 | LG+R8 | 15 of 21; **all 11 Lig4 references correct**; failures are 5 LIG1 references (animal and fungal LIG1 do not resolve as one supported clade) and *S. pombe* adl1 (LIG3, placed with LIG6) |
| Pol X, annotated | 87 | 205 | LG+I+R5 | 13 of 15; the two Pol4 references are ambiguous |
| RecA/Rad51, annotated (first build) | 362 | **41** | Q.PFAM+R5 | 13 of 31 (**not usable**, too few columns) |
| ligase, rescue layer | 327 | **13** | Q.PFAM+G4 | 7 of 21 (**not usable**) |
| Pol X, rescue layer | 178 | **65** | LG+G4 | 11 of 15 |

### Ligases: is it Lig4? (annotated layer; `results/annotated/lig_by_group.tsv`, `lig_classification.tsv`, `lig_architecture.tsv`, `lig4_vs_tblastn.tsv`)

Observed:

| Group (genomes) | LIG4 | LIG1 | LIG3 | Annotated ligases with the Lig4 domain PF11411 / any BRCT |
|---|---|---|---|---|
| Mucoromycota (21) | 20 | uncertain 21 | 12 | 8 of 81 / 22 |
| Root taxa (15) | 14 | uncertain 14 | 8 | 3 of 56 / 8 |
| Zoopagales (8) | 4 | 1 assigned, 6 uncertain | 0 | 0 of 14 / 3 |
| Entomophthorales + *Basidiobolus* (4) | 1 | 3 | 0 | 1 of 6 / 1 |
| Harpellales (4) | **0** | 4 assigned | 0 | 0 of 4 / 0 |
| Dimargaritales (3) | **0** | uncertain 3 | 0 | 0 of 4 / 0 |
| *Ramicandelaber* (1) | **0** | 0 (no ligase placed) | 0 | 0 of 1 / 0 |
| Core Kickxellales (14; 10 are *Coemansia* strains) | **0** | uncertain 14 | 0 | 0 of 14 / 0 |
| *Mycoemilia* + *Spiromyces* (2) | **0** | uncertain 1 | 0 | 0 of 1 / 0 |

- **No Kickxellomycota genome has a protein in the Lig4 clade, and no Kickxellomycota ligase carries the Lig4-specific domain or a BRCT.** Lig4 is present in the
  Mucoromycota, root taxa and Zoopagales, and in one *Basidiobolus* protein.
- **The strong tblastn "Lig4" hits in Kickxellomycota were not Lig4**: the Kickxellomycota genomes with a strong round-2 tblastn Lig4 hit other than *Ramicandelaber* (*Dispira*, *Mycoemilia*,
  *Zancudomyces*, *Capniomyces*) have an annotated ligase that places with LIG1 (assigned in the Harpellales, uncertain elsewhere) and none in the Lig4 clade.
  *Ramicandelaber* is the exception: it has no annotated ligase placed anywhere, though tblastn is strong there (a missed or short gene; the rescue could not settle it, below).
- The LIG1 and LIG3 calls are weaker than the Lig4 call (LIG1 leave-one-out 2 of 7); read "LIG1 uncertain" as "a ligase of the non-Lig4 kind".
- Interpretation (not shown by the data alone): Lig4 was lost on the Kickxellomycota stem, before the Ku loss on the Kickxellales + Harpellales stem. This **supports** the Q1 Pfam
  result ("Lig4 lost before Ku") and removes the doubt that the tblastn hits had raised.

Annotated layer only; Lig4 absence is from annotated proteins, so a missed gene model in a genome would not be seen except through the rescue layer (below).

### Pol X (annotated layer)

Only Pol beta is assigned (14 members). All other Pol X members (36 Pol4-like, 20 Pol lambda-like) are low-support, and the two Pol4 references fail leave-one-out, so
this tree **cannot tell Pol4 from Pol lambda or mu**. No Kickxellomycota genome has an annotated Pol X protein at all (none in the Dimargaritales, Harpellales, Kickxellales,
*Mycoemilia*/*Spiromyces* or *Ramicandelaber* genomes). That is consistent with the tblastn Pol X absence but is not a classification.

### What the rescue layer added

The rescue added 124 ligase and 92 Pol X candidate models that are not re-finds of annotated genes (62 and 38 of them in the core Kickxellales). **None could be placed**: ligase
113 ambiguous and 11 low-support LIG6-like, Pol X 89 ambiguous and one Pol beta. **This is not evidence about the models**: the rescue-layer ligase and Pol X trees were built on only
13 and 65 alignment columns (trimAl removed most of the divergent rescued models' columns), so they could not place anything and their leave-one-out fell (7 of 21 for ligases). The rescue layer is therefore
**untested at tree level**, not shown to be wrong or uninformative. It does not change the Lig4 answer above, which rests on the annotated layer, but it also means the rescue has **not yet excluded**
a Lig4 (or Pol4) that annotation missed. Many of the rescued models may be other nucleotidyltransferases (the ligase domain PF01068 is shared with mRNA capping enzyme); that was not checked.

## Pending (RecA / Dmc1 / Rad55 / Rad57 is NOT classified yet)

The RecA/Rad51 tree above has 41 columns and is not usable. A re-trimmed tree is queued as job **32026340** (`06b_retrim_tree.sh`, LAYER=annotated FAM=reca; it keeps the
old tree as `reca.treefile.automated1`). When it has finished, on the HPCC:

    cd ~/projects/BFD_dashboard/bench/kx2_ku_family
    sacct -X -j 32026340 -o JobID,State,Elapsed                     # expect COMPLETED
    grep -E "Input data|Model of substitution" trees_annotated/reca.iqtree   # the alignment should have a few hundred columns, not 41
    sbatch scripts/07_classify.sh annotated                          # re-classifies all three families; reads trees_annotated/*.treefile
    # then read: trees_annotated/reca_loo.tsv (leave-one-out), reca_by_group.tsv (Rad51 / Dmc1 / Rad55 / Rad57 per group), reca_classification.tsv

Check the leave-one-out first: if Dmc1, Rad55 and Rad57 references are not re-assigned correctly, their absence calls are not usable. Optional (rescue layer): `LAYER=rescue FAM=lig sbatch --export=ALL scripts/06b_retrim_tree.sh`
(and `FAM=polx`, `FAM=reca` once its tree exists), then `sbatch scripts/07_classify.sh rescue`. Copy the small tables and the Newick tree into `results/annotated/` afterwards.

## Limits

- **First pass**: single maximum-likelihood trees; no alternative topologies or AU tests; support is UFBoot/SH-aLRT with conventional cut-offs.
- The ligase result is for **annotated** proteins. Absence of a Lig4-clade protein in a genome is only as good as the annotation, and the rescue that was meant to cover that gap is untested at tree level.
- Reference sets are small for several clades (Rad55 and Rad57: two each; Pol4: two); LIG1 does not resolve as one clade.
- Genetic code not checked per genome; one genome per genus for most backbone taxa, 10 *Coemansia* strains, so "core Kickxellales 0 of 14" is 14 genomes of 5 genera, mostly one.
- Pol X validation of the rescue rests on 2 loci. PF01068 and PF08423 also occur in non-family proteins (capping enzyme, other RecA-like ATPases), which sit outside the reference clades.
- No literature is cited; Pfam names and accessions were read from the data (bfd.swissprot_subject / bfd.pfam), each reference length matched UniProt.
