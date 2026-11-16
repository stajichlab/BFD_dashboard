# T5: where did fungal GH65 acid trehalase come from, and is "GH65 presence" the same as acid trehalase?

Date: 2026-11-15. Status: **design and kill rules fixed before any tree was built; results are added below the line "Results".**
Part of the trehalose analysis ([`../README.md`](../README.md), section 7, hypothesis T5). This folder does not edit that README.

## The hypothesis, and how it is read here

T5 as written in the parent README: *"GH65 acid trehalase is a Dikarya-stem innovation with repeated losses (0 % in Mucoromycetes, 3 % in
Agaricomycetes), not a single gain"*, to be tested by *"a gene tree from a curated seed alignment over all genera + a genus tree; also test the HGT alternative
(bacterial GH65 sisters, near-identical GC to host)"*, and killed by *"GH65 proteins that are arbutinases/alpha-xylosidases rather than trehalases"*.
The phrase "not a single gain" is ambiguous next to "Dikarya-stem innovation". **It is read here as: one origin of the fungal acid trehalase on or before the Dikarya stem, followed
by repeated losses, rather than patchiness produced by many independent gains (including horizontal transfer from bacteria).** The competing readings are tested, not assumed:

| | Reading | What it predicts in the gene tree and the data |
|---|---|---|
| H-T5 | one origin at/near the Dikarya stem, many losses | fungal trehalases in one clade, only Dikarya genera inside it, several independent loss events on the genus tree |
| H-old | origin before Dikarya (Mortierellomycetes carry GH65, README 4.2) | the clade also contains Mucoromycota; losses in Mucorales, Tremellomycetes |
| H-multi | several independent gains (e.g. HGT) | fungal GH65 proteins fall into more than one clade separated by non-fungal sequences |
| H-HGT | one transfer from bacteria | the fungal clade nests among bacterial sequences only; gene GC and intron content differ from the host genome |
| H-mis | the Pfam column "GH65" is mostly not acid trehalase | most fungal PF03632 proteins group with phosphorylase, hydrolase or animal glucosidase references |

## Design (fixed in advance)

| Choice | Value | Why |
|---|---|---|
| Genomes | the 1,654 atlas genus representatives (`atlas_data/phylo/genus_representatives.tsv`; best genome per genus) | one genome per genus, the same genomes as the atlas genus tree used for the loss count |
| Candidates | every protein with the GH65 central domain [PF03632](https://www.ebi.ac.uk/interpro/entry/pfam/PF03632/) (Glyco_hydro_65m), `domain_i_evalue <= 1e-5`, from the local `BFD.duckdb`; architecture noted with [PF03636](https://www.ebi.ac.uk/interpro/entry/pfam/PF03636/) (N-terminal) and [PF03633](https://www.ebi.ac.uk/interpro/entry/pfam/PF03633/) (C-terminal) | same domain rule as the trehalose survey (CAZy family [GH65](https://cazy.org/famGH65)) |
| References (the gate) | every reviewed UniProtKB entry with a CAZy GH65 or Pfam PF03632 cross-reference, length 500-1500 aa; near-identical same-species entries collapsed (`01_fetch_references.py`). Clade label from the UniProt EC number: `TREHALASE_FUNGI` (EC 3.2.1.28, Fungi), `ATHL1_METAZOA` (EC 3.2.1.107, animal glucosylgalactosylhydroxylysine glucosidase "acid trehalase-like"), kojibiose hydrolase / phosphorylase, trehalose phosphorylase, trehalose-6-phosphate phosphorylase, maltose phosphorylase, nigerose, glucosylglycerol, glucosylrhamnose phosphorylase, other | the non-trehalase groups are the kill check; accessions and EC numbers come from UniProt, not from memory |
| Context taxa (not in the gate) | unreviewed PF03632 entries, one per genus, sampled per bacterial phylum (14 genera each), Archaea (20), non-fungal Eukaryota in three pools (animals 15, plants 10, others 25), seed 20261113; function not used | to see what the fungal clade is closest to; UniProt returns at most 4 x 500 entries per query in its own order, so this is a spread, not a random sample |
| Alignment, tree | mafft `--auto`, trimAl `-gt 0.4 -resoverlap 0.5 -seqoverlap 50` (columns present in >= 40% of sequences; the same setting as the RecA tree), IQ-TREE 3.0.1, a fixed model (Q.PFAM+R8, as the RecA tree) with 1,000 UFBoot and 1,000 SH-aLRT | the family is diverged; trimming and model are stated before the run and are not tuned to the result |
| Classification | `fam_lib.classify`: a candidate is assigned to a reference group when the nearest ancestor containing references holds one group and has UFBoot >= 90 and SH-aLRT >= 80 | the rule used for the Ku, ligase and RecA trees |
| Species tree | the atlas genus tree (`atlas_data/phylo/results/protein/tree/fungi_odb10/final_tree.nw`, 1,654 genera), read-only | already built for exactly these genomes |
| Quality filter for absences | genus representative with BUSCO complete >= 85 and N50 >= 50 kb; others are "unknown", not "absent" (an unfiltered count is also shown) | an absence in a fragmented or incomplete genome is weak |

## Kill rules (fixed before the run)

**Gate.** The tree's classification of fungal proteins is used only if all `TREHALASE_FUNGI` references are re-assigned to their own group in leave-one-out and the smallest clade
holding them has UFBoot >= 90. If the gate fails, no test below is evaluated and the result is "tree cannot decide".

1. **R1, mis-assignment (H-mis).** Population: fungal PF03632 proteins in the genus representatives. Killed if >= 30% are *assigned with support* to a group other than `TREHALASE_FUNGI`
   (ATHL1, phosphorylases, hydrolases, other). Ambiguous or unsupported proteins are "unresolved", reported separately, and do not count as mis-assigned.
2. **R2, one origin (against H-multi and H-HGT).** Take the smallest clade holding every protein assigned to `TREHALASE_FUNGI`. One origin is supported if that clade contains no non-fungal tip
   and its root has UFBoot >= 90; rejected if it contains any bacterial, archaeal or other non-fungal tip, or if assigned fungal proteins sit in two supported (UFBoot >= 90) clades separated by non-fungal tips.
3. **R3, Dikarya stem (H-T5 against H-old).** If any genus in the `TREHALASE_FUNGI` clade is outside Ascomycota and Basidiomycota, in a genome with BUSCO >= 85 and N50 >= 50 kb and with the full domain
   architecture (PF03636 + PF03632 + PF03633), the origin is older than the Dikarya stem (H-old) and the "Dikarya-stem" part of T5 is rejected. If only Dikarya genera are inside, it is supported.
4. **R4, repeated losses.** Dollo parsimony on the genus tree: presence = genus has a protein in the `TREHALASE_FUNGI` clade; a genus counts as absent only if its representative passes the quality
   filter; the gain is placed at the MRCA of the presence tips; losses = maximal clades inside that subtree in which every scored genus is absent (so they are not nested). Supported if >= 3 losses
   in >= 2 classes; rejected if <= 2 losses. Reported for the filtered and the unfiltered set, plus the number of losses expected when the presence labels are permuted over the scored tips (1,000
   permutations, same prevalence) to show whether absences are clustered.
5. **R5, horizontal transfer (H-HGT).** Two lines, both needed. (a) **Topology, rooting-free:** climb from the root of the `TREHALASE_FUNGI` clade to the first ancestor that contains a non-fungal tip; if every
   non-fungal tip there is bacterial or archaeal (UFBoot >= 90 at that node) the topology is HGT-consistent, if it includes other eukaryotes (animal, amoebozoan, plant ...) it is vertical-consistent.
   (b) **Composition, from the local database:** for each genus with a protein in the clade, the length-weighted GC of its coding sequence (`gene_CDS.GC_content`) against all genes of the same genome (z-score), and the
   number of introns against the genome. HGT is supported only if (a) is bacterial and in more than 25% of genomes (five times the 5% expected by chance) the GC z-score is beyond +-2; a vertical call needs (a)
   eukaryotic and the fraction beyond +-2 at or below 10%. Anything else is "inconclusive".

## Scripts (numbered; run on the HPCC from a purgeable bench folder `bench/trehalose_T5/`)

1. `01_fetch_references.py` UniProt references and context taxa, with the length of every sequence verified against UniProt.
2. `02_candidates.py` PF03632 candidates and the genome table from the local `BFD.duckdb`.
3. `03_extract_proteins.py`, `04_build_inputs.py` protein sequences from the BFD proteomes, combined with the references.
4. `05_tree.sh` mafft, trimAl, IQ-TREE (SLURM, `$SCRATCH`).
5. `06_classify.py` leave-one-out gate, classification, R1-R3.
6. `07_dollo.py` R4. 7. `08_hgt_gc.py` R5(b). Fixture tests: `test_t5.py`.

## Limits stated now

- Reviewed GH65 references exist for only five fungal trehalases, all Ascomycota, so basidiomycete and mucoromycete proteins are placed by nesting with them, not by their own references.
- The context taxa are a spread of UniProt sequences with the Pfam domain, not a complete sample of bacterial GH65.
- A gene tree places proteins; it does not measure enzyme activity. Statements about function are citations to verify (none are cited here).

---
## Results

(added after the runs)
