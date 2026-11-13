# Ku loss in Kickxellomycota: the interacting machinery and possible compensation

Date: 2026-11-11. Status: **first pass, exploratory** (MotherDuck SQL on `bfd`; scratch tables `bfd_scratch.kx2_ku_*` and
`kx2_kc_*`; nothing here has been confirmed by tblastn on the assemblies). Ku70/Ku80 is absent from Kickxellales and Harpellales
in BFD (Analyses A, B and E). Two questions follow: do the proteins that work with Ku go with it, and did anything change that could
compensate?

| Question | Directory | What it holds |
|---|---|---|
| Q1: are Ku-interacting proteins also lost, and where on the tree? | [`sql_coloss/`](sql_coloss/README.md) | 46 gene rules, presence matrix, Dollo parsimony on the first-pass backbone tree, control-based false-absence estimate, tblastn worklist |
| Q2: alternative end joining, recombination and telomere changes? | [`sql_compensation/`](sql_compensation/README.md) | gene calls and copy numbers, telomere tables, replication test across Ku-negative lineages |
| Follow-up: do other meiosis genes go with Dmc1? | [`meiosis_toolkit/`](meiosis_toolkit/README.md) | annotation rules, tblastn, MutS/MutL/kleisin/HORMA trees; Pch2 lost with Ku, Spo11/Hop2/Mnd1/Msh4 retained, Msh5 partial, Mlh3 not seen in Kickxellomycota |

Each directory has numbered SQL files in run order, a README with the gene list and identification rules, and a verification record.
The tree is the first-pass FastTree backbone in `workflows/phylo/README.md` (63 genus representatives, tips are LOCUSTAGs).

## What the data show (observed)

1. **Ku is lost once, on the stem of Kickxellales + Harpellales, in this tree.** Dimargaritales and *Ramicandelaber* keep it,
   and it is absent from all 107 *Coemansia* strains and all 11 stem genera. The relaxed single-domain rule gives 2 losses on the tree and the strict rule 5; the
   extra ones are single genera that are probably annotation misses (*Tieghemiomyces* is a 55%-duplicated genome).
   **tblastn confirms this** (see [`tblastn/`](tblastn/README.md), round 2): Ku is found in *Tieghemiomyces*, *Caulochytrium* and
   *Conidiobolus* (annotation misses), and not found in any of the 20 genomes of the 11 stem genera, with all four Ku positive
   controls recovered. *Piptocephalis* and *Neoconidiobolus* are also none and *Capillidium* weak, so further independent Ku losses
   or strong divergence are not excluded.
2. **Much of the rest of classical end joining looks gone before Ku; tblastn and a gene-family tree confirm Lig4 and Pol4.** By Pfam,
   on the stem of Entomophthorales + Kickxellomycota (18 genera) Lig4, XLF/Nej1, Pol4 and aprataxin are lost, and DNA-PKcs earlier
   still (27 genera), separately from the Dikarya-stem loss because Mucoromycota keep it.
   - **Pol4:** tblastn is supported (strong in *Rhizopus* and Zoopagales; none in all 107 *Coemansia* strains and the other
     Kickxellomycota genomes, controls passing; few queries). The Pol X tree cannot separate Pol4 from Pol lambda/mu, and no
     Kickxellomycota genome has an annotated Pol X protein.
   - **Lig4:** the "strong" tblastn hits only reflected the LIG1-like ligase core. A ligase tree ([`family_trees/`](family_trees/README.md);
     203 sequences, LG+R8, all 11 Lig4 references re-assigned correctly in leave-one-out) places **no annotated protein of any
     Kickxellomycota genome in the Lig4 clade** (Dimargaritales 0/3, Harpellales 0/4, core Kickxellales 0/14, *Mycoemilia* and
     *Spiromyces* 0/2, *Ramicandelaber* 0/1), and no Kickxellomycota ligase has the Lig4 domain or a BRCT. Annotated ligases in
     *Dispira*, *Mycoemilia*, *Zancudomyces* and *Capniomyces* place with LIG1. Lig4 is present in Mucoromycota (20 of 21), the root
     taxa (14 of 15), Zoopagales (4 of 8) and *Basidiobolus*. So Lig4 was lost on the Kickxellomycota stem, before Ku was lost on the
     stem of Kickxellales + Harpellales. The C-terminal (BRCT-region) tblastn alignments seen in non-core genomes are probably other
     BRCT-containing proteins (not checked). **Limits:** *Ramicandelaber* has no annotated ligase placed although tblastn is strong
     there (unresolved); Entomophthorales are not covered by the classification; LIG1 does not resolve as one supported clade, so
     LIG1/LIG3 calls are weak; the rescue layer (miniprot models) could not be placed (too few alignment columns), so it has not
     excluded a Lig4 that annotation missed.
   - XLF/Nej1, aprataxin and DNA-PKcs are uninformative by tblastn (controls fail, or only PIKK-like weak hits). ATR is present in
     every genome. Nothing else is co-lost with Ku itself.
3. **The other arms are retained:** MRN (Mre11, Rad50), Dna2, TERT, Est1, Ten1, PNKP, TDP1, and Rad51, Rad52, Mus81 and FEN1.
   Sae2 is lost in Kickxellales but kept in Harpellales; Rif1 is retained in core Kickxellales but absent from Harpellales.
4. **No compensation was found.** POLQ, LIG3 and PARP1-like are absent in the Ku-retaining Dimargaritales and *Ramicandelaber* as well,
   so they were not available to compensate; only one LIG1-type ligase remains. HR genes are not expanded; the Rad51 family and
   paralogs are slightly smaller. **RecA-family tree** ([`family_trees/RECA_RESULT.md`](family_trees/RECA_RESULT.md); 327 sites, leave-one-out
   19 of 31, Dmc1 5/5): the meiotic recombinase **Dmc1 is absent from the five core Kickxellales genera** (all 10 *Coemansia*
   strains in the tree, *Kickxella*, *Linderina*, *Dipsacomyces*, *Martensiomyces*; also *Spiromyces* and *Ramicandelaber* with weaker
   evidence) but present in *Mycoemilia*, Harpellales (4/4), Dimargaritales (3/3), Entomophthorales (4/4) and most Mucoromycota (18/21),
   so it is a later event than Ku loss, on the stem of core Kickxellales (the clade with the strongest intron and size reduction).
   A tblastn locus count supports one Rad51-like locus in 100 of 111 core Kickxellales genomes (the six with two are duplicated Rad51
   loci), with limited sensitivity. The earlier hypothesis of a Rad55/Rad57 loss in Kickxellales is **not supported**: they are
   absent from Mucoromycota too (0 of 21), and Rad55 is unrecoverable by the tree. Rad51 itself is present but not callable by the tree.
5. **Telomeres: nothing Ku-specific.** Harpellales show canonical CCCTAA repeats; *Coemansia* shows no detected tract, but so do most
   Mortierellomycota and Glomeromycota, and the detector only scans canonical repeats.
6. **Replication across fungi:** 6 independent Ku-negative lineages (5 usable); Lig4 is lower in 3 of 5 and never higher; no
   compensatory feature replicates. With 5 lineages, a sign test cannot do better than p about 0.06.

## Interpretation (not shown by the data)

Ku loss may be a late step in a lineage whose classical NHEJ was already dismantled, not the trigger. The resection and telomerase
arms were not lost with it. If the cell relies on the one LIG1-type ligase plus Rad51/Rad52, **testable predictions** are higher
gene-targeting efficiency and fewer random integrations in Harpellales and *Coemansia*, and non-canonical telomere repeats in
*Coemansia*. There is no transformation or telomere-sequence data in BFD to test them.

## Caveats

- **tblastn status (2026-11-11):** done for Ku, Lig4, Pol4, XLF/Nej1, aprataxin, DNA-PKcs and ATR in 31 genomes, in two rounds (round 2
  used closer-relative queries and was a post-hoc change after a control failed); see [`tblastn/`](tblastn/README.md). Sae2 and the
  low-confidence families were not checked, and weak or strong hits were not tested by reciprocal best hit. The caveat below
  still applies to everything not covered.
- **Every absence needs tblastn.** BFD annotation under-calls genes and BUSCO is 61-93% in these groups. Control genes that must be
  present are missed in 12.5-14% of calls with the multi-domain rules (2-6% with a single key domain), so single-gene absences are weak.
  Ku absence in all 107 strains and 11 stem genera, against presence in sister lineages by the same rule, is far outside random misses
  but would not show a systematic divergence. A 14-item worklist (211 gene-by-genome rows; 152 at priority 1) is in
  `sql_coloss/14_tblastn_worklist.sql`.
- **The Q1 control rules were revised after a first look** (Rad51 and Pol delta switched to family-level domain rules, Pol alpha,
  Rad54 and Rad52 demoted), so the control false-absence rate is partly tuned; this is stated in the Q1 README.
- **Low-confidence families:** XRCC4, PAXX, TDP2, Nbs1/Xrs2, ATRIP, 53BP1, CST, Pot1, Cdc13, Rap1, Taz1 and TPP1 are rare even in
  outgroups and are uninformative. ATM/ATR in Harpellales are called absent by the strict rules but a relaxed census finds
  diverged ATM-like and ATR-like proteins in 4 of 4, so they are diverged, not clearly lost. POLQ is called from architecture only.
- **The tree is first-pass.** Entomophthorales are sister to Kickxellomycota here, which affects the placement of the stem losses;
  *Linderina* and *Basidiobolus* are unstable. Parsimony is not a likelihood reconstruction.
- **Verification re-ran the same inputs** (deterministic and matching the saved SQL), not an independent re-derivation. Q1 steps
  01-04b, 07-12 and 14, and Q2's name-based classes, genus, effect-size and replication tables were not re-run.
- Telomere tables list terminal tracts only (the internal count is 0 for every genome), so interstitial repeats cannot be assessed.
- Literature statements (for example the function of Rad54's N-terminal domain, or typical telomere effects of Ku deletion) are
  citations to verify; none are cited here.

## Next

1. tblastn is done for 128 genomes (all 107 *Coemansia* strains) and the Lig4 identity question is answered by the ligase tree
   (above). Still open: the RecA-family classification of Dmc1, Rad55 and Rad57 (the first RecA tree had only 41 alignment columns
   and was not used; a re-trim is queued, job 32026340), a re-built rescue layer, Sae2 (uninformative by tblastn), a *Ramicandelaber*
   ligase check, and a long-read terminal-repeat scan for *Coemansia* telomeres.
2. Re-do the parsimony as a likelihood or Dollo reconstruction on the final tree (IQ-TREE, ASTRAL, rooted, strains dereplicated).
3. Link to the reduction analyses: the tree puts Ku loss on the stem of Kickxellales + Harpellales, before the intron and size
   reduction of core Kickxellales, but this is a placement on a first-pass tree, not a causal test.
