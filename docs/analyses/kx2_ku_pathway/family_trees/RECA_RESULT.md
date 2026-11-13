# How the RecA/Rad51 family is organized in Kickxellomycota (re-trimmed tree)

Date: 2026-11-11. Status: **first pass, annotated layer only.** Tree: 362 sequences (31 references, the rest annotated BFD proteins
with RecA-fold domains from 73 genomes), **327 alignment sites** after trimming, Q.PFAM+R8, IQ-TREE; the first RecA tree (41 sites)
was discarded. Files: `results/annotated/reca.*` (`reca_loo.tsv`, `reca_by_group.tsv`, `reca_by_genome.tsv`, `reca.treefile`).
References: Rad51 (7), Dmc1 (5), Rad55 (2), Rad57 (2), Rad51B/C/D, XRCC2/3, and archaeal RadA as outgroup (`../references.tsv`).

## Which clades the tree can recover (leave-one-out, 19 of 31 references)

| Clade | references re-assigned correctly | Use |
|---|---|---|
| Dmc1 | 5 of 5 | **usable**: an absence call is meaningful |
| Rad57 | 2 of 2 | usable (yeast references only) |
| Rad51B, Rad51C, Rad51D, XRCC2, RadA | all, except one Rad51D at low support | usable, but see caveat on fungal paralogs |
| XRCC3 | 0 of 2 correct (both assigned to XRCC3 at low support) | weak |
| Rad55 | 0 of 2 | **not recoverable**: no Rad55 call can be used |
| Rad51 | 0 of 7 pass the strict rule (2 low support, 5 ambiguous) | **the tree cannot call Rad51**; Rad51 is known present from annotation and tblastn (Rad51 control strong in 126 of 128 genomes) |

## Organization (observed)

1. **Number of RecA-fold proteins per genome** (annotated, tree-independent): *Coemansia* 4 (3-5), other Kickxellales 4, Harpellales 4
   (3-5), Dimargaritales + *Ramicandelaber* 5.5, Zoopagales 4, Entomophthorales + *Basidiobolus* 4.5, Mucoromycota 6 (3-8), root taxa 5.
2. **Dmc1** (the meiosis-specific recombinase) is present in Dimargaritales 3/3, Harpellales 4/4, Entomophthorales + *Basidiobolus*
   4/4, Mucoromycota 18/21 and root taxa 12/15, and in *Mycoemilia*. It is **absent from every genome of the five core
   Kickxellales genera** (*Coemansia* all 10 strains in the tree, *Kickxella*, *Linderina*, *Dipsacomyces*, *Martensiomyces*) and from *Spiromyces*
   (BUSCO 44%) and *Ramicandelaber* (one genome, weaker evidence).
3. **Rad55 and Rad57 are not found in any group**, including Mucoromycota (0 of 21) and root taxa, so their absence in Kickxellales
   is **not** a Kickxellales-specific loss. The references are only yeast proteins and Rad55 is unrecoverable in leave-one-out, so
   they may simply not be present as such outside Ascomycota; the earlier hypothesis of a Rad55/Rad57 loss in Kickxellales
   (Q2 report) is not supported.
4. **Rad51B/C/D, XRCC2/3-like proteins:** many calls are "uncertain" (for example Rad51C-like 13 of 14 core Kickxellales and 9 of 21
   Mucoromycota). Animal and plant paralog references do not map cleanly onto fungal paralogs, so these are better read as
   "other RecA-fold paralogs" than as named members.

## Independent check with tblastn loci (does not use the annotation)

`../../tblastn/reca_loci.py`, output `../../tblastn/results/round3_reca_loci.txt`: distinct loci hit by Rad51 or Dmc1 queries
(e-value <= 1e-25, >= 150 aa, HSPs merged within 20 kb). 100 of 111 core Kickxellales genomes have exactly one locus; Dimargaritales 3 of 3
have two. In genomes that really have a second, Dmc1-like gene, that locus scores higher with the Dmc1 queries than with the Rad51 queries
(*Furculomyces* 238 vs 291 bits, *Capillidium* 294 vs 337). The 6 core Kickxellales genomes with two loci have **identical Rad51-like scores at both
loci** (a duplicated Rad51 locus in duplicated assemblies, as in *Tieghemiomyces*), and none has a Dmc1-like second locus.
**Limits:** the count under-detects (2 of 4 Harpellales show one locus although the tree places a Dmc1 in all four), so one locus is
consistent with Dmc1 loss but does not prove it.

## Interpretation (not shown by the data)

Dmc1 appears to have been lost once, on the stem of the five core Kickxellales genera, after Ku was lost on the stem of
Kickxellales + Harpellales (Harpellales keep Dmc1; *Mycoemilia* keeps it). This is the same clade that shows the strongest intron
and genome-size reduction; Sae2 is also called lost in Kickxellales by Pfam (tblastn is uninformative for Sae2, so that call is
unconfirmed). A loss of the meiotic recombinase would predict an eroded meiotic program in core Kickxellales (the order is
described as reproducing sexually by zygospores, citation to verify; any link is a hypothesis).

## Caveats

- Annotated layer only; the rescued models were not usable at tree level. A Dmc1 missed by annotation in every core genome is
  possible but unlikely given the tblastn locus pattern and the earlier name-based result (1 of 101 genomes).
- Only 10 of the 107 *Coemansia* strains are in the tree (the 97 others were not run); the Dmc1 call for them rests on tblastn
  loci (one locus in nearly all).
- Single ML tree; Rad51 is not callable, and Rad55 and XRCC3 are weak.

## Next

Done 2026-11-12: other meiosis genes were tested in [`../meiosis_toolkit/`](../meiosis_toolkit/README.md) (Pch2 lost earlier with Ku;
Spo11, Hop2, Mnd1, Msh4 retained; Msh5 partial). Still open: run the RecA classification for all 107 *Coemansia* strains, and check
*Ramicandelaber* and *Spiromyces* with higher-quality assemblies.
