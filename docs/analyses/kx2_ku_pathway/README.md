# Ku loss in Kickxellomycota: the interacting machinery and possible compensation

Date: 2026-11-11. Status: **first pass, exploratory** (MotherDuck SQL on `bfd`; scratch tables `bfd_scratch.kx2_ku_*` and
`kx2_kc_*`; nothing here has been confirmed by tblastn on the assemblies). Ku70/Ku80 is absent from Kickxellales and Harpellales
in BFD (Analyses A, B and E). Two questions follow: do the proteins that work with Ku go with it, and did anything change that could
compensate?

| Question | Directory | What it holds |
|---|---|---|
| Q1: are Ku-interacting proteins also lost, and where on the tree? | [`sql_coloss/`](sql_coloss/README.md) | 46 gene rules, presence matrix, Dollo parsimony on the first-pass backbone tree, control-based false-absence estimate, tblastn worklist |
| Q2: alternative end joining, recombination and telomere changes? | [`sql_compensation/`](sql_compensation/README.md) | gene calls and copy numbers, telomere tables, replication test across Ku-negative lineages |

Each directory has numbered SQL files in run order, a README with the gene list and identification rules, and a verification record.
The tree is the first-pass FastTree backbone in `workflows/phylo/README.md` (63 genus representatives, tips are LOCUSTAGs).

## What the data show (observed)

1. **Ku is lost once, on the stem of Kickxellales + Harpellales, in this tree.** Dimargaritales and *Ramicandelaber* keep it,
   and it is absent from all 107 *Coemansia* strains and all 11 stem genera. The relaxed single-domain rule gives 2 losses on the tree and the strict rule 5; the
   extra ones are single genera that are probably annotation misses (*Tieghemiomyces* is a 55%-duplicated genome).
2. **The rest of classical end joining was already gone before Ku.** On the stem of Entomophthorales + Kickxellomycota (18 genera)
   Lig4, XLF/Nej1, Pol4 and aprataxin are lost; DNA-PKcs is lost earlier still (27 genera, with Basidiobolus and Zoopagales), a loss
   separate from the Dikarya-stem loss because Mucoromycota keep it. Nothing else is co-lost with Ku itself.
3. **The other arms are retained:** MRN (Mre11, Rad50), Dna2, TERT, Est1, Ten1, PNKP, TDP1, and Rad51, Rad52, Mus81 and FEN1.
   Sae2 is lost in Kickxellales but kept in Harpellales; Rif1 is retained in core Kickxellales but absent from Harpellales.
4. **No compensation was found.** POLQ, LIG3 and PARP1-like are absent in the Ku-retaining Dimargaritales and *Ramicandelaber* as well,
   so they were not available to compensate; only one LIG1-type ligase remains. HR genes are not expanded; the Rad51 family and
   paralogs are slightly smaller. A possible Kickxellales-specific loss of Dmc1 (and Rad55/Rad57) is a later event than Ku loss.
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

1. Run the tblastn worklist on the HPCC (priority 1 first; Ku against all stem genera and *Coemansia* strains, with Ku in
   *Ramicandelaber*, *Dimargaris* and *Rhizopus* as positive controls). A long-read terminal-repeat scan for *Coemansia* telomeres.
2. Re-do the parsimony as a likelihood or Dollo reconstruction on the final tree (IQ-TREE, ASTRAL, rooted, strains dereplicated).
3. Link to the reduction analyses: the tree puts Ku loss on the stem of Kickxellales + Harpellales, before the intron and size
   reduction of core Kickxellales, but this is a placement on a first-pass tree, not a causal test.
