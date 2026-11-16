# T1 + T2: are TPS and TPP really absent from lichen-forming fungi and *Scedosporium*?

Status: **design and decision rules, written before any tblastn result exists** (commit order: this file first, results later).
Parent analysis: [`../README.md`](../README.md) sections 4.5 and 7 (T1, T2).

## Question

BFD has no TPS (trehalose-6-phosphate synthase, [PF00982](https://www.ebi.ac.uk/interpro/entry/pfam/PF00982/)) or TPP
(trehalose-6-phosphate phosphatase, [PF02358](https://www.ebi.ac.uk/interpro/entry/pfam/PF02358/)) call for about 40% of
Lecanoromycetes genera and for *Scedosporium*/*Lomentospora*. Pfam and annotation can both miss a gene, so the assemblies are
searched directly. Three explanations were listed in the parent README: (a) the gene is in the assembly but was never predicted,
(b) the gene is present but too diverged for Pfam and the SwissProt-name check, (c) the pathway is really lost.

- **T1:** TPS and TPP are absent from the assemblies of Lecanoromycetes and of the *Scedosporium* clade (two clades, so at least
  two losses).
- **T2:** if so, losses should be tolerated only where a partner supplies carbon (lichen-forming, lichenicolous, mycoparasitic),
  not in free-living saprotrophs.

## Design

- **Search:** `tblastn` of protein queries against the unmasked assembly (`input_clean_genomes/<ASMID>.fa.gz`), BLAST+ 2.17.0,
  e-value 1e-3, up to 50 targets, default masking. `tblastn` already searches the six translated frames of the contigs, so no
  separate six-frame search is needed; a hit is a protein-to-contig alignment. Same machinery and thresholds as
  [`../../kx2_ku_pathway/tblastn/`](../../kx2_ku_pathway/tblastn/README.md).
- **Queries** (`queries.tsv`, fixed before the run): reviewed fungal SwissProt TPS1 proteins (`TPS1A_ASPFU`, `TPS1_BOTFB`,
  `TPS1_EMENI`, `TPS1_CANAL`, `TPS1_YEAST`, `TPS1_SCHPO`), each length-checked against `bfd.swissprot_subject`; for closer relatives,
  the best PF00982 and the best PF02358 protein of five Pfam-positive control genomes (2 Lecanoromycetes, 2 Sordariomycetes,
  1 Eurotiomycetes), **trimmed to the Pfam domain envelope (+-5 aa)**; Rad51 as a does-the-search-work control. Calls: **TPS** =
  the TPS1 references plus the PF00982-trimmed BFD queries; **TPP** = the PF02358-trimmed BFD queries only. The TPS2-type
  references (`TPS2_YEAST`, `TPS2_SCHPO`, `TPS2_ZYGRO`) carry both domains, so a true TPS gene would hit them over half their
  length and look like a TPP hit; they are searched under the key `TPS2like` for context and are **not** used in any call.
  *Amendment, made before any tblastn was run:* the first draft of this README used the TPS2-type proteins as TPP queries and
  untrimmed BFD proteins; that was changed once the cross-hit problem was seen. A query that comes from the genome being
  searched is dropped (leave-one-out).
- **Genome strata** (`genomes.tsv`, chosen by the rules in `01_select_genomes.py`, not by looking at hits):
  1. `lecan_primary`: 10 Lecanoromycetes with no Pfam TPS and no Pfam TPP, N50 >= 50 kb, best BUSCO, at most 2 per genus.
  2. `lecan_fragmented`: up to 4 more with N50 < 50 kb (weak evidence by design; reported separately).
  3. `scedo`: 5 *Scedosporium*/*Lomentospora*, one per species, best BUSCO then N50.
  4. `other_loss`: up to 4 other genomes of the credible-loss list outside Lecanoromycetes and *Scedosporium*, one per genus,
     best BUSCO (the rule picked relatives of *Scedosporium* in Microascales, which sit in the same clade as the `scedo` stratum;
     these are reported separately and are not counted in the clade verdict).
  5. `control_pos`: 10 Pfam-positive genomes (5 Lecanoromycetes of different genera, 3 Sordariomycetes close to *Scedosporium*,
     2 Eurotiomycetes), N50 >= 50 kb.

## Calls (thresholds as in the Ku check, fixed in advance)

Per genome and gene key (TPS, TPP, Rad51), the best query is the one with the highest merged query coverage x bit score.

- **strong:** best e-value <= 1e-10 and merged HSPs cover >= 30% of the query on one scaffold. **full-length:** strong and
  coverage >= 70%. **weak:** e-value <= 1e-5 but not strong. **none:** otherwise.

An **informative** genome has Rad51 strong in that genome and N50 >= 50 kb. The **assay passes** for a gene only if that gene is
strong in at least 8 of the 10 `control_pos` genomes; otherwise every negative for that gene is uninformative and is reported
as such.

Per genome (both TPS and TPP absent by Pfam):

- **annotation gap:** TPS or TPP is full-length strong. The pathway is in the assembly and was not annotated.
- **candidate loss:** informative genome, assay passes, and both TPS and TPP are none.
- **ambiguous:** anything else (weak hits, strong with 30 to 70% coverage, non-informative genome).

## Pre-registered verdicts

Clade-level, using informative genomes of `lecan_primary` and of `scedo` separately:

- **T1 supported for a clade** if at least half of its informative genomes are candidate losses and none of the rest is a
  full-length strong hit. **T1 killed for that clade** if at least half show an annotation gap. Otherwise **mixed**.
- **T1 overall** is supported only if both clades are supported. If one clade is killed, T1 as written ("lost at least twice") is
  not supported. A tblastn "candidate loss" is not proof: a diverged gene or a missing assembly region could still hide it.
- **T2** is evaluated only on genomes called candidate loss. It is **killed** if a candidate loss falls in a genome whose
  FungalTraits primary lifestyle is a free-living saprotroph; the lifestyle label and its source are written next to each genome.
  The number of independent candidate-loss lineages is small, so T2 is at best consistent, never confirmed, by this test.

## Not done here

- No profile HMM from six-frame lichen translations (the parent README's second step); only if the result is mixed.
- No phylogeny of independent losses; clades are counted by taxonomy (two clades: Lecanoromycetes, *Scedosporium*).
- No wet-lab validation. Any statement about the biology of trehalose use is a citation to verify.

## Files

`01_select_genomes.py` (HPCC job, writes `genomes.tsv`, `queries.tsv`, `selection_log.txt`), `02_fetch_queries.sh`,
`03_run_tblastn.sh`, `04_summarize.py`, `test_t1t2.py` (fixture test), `results/`.
