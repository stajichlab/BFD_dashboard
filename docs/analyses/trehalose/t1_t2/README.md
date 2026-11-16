# T1 + T2: are TPS and TPP really absent from lichen-forming fungi and *Scedosporium*?

Status: **done 2026-11-16: T1 supported (candidate loss), T2 killed by its pre-registered rule**; see Results at the end. The design and
decision rules below were written, and committed, before any tblastn result existed.
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
`03_run_tblastn.sh`, `04_summarize.py`, `05_lifestyle.py`, `test_t1t2.py` (fixture test, 7 tests pass on the HPCC), `results/`.
The raw hit tables stay in the bench folder (`bench/trehalose_T1T2/hits/`, one file per genome).

---
## Results (2026-11-16)

Order of events: this design and the selection rules were committed first; one amendment (domain-trimmed queries, TPS2-type
proteins excluded from calls) was committed before any search ran; the searches ran after. Jobs: selection 32116944, tblastn array
32117116 (33 genomes, all completed), summary and lifestyle join on the login-node queue.

### The search works

| Check | Result |
|---|---|
| TPS strong in positive controls | 10 of 10 (all full-length, 94 to 100% of the query) |
| TPP strong in positive controls | 10 of 10 (all full-length, 99 to 100% of the domain query) |
| Rad51 strong | 33 of 33 genomes |
| All 10 UniProt query lengths | equal to `bfd.swissprot_subject.s_len` |

So both assays pass and every assembly is searchable.

### Outcomes by stratum (`results/t1t2_genomes.tsv`)

| Stratum | Genomes | Informative | Candidate loss | Annotation gap | Ambiguous |
|---|---|---|---|---|---|
| `lecan_primary` | 10 | 10 | 10 | 0 | 0 |
| `scedo` | 5 | 5 | 5 | 0 | 0 |
| `other_loss` (Microascales relatives, *Pycnora*) | 4 | 4 | 4 | 0 | 0 |
| `lecan_fragmented` (N50 < 50 kb) | 4 | 0 | 0 | 1 | 3 |

In all 19 candidate-loss genomes the TPS and TPP queries produce **no alignment at all** at e-value 1e-3 (the best-hit columns are
empty), so these are not near-misses. The one annotation gap is *Flavoplaca oasis* (N50 28.6 kb): TPS is absent but TPP is
found full-length (99% of the domain query, e-value 1.8e-29) although Pfam called it absent. That fits the parent README's concern
that fragmented lichen assemblies are gap-prone, and is why the fragmented stratum is kept out of the verdict.

### Verdicts, by the rules fixed above

- **T1: supported, as a candidate loss.** Both clades pass: `lecan_primary` 10 of 10 informative genomes are candidate losses and
  none is an annotation gap; `scedo` 5 of 5. The 10 lichen genomes span three orders (Caliciales 4, Teloschistales 3,
  Peltigerales 3) and the Microascales genomes are a separate clade, so by taxonomy this is at least two clades, consistent with
  several independent losses. Independence was **not** tested (no tree).
- **T2: killed by its pre-registered rule.** 8 of the 19 candidate-loss genomes are in genera whose FungalTraits 1.2 primary
  lifestyle is a free-living saprotroph (*Scedosporium* x4, *Lomentospora*, *Scopulariopsis*, *Cephalotrichum* and
  *Basipetospora*, all Microascales; `results/t2_lifestyle_by_genus.tsv`). The other 11 are lichenized (9 genera) or unlabelled
  (*Gallowayella*, no FungalTraits entry; not counted either way).

### What this does and does not show

- A candidate loss means TPS and TPP are not recoverable by protein-to-contig search in a contiguous assembly whose controls
  pass. It does not rule out a very diverged gene or a gene in an unassembled region.
- The 10 lichen and 5 *Scedosporium* genomes were selected because Pfam found no TPS and no TPP, so these results validate those
  Pfam absences; they do not estimate how common the loss is in the class (that remains the earlier 40% of genera, now with
  a confirmed absence in the contiguous genomes tested).
- The killing of T2 rests on one clade: all eight saprotroph genomes are Microascales, so it is one lineage (probably one loss),
  not eight events. The rule used the FungalTraits primary lifestyle at genus level; whether such genera have a sugar-rich niche
  was not judged. The genome-level labels were not independently checked (citation to verify).
- Absence of TPS and TPP says nothing about whether these fungi contain trehalose or obtain it another way; alternative routes
  were not searched.

### Not done

No profile HMM from six-frame translations (no mixed result, so not triggered); no tree of independent losses; the four
fragmented lichen genomes are not interpreted beyond the single annotation gap.
