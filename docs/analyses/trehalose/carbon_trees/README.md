# Carbon-storage family trees (one genome per genus)

Status: **design and decision rules, written before any tree was built** (committed first; results are added after).
Parent: [`../README.md`](../README.md) (module definitions), [`../t5/`](../t5/) (the GH65 tree, same pipeline and settings).

## Why

The trehalose survey found that copy number, not presence, separates many lineages (for example a median of 7 TPS and 7 TPP copies in
every genus of Mucoromycetes and 5 in Mortierellomycetes, against 3 and 2 in most classes). A count cannot say whether extra copies are
**ancient paralogs kept by every genus** or **recent lineage-specific duplications**. A gene tree can.

## Families and sequences

One representative genome per genus (`is_genus_rep`, annotated; 1,669 genomes). All proteins carrying the family's Pfam domain
(`domain_i_evalue <= 1e-5`, version stripped), trimmed to the best-scoring domain envelope +-10 aa. A protein shorter than 60% of the
family's median trimmed length is dropped.

| Family | Pfam | Sequences (counted before design) | Method |
|---|---|---|---|
| TPS | PF00982 | 5,226 | FastTree first pass |
| TPP | PF02358 | 3,639 | FastTree first pass |
| NTH (neutral trehalase) | PF01204 | 3,885 | FastTree first pass |
| UGP (UDP-glucose pyrophosphorylase) | PF01704 | 3,866 | FastTree first pass |
| GSY (glycogen synthase) | PF05693 | 1,839 | FastTree + IQ-TREE |
| GPH (glycogen phosphorylase) | PF00343 | 1,742 | FastTree + IQ-TREE |
| GDB (glycogen debranching) | PF06202 | 1,897 | FastTree + IQ-TREE |
| GPD (glycerol-3-phosphate dehydrogenase) | PF01210 | 2,262 | FastTree + IQ-TREE |

Not built: **PGM** ([PF02878](https://www.ebi.ac.uk/interpro/entry/pfam/PF02878/), 4,564 sequences, but the domain is only about 75 aa, too short for a
reliable tree) and **GH65** (already built in `../t5/`). Accessions are linked to InterPro on first mention in the parent README.

## Pipeline (fixed before the run, as in T5)

mafft 7.505 `--auto`; trimAl 1.5.1 `-gt 0.4 -resoverlap 0.5 -seqoverlap 50`; FastTree 2.1.11 `-lg -gamma` for all families;
IQ-TREE 3.0.1 `-m Q.PFAM+R8` with 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates for the four families of at most 2,300
sequences. Large families get only the first-pass tree (no bootstrap; FastTree local support values are reported but are not
bootstrap values). Trees are unrooted: no outgroup was added, so root-dependent statements are not made.

## Question and decision rule (one pre-registered test)

**TP1.** In classes with elevated copy number (**Mucoromycetes, Mortierellomycetes**) the extra TPS and TPP copies are **ancient
paralogs retained by every genus** rather than recent duplications.

- Measure: for each genus-representative genome with at least 2 copies, the median patristic distance between its own copies; the
  reference is the median patristic distance between copies drawn from **different genera of the same class** (at most 2,000
  random pairs per class, fixed seed). **Paralog index PI** = median over genomes of (within-genome distance / between-genus
  distance). Computed per family and class.
- **Supported** (ancient paralogs): PI >= 0.8 in both classes and both families. **Killed** (recent duplications): PI < 0.5 in
  either class for either family. Otherwise **mixed**.
- Baseline for context, not part of the rule: PI for Sordariomycetes, Agaricomycetes and Eurotiomycetes (3 or 4 copies).
- A within-genome distance below the between-genus distance means copies are more alike than different genera are, the
  signature of lineage-specific duplication; a ratio near or above 1 means copies diverged before the genera did.

The same index is reported for the other families and classes with at least 8 genera, as description only.

## Limits stated in advance

- First-pass FastTree trees are approximate and unrooted; branch lengths are model-based estimates, so PI carries error that is not
  quantified here (a bootstrap of PI over genomes is reported, not over trees).
- Gene-model errors (split or merged genes) can inflate or deflate within-genome distances.
- One genome per genus does not make genera independent; the classes are not independent either.
- No outgroup, so paralog clades are described by tree topology within the fungal set only; claims about the origin of TPS and TPP
  paralogs (or about function) are citations to verify and none is made here.

## Files (planned)

`01_extract_domains.py` (HPCC job), `02_trees.sh` (SLURM array, one task per family), `03_paralog_index.py`, `test_carbon_trees.py`, `results/`.
