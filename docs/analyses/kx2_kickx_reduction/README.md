# Kickxellomycota reduction: four hypotheses against the yeast lineages (Analysis E)

Date: 2026-11-10. Status: **first pass, exploratory** (MotherDuck SQL on `bfd`; scratch tables `bfd_scratch.kx2_h1_*` to
`kx2_h4_*`). The summary, findings and caveats are in section E of [../../FUNCTIONAL_ANALYSES.md](../../FUNCTIONAL_ANALYSES.md);
this directory holds the SQL that produced them. A narrative report with figures has not been written yet.

| Directory | Hypothesis | Outcome (see section E) |
|---|---|---|
| [`sql_h1/`](sql_h1/README.md) | Multidomain / architecture simplification (domains, bigrams, fusion enzymes, protein-length decomposition) | not supported |
| [`sql_h2/`](sql_h2/README.md) | Genome-size accounting: gene number vs intron vs intergenic compaction (UTR and overlap not measurable) | not supported |
| [`sql_h3/`](sql_h3/README.md) | Selection vs drift proxies (ENC vs GC3, tRNA-codon match, amino acids, ortholog identity) | weak lean |
| [`sql_h4/`](sql_h4/README.md) | Loss composition (complete loss vs contraction), scaling residuals, expansions | mixed |
| [`intron_crosscheck/`](intron_crosscheck/README.md) | Follow-up: are the low intron counts a gene-model artifact? miniprot alignments of 1,122 BUSCO proteins on 83 genomes | deficit confirmed (pre-registered rule met) |

Each directory has numbered SQL files in run order and a README with the dependency graph and a verification record. The verification
re-ran later steps from the existing scratch tables and matched the reported numbers; the earlier steps are saved but were **not**
re-run, and no end-to-end rebuild from empty scratch has been done (see "Reproducibility" in section E). H2 and H3 start from the
Analysis A scratch table `kx_genomes`, built by SQL in [../kickxellomycota/](../kickxellomycota/README.md).

The follow-up on Ku loss and the genes that work with it is in [../kx2_ku_pathway/](../kx2_ku_pathway/README.md).
