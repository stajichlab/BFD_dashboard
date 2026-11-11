# kx2_h1: H1 test (multidomain proteins / protein length, Kickxellomycota vs yeasts)

SQL used for the H1 analysis, run in MotherDuck with default database `bfd_scratch` (reads `bfd.*`, writes only `bfd_scratch.kx2_h1_*`). Run in numeric order; each file was copied from the final working query that was run, with only comment headers added.

## Run order and dependencies

```
kw_rep_genomes
  01 kx2_h1_genomes
       02 kx2_h1_prot   (also reads bfd.gene_proteins, pfam, gene_transcripts, idp_summary)
            03 kx2_h1_genome_metrics
            04 kx2_h1_arch_metrics
                 05 kx2_h1_genome_all  (genomes + 03 + 04)
                      06 kx2_h1_genus
                           07 kx2_h1_cliff      (vs Mucorales)
                           08 kx2_h1_cliff_ref  (core vs unreduced Zoopagomycota/Kickxellomycota)
            09 kx2_h1_fam_gen (+ bfd.swissprot_best)
                 10 kx2_h1_decomp
            11 kx2_h1_fusion
            12 kx2_h1_longfam
            13 artefact checks (result queries)
            14 fusion marker discovery (exploratory result query)
```

## Manual or non-SQL steps

- Everything was run through the MotherDuck MCP query tool; no Python. Results in the report were read from the printed output of the result queries in each file; they were not saved to files. The Pfam name lookup in the comment of file 14 was run by hand and is not saved as its own file.
- Medians over genera per group (report table) are the second query in 06; the Cliff table is the second query in 07.

## Verification (2026-11-11)

Steps 03 to 07 were re-run into `kx2_h1_chk_*` tables from the existing `kx2_h1_prot`.
- Row counts identical: genome_metrics 740, arch_metrics 727, genome_all 740, genus 251, cliff 300.
- 0 differing rows in cliff (all 300 deltas identical); 0 differing genus rows (mean_len, f_multi, bigram_per1000, promisc).
- Core Kickxellales filtered medians: mean_len 455.1, multi 0.2818, bigrams 591.5, archs 818.1; Mucorales 447.1, 0.332, 633.5, 794.1. Cliff vs Mucorales for core: mean_len 0.09, f_multi -0.58, bigram -0.75, promisc -0.87. All identical.
- The chk tables were dropped; none remain.
- Rarefaction is deterministic (`ORDER BY hash(protein_id)`, no random seed), which is why it reproduces.

Not re-run: steps 01, 02 and 09 to 12 (they were not re-executed to avoid rebuilding the 5M-row table), so only steps 03 to 07 are verified; 02 and 09 to 12 stand on their original single run.

## Caveats and things that are not exactly as described in the report

- `arch_per1000` and `bigram_per1000` are per 1000 **Pfam-annotated** proteins in the 2000-protein sample (ns 53 to 1867, median 1571), not per 1000 sampled proteins. The report's wording ("per 1000 proteins, rarefied to 2000") should read this way. Groups with lower Pfam hit rates (Microsporidia) have lower ns, so their diversity is computed on fewer proteins.
- Promiscuity counts only domains that have at least one adjacent partner.
- The `IDP_residues` column in `kx2_h1_prot` is all NULL (id suffix mismatch); `idp_res` is the working column. IDP coverage is 2.56M of 5.0M proteins, so disorder medians are indicative only.
- In `kx2_h1_fam_gen`, shares for Microsporidia sum to 0.913 because some q genomes have no rows in `kx2_h1_prot`; the decomposition for Microsporidia is therefore approximate. Lost/gained families depend on SwissProt hit thresholds (e-value 1e-10, 30% identity) and are indicative only.
- No Pfam cutoff sensitivity (1e-3, 1e-10) was run.
- Pfam accessions used as fusion markers come from the architectures of SwissProt-anchored proteins; names were looked up with the one-off query described in file 14.
