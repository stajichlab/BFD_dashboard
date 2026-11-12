# Is the funannotate annotation of *Coemansia* complete? Genome-mode vs protein-mode BUSCO

Date: 2026-11-11. Status: **first pass, exploratory** (MotherDuck SQL; one HPCC read of existing BUSCO summaries).

**Question.** The tblastn checks (see [../../kx2_ku_pathway/](../../kx2_ku_pathway/README.md)) agreed with the annotation for every
*Coemansia* gene where there was a definite answer, but the assemblies are fragmented, so is the annotation losing genes the
assemblies contain? BFD already has the data to test this at the level of conserved genes: BUSCO was run in genome mode (miniprot on
the assembly) and in protein mode (on the annotated proteome), both with `fungi_odb12` (BUSCO 6.1.0).

**Data.** Protein-mode summaries are `Fungi_BFD_runs/results/genome_stats/BUSCO_protein/<shard>/<LOCUSTAG>.BUSCO_summary.fungi_odb12.txt`
(read only; `protein_busco.py` collects them for the 63 backbone genera plus the 107 *Coemansia* strains, 169 genomes, none missing;
`protein_busco_169.tsv`). Genome-mode values are `bfd.busco_genome.complete_pct`. The query is `busco_gap.sql`.

**Result (observed).** Gap = genome-mode minus protein-mode complete %.

| Group | n | genome BUSCO | protein BUSCO | median gap | range | n with gap > 10 |
|---|---|---|---|---|---|---|
| *Coemansia* (strains) | 107 | 80.4 | 82.4 | -0.2 | -3.2 to +3.9 | 0 |
| other Kickxellales | 6 | 87.6 | 77.5 | +1.2 | -5.7 to +16.9 | 1 |
| Harpellales | 4 | 79.8 | 85.0 | -5.5 | -6.0 to -4.0 | 0 |
| Dimargaritales + *Ramicandelaber* | 4 | 82.1 | 85.7 | -3.7 | -4.9 to -1.7 | 0 |
| Zoopagales | 8 | 74.1 | 79.3 | -4.7 | -6.3 to +1.9 | 0 |
| Entomophthorales + *Basidiobolus* | 4 | 78.6 | 90.2 | -9.8 | -12.9 to -0.7 | 0 |
| Mucoromycota group | 21 | 95.9 | 95.4 | +0.7 | -4.9 to +27.2 | 1 |
| root (chytrids etc.) | 15 | 83.8 | 87.5 | -1.2 | -7.7 to +14.0 | 1 |

For *Coemansia* the annotation recovers essentially every BUSCO gene the assembly contains (median gap -0.2 points; no strain
above 10). Negative gaps mean protein-mode BUSCO is higher than genome-mode, which happens when genome-mode alignment misses genes
that the gene finder recovers (largest in the large Entomophthorales genomes).

**What this does and does not show (interpretation).**
- It argues against a systematic funannotate failure in *Coemansia* for conserved single-copy genes.
- It cannot see genes absent from the assembly (fragmentation), species-specific or fast-evolving genes, or intron mis-calls:
  BUSCO checks that a gene is found, not that its exon structure is right. The intron question is tested separately in
  [../intron_crosscheck/](../intron_crosscheck/README.md).
- Both BUSCO modes use the same marker set, so a marker lost to assembly fragmentation is missing from both and does not show as a gap.
- Context: the BFD team's own "BUSCO gap" set (52 genomes with genome >= 90% but protein < 50%, mostly TE-rich ectomycorrhizal
  Agaricomycetes) is documented in `Fungi_BFD/todo/tf_zero_hit_genomes.md`. No *Coemansia* strain looks like those.
