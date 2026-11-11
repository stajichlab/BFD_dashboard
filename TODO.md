[ ] - Setup example dataset to demonstrate pipeline run
[ ] - Setup example access to BFD database
[ ] - When a BFD db is accessible remotely setup access to it for sub-set query

## Upstream follow-ups in BFD/Fungi_BFD_runs (do when HPCC is back online)
Noted 2026-11-10 while HPCC was in maintenance (`20261110-Shutdown`). Details: docs/FUNCTIONAL_ANALYSES.md, "Upstream follow-ups".
[ ] - Gene annotation under-calls in MAGs: 33 of 44 near-universal Pfam families missing from all 6 annotated Antarctolithica obscura proteomes are present in the assemblies by blastx (incl. ribosomal proteins L22e, L34e, S21e). Find out why the annotation pipeline missed them (gene-model prediction, training, filtering, contig-edge or short-exon genes?), and whether the same happens in other small or fragmented genomes. Also find out why Mars-9.1 and Mars-19.38 were skipped as too small.
[ ] - Fix gene_introns.codon_position / codon_frame (off by one exon: they describe the next CDS junction) in the table builder, then rebuild gene_introns.
[ ] - Check the phylum of "Saccharomycotina sp. YJ-2025a/b" (FF90B97F, F6EB82EF); they group with Basidiomycota in phylogenomics. Excluded from the current analyses.
[ ] - Write a NOTE file with these items into /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/ once HPCC is reachable.
