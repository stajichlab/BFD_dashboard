# BFD local DuckDB: structure review, MotherDuck compatibility, proposed views

Date: 2026-11-11. Author: Claude agent (read-only review).
Folder: `/rhome/jstajich/.claude/jobs/3682abe7/tmp/viewreview/`

Files in this folder:

| File | Content |
|---|---|
| `compat_views.sql` | 9 views for the missing MotherDuck names, 42 `COMMENT ON` statements, 1 commented-out decision block (`asm_stats`) |
| `extra_views.sql` | 10 helper views and 1 macro (`in_genome`). Needs `compat_views.sql` first (`samples`) |
| `extra_views_atlas.sql` | Optional. 4 views over `BFD_dashboard/atlas_data` parquet. Decision needed (section 7) |
| `REPORT.md` | This report |
| `gen_compat.py` | Regenerates `compat_views.sql`; it reads the descriptions from `md_load_bfd.py` |
| `test/` | Test scripts and logs |
| `work/` | Helper scripts and the raw structure dump (`struct.txt`) |

Scope of the checks. I read the published database `Fungi_BFD_runs/BFD_function/db/BFD.duckdb` (2026-11-08) read-only. I wrote only in this folder. I submitted no job and changed nothing under `/bigdata`. The MotherDuck MCP server did not connect in this session, so I did not query MotherDuck. The MotherDuck table list comes from the load log `BFD_dashboard/md_load_bfd.31981883.log` (2026-11-09, 36 tables). The MotherDuck columns are taken to be the parquet columns, because the loader runs `CREATE OR REPLACE TABLE ... AS SELECT * FROM read_parquet(...)` (`bin/md_load_bfd.py:190`). That equality is **not verified** against MotherDuck itself.

No file in `tables/` had a modification time after 2026-11-08 22:13 when I listed it (2026-11-11 about 19:55). I saw no sign of the merge job in the files that I read.

---

## 1. Structure

### 1.1 Database objects (published BFD.duckdb, 2026-11-08)

38 views, 1 base table (`swissprot_subject`), 3 macros (`swissprot_hits_file`, `swissprot_hits_for`, `swissprot_model_fungi_for`). No view had a comment (`duckdb_views().comment` is NULL for all).

All view paths point to `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/BFD_function/tables/...`. `BFD_function/tables` is a symlink to `Fungi_BFD_runs/tables`. The dashboard wrapper `workflows/build_bfd_db.sh` uses `Fungi_BFD_runs/tables` directly. My SQL files use `Fungi_BFD_runs/tables`.

Row counts below are from `docs/BFD_SCHEMA.md` (2026-11-08). I recounted only the tables used in the tests; those counts agree.

| Object | Type | Source | Rows | Key columns | Columns (type) |
|---|---|---|---|---|---|
| `_table_manifest` | VIEW | `_manifest/*.json` | not counted |  | table (VARCHAR), path (VARCHAR), row_count (BIGINT), built_at (TIMESTAMP), built_by (VARCHAR), merge_run_id (VARCHAR), schema_version (BIGINT), columns (VARCHAR[]), parquet_sha256 (VARCHAR), genomes (BIGINT), interim (BOOLEAN), note (VARCHAR), estimator (VARCHAR), code (VARCHAR), description (VARCHAR) |
| `aa_frequency` | VIEW | `aa_freq.parquet` | 477,535 | species_prefix | species_prefix (VARCHAR), amino_acid (VARCHAR), frequency (DOUBLE) |
| `asm_stats` | VIEW | `asm_stats.parquet` (joined to species.parquet for LOCUSTAG; columns cast/renamed) | 23,584 | LOCUSTAG, ASMID | LOCUSTAG (VARCHAR), ASMID (VARCHAR), contig_count (BIGINT), TOTAL_LENGTH (BIGINT), min_contig_bp (BIGINT), max_contig_bp (BIGINT), median_contig_bp (DOUBLE), mean_contig_bp (DOUBLE), L50 (BIGINT), N50_bp (BIGINT), L90 (BIGINT), N90_bp (BIGINT), GC_PERCENT (DOUBLE), n_gap_count (BIGINT), total_n_bases (BIGINT), masked_bases (BIGINT), masked_pct (DOUBLE), t2t_scaffolds (BIGINT), telomere_fwd (BIGINT), telomere_rev (BIGINT), annotation_status (VARCHAR) |
| `busco_genome` | VIEW | `busco_genome.parquet` (joined to species.parquet for LOCUSTAG; columns cast/renamed) | 23,582 | LOCUSTAG, ASMID | LOCUSTAG (VARCHAR), ASMID (VARCHAR), complete_pct (DOUBLE), single_pct (DOUBLE), duplicated_pct (DOUBLE), fragmented_pct (DOUBLE), missing_pct (DOUBLE), n_markers (BIGINT), lineage (VARCHAR) |
| `cazy` | VIEW | `cazy.cazymes_hmm.parquet` | 9,164,068 | species_prefix, protein_id | species_prefix (VARCHAR), HMM_id (VARCHAR), profile_length (BIGINT), protein_id (VARCHAR), protein_length (BIGINT), evalue (DOUBLE), q_start (BIGINT), q_end (BIGINT), s_start (BIGINT), s_end (BIGINT), coverage (DOUBLE) |
| `cazy_overview` | VIEW | `cazy.overview.parquet` | 9,310,657 | species_prefix, protein_id | species_prefix (VARCHAR), protein_id (VARCHAR), EC (VARCHAR), cazyme_fam (VARCHAR), sub_fam (VARCHAR), diamond_fam (VARCHAR), substrate (VARCHAR), toolcount (BIGINT) |
| `codon_frequency` | VIEW | `codon_freq.parquet` | 1,500,925 | species_prefix | species_prefix (VARCHAR), codon (VARCHAR), frequency (DOUBLE) |
| `gene_CDS` | VIEW | `gene_CDS.parquet` (+ computed species_prefix) | 753,622,745 | species_prefix, transcript_id | cds_id (VARCHAR), transcript_id (VARCHAR), order (BIGINT), chrom (VARCHAR), start (BIGINT), end (BIGINT), strand (BIGINT), GC_content (DOUBLE), species_prefix (VARCHAR) |
| `gene_exons` | VIEW | `gene_exons.parquet` (+ computed species_prefix) | 758,181,061 | species_prefix, transcript_id | exon_id (VARCHAR), transcript_id (VARCHAR), order (BIGINT), chrom (VARCHAR), start (BIGINT), end (BIGINT), strand (BIGINT), GC_content (DOUBLE), species_prefix (VARCHAR) |
| `gene_info` | VIEW | `gene_info.parquet` | 249,546,610 | LOCUSTAG, gene_id | gene_id (VARCHAR), LOCUSTAG (VARCHAR), chrom (VARCHAR), start (BIGINT), end (BIGINT), strand (BIGINT), gene_type (VARCHAR) |
| `gene_intergenic_distances` | VIEW | `gene_intergenic_distances.parquet` | 233,406,624 | species_prefix | species_prefix (VARCHAR), left_gene (VARCHAR), right_gene (VARCHAR), distance (BIGINT) |
| `gene_introns` | VIEW | `gene_introns.parquet` (+ computed species_prefix) | 508,634,451 | species_prefix, transcript_id | intron_id (VARCHAR), transcript_id (VARCHAR), intron_number (BIGINT), chrom (VARCHAR), start (BIGINT), end (BIGINT), strand (BIGINT), splice_5 (VARCHAR), splice_3 (VARCHAR), GC_content (DOUBLE), codon_position (BIGINT), codon_frame (BIGINT), species_prefix (VARCHAR) |
| `gene_proteins` | VIEW | `gene_proteins.parquet` (+ computed species_prefix) | 244,988,294 | species_prefix, gene_id, transcript_id, protein_id | gene_id (VARCHAR), protein_id (VARCHAR), transcript_id (VARCHAR), length (BIGINT), md5checksum (VARCHAR), species_prefix (VARCHAR) |
| `gene_transcripts` | VIEW | `gene_transcripts.parquet` (+ computed species_prefix) | 249,546,610 | species_prefix, gene_id, transcript_id | gene_id (VARCHAR), transcript_id (VARCHAR), chrom (VARCHAR), start (BIGINT), end (BIGINT), strand (BIGINT), is_partial (BOOLEAN), has_start_codon (VARCHAR), has_stop_codon (VARCHAR), CDS_length (BIGINT), md5checksum (VARCHAR), species_prefix (VARCHAR) |
| `gene_trna` | VIEW | `gene_trnas.parquet` (+ computed species_prefix) | 4,558,316 | species_prefix, gene_id | gene_id (VARCHAR), amino_acid (VARCHAR), codon (VARCHAR), species_prefix (VARCHAR) |
| `idp` | VIEW | `idp.parquet` | 233,329,228 | species_prefix, protein_id | species_prefix (VARCHAR), protein_id (VARCHAR), IDP_start (BIGINT), IDP_end (BIGINT), IDP_length (BIGINT), mean_score (DOUBLE) |
| `idp_summary` | VIEW | `idp_summary.parquet` | 144,924,594 | species_prefix, protein_id | species_prefix (VARCHAR), protein_id (VARCHAR), IDP_residues (BIGINT), IDP_fraction (DOUBLE), length (BIGINT) |
| `merops` | VIEW | `merops.parquet` | 43,462,876 | species_prefix, protein_id | species_prefix (VARCHAR), protein_id (VARCHAR), merops_id (VARCHAR), percent_identity (DOUBLE), aln_length (BIGINT), mismatches (BIGINT), gap_openings (BIGINT), q_start (BIGINT), q_end (BIGINT), s_start (BIGINT), s_end (BIGINT), evalue (DOUBLE), bitscore (DOUBLE) |
| `pfam` | VIEW | `pfam.parquet` (+ computed species_prefix) | 368,209,721 | species_prefix, protein_id | protein_id (VARCHAR), pfam_id (VARCHAR), pfam_acc (VARCHAR), pfam_len (BIGINT), full_seq_e_value (DOUBLE), full_seq_score (DOUBLE), full_seq_bias (DOUBLE), domain_num (BIGINT), domain_num_of (BIGINT), domain_c_evalue (DOUBLE), domain_i_evalue (DOUBLE), domain_score (DOUBLE), domain_bias (DOUBLE), hmm_from (BIGINT), hmm_to (BIGINT), ali_from (BIGINT), ali_to (BIGINT), env_from (BIGINT), env_to (BIGINT), species_prefix (VARCHAR) |
| `predgpi` | VIEW | `predgpi.parquet` | 227,238,177 | species_prefix, protein_id | species_prefix (VARCHAR), protein_id (VARCHAR), source (VARCHAR), feature (VARCHAR), start (BIGINT), end (BIGINT), score (DOUBLE), strand (VARCHAR), phase (VARCHAR), attributes (VARCHAR) |
| `signalp` | VIEW | `signalp.signal_peptide.parquet` | 18,767,284 | species_prefix, protein_id | species_prefix (VARCHAR), protein_id (VARCHAR), peptide_start (BIGINT), peptide_end (BIGINT), probability (DOUBLE) |
| `species` | VIEW | `species.parquet` | 23,680 | LOCUSTAG, ASMID | ASMID (VARCHAR), LOCUSTAG (VARCHAR), PHYLUM (VARCHAR), SUBPHYLUM (VARCHAR), CLASS (VARCHAR), ORDER (VARCHAR), FAMILY (VARCHAR), GENUS (VARCHAR), SPECIES (VARCHAR), BUSCO_LINEAGE (VARCHAR), SUBCLASS (VARCHAR), NCBI_TAXONID (BIGINT), STRAIN (VARCHAR) |
| `swissprot_best` | VIEW | `swissprot_best.parquet` | 157,614,964 | species_prefix, protein_id, sp_id | protein_id (VARCHAR), species_prefix (VARCHAR), sp_id (INTEGER), accession (VARCHAR), entry_name (VARCHAR), org_code (VARCHAR), protein_name (VARCHAR), pident_x10 (USMALLINT), qcov_x10 (USMALLINT), scov_x10 (USMALLINT), aln_length (INTEGER), q_start (INTEGER), q_end (INTEGER), s_start (INTEGER), s_end (INTEGER), q_len (INTEGER), s_len (INTEGER), bitscore (FLOAT), neg_log10_evalue (FLOAT) |
| `swissprot_hits` | VIEW | `swissprot_hits/swissprot_hits.*.parquet` | 2,192,440,076 | species_prefix, protein_id, sp_id | protein_id (VARCHAR), species_prefix (VARCHAR), sp_id (INTEGER), rk (UTINYINT), pident_x10 (USMALLINT), qcov_x10 (USMALLINT), scov_x10 (USMALLINT), aln_length (INTEGER), q_start (INTEGER), q_end (INTEGER), s_start (INTEGER), s_end (INTEGER), q_len (INTEGER), bitscore (FLOAT), neg_log10_evalue (FLOAT) |
| `swissprot_subject` | TABLE | `swissprot_subject.parquet (copied)` | 575,503 | sp_id | sp_id (INTEGER), accession (VARCHAR), entry_name (VARCHAR), org_code (VARCHAR), s_len (INTEGER), protein_name (VARCHAR), gene_name (VARCHAR), gene_synonyms (VARCHAR), organism (VARCHAR), ncbi_taxid (BIGINT), go_terms (VARCHAR), ec_numbers (VARCHAR), interpro (VARCHAR), pfam (VARCHAR) |
| `swissprot_subject_summary` | VIEW | `swissprot_subject_summary.parquet` | 333,632 | sp_id | sp_id (INTEGER), n_hits (BIGINT), n_genomes (BIGINT), n_best_hits (DOUBLE), max_bitscore (FLOAT), max_pident_x10 (USMALLINT) |
| `targetp` | VIEW | `targetP.parquet` | 178,535 | species_prefix, protein_id | species_prefix (VARCHAR), protein_id (VARCHAR), prediction (VARCHAR), probability (DOUBLE), cleavage_position_start (BIGINT), cleavage_position_end (BIGINT), cleavage_probability (DOUBLE), motif (VARCHAR) |
| `telomere_summary` | VIEW | `telomere_summary.parquet` | 23,582 | LOCUSTAG, ASMID | ASMID (VARCHAR), LOCUSTAG (VARCHAR), SPECIES (VARCHAR), STRAIN (VARCHAR), telomere_scaffolds (BIGINT), telomere_scaffolds_both_ends (BIGINT), telomere_tracts (BIGINT), telomere_total_length_bp (BIGINT), telomere_total_repeats (BIGINT), telomere_5prime_count (BIGINT), telomere_3prime_count (BIGINT), telomere_plus_count (BIGINT), telomere_minus_count (BIGINT), telomere_terminal_count (BIGINT), telomere_internal_count (BIGINT), telomere_top_monomer (VARCHAR), telomere_monomers (VARCHAR) |
| `telomere_tracts` | VIEW | `telomere_tracts.parquet` | 73,242 | LOCUSTAG, ASMID | tract_id (VARCHAR), ASMID (VARCHAR), LOCUSTAG (VARCHAR), scaffold (VARCHAR), end_type (VARCHAR), strand (VARCHAR), monomer (VARCHAR), repeat_count (BIGINT), tract_length (BIGINT), start (BIGINT), end_coord (BIGINT), terminal (BOOLEAN), distance_to_end (BIGINT), tract_seq (VARCHAR), flank_seq (VARCHAR) |
| `tf_inventory` | VIEW | `tf_inventory.parquet` | 10,946,927 | species_prefix, protein_id | species_prefix (VARCHAR), protein_id (VARCHAR), pfam_acc (VARCHAR), pfam_name (VARCHAR), tf_family (VARCHAR), domain_i_evalue (DOUBLE), domain_score (DOUBLE) |
| `tmhmm` | VIEW | `tmhmm.parquet` | 45,707,260 | species_prefix, protein_id | species_prefix (VARCHAR), protein_id (VARCHAR), len (BIGINT), ExpAA (DOUBLE), First60 (DOUBLE), PredHel (BIGINT), Topology (VARCHAR) |
| `v_protein_annotation` | VIEW | `other views` | not counted | LOCUSTAG, protein_id | protein_id (VARCHAR), LOCUSTAG (VARCHAR), length (BIGINT), pfam_domain_count (BIGINT), has_cazyme (BOOLEAN), has_merops (BOOLEAN), has_signal_peptide (BOOLEAN), targetp_class (VARCHAR), wolfpsort_loc (VARCHAR), tmhmm_helices (BIGINT), is_gpi_anchored (BOOLEAN), IDP_fraction (DOUBLE) |
| `v_species_summary` | VIEW | `other views` | not counted | LOCUSTAG | LOCUSTAG (VARCHAR), SPECIES (VARCHAR), GENUS (VARCHAR), FAMILY (VARCHAR), ORDER (VARCHAR), CLASS (VARCHAR), SUBPHYLUM (VARCHAR), PHYLUM (VARCHAR), BUSCO_LINEAGE (VARCHAR), TOTAL_LENGTH (BIGINT), GC_PERCENT (DOUBLE), contig_count (BIGINT), N50_bp (BIGINT), masked_pct (DOUBLE), gene_count (BIGINT), mean_protein_length (DOUBLE), pfam_proteins (BIGINT), cazy_proteins (BIGINT), merops_proteins (BIGINT), signalp_count (BIGINT), tmhmm_count (BIGINT), idp_high_disorder (BIGINT) |
| `v_swissprot_80_80` | VIEW | `other views` | not counted | species_prefix, protein_id, sp_id | protein_id (VARCHAR), species_prefix (VARCHAR), sp_id (INTEGER), rk (UTINYINT), pident_x10 (USMALLINT), qcov_x10 (USMALLINT), scov_x10 (USMALLINT), aln_length (INTEGER), q_start (INTEGER), q_end (INTEGER), s_start (INTEGER), s_end (INTEGER), q_len (INTEGER), bitscore (FLOAT), neg_log10_evalue (FLOAT) |
| `v_swissprot_model_fungi` | VIEW | `other views` | not counted | species_prefix, protein_id, sp_id | protein_id (VARCHAR), species_prefix (VARCHAR), sp_id (INTEGER), rk (UTINYINT), pident_x10 (USMALLINT), qcov_x10 (USMALLINT), scov_x10 (USMALLINT), aln_length (INTEGER), q_start (INTEGER), q_end (INTEGER), s_start (INTEGER), s_end (INTEGER), q_len (INTEGER), bitscore (FLOAT), neg_log10_evalue (FLOAT), org_code (VARCHAR), accession (VARCHAR), entry_name (VARCHAR), protein_name (VARCHAR) |
| `v_telomere_both_ends` | VIEW | `other views` | not counted | ASMID | ASMID (VARCHAR), scaffold (VARCHAR) |
| `wgd_ks` | VIEW | `wgd.ks.parquet` | 25,831,001 | species_prefix, genome | genome (VARCHAR), species_prefix (VARCHAR), pair (VARCHAR), family (VARCHAR), g1 (VARCHAR), g2 (VARCHAR), gene1 (VARCHAR), gene2 (VARCHAR), node (VARCHAR), status (VARCHAR), N (DOUBLE), S (DOUBLE), dN (DOUBLE), dN/dS (DOUBLE), dS (DOUBLE), alignmentlength (DOUBLE), strippedalignmentlength (DOUBLE), alignmentcoverage (DOUBLE), alignmentidentity (DOUBLE), t (DOUBLE), pair_codons (DOUBLE), ng_dS (DOUBLE), ng_codons (DOUBLE), copies1 (DOUBLE), copies2 (DOUBLE), npairs (DOUBLE) |
| `wgd_ks_summary` | VIEW | `wgd.ks.summary.parquet` | 7,384 | species_prefix, genome | genome (VARCHAR), species_prefix (VARCHAR), n_pairs_raw (BIGINT), n_pairs (BIGINT), n_pairs_with_ds (BIGINT), n_families (BIGINT), n_nodes (BIGINT), n_nodes_with_ds (BIGINT), hybrid (BIGINT), n_identical_groups (BIGINT), n_genes_in_identical_groups (BIGINT), largest_identical_group (BIGINT), identical_pair_share (DOUBLE) |
| `wolfpsort` | VIEW | `wolfpsort.parquet` | 244,790,951 | species_prefix, protein_id | species_prefix (VARCHAR), protein_id (VARCHAR), localization (VARCHAR), score (DOUBLE) |


### 1.2 Parquet files in `tables/`

Rows are from the parquet footers (no data scan). `swissprot_hits/` holds 256 bucket files.

| Parquet in tables/ | Size (MB) | Rows (footer) | Local DB object(s) on 2026-11-08 | MotherDuck table | Columns (type) |
|---|---|---|---|---|---|
| `aa_freq.parquet` | 1.1 | 477,535 | `aa_frequency` | `aa_freq` | species_prefix (VARCHAR), amino_acid (VARCHAR), frequency (DOUBLE) |
| `asm_stats.parquet` | 2.3 | 23,584 | `asm_stats` | `asm_stats` | ASMID (VARCHAR), SPECIES (VARCHAR), STRAIN (VARCHAR), contig_count (BIGINT), total_length_bp (BIGINT), min_contig_bp (BIGINT), max_contig_bp (BIGINT), median_contig_bp (BIGINT), mean_contig_bp (DOUBLE), L50 (BIGINT), N50_bp (BIGINT), L90 (BIGINT), N90_bp (BIGINT), gc_pct (DOUBLE), n_gap_count (BIGINT), total_n_bases (BIGINT), masked_bases (BIGINT), masked_pct (DOUBLE), t2t_scaffolds (BIGINT), telomere_fwd (BIGINT), telomere_rev (BIGINT), annotation_status (VARCHAR) |
| `busco_genome.parquet` | 0.7 | 23,582 | `busco_genome` | `busco_genome` | ASMID (VARCHAR), complete_pct (DOUBLE), single_pct (DOUBLE), duplicated_pct (DOUBLE), fragmented_pct (DOUBLE), missing_pct (DOUBLE), n_markers (BIGINT), lineage (VARCHAR) |
| `cazy.cazymes_hmm.parquet` | 256.9 | 9,164,068 | `cazy` | `cazy_cazymes_hmm` | species_prefix (VARCHAR), HMM_id (VARCHAR), profile_length (BIGINT), protein_id (VARCHAR), protein_length (BIGINT), evalue (DOUBLE), q_start (BIGINT), q_end (BIGINT), s_start (BIGINT), s_end (BIGINT), coverage (DOUBLE) |
| `cazy.overview.parquet` | 127.8 | 9,310,657 | `cazy_overview` | `cazy_overview` | species_prefix (VARCHAR), protein_id (VARCHAR), EC (VARCHAR), cazyme_fam (VARCHAR), sub_fam (VARCHAR), diamond_fam (VARCHAR), substrate (VARCHAR), toolcount (BIGINT) |
| `codon_freq.parquet` | 2.6 | 1,500,925 | `codon_frequency` | `codon_freq` | species_prefix (VARCHAR), codon (VARCHAR), frequency (DOUBLE) |
| `gene_CDS.parquet` | 7,983.0 | 753,622,745 | `gene_CDS` | `gene_CDS` | cds_id (VARCHAR), transcript_id (VARCHAR), order (BIGINT), chrom (VARCHAR), start (BIGINT), end (BIGINT), strand (BIGINT), GC_content (DOUBLE) |
| `gene_exons.parquet` | 7,987.6 | 758,181,061 | `gene_exons` | `gene_exons` | exon_id (VARCHAR), transcript_id (VARCHAR), order (BIGINT), chrom (VARCHAR), start (BIGINT), end (BIGINT), strand (BIGINT), GC_content (DOUBLE) |
| `gene_info.parquet` | 1,585.7 | 249,546,610 | `gene_info` | `gene_info` | gene_id (VARCHAR), LOCUSTAG (VARCHAR), chrom (VARCHAR), start (BIGINT), end (BIGINT), strand (BIGINT), gene_type (VARCHAR) |
| `gene_intergenic_distances.parquet` | 937.5 | 233,406,624 | `gene_intergenic_distances` | `gene_intergenic_distances` | species_prefix (VARCHAR), left_gene (VARCHAR), right_gene (VARCHAR), distance (BIGINT) |
| `gene_introns.parquet` | 6,819.7 | 508,634,451 | `gene_introns` | `gene_introns` | intron_id (VARCHAR), transcript_id (VARCHAR), intron_number (BIGINT), chrom (VARCHAR), start (BIGINT), end (BIGINT), strand (BIGINT), splice_5 (VARCHAR), splice_3 (VARCHAR), GC_content (DOUBLE), codon_position (BIGINT), codon_frame (BIGINT) |
| `gene_proteins.parquet` | 5,532.3 | 244,988,294 | `gene_proteins` | `gene_proteins` | gene_id (VARCHAR), protein_id (VARCHAR), transcript_id (VARCHAR), length (BIGINT), md5checksum (VARCHAR) |
| `gene_transcripts.parquet` | 6,768.9 | 249,546,610 | `gene_transcripts` | `gene_transcripts` | gene_id (VARCHAR), transcript_id (VARCHAR), chrom (VARCHAR), start (BIGINT), end (BIGINT), strand (BIGINT), is_partial (BOOLEAN), has_start_codon (VARCHAR), has_stop_codon (VARCHAR), CDS_length (BIGINT), md5checksum (VARCHAR) |
| `gene_trnas.parquet` | 21.3 | 4,558,316 | `gene_trna` | `gene_trnas` | gene_id (VARCHAR), amino_acid (VARCHAR), codon (VARCHAR) |
| `idp.parquet` | 4,197.4 | 233,329,228 | `idp` | `idp` | species_prefix (VARCHAR), protein_id (VARCHAR), IDP_start (BIGINT), IDP_end (BIGINT), IDP_length (BIGINT), mean_score (DOUBLE) |
| `idp_summary.parquet` | 2,203.7 | 144,924,594 | `idp_summary` | `idp_summary` | species_prefix (VARCHAR), protein_id (VARCHAR), IDP_residues (BIGINT), IDP_fraction (DOUBLE), length (BIGINT) |
| `infernal_rfam.parquet` | 0.1 | 3,925 | **none** | `infernal_rfam` | species_prefix (VARCHAR), target_name (BIGINT), target_acc (VARCHAR), query_name (VARCHAR), contig (VARCHAR), mdl_from (VARCHAR), mdl_to (VARCHAR), seq_from (BIGINT), seq_to (BIGINT), strand (BIGINT), evalue (DOUBLE), score (DOUBLE) |
| `merops.parquet` | 1,516.6 | 43,462,876 | `merops` | `merops` | species_prefix (VARCHAR), protein_id (VARCHAR), merops_id (VARCHAR), percent_identity (DOUBLE), aln_length (BIGINT), mismatches (BIGINT), gap_openings (BIGINT), q_start (BIGINT), q_end (BIGINT), s_start (BIGINT), s_end (BIGINT), evalue (DOUBLE), bitscore (DOUBLE) |
| `pfam.parquet` | 12,049.7 | 368,209,721 | `pfam` | `pfam` | protein_id (VARCHAR), pfam_id (VARCHAR), pfam_acc (VARCHAR), pfam_len (BIGINT), full_seq_e_value (DOUBLE), full_seq_score (DOUBLE), full_seq_bias (DOUBLE), domain_num (BIGINT), domain_num_of (BIGINT), domain_c_evalue (DOUBLE), domain_i_evalue (DOUBLE), domain_score (DOUBLE), domain_bias (DOUBLE), hmm_from (BIGINT), hmm_to (BIGINT), ali_from (BIGINT), ali_to (BIGINT), env_from (BIGINT), env_to (BIGINT) |
| `predgpi.parquet` | 1,994.0 | 227,238,177 | `predgpi` | `predgpi` | species_prefix (VARCHAR), protein_id (VARCHAR), source (VARCHAR), feature (VARCHAR), start (BIGINT), end (BIGINT), score (DOUBLE), strand (VARCHAR), phase (VARCHAR), attributes (VARCHAR) |
| `samples.parquet` | 1.1 | 23,680 | **none** | `samples` | ASMID (VARCHAR), SPECIES_IN (VARCHAR), STRAIN (VARCHAR), BIOPROJECT (VARCHAR), NCBI_TAXONID (BIGINT), BUSCO_LINEAGE (VARCHAR), PHYLUM (VARCHAR), SUBPHYLUM (VARCHAR), CLASS (VARCHAR), SUBCLASS (VARCHAR), ORDER (VARCHAR), FAMILY (VARCHAR), GENUS (VARCHAR), SPECIES (VARCHAR), TRANSL_TABLE (BIGINT), LOCUSTAG (VARCHAR) |
| `signalp.signal_peptide.parquet` | 277.3 | 18,767,284 | `signalp` | `signalp_signal_peptide` | species_prefix (VARCHAR), protein_id (VARCHAR), peptide_start (BIGINT), peptide_end (BIGINT), probability (DOUBLE) |
| `species.parquet` | 0.9 | 23,680 | `asm_stats`, `busco_genome`, `species` | `species` | ASMID (VARCHAR), LOCUSTAG (VARCHAR), PHYLUM (VARCHAR), SUBPHYLUM (VARCHAR), CLASS (VARCHAR), ORDER (VARCHAR), FAMILY (VARCHAR), GENUS (VARCHAR), SPECIES (VARCHAR), BUSCO_LINEAGE (VARCHAR), SUBCLASS (VARCHAR), NCBI_TAXONID (BIGINT), STRAIN (VARCHAR) |
| `swissprot_annot.parquet` | 16.8 | 575,503 | **none** | `swissprot_annot` | accession (VARCHAR), entry_name (VARCHAR), protein_name (VARCHAR), gene_name (VARCHAR), gene_synonyms (VARCHAR), organism (VARCHAR), go_terms (VARCHAR), ec_numbers (VARCHAR), interpro (VARCHAR), pfam (VARCHAR), length (BIGINT), ncbi_taxid (BIGINT) |
| `swissprot_best.parquet` | 5,111.6 | 157,614,964 | `swissprot_best` | `swissprot_best` | protein_id (VARCHAR), species_prefix (VARCHAR), sp_id (INTEGER), accession (VARCHAR), entry_name (VARCHAR), org_code (VARCHAR), protein_name (VARCHAR), pident_x10 (USMALLINT), qcov_x10 (USMALLINT), scov_x10 (USMALLINT), aln_length (INTEGER), q_start (INTEGER), q_end (INTEGER), s_start (INTEGER), s_end (INTEGER), q_len (INTEGER), s_len (INTEGER), bitscore (FLOAT), neg_log10_evalue (FLOAT) |
| `swissprot_hits/ (256 files)` | 39,635.8 | 2,192,440,076 (BFD_SCHEMA.md; not recounted) | `swissprot_hits` | not loaded (excluded by default) | protein_id (VARCHAR), species_prefix (VARCHAR), sp_id (INTEGER), rk (UTINYINT), pident_x10 (USMALLINT), qcov_x10 (USMALLINT), scov_x10 (USMALLINT), aln_length (INTEGER), q_start (INTEGER), q_end (INTEGER), s_start (INTEGER), s_end (INTEGER), q_len (INTEGER), bitscore (FLOAT), neg_log10_evalue (FLOAT) |
| `swissprot_subject.parquet` | 26.4 | 575,503 | `swissprot_subject (table)` | `swissprot_subject` | sp_id (INTEGER), accession (VARCHAR), entry_name (VARCHAR), org_code (VARCHAR), s_len (INTEGER), protein_name (VARCHAR), gene_name (VARCHAR), gene_synonyms (VARCHAR), organism (VARCHAR), ncbi_taxid (BIGINT), go_terms (VARCHAR), ec_numbers (VARCHAR), interpro (VARCHAR), pfam (VARCHAR) |
| `swissprot_subject_summary.parquet` | 3.7 | 333,632 | `swissprot_subject_summary` | `swissprot_subject_summary` | sp_id (INTEGER), n_hits (BIGINT), n_genomes (BIGINT), n_best_hits (DOUBLE), max_bitscore (FLOAT), max_pident_x10 (USMALLINT) |
| `targetP.parquet` | 4.8 | 178,535 | `targetp` | `targetP` | species_prefix (VARCHAR), protein_id (VARCHAR), prediction (VARCHAR), probability (DOUBLE), cleavage_position_start (BIGINT), cleavage_position_end (BIGINT), cleavage_probability (DOUBLE), motif (VARCHAR) |
| `telomere_summary.parquet` | 1.2 | 23,582 | `telomere_summary` | `telomere_summary` | ASMID (VARCHAR), LOCUSTAG (VARCHAR), SPECIES (VARCHAR), STRAIN (VARCHAR), telomere_scaffolds (BIGINT), telomere_scaffolds_both_ends (BIGINT), telomere_tracts (BIGINT), telomere_total_length_bp (BIGINT), telomere_total_repeats (BIGINT), telomere_5prime_count (BIGINT), telomere_3prime_count (BIGINT), telomere_plus_count (BIGINT), telomere_minus_count (BIGINT), telomere_terminal_count (BIGINT), telomere_internal_count (BIGINT), telomere_top_monomer (VARCHAR), telomere_monomers (VARCHAR) |
| `telomere_tracts.parquet` | 17.9 | 73,242 | `telomere_tracts` | `telomere_tracts` | tract_id (VARCHAR), ASMID (VARCHAR), LOCUSTAG (VARCHAR), scaffold (VARCHAR), end_type (VARCHAR), strand (VARCHAR), monomer (VARCHAR), repeat_count (BIGINT), tract_length (BIGINT), start (BIGINT), end_coord (BIGINT), terminal (BOOLEAN), distance_to_end (BIGINT), tract_seq (VARCHAR), flank_seq (VARCHAR) |
| `tf_inventory.parquet` | 174.6 | 10,946,927 | `tf_inventory` | `tf_inventory` | species_prefix (VARCHAR), protein_id (VARCHAR), pfam_acc (VARCHAR), pfam_name (VARCHAR), tf_family (VARCHAR), domain_i_evalue (DOUBLE), domain_score (DOUBLE) |
| `tmhmm.parquet` | 1,897.7 | 45,707,260 | `tmhmm` | `tmhmm` | species_prefix (VARCHAR), protein_id (VARCHAR), len (BIGINT), ExpAA (DOUBLE), First60 (DOUBLE), PredHel (BIGINT), Topology (VARCHAR) |
| `wgd.ks.nodes.parquet` | 596.8 | 28,982,469 | **none** | `wgd_ks_nodes` | genome (VARCHAR), species_prefix (VARCHAR), family (VARCHAR), node (VARCHAR), kind (VARCHAR), node_pairs (BIGINT), node_rep_pairs (BIGINT), node_pass (BIGINT), node_sampled (BIGINT), node_gated (BIGINT), n_valid (BIGINT), dS_mean (DOUBLE), dS_min (DOUBLE), dS_max (DOUBLE), identity_mean (DOUBLE) |
| `wgd.ks.parquet` | 2,473.5 | 25,831,001 | `wgd_ks` | `wgd_ks` | genome (VARCHAR), species_prefix (VARCHAR), pair (VARCHAR), family (VARCHAR), g1 (VARCHAR), g2 (VARCHAR), gene1 (VARCHAR), gene2 (VARCHAR), node (VARCHAR), status (VARCHAR), N (DOUBLE), S (DOUBLE), dN (DOUBLE), dN/dS (DOUBLE), dS (DOUBLE), alignmentlength (DOUBLE), strippedalignmentlength (DOUBLE), alignmentcoverage (DOUBLE), alignmentidentity (DOUBLE), t (DOUBLE), pair_codons (DOUBLE), ng_dS (DOUBLE), ng_codons (DOUBLE), copies1 (DOUBLE), copies2 (DOUBLE), npairs (DOUBLE) |
| `wgd.ks.summary.parquet` | 0.3 | 7,384 | `wgd_ks_summary` | `wgd_ks_summary` | genome (VARCHAR), species_prefix (VARCHAR), n_pairs_raw (BIGINT), n_pairs (BIGINT), n_pairs_with_ds (BIGINT), n_families (BIGINT), n_nodes (BIGINT), n_nodes_with_ds (BIGINT), hybrid (BIGINT), n_identical_groups (BIGINT), n_genes_in_identical_groups (BIGINT), largest_identical_group (BIGINT), identical_pair_share (DOUBLE) |
| `wolfpsort.parquet` | 1,890.3 | 244,790,951 | `wolfpsort` | `wolfpsort` | species_prefix (VARCHAR), protein_id (VARCHAR), localization (VARCHAR), score (DOUBLE) |


Other files in `tables/`: `BUSCO.csv.gz`, `_manifest/` (JSON per table, read by `_table_manifest`), `wgd.ks.README.md`, `telomeres.parquet.bak.202608051217`, six empty `.lock.*` files (2026-11-06).

Parquet files with **no** object in the published database: `samples.parquet`, `swissprot_annot.parquet`, `wgd.ks.nodes.parquet`, `infernal_rfam.parquet`. This confirms the comparison you gave me. `cazy.cazymes_hmm.parquet` has a view, named `cazy`.

### 1.3 Join keys

| Key | Where it is | Notes |
|---|---|---|
| `LOCUSTAG` (8 hex characters) | `species`, `asm_stats`, `busco_genome`, `gene_info`, `telomere_summary`, `telomere_tracts`, `samples.parquet` | The genome key. `species` has 23,680 rows, all with a distinct LOCUSTAG |
| `ASMID` | `species`, `asm_stats`, `busco_genome`, `telomere_*`, `samples.parquet` | Assembly key. `asm_stats.parquet` and `busco_genome.parquet` hold only ASMID; the views add LOCUSTAG by a join to `species.parquet` |
| `species_prefix` | every per-protein result table, `gene_intergenic_distances`, `aa_frequency`, `codon_frequency`, `wgd_*` | Equal to LOCUSTAG. Stored in the parquet for most tables. **Computed by the view** (`string_split(id,'_')[1]`) for `pfam`, `gene_proteins`, `gene_transcripts`, `gene_exons`, `gene_CDS`, `gene_introns`, `gene_trna` |
| `genome` | `wgd_ks`, `wgd_ks_summary`, `wgd.ks.nodes.parquet` | The **sample name** (for example `Aaosphaeria_arxii_CBS_175.79`), not the LOCUSTAG. These tables also have `species_prefix`. In `wgd_ks_summary`, `genome = species_prefix` in 0 of 7,384 rows |
| `gene_id` | `gene_info`, `gene_transcripts`, `gene_proteins`, `gene_trna`; `left_gene`/`right_gene` in `gene_intergenic_distances` | Form `F0000039_000001` |
| `transcript_id` | `gene_transcripts`, `gene_proteins`, `gene_exons`, `gene_CDS`, `gene_introns` | Form `F0000039_000001-T1` |
| `protein_id` | `gene_proteins` and the result tables | Form `F0000039_000001-T1`. Exception: `idp` and `idp_summary` (below) |
| `sp_id` | `swissprot_hits`, `swissprot_best`, `swissprot_subject`, `swissprot_subject_summary` | SwissProt subject key |

### 1.4 Inconsistencies found

Each item was checked on the data unless marked.

1. **Three names for the genome key**: `LOCUSTAG`, `species_prefix`, and (in the atlas files) `genome`. In the wgd tables `genome` means the sample name instead. A join between an atlas file and a wgd table on `genome` returns nothing.
2. **`idp` and `idp_summary` use the gene-id form in `protein_id`** (`F77FFF1D_000001`, no `-T1`). For genome F77FFF1D, 0 of 6,252 `idp_summary` rows join to `gene_proteins.protein_id`; 6,252 of 6,252 join to `gene_proteins.gene_id`. Fungi_BFD commit `f27fc8e` ("gather_AIUPred: use protein id as protein_id") changes the builder; the parquet (2026-11-06) is older than that commit.
3. **`infernal_rfam.parquet` is malformed.** `species_prefix` holds the sample name (`Aaosphaeria_arxii_CBS_175`), not the LOCUSTAG. The other columns are shifted by one: `target_name` (BIGINT) is a row number, `target_acc` is the model name (`U6`), `query_name` is the Rfam accession. 3,925 rows, 14 genomes.
4. **`asm_stats` column names differ between the view and the parquet.** View: `TOTAL_LENGTH`, `GC_PERCENT`, plus `LOCUSTAG`; no `SPECIES`, `STRAIN`. Parquet (and MotherDuck): `total_length_bp`, `gc_pct`, `SPECIES`, `STRAIN`. `median_contig_bp` is DOUBLE in the view and BIGINT in the parquet. The atlas label files use the parquet names.
5. **`gene_introns.codon_position` and `codon_frame` are off by one exon** in the 2026-11-08 parquet (section 4.3).
6. **Text booleans.** `gene_transcripts.has_start_codon` and `has_stop_codon` are VARCHAR with the values `'TRUE'`, `'FALSE'` and `'NULL'` (the text). `is_partial` in the same table is BOOLEAN.
7. **Positive-only tables.** `signalp`, `tmhmm` (minimum `PredHel` = 1 in the genome checked), `targetp`, `cazy_overview`, `merops`, `tf_inventory` list only proteins with a hit. `predgpi` is the opposite: it has a `Chain` row for every protein (10,275 of 10,275 in F77FFF1D) and a `GPI-anchor` row for 165.
8. **Partial coverage.** `targetp`: 125 genomes. `infernal_rfam`: 14. `wgd_ks_summary`: 7,384. `asm_stats`: 23,584 rows, 2 with no LOCUSTAG; 98 `species` rows have no `asm_stats` row. One prefix, `F0B62757`, has rows in `aa_frequency` and `cazy_overview` but no row in `species`.
9. **Empty taxonomy columns.** `species.SUBPHYLUM` and `species.SUBCLASS` are NULL in all 23,680 rows. `CLASS` is filled in 23,462.
10. **Reserved words and quoting.** `species."ORDER"` (the atlas files use `ORDER_`), `"order"`, `"start"`, `"end"` in the gene tables, `wgd_ks."dN/dS"`. `telomere_tracts` uses `end_type` and `end_coord` to avoid `end`; the gene tables do not.
11. **Different names and types for the same thing.**
    - Protein length: `gene_proteins.length`, `idp_summary.length`, `tmhmm.len`, `cazy.protein_length`, `swissprot_*.q_len`.
    - E-value: `evalue` (merops, cazy), `domain_i_evalue` / `full_seq_e_value` (pfam), `neg_log10_evalue` (swissprot).
    - Identity: `merops.percent_identity` (percent), `swissprot_*.pident_x10` (percent x 10), `wgd_ks.alignmentidentity`.
    - Pfam name: `pfam.pfam_id` is the name and `pfam.pfam_acc` the accession with a version suffix; `tf_inventory` calls the name `pfam_name`.
    - Strand: BIGINT (-1/1) in the gene tables, VARCHAR in `telomere_tracts` and `predgpi`.
    - `gene_trna.codon`: the MotherDuck comment calls it the anticodon. Not checked.
    - Letter case is mixed: `IDP_fraction`, `PredHel`, `ExpAA`, `GC_content`, `CDS_length`, `TOTAL_LENGTH`, `N50_bp`.
12. **Duplicate content.** `swissprot_annot.parquet` and `swissprot_subject.parquet` have the same 575,503 accessions. `swissprot_subject` adds `sp_id`, `org_code`, `s_len`; `swissprot_annot` has `length`.
13. **Sort order.** From the parquet footers (strict overlap of row-group min/max):
    - Sorted, 0 or 1 overlapping row groups: `pfam` (protein_id); `signalp`, `tmhmm`, `predgpi`, `wolfpsort`, `merops`, `idp`, `idp_summary`, `swissprot_best` (species_prefix and protein_id); `wgd.ks`, `wgd.ks.nodes`.
    - Nearly sorted, neighbours overlap: `gene_info` (1,647 of 2,015 row groups), `gene_CDS` (4,809 of 6,086), `gene_introns` (3,392 of 4,108), `gene_intergenic_distances` (1,588 of 1,885). Lookups are still fast (0.3-0.9 s measured), because only neighbours overlap. This differs from the statement in `docs/BFD_Introduction.md` ("0 row groups whose id range overlaps"). I did not check `gene_proteins`, `gene_transcripts`, `gene_exons`.
    - Not sorted: `cazy.overview`, `cazy.cazymes_hmm`, `tf_inventory`, `targetP`, `gene_trnas`. All are small (5-257 MB).
14. **Defects in two existing pipeline views** (not changed by my files; they need a fix in `build_BFD_duckDB.sh`):
    - `v_protein_annotation.is_gpi_anchored` is built from `SELECT DISTINCT protein_id FROM predgpi` with no filter on `feature`. `predgpi` has a `Chain` row for every protein, so the flag is true for every protein with a PredGPI row. Checked on F77FFF1D: 10,275 proteins with a row, 165 with `feature = 'GPI-anchor'`. I did not run the view itself (it scans `pfam`).
    - `v_protein_annotation.IDP_fraction` joins `idp_summary` on `protein_id`. Because of item 2 the join finds nothing (0 of 6,252 for F77FFF1D).
    - `v_protein_annotation.targetp_class` is `'OTHER'` for every protein of the genomes without TargetP data (all but 125).
    - `v_species_summary` uses inner joins to `asm_stats` and `gene_proteins`, so genomes without proteins are dropped. Not a defect, but it is not "one row per genome".

---

## 2. MotherDuck compatibility

Rule tested: SQL written for MotherDuck uses `bfd.<table>` or the bare table name. With `lib/bfd_connect.py` the local file is attached as `bfd`. With `duckdb BFD.duckdb` the catalog is named `BFD`; `bfd.species` and `bfd.main.species` both resolve (tested with CLI 1.1.3).

"Runs" means: the name resolves and the parquet columns are present with the same names.

| MotherDuck table | Local object (2026-11-08) | Runs before the fix? | Fix in `compat_views.sql` | Runs after? |
|---|---|---|---|---|
| `aa_freq` | `aa_frequency` | No (name) | `CREATE OR REPLACE VIEW aa_freq AS SELECT * FROM aa_frequency` | Yes |
| `codon_freq` | `codon_frequency` | No (name) | alias of `codon_frequency` | Yes |
| `gene_trnas` | `gene_trna` (+ `species_prefix`) | No (name) | `SELECT * EXCLUDE (species_prefix) FROM gene_trna` | Yes |
| `signalp_signal_peptide` | `signalp` | No (name) | alias of `signalp` | Yes |
| `cazy_cazymes_hmm` | `cazy` | No (name) | alias of `cazy` | Yes |
| `samples` | none | No | view over `samples.parquet` | Yes |
| `swissprot_annot` | none | No | view over `swissprot_annot.parquet` | Yes |
| `wgd_ks_nodes` | none | No | view over `wgd.ks.nodes.parquet` | Yes |
| `infernal_rfam` | none | No | view over `infernal_rfam.parquet` (data defect, section 1.4 item 3) | Yes |
| `targetP` | `targetp` | Yes (names are case-insensitive; tested) | none | Yes |
| `asm_stats` | `asm_stats` | **Partly.** `ASMID`, `N50_bp`, `annotation_status`, `contig_count`, `masked_pct` and the other shared names work. `total_length_bp`, `gc_pct`, `SPECIES`, `STRAIN` fail (Binder Error, tested) | Commented-out block in `compat_views.sql` section 3 adds `total_length_bp` and `gc_pct`. **Decision needed** | Partly until decided |
| `busco_genome` | `busco_genome` | Yes. The view adds `LOCUSTAG` | none | Yes |
| `species` | `species` | Yes (same 13 columns) | none | Yes |
| `gene_info`, `gene_intergenic_distances` | same name | Yes (same columns) | none | Yes |
| `gene_transcripts`, `gene_proteins`, `gene_exons`, `gene_CDS`, `gene_introns`, `pfam` | same name | Yes. The views add `species_prefix` | none | Yes |
| `cazy_overview`, `merops`, `tf_inventory`, `wolfpsort`, `predgpi`, `tmhmm`, `idp`, `idp_summary` | same name | Yes (same columns) | none | Yes |
| `telomere_summary`, `telomere_tracts`, `wgd_ks`, `wgd_ks_summary` | same name | Yes (same columns) | none | Yes |
| `swissprot_best`, `swissprot_subject`, `swissprot_subject_summary` | same name (`swissprot_subject` is a base table) | Yes (same columns) | none | Yes |

Total: 36 MotherDuck tables. 26 ran unchanged, 9 needed a view (all 9 now run), 1 (`asm_stats`) runs only for the shared column names.

`swissprot_hits` was not in MotherDuck (excluded by default in `md_load_bfd.py:35`). It exists locally as a view and through the lookup macros.

One real MotherDuck query was run as a test: the representative query of `docs/analyses/kingdom_domains/README.md:37-42` (as a SELECT). It returned 7,435 genomes, the number the README reports.

Things that still differ after the fix:

- `SELECT *` on the gene tables, `pfam` and `busco_genome` returns one more column locally (`species_prefix` or `LOCUSTAG`). A query that depends on column position or uses `UNION` with `SELECT *` can break. I found no such query in the files I read.
- MotherDuck had `COMMENT ON TABLE` for every table. The local views had none. `compat_views.sql` adds the same texts with `COMMENT ON VIEW` (42 statements). The texts are copied from `DEFAULT_COMMENTS`; I added a short note where the local view differs (`asm_stats`, `busco_genome`, `gene_trna`, the renamed views, `infernal_rfam`).
- `bfd_scratch` is a separate matter (`lib/bfd_connect.py`); not reviewed.
- `lib/bfd_connect.py` defaults to `BFD_dashboard/db/BFD.duckdb`. That folder does not exist yet (checked 2026-11-11).

---

## 3. Query patterns

Sources read: the three analysis folders under `BFD_dashboard/docs/analyses/` (READMEs and scripts), `docs/FUNCTIONAL_ANALYSES.md`, `bin/atlas/sql/*.sql` (9 files), `bin/atlas/*.py`, `web/chat.html`, `lib/bfd_data.py`, `lib/subset.py`, `lib/bfd_connect.py`. Paths below are relative to `BFD_dashboard/`.

`Fungi_BFD/analysis/` holds pipeline-performance analyses. I found no SQL there against the BFD tables (2 files mention duckdb; neither queries these tables). I found no skill or MCP file that describes the BFD tables. `~/.claude.json` contains a `motherduck` entry; I did not open it. `Fungi_BFD/.claude/settings.local.json` mentions duckdb only in a permission rule.

### 3.1 How the tables are reached

1. **Most analysis code does not query the database tables. It reads the atlas intermediates** in `atlas_data/` (parquet built once by jobs J0-J9):
   - `labels/genome_labels.parquet` and `genome_labels_ext.parquet`: `docs/analyses/kickxellomycota/scripts/k1.py:4-5`, `docs/analyses/antarctolithica/scripts/01_build_set.py:6`, `bin/atlas/sql/J3_protein_core.sql:35,50`.
   - `A2/genome_gene_stats.parquet`: `k1.py:6`, `antarctolithica/scripts/04_build_features.py:12`.
   - `A5/secretome`, `cazy_genome`, `merops_genome`, `tf_genome`, `idp_genome`, `cazy_class_counts`, `merops_family_counts`, `tf_counts`, `wolfpsort_counts`: `k1.py:7-11`, `04_build_features.py:46-56`, `kickxellomycota/scripts/k6_pf.py:6`.
   - `A7/pfam_genome_matrix.parquet`: `k6_pf.py:5`, `04_build_features.py:42`.
2. **The atlas jobs read the parquet files directly** (`FROM '{T}/pfam.parquet'`), not the database: all of `bin/atlas/sql/*.sql`.
3. **MotherDuck SQL** (`bfd.<table>`) appears in `docs/analyses/kingdom_domains/README.md:37-49` and `docs/analyses/kickxellomycota/README.md:447-457`. Tables: `species`, `asm_stats`, `busco_genome`, `pfam`, `gene_introns`, `gene_transcripts`, `gene_CDS`, `swissprot_best`.
4. **`lib/bfd_data.py` uses the pipeline names** (`aa_frequency`, `codon_frequency`, `signalp`, `targetp`, `v_species_summary`): lines 61-63, 160, 252, 258, 445-466.
5. **`web/chat.html` describes a different schema** (the WASM dashboard files `genome`, `protein`, `pfam_by_domain`, `pfam_by_genome`, `genome_pfam_counts`, `genome_features`, `pfam_dict`; lines 81-124). It is not the BFD.duckdb schema. It shows which small tables were built to keep the chat model away from the large ones.

### 3.2 Patterns that repeat

| # | Pattern | Evidence |
|---|---|---|
| P1 | Join `species` + `asm_stats` + `busco_genome`, then keep one genome per species (best BUSCO, then N50) | `kingdom_domains/README.md:37-42`; `kickxellomycota/README.md:447-448`; `bin/atlas/build_labels.py:102-112,152,156-157`; `lib/subset.py:106-113`; `bin/atlas/sql/J3_protein_core.sql:35-36,50-51` |
| P2 | Quality flags `q_ann`, `q_busco_rel`, `q_dup`, `q_n50` used as filters | `k1.py:15,21`; `kickxellomycota/scripts/k3_groups.py:45`; `bin/atlas/hs_labels.py:105`; `docs/FUNCTIONAL_ANALYSES.md:18-21` |
| P3 | Genome from an id: `split_part(id,'_',1)`, `left(id,8)`, `string_split` | `kingdom_domains/README.md:45,48`; `kickxellomycota/README.md:450,455`; `kickxellomycota/scripts/k4_scan.py:11,14`; `bin/atlas/sql/J0_coverage.sql:3,9,11,13`; `J1_gene_stats.sql:16,32,51,54`; `J2_pfam.sql:4` |
| P4 | Pfam gate `domain_i_evalue <= 1e-5`, version suffix removed, count distinct proteins per genome and Pfam | `bin/atlas/sql/J2_pfam.sql:4-10`; `kingdom_domains/README.md:44-49`; `kickxellomycota/README.md:455-456` |
| P5 | Complete-model flag `NOT is_partial AND has_start_codon = 'TRUE' AND has_stop_codon = 'TRUE'` | `bin/atlas/sql/J1_gene_stats.sql:21-26`; `J3_protein_core.sql:10-12`; `kickxellomycota/README.md:451` |
| P6 | Secretome: SignalP joined to TMHMM and PredGPI (`feature = 'GPI-anchor'`); strict = probability >= 0.9, <= 1 helix, no GPI | `bin/atlas/sql/J4_function_counts.sql:5-34` and again `:39-46`; counts in `lib/bfd_data.py:445-458` |
| P7 | Per-genome family counts: CAZy (`toolcount >= 2`, split `cazyme_fam` on `+`), MEROPS (best hit, `evalue <= 1e-10`), TF (`domain_i_evalue <= 1e-5`) | `bin/atlas/sql/J4_function_counts.sql:55-100`; `lib/bfd_data.py:338-359` |
| P8 | SwissProt gates on the x10 columns: `qcov_x10 >= 800 AND scov_x10 >= 800 AND pident_x10 >= 400` (or 300) | `bin/atlas/sql/J3_protein_core.sql:37,70`; `kickxellomycota/README.md:459`; `k4_scan.py:51` |
| P9 | Intron position recomputed from cumulative `gene_CDS` lengths, because the table columns are wrong | `kickxellomycota/README.md:16,39,452-454,562`; `docs/FUNCTIONAL_ANALYSES.md:32,217` |
| P10 | Count of genomes with data per table, to tell "zero" from "not run" | `lib/bfd_data.py:61-75,127`; `bin/atlas/sql/J0_coverage.sql` (all); `kickxellomycota/README.md:59` |

### 3.3 Traps named in the sources

- "No label table in the queryable database ... In MotherDuck I had to rebuild representatives in SQL" (`kickxellomycota/README.md:55`). The rebuilt set had 1,919 genomes against 1,944 from the atlas (`:415`).
- "No genome x Pfam matrix in the release. Any presence/absence question scans 368 M rows" (`:56`).
- "`SUBPHYLUM` is empty" (`:57`; confirmed, section 1.4 item 9).
- "Key-form traps (`protein_id` vs `gene_id` forms, `idp_summary` keys, strings 'TRUE'/'FALSE')" (`:58`; confirmed, items 2 and 6).
- "never infer absence from tables that list only positives (tmhmm, signalp, targetP)" (`:59`).
- `web/chat.html:88` tells the model to use `pfam_acc LIKE 'PF02364%'`. On the full `pfam` table a LIKE cannot use the row-group statistics.
- `species_prefix = 'X'` on `pfam` or a gene table is pushed into the scan as `array_extract(string_split(protein_id,'_'),1) = 'X'` (EXPLAIN output, section 6). That is a row filter, not a range, so every row group is read. `build_BFD_duckDB.sh:69-71` reports 8.5 s for one genome in `gene_proteins` through the view. A range on the id returns in 0.1-0.4 s (measured here).

---

## 4. Proposed views

### 4.1 `compat_views.sql` (9 views + comments)

See section 2. All 9 views have the same column names and types as the parquet file, and the same row count (checked for each).

### 4.2 `extra_views.sql` (10 views + 1 macro)

| Object | Pattern | What it gives | Validation |
|---|---|---|---|
| `in_genome(id, locustag)` (macro) | P3, section 3.3 | `id >= locustag\|\|'_' AND id < locustag\|\|'`'`. A genome filter that DuckDB pushes into the parquet scan as a range | EXPLAIN shows the range in the READ_PARQUET filter. `pfam` for F07D2B4E: 17,166 rows in 0.2-0.4 s (the count in `BFD_Introduction.md` is also 17,166) |
| `v_genome_labels` | P1, P2 | One row per genome: taxonomy, assembly, BUSCO, `q_ann`, `q_busco_rel`, `q_busco_rel50`, `q_dup`, `q_n50`, `is_representative`, `is_genus_rep`. Rules copied from `bin/atlas/build_labels.py`. Reads only `species`, `samples`, `asm_stats`, `busco_genome` | Compared with `atlas_data/labels/genome_labels.parquet` on 31 shared columns for all 23,680 genomes: 0 differing cells. 7,435 species representatives, 1,669 genus representatives. 0.1 s |
| `v_genome_coverage` | P10 | One row per genome, `has_*` flags for 13 tables | 23,680 rows. has_targetp = 125, has_wgd_ks = 7,384, has_signalp = 23,410, has_tmhmm = 23,409. 6-13 s (python 1.5.5), 13-26 s (CLI 1.1.3) |
| `v_secretome_proteins` | P6 | One row per SignalP protein: `signalp_prob`, `tm_helices`, `is_gpi`, `secreted_noTM`, `secreted_strict` | F77FFF1D: 921 rows, 694 strict; equal to the atlas row. One genome 0.1 s (filter reaches all three scans in both DuckDB versions) |
| `v_genome_secretome` | P6 | Per-genome counts | Equal to `atlas_data/A5/secretome.parquet` in 23,410 of 23,410 genomes (5 columns). 5-9 s |
| `v_cazy_family_counts` | P7 | Proteins per genome and CAZy family | Equal to `A5/cazy_class_counts.parquet`: 2,242,163 of 2,242,163 rows. 1.3 s (1.5.5), 11-16 s (1.1.3) |
| `v_tf_family_counts` | P7 | Proteins per genome and TF family | Equal to `A5/tf_counts.parquet`: 413,480 of 413,480 rows. 0.3-2 s |
| `v_merops_best` | P7 | Best MEROPS hit per protein (`evalue <= 1e-10`) | 6,750,232 rows, the row count of `A5/merops_best.parquet`. F77FFF1D: 363 = atlas `n_merops`. 1-3 s |
| `v_genome_function_counts` | P6, P7 (`k1.py:4-11` join) | One row per genome: secretome counts, `n_cazy`, `n_cazy_secreted`, `n_merops`, `n_merops_secreted`, `n_tf_proteins` | Equal to `A5/cazy_genome` (23,392 of 23,392), `merops_genome` (23,404 of 23,404), `tf_genome` (23,394 of 23,394), `secretome` (23,410 of 23,410). 9-12 s (1.5.5), 13-20 s (1.1.3) |
| `v_gene_transcripts` | P5 | `gene_transcripts` plus BOOLEAN `has_start`, `has_stop`, `q_model` | F77FFF1D: 8,783 `q_model` of 10,377; the text `'NULL'` rows (102) give NULL flags and `q_model` = false. Tested with a genome filter and `LIMIT 5` only |
| `v_gene_introns_phase` | P9 | `gene_introns` with `cds_nt_upstream`, corrected `codon_position` and `codon_frame`, `in_cds`; raw values kept as `*_raw` | Section 4.3 |

### 4.3 `v_gene_introns_phase`: is it still needed?

Yes, on 2026-11-11.

- `gene_introns.parquet` was written 2026-11-08 22:13. The builder fix is Fungi_BFD commit `34c8d14` (2026-11-09 21:24, "compute intron position correctly", on `main`). The parquet is older than the fix.
- The bug is in the data. Transcript `F77FFF1D_000014-T1` has CDS segments of 146, 82, 103, 397, 192, 688 nt. Intron 0 follows 146 coding bases, so the correct values are position 48, phase 2. The table holds 76 and 0, which belong to the next junction (228 nt). The last intron holds 536 and 0 (1,608 / 3: the whole CDS).
- The view gives 48/2, 76/0, 110/1, 242/2, 306/2 for the five introns.
- Check on three whole genomes (F77FFF1D, F07D2B4E, F2D1C788; 56,585 introns): the corrected value of intron n equals the raw value of intron n-1 in 31,199 of 31,201 later introns. 0 introns have corrected = raw. All introns were matched to a CDS junction (`in_cds` true).
- The 2 mismatches are in two 5'-partial transcripts (`has_start_codon = 'FALSE'`). For these the reading frame does not start at the first CDS base and the tables do not hold the GFF phase. The view comment states this limit. I did not check what the fixed builder returns for such models.
- The `in_cds = false` branch (an intron outside the CDS) did not occur in the genomes checked. It is **not exercised** by the tests. `gene_exons` minus `gene_CDS` rows = 4,558,316 = the tRNA row count, so UTR exons are probably absent everywhere; not verified.
- With a filter on `transcript_id` the view takes 0.5-0.9 s (one transcript or one genome, both DuckDB versions). **Without a filter it windows over all 754 M `gene_CDS` rows. I did not run that.**
- Drop the view when `gene_introns.parquet` is rebuilt with the fix. The definitions are the same as in the fixed builder (`cds_len_upstream // 3`, `% 3`).

### 4.4 `extra_views_atlas.sql` (4 views, optional)

| View | File | Rows | Pattern |
|---|---|---|---|
| `genome_labels` | `atlas_data/labels/genome_labels.parquet` | 23,680 | `kickxellomycota/README.md:55`; adds FungalTraits `growth_form`, `primary_lifestyle`, `yeast_origin` to what `v_genome_labels` gives |
| `genome_labels_ext` | `atlas_data/labels/genome_labels_ext.parquet` | 23,680 | `k1.py:4-5` (`hs_group`, `septation`, `obligate_mode`) |
| `pfam_genome_matrix` | `atlas_data/A7/pfam_genome_matrix.parquet` (355 MB) | 104,994,592 | P4; `kickxellomycota/README.md:56`. Sorted by `genome`. F07D2B4E: `sum(n_domains)` = 15,686 = count of `pfam` rows with `domain_i_evalue <= 1e-5` for that genome |
| `genome_gene_stats` | `atlas_data/A2/genome_gene_stats.parquet` | 23,452 | `k1.py:6`; `kickxellomycota/README.md:546` |

Answer to your question on the label file: representatives, quality flags, growth form and lifestyle are in `BFD_dashboard/atlas_data/labels/genome_labels.parquet` (42 columns, 23,680 rows, written 2026-11-08 by `bin/atlas/build_labels.py`). It **can** be exposed as a view over parquet, and the view works (tested). Two facts limit it. The file is outside `tables/`, so the build script needs a second path. The file is built by a dashboard job from `tables/` plus the FungalTraits CSV, so it can be older than `tables/`. `v_genome_labels` avoids both limits for every column except the FungalTraits ones, because it is computed from `tables/` at query time and it matches the file exactly today.

---

## 5. Rejected or not done

| Idea | Decision | Reason |
|---|---|---|
| Live per-genome Pfam count view over `pfam` | Rejected | Needs a scan of all 368 M rows (12 GB) for any genome set. Use `pfam_genome_matrix` (atlas), or ask the pipeline to publish that matrix into `tables/` |
| "Genome key on every protein-level table" | Not needed as a view | Every protein-level view already has `species_prefix`. The real problem is that the computed `species_prefix` on `pfam` and the gene tables cannot prune row groups. `in_genome()` addresses that |
| Rename `species_prefix` to `LOCUSTAG` (or add `genome`) on all views | Rejected | It changes the pipeline names that `lib/bfd_data.py` and the docs use. The table comments state the equality instead |
| TargetP or DeepLoc flags in the secretome view | Rejected | `targetp` covers 125 genomes. No DeepLoc parquet is in `tables/`. The atlas secretome rule (J4) does not use them. `v_genome_coverage.has_targetp` shows where TargetP exists |
| Any view over `swissprot_hits` | Rejected | Your constraint. The existing macros cover lookups |
| `secretome_proteins_for(locustag)` table macro | Written, tested, then removed | A `species_prefix` filter on `v_secretome_proteins` already reaches all three scans in DuckDB 1.1.3 and 1.5.5 (0.1 s) |
| Replace `v_protein_annotation` with a fixed version | Not done | It is a pipeline view, and a test needs a scan of `pfam`. The two defects (section 1.4 item 14) are one-line fixes in `build_BFD_duckDB.sh`: add `WHERE feature = 'GPI-anchor'`; join `idp_summary` on `gene_id` until the idp tables are rebuilt |
| Normalize the `idp`/`idp_summary` key in a view | Not done | The builder fix exists upstream (`f27fc8e`). A view that appends `-T1` would be wrong after the rebuild |
| Taxon synonym table, order-level groups, orthogroups, repeat classes | Not done | Asked for in `kickxellomycota/README.md:57,546`, but there is no source data in `tables/` |
| Corrected `infernal_rfam` view (shift the columns back) | Not done | 14 genomes; the shift should be fixed in the merge step. The comment on the view states the defect |
| `swissprot_annot` as an alias of `swissprot_subject` | Rejected | The column sets differ (`length` vs `s_len`, no `sp_id`). A view over the parquet keeps the MotherDuck columns |
| Per-genome coverage for `pfam`, the gene tables and `idp` in `v_genome_coverage` | Left out | Needs a full scan of each. `atlas_data/coverage/*.parquet` holds only totals (for example pfam: 23,410 genomes) |

---

## 6. What was tested

Two test set-ups. Both ran inside SLURM job 31992150 (node r11, 4 CPUs), with `threads=4` and `memory_limit=12GB`.

**A. In-memory, python duckdb 1.5.5** (`test/test_inmemory.py`, log `test/test_inmemory.log`). In-memory connection; published `BFD.duckdb` attached `READ_ONLY` as `bfd`; `SET search_path='memory.main,bfd.main'`; then every statement of the three SQL files.
- 47 statements ran without error: 23 `CREATE OR REPLACE VIEW`, 1 macro, 23 `COMMENT ON VIEW` for the new views.
- 33 `COMMENT ON` statements target objects of the read-only `bfd` catalog. They cannot run in this set-up and were skipped. They were tested in set-up B.
- `SELECT * ... LIMIT 5` ran on all 23 views (`v_gene_introns_phase` with a transcript filter).
- `count(*)` ran on 20 views. Not counted: `v_secretome_proteins` (its total is `v_genome_secretome.sum(n_signalp)` = 18,767,284), `v_gene_transcripts`, `v_gene_introns_phase` (gene tables; filtered tests only).
- The comparisons with the atlas files in section 4.2 ran here.

**B. DuckDB CLI 1.1.3** (`module load duckdb`, the version the build script uses). I built a replica of the published catalog in this folder from the view and macro definitions stored in the published file (`test/replicate_pipeline_db.sql`; 38 views, `swissprot_subject` as a table, 3 macros). Then:
- `duckdb BFD.duckdb < compat_views.sql`, `< extra_views.sql`, `< extra_views_atlas.sql`, each **twice**: exit code 0 and 0 errors in all six runs (`test/run_*.log` are empty). After the runs: 61 views, 55 with a comment; `swissprot_subject` has a comment.
- `test/test_cli.sql` and `test/test_cli_transcripts.sql`: 56 queries, 0 errors (`LIMIT 5` and counts for every view, the `bfd.` and `bfd.main.` prefixes, `targetP`, the kingdom_domains representative query).
- The replica file (334 MB) was deleted after the tests. The SQL to rebuild it is kept.

**The commented-out `asm_stats` block** was run in set-up A only (uncommented, created in memory, `test/test_extra_checks.py`): 23,584 rows; the first 21 columns are identical to the pipeline view; `total_length_bp = TOTAL_LENGTH` and `gc_pct = GC_PERCENT` in all rows. It was **not** run with CLI 1.1.3.

Timings (seconds; one run each, shared node, file cache state not controlled):

| Query | CLI 1.1.3 | python 1.5.5 |
|---|---|---|
| `pfam WHERE in_genome(protein_id, 'F07D2B4E')` (17,166 rows) | 0.22 | 0.33 |
| `v_genome_labels` count | 0.08 | 0.10 |
| `v_genome_coverage` count / `LIMIT 5` | 22.5 / 13.4 | 6.2 / 6.8 |
| `v_secretome_proteins WHERE species_prefix = 'F77FFF1D'` | 0.10 | 0.08 |
| `v_genome_secretome` count | 8.0 | 4.9 |
| `v_cazy_family_counts` count / `LIMIT 5` | 11.3 / 16.4 | 1.3 / 1.9 |
| `v_tf_family_counts` count | 0.34 | 0.28 |
| `v_merops_best` count | 1.4 | 1.0 |
| `v_genome_function_counts` count | 12.9 | 8.9 |
| `v_gene_introns_phase`, one transcript | 0.53-0.73 | 0.47-0.74 |
| `v_gene_introns_phase`, one genome (`in_genome`) | 0.51 | 0.44-0.57 |
| `v_gene_transcripts`, one genome | 0.43 | not timed |
| `pfam_genome_matrix WHERE genome = ...` | 0.02 | 0.06 |

**Not tested:**
- MotherDuck itself (the MCP server did not connect). Column equality with MotherDuck is inferred from the loader.
- The files against a database built by a fresh run of `build_BFD_duckDB.sh`. The replica has the same definitions as the 2026-11-08 file; a rebuild after the current merge can differ.
- Unfiltered queries on `v_gene_introns_phase` and `v_gene_transcripts`; full scans of `pfam`, the gene tables and `idp`.
- `v_protein_annotation` and `v_species_summary` were not run.
- The `in_cds = false` branch of `v_gene_introns_phase`.
- Substitution of the `TABLES_DIR` / `ATLAS_DIR` path by a build script. The files hold literal paths; a `sed` over the file is the intended method.
- Behaviour when an optional pipeline view is missing: the matching `COMMENT ON` line then fails. In my runs all were present.

---

## 7. Decisions for a human

1. **`asm_stats` column names.** Add `total_length_bp` and `gc_pct` (and maybe `SPECIES`, `STRAIN`) to the schema-contract view, or leave MotherDuck-era SQL to be edited. In the files I read, no SQL names these columns on `bfd.asm_stats`; the uses of `total_length_bp` and `gc_pct` (62 lines in 17 files) read them from the atlas label file. The tested block is in `compat_views.sql` section 3, commented out.
2. **Atlas views.** Decide whether the database may read `BFD_dashboard/atlas_data` (run `extra_views_atlas.sql`), or whether the pipeline should publish `genome_labels`, `pfam_genome_matrix` and `genome_gene_stats` into `tables/`. Without them, growth form, lifestyle and a cheap Pfam presence matrix are not in the database.
3. **Where the two (or three) SQL files run.** They are not wired into `build_BFD_duckDB.sh` or `workflows/build_bfd_db.sh`. The wrapper's check "views lost" would not notice if they were skipped on a rebuild.
4. **`gene_trnas` alias.** I excluded `species_prefix` to match the parquet exactly. If you prefer the genome key there, change it to `SELECT * FROM gene_trna`.
5. **Pipeline fixes found on the way** (not in my files): `v_protein_annotation` GPI flag and IDP join; `infernal_rfam.parquet` column shift; `idp`/`idp_summary` key form; the sort-order statement in `BFD_Introduction.md` for the gene tables.
6. **`v_gene_introns_phase`** must be dropped after `gene_introns.parquet` is rebuilt with commit `34c8d14`. I do not know when that rebuild is planned.
