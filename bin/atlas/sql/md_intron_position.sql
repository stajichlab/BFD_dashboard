-- Intron position queries for the MotherDuck copy of BFD (database `bfd`).
-- Results and hypotheses: docs/summary_atlas/HS_I_intron_position.md
--
-- Run each query on its own (MotherDuck MCP `query`, the web UI, or
-- `duckdb md:bfd`). The MCP has a 55 s limit: run them one at a time.
--
-- WHY gene_CDS AND NOT gene_introns.codon_position
--   In the 2026-11-06 snapshot gene_introns.codon_position / codon_frame are
--   off by one exon: they hold the CDS offset at the END of the exon that
--   follows the intron (so the last intron of every transcript reports the
--   full CDS length and phase 0). Fixed in Fungi_BFD
--   nextflow/bin/build_genestats_table.py; until gene_introns is rebuilt,
--   positions are recomputed here from gene_CDS:
--     cum      = coding bases 5' of the intron (running sum of CDS segment lengths)
--     phase    = cum % 3            (0 between codons, 1 after base 1, 2 after base 2)
--     rel_pos  = cum / CDS_length   (0 = start codon, 1 = stop codon)
--     aa_pos   = cum // 3 + 1       (residue the intron falls in or before)
--   gene_CDS."order" is 0-based in transcript (5' to 3') order on both strands.
--
-- SHARED DEFINITIONS (repeated in each query; MotherDuck read-only cannot keep views)
--   rep : one genome per species (best BUSCO complete, then N50), annotated,
--         BUSCO complete >= 25th percentile of its CLASS (atlas Q-ann, Q-rep, Q-busco-rel).
--   tx  : complete models only (NOT is_partial, CDS_length % 3 = 0) (atlas Q-model).
--   Unit of replication = genus (median of representative genomes), as in the atlas.
--   Growth-form labels (FungalTraits) are not in MotherDuck; groups here are taxonomic.

-- =====================================================================
-- I1a. Intron position, phase and length by class
--      frac_first10pct under a uniform distribution along the CDS = 0.10
-- =====================================================================
WITH q AS (
  SELECT s.LOCUSTAG AS tag, s.PHYLUM, s.CLASS, s.GENUS, s.SPECIES, b.complete_pct, a.total_length_bp,
         quantile_cont(b.complete_pct, 0.25) OVER (PARTITION BY s.CLASS) AS class_p25,
         row_number() OVER (PARTITION BY s.SPECIES ORDER BY b.complete_pct DESC NULLS LAST, a.N50_bp DESC NULLS LAST) AS rk
  FROM species s JOIN busco_genome b USING (ASMID) JOIN asm_stats a USING (ASMID)
  WHERE a.annotation_status = 'annotated'),
rep AS (SELECT * FROM q WHERE rk = 1 AND complete_pct >= class_p25),
tx AS (
  SELECT transcript_id, CDS_length FROM gene_transcripts
  WHERE NOT is_partial AND CDS_length % 3 = 0 AND split_part(transcript_id, '_', 1) IN (SELECT tag FROM rep)),
c AS (
  SELECT transcript_id, "order" AS ord,
         sum("end" - start + 1) OVER w AS cum,
         count(*) OVER (PARTITION BY transcript_id) AS n_cds,
         greatest(lead(start) OVER w - "end", start - lead("end") OVER w) - 1 AS intron_len
  FROM gene_CDS WHERE split_part(transcript_id, '_', 1) IN (SELECT tag FROM rep)
  WINDOW w AS (PARTITION BY transcript_id ORDER BY "order")),
ip AS (  -- one row per CDS intron
  SELECT split_part(c.transcript_id, '_', 1) AS tag, cum, cum % 3 AS phase,
         cum / tx.CDS_length AS rel_pos, intron_len
  FROM c JOIN tx USING (transcript_id) WHERE ord < n_cds - 1),
gi AS (
  SELECT tag, count(*) AS n_introns,
         avg((phase = 0)::INT) AS ph0, avg((phase = 1)::INT) AS ph1,
         avg((rel_pos < 0.1)::INT) AS first10, avg((rel_pos >= 0.9)::INT) AS last10,
         avg((cum <= 30)::INT) AS within10codons, median(rel_pos) AS med_rel,
         median(intron_len) FILTER (rel_pos < 0.1) AS len_first10,
         median(intron_len) FILTER (rel_pos >= 0.1) AS len_rest
  FROM ip GROUP BY 1),
gt AS (SELECT split_part(transcript_id, '_', 1) AS tag, count(*) AS n_tx FROM tx GROUP BY 1),
gs AS (  -- one row per representative genome
  SELECT rep.*, gt.n_tx, coalesce(gi.n_introns, 0) / gt.n_tx AS ipt, gi.* EXCLUDE (tag, n_introns)
  FROM rep JOIN gt USING (tag) LEFT JOIN gi USING (tag)),
gen AS (  -- genus medians
  SELECT PHYLUM, CLASS, GENUS, median(ipt) AS ipt, median(total_length_bp) / 1e6 AS mb,
         median(first10) AS first10, median(last10) AS last10, median(within10codons) AS w10,
         median(med_rel) AS med_rel, median(ph0) AS ph0, median(ph1) AS ph1,
         median(len_first10 / len_rest) AS len_ratio, median(len_first10) AS len_first10, median(len_rest) AS len_rest
  FROM gs GROUP BY ALL)
SELECT PHYLUM, CLASS, count(*) AS n_genera,
       round(median(mb), 1) AS genome_mb,
       round(median(ipt), 3) AS introns_per_tx,
       round(median(first10), 3) AS frac_first10pct,
       round(median(last10), 3) AS frac_last10pct,
       round(median(w10), 3) AS frac_within_10_codons,
       round(median(med_rel), 3) AS median_rel_pos,
       round(median(ph0), 3) AS phase0, round(median(ph1), 3) AS phase1,
       round(median(len_first10), 0) AS len_5p, round(median(len_rest), 0) AS len_rest,
       round(median(len_ratio), 2) AS len_ratio_5p_rest
FROM gen GROUP BY ALL HAVING count(*) >= 3
ORDER BY introns_per_tx;

-- =====================================================================
-- I1b. Within-clade Spearman: does 5' bias rise as intron density (or
--      genome size) falls? Genus medians; genomes with >= 50 CDS introns.
-- =====================================================================
WITH q AS (
  SELECT s.LOCUSTAG AS tag, s.PHYLUM, s.CLASS, s.GENUS, s.SPECIES, b.complete_pct, a.total_length_bp,
         quantile_cont(b.complete_pct, 0.25) OVER (PARTITION BY s.CLASS) AS class_p25,
         row_number() OVER (PARTITION BY s.SPECIES ORDER BY b.complete_pct DESC NULLS LAST, a.N50_bp DESC NULLS LAST) AS rk
  FROM species s JOIN busco_genome b USING (ASMID) JOIN asm_stats a USING (ASMID)
  WHERE a.annotation_status = 'annotated'),
rep AS (
  SELECT *, CASE
      WHEN PHYLUM = 'Microsporidia' THEN 'Microsporidia'
      WHEN CLASS IN ('Saccharomycetes','Pichiomycetes','Dipodascomycetes','Trigonopsidomycetes','Lipomycetes',
                     'Alloascoideomycetes','Sporopachydermiomycetes','Phaffomycetes','Alaninomycetes') THEN 'Saccharomycotina'
      WHEN CLASS IN ('Schizosaccharomycetes','Pneumocystomycetes','Taphrinomycetes','Neolectomycetes',
                     'Archaeorhizomycetes','Novakomycetes') THEN 'Taphrinomycotina'
      WHEN PHYLUM = 'Ascomycota' THEN 'Pezizomycotina'
      WHEN CLASS IN ('Ustilaginomycetes','Exobasidiomycetes','Malasseziomycetes','Moniliellomycetes') THEN 'Ustilaginomycotina'
      WHEN PHYLUM = 'Basidiomycota' THEN 'Basidiomycota_other'
      ELSE 'early_diverging' END AS grp
  FROM q WHERE rk = 1 AND complete_pct >= class_p25),
tx AS (
  SELECT transcript_id, CDS_length FROM gene_transcripts
  WHERE NOT is_partial AND CDS_length % 3 = 0 AND split_part(transcript_id, '_', 1) IN (SELECT tag FROM rep)),
c AS (
  SELECT transcript_id, "order" AS ord, sum("end" - start + 1) OVER w AS cum,
         count(*) OVER (PARTITION BY transcript_id) AS n_cds
  FROM gene_CDS WHERE split_part(transcript_id, '_', 1) IN (SELECT tag FROM rep)
  WINDOW w AS (PARTITION BY transcript_id ORDER BY "order")),
gi AS (
  SELECT split_part(c.transcript_id, '_', 1) AS tag, count(*) AS n_introns,
         avg((cum / tx.CDS_length < 0.1)::INT) AS first10, avg((cum % 3 = 0)::INT) AS ph0
  FROM c JOIN tx USING (transcript_id) WHERE ord < n_cds - 1 GROUP BY 1),
gt AS (SELECT split_part(transcript_id, '_', 1) AS tag, count(*) AS n_tx FROM tx GROUP BY 1),
gs AS (
  SELECT rep.grp, rep.GENUS, rep.total_length_bp, gi.n_introns / gt.n_tx AS ipt, gi.first10, gi.ph0
  FROM rep JOIN gt USING (tag) JOIN gi USING (tag) WHERE gi.n_introns >= 50),
gen AS (
  SELECT grp, GENUS, median(ipt) AS ipt, median(first10) AS first10, median(ph0) AS ph0, median(total_length_bp) AS len
  FROM gs GROUP BY ALL),
r AS (
  SELECT *, rank() OVER (PARTITION BY grp ORDER BY ipt) AS r_ipt, rank() OVER (PARTITION BY grp ORDER BY first10) AS r_f,
         rank() OVER (PARTITION BY grp ORDER BY len) AS r_len, rank() OVER (PARTITION BY grp ORDER BY ph0) AS r_ph
  FROM gen)
SELECT grp, count(*) AS n_genera,
       round(median(ipt), 2) AS med_introns_per_tx,
       round(median(first10), 3) AS med_first10,
       round(quantile_cont(first10, 0.1), 3) AS p10_first10,
       round(quantile_cont(first10, 0.9), 3) AS p90_first10,
       round(corr(r_ipt, r_f), 3) AS rho_density_vs_first10,
       round(corr(r_len, r_f), 3) AS rho_genomesize_vs_first10,
       round(corr(r_len, r_ipt), 3) AS rho_genomesize_vs_density,
       round(corr(r_ipt, r_ph), 3) AS rho_density_vs_phase0
FROM r GROUP BY 1 ORDER BY med_introns_per_tx;

-- =====================================================================
-- I4. Which genes keep introns: cytosolic ribosomal-protein (RP) genes
--     vs all other genes, by class. RP = best Swiss-Prot hit named
--     "Large|Small ribosomal subunit protein [eu][LS]n" (mitochondrial
--     "...m" and bacterial-type "bL/bS" names are excluded by the regex).
--     Genomes with < 40 RP genes found are dropped.
-- =====================================================================
WITH q AS (
  SELECT s.LOCUSTAG AS tag, s.PHYLUM, s.CLASS, s.GENUS, s.SPECIES, b.complete_pct,
         quantile_cont(b.complete_pct, 0.25) OVER (PARTITION BY s.CLASS) AS class_p25,
         row_number() OVER (PARTITION BY s.SPECIES ORDER BY b.complete_pct DESC NULLS LAST, a.N50_bp DESC NULLS LAST) AS rk
  FROM species s JOIN busco_genome b USING (ASMID) JOIN asm_stats a USING (ASMID)
  WHERE a.annotation_status = 'annotated'),
rep AS (SELECT * FROM q WHERE rk = 1 AND complete_pct >= class_p25),
tx AS (
  SELECT transcript_id FROM gene_transcripts
  WHERE NOT is_partial AND CDS_length % 3 = 0 AND split_part(transcript_id, '_', 1) IN (SELECT tag FROM rep)),
nc AS (
  SELECT transcript_id, count(*) - 1 AS n_introns FROM gene_CDS
  WHERE split_part(transcript_id, '_', 1) IN (SELECT tag FROM rep) GROUP BY 1),
rp AS (
  SELECT DISTINCT protein_id AS transcript_id FROM swissprot_best   -- protein_id = transcript_id
  WHERE species_prefix IN (SELECT tag FROM rep)
    AND regexp_matches(protein_name, '^(Large|Small) ribosomal subunit protein [eu][LS][0-9]+[A-Z]?$')),
tt AS (
  SELECT split_part(tx.transcript_id, '_', 1) AS tag, (rp.transcript_id IS NOT NULL) AS is_rp, nc.n_introns
  FROM tx JOIN nc USING (transcript_id) LEFT JOIN rp USING (transcript_id)),
gs AS (
  SELECT tag, count(*) FILTER (is_rp) AS n_rp,
         avg((n_introns > 0)::INT) FILTER (is_rp) AS rp_w,
         avg((n_introns > 0)::INT) FILTER (NOT is_rp) AS other_w,
         sum(n_introns) FILTER (is_rp) / nullif(sum(n_introns), 0) AS share_rp,
         avg(n_introns) AS ipt
  FROM tt GROUP BY 1 HAVING count(*) FILTER (is_rp) >= 40),
gen AS (
  SELECT rep.PHYLUM, rep.CLASS, rep.GENUS, median(ipt) AS ipt, median(rp_w) AS rp_w, median(other_w) AS other_w,
         median(share_rp) AS share_rp, median(n_rp) AS n_rp
  FROM gs JOIN rep USING (tag) GROUP BY ALL)
SELECT PHYLUM, CLASS, count(*) AS n_genera,
       round(median(ipt), 2) AS introns_per_tx,
       round(median(n_rp), 0) AS n_rp_genes,
       round(median(rp_w), 3) AS rp_with_intron,
       round(median(other_w), 3) AS other_with_intron,
       round(median(rp_w / nullif(other_w, 0)), 2) AS rp_enrichment,
       round(median(share_rp), 4) AS share_of_introns_in_rp
FROM gen GROUP BY ALL HAVING count(*) >= 5
ORDER BY introns_per_tx;

-- =====================================================================
-- I5 (pilot, 10 genomes). Intron position relative to protein features:
--     Pfam domain envelopes and disordered regions (idp), observed vs
--     expected under a uniform position within the same protein; RP
--     intron length and 5' position. Edit the VALUES list to change genomes.
--     Note: idp.protein_id has no '-T1' suffix; pfam / swissprot_best do.
--     To scale: replace the per-residue unnest (domcov) with interval
--     merging, and run on the HPCC parquet (pfam is 368 M rows).
-- =====================================================================
WITH g(tag, label) AS (VALUES
  ('F5A7DA31','E.cuniculi'), ('FB2702AD','Y.lipolytica'), ('F591A93A','B.dendrobatidis'),
  ('F0200A73','S.cerevisiae'), ('F0371FFB','P.jirovecii'), ('F56E78F8','C.neoformans'),
  ('FC69C3D3','N.crassa'), ('F4E3794F','R.delemar'), ('FFECC5E4','S.pombe'), ('F1F5452B','C.cinerea')),
tx AS (
  SELECT transcript_id, CDS_length, CDS_length / 3 - 1 AS aa_len FROM gene_transcripts
  WHERE split_part(transcript_id, '_', 1) IN (SELECT tag FROM g) AND NOT is_partial AND CDS_length % 3 = 0),
c AS (
  SELECT transcript_id, "order" AS ord, sum("end" - start + 1) OVER w AS cum,
         count(*) OVER (PARTITION BY transcript_id) AS n_cds,
         greatest(lead(start) OVER w - "end", start - lead("end") OVER w) - 1 AS intron_len
  FROM gene_CDS WHERE split_part(transcript_id, '_', 1) IN (SELECT tag FROM g)
  WINDOW w AS (PARTITION BY transcript_id ORDER BY "order")),
rp AS (
  SELECT DISTINCT protein_id AS transcript_id FROM swissprot_best
  WHERE species_prefix IN (SELECT tag FROM g)
    AND regexp_matches(protein_name, '^(Large|Small) ribosomal subunit protein [eu][LS][0-9]+[A-Z]?$')),
tt AS (
  SELECT split_part(tx.transcript_id, '_', 1) AS tag, tx.transcript_id, (rp.transcript_id IS NOT NULL) AS is_rp,
         max(c.n_cds) - 1 AS n_introns
  FROM tx JOIN c USING (transcript_id) LEFT JOIN rp USING (transcript_id) GROUP BY ALL),
ip AS (
  SELECT tt.tag, tt.is_rp, c.cum % 3 AS phase, c.cum / tx.CDS_length AS rel_pos, c.intron_len,
         c.cum // 3 + 1 AS aa_pos, tx.aa_len, tx.transcript_id
  FROM c JOIN tx USING (transcript_id) JOIN tt USING (transcript_id) WHERE c.ord < c.n_cds - 1),
dom AS (
  SELECT protein_id AS transcript_id, env_from AS s, env_to AS e FROM pfam
  WHERE split_part(protein_id, '_', 1) IN (SELECT tag FROM g) AND domain_i_evalue < 1e-3),
domcov AS (  -- residues covered by >= 1 Pfam envelope
  SELECT transcript_id, count(DISTINCT r) AS cov
  FROM (SELECT transcript_id, unnest(range(s, e + 1)) AS r FROM dom) GROUP BY 1),
idr AS (
  SELECT protein_id || '-T1' AS transcript_id, IDP_start AS s, IDP_end AS e, IDP_length AS l FROM idp
  WHERE species_prefix IN (SELECT tag FROM g)),
idrcov AS (SELECT transcript_id, sum(l) AS cov FROM idr GROUP BY 1),
x AS (
  SELECT ip.*,
         EXISTS (SELECT 1 FROM dom d WHERE d.transcript_id = ip.transcript_id AND ip.aa_pos BETWEEN d.s AND d.e) AS in_dom,
         EXISTS (SELECT 1 FROM idr d WHERE d.transcript_id = ip.transcript_id AND ip.aa_pos BETWEEN d.s AND d.e) AS in_idr,
         coalesce(dc.cov, 0) / ip.aa_len AS dom_frac, coalesce(ic.cov, 0) / ip.aa_len AS idr_frac
  FROM ip LEFT JOIN domcov dc USING (transcript_id) LEFT JOIN idrcov ic USING (transcript_id)),
a AS (
  SELECT tag, count(*) FILTER (is_rp) AS n_rp,
         round(avg((n_introns > 0)::INT) FILTER (is_rp), 3) AS rp_with_intron,
         round(avg((n_introns > 0)::INT) FILTER (NOT is_rp), 3) AS other_with_intron
  FROM tt GROUP BY 1),
b AS (
  SELECT tag, count(*) AS n_introns,
         round(avg(in_dom::INT), 3) AS obs_in_pfam, round(avg(dom_frac), 3) AS exp_in_pfam,
         round(avg(in_dom::INT) / nullif(avg(dom_frac), 0), 3) AS pfam_obs_exp,
         round(avg(in_idr::INT), 3) AS obs_in_idr, round(avg(idr_frac), 3) AS exp_in_idr,
         round(avg(in_idr::INT) / nullif(avg(idr_frac), 0), 3) AS idr_obs_exp,
         round(avg((phase = 1)::INT) FILTER (in_dom), 3) AS ph1_in_dom,
         round(avg((phase = 1)::INT) FILTER (NOT in_dom), 3) AS ph1_out_dom,
         round(avg(is_rp::INT), 3) AS share_introns_in_rp,
         median(intron_len) FILTER (is_rp) AS rp_intron_len,
         median(intron_len) FILTER (NOT is_rp) AS other_intron_len,
         round(avg((rel_pos < 0.1)::INT) FILTER (is_rp), 3) AS rp_first10,
         round(avg((rel_pos < 0.1)::INT) FILTER (NOT is_rp), 3) AS other_first10
  FROM x GROUP BY 1)
SELECT g.label, a.* EXCLUDE (tag), b.* EXCLUDE (tag)
FROM g JOIN a USING (tag) JOIN b USING (tag)
ORDER BY other_with_intron;

-- =====================================================================
-- I0. Check of the stored columns against the recomputed position
--     (one genome). After gene_introns is rebuilt with the fix, n_match
--     should equal n; in the 2026-11-06 snapshot it does not.
-- =====================================================================
WITH c AS (
  SELECT transcript_id, "order" AS ord, sum("end" - start + 1) OVER (PARTITION BY transcript_id ORDER BY "order") AS cum
  FROM gene_CDS WHERE transcript_id LIKE 'F0200A73\_%' ESCAPE '\')
SELECT count(*) AS n,
       count(*) FILTER (i.codon_position = c.cum // 3 AND i.codon_frame = c.cum % 3) AS n_match
FROM gene_introns i JOIN c ON c.transcript_id = i.transcript_id AND c.ord = i.intron_number
WHERE i.transcript_id LIKE 'F0200A73\_%' ESCAPE '\';
