-- J9: intermediates for the size-transition views (HS_*: G1, G3, P2, S1, C1, C3). Output dir A11.
-- Representative genomes only (is_representative); the HS filter (hs_keep) is applied at render time.
CREATE OR REPLACE TABLE sets AS SELECT "set" AS set_name, pfam_acc FROM read_csv('{A}/labels/pfam_trait_sets.tsv', delim='\t', header=true, all_varchar=true);
CREATE OR REPLACE TABLE rep AS SELECT LOCUSTAG AS genome FROM '{A}/labels/genome_labels.parquet' WHERE is_representative;
-- C3 gene set: SwissProt sp_id -> gene (module). Presence = best SwissProt hit (A3 protein_core.sp_id) in the set.
CREATE OR REPLACE TABLE c3genes (sp_id INTEGER, gene VARCHAR, module VARCHAR, entry VARCHAR);
INSERT INTO c3genes VALUES
 (297270,'HEX1','woronin','HEX1_NEUCR'),(556249,'HEX1','woronin','HEX1_EMENI'),(398704,'HEX1','woronin','HEXA_ASPFU'),
 (575006,'WSC','woronin','WSC_NEUCR'),(398868,'WSC','woronin','WSCA_ASPFU'),
 (471618,'PRM1','fusion','PRM1_NEUCR'),(409790,'PRM1','fusion','PRM1_EMENI'),(398559,'PRM1','fusion','PRM1_ASPFU'),(25452,'PRM1','fusion','PRM1_ASPNC'),(364706,'PRM1','fusion','PRM1_ASPOR'),
 (56372,'PRM1','fusion','PRM1_BOTFB'),(408903,'PRM1','fusion','PRM1_CANAL'),(78260,'PRM1','fusion','PRM1_COPC7'),(62305,'PRM1','fusion','PRM1_SCLS1'),(273122,'PRM1','fusion','PRM1_YEAST'),(556001,'PRM1','fusion','PRM1_SCHPO'),
 (258977,'FIG1','fusion','FIG1_YEAST'),(575028,'HAM7','fusion','X325_NEUCR'),
 (240526,'FUS3_MAPK','fusion_mapk','FUS3_YEAST'),(238647,'FUS3_MAPK','fusion_mapk','KSS1_YEAST'),(187244,'FUS3_MAPK','fusion_mapk','PMK1_PYRO7'),(249837,'FUS3_MAPK','fusion_mapk','SPK1_SCHPO'),
 (247828,'STE50','fusion_mapk','STE50_YEAST'),(187270,'STE50','fusion_mapk','MST50_PYRO7'),
 (246276,'STE11','fusion_mapk','STE11_YEAST'),(187279,'STE11','fusion_mapk','MST11_PYRO7'),
 (216590,'STE7','fusion_mapk','STE7_YEAST'),(187304,'STE7','fusion_mapk','MST7_PYRO7'),(327839,'STE7','fusion_mapk','STE7_SCHPO'),
 (237703,'STE12','fusion_mapk','STE12_YEAST'),(257630,'STE12','fusion_mapk','STE11_SCHPO'),
 (304112,'SLT2_MAPK','cwi_control','SLT2_YEAST'),(398583,'SLT2_MAPK','cwi_control','MPKA_ASPFU'),(187263,'SLT2_MAPK','cwi_control','MPS1_PYRO7'),(409179,'SLT2_MAPK','cwi_control','MKC1_CANAL'),(524143,'SLT2_MAPK','cwi_control','SPM1_SCHPO'),
 (254040,'HOG1_MAPK','hog_control','HOG1_YEAST'),(530151,'HOG1_MAPK','hog_control','HOG1_NEUCR'),(555829,'HOG1_MAPK','hog_control','HOG1_EMENI'),(398648,'HOG1_MAPK','hog_control','HOG1_ASPFU'),(565099,'HOG1_MAPK','hog_control','HOG1_PYRO7'),
 (225933,'HOG1_MAPK','hog_control','HOG1_GIBZE'),(15426,'HOG1_MAPK','hog_control','HOG1_BOTFB'),(524007,'HOG1_MAPK','hog_control','HOG1_CANAL'),(316046,'HOG1_MAPK','hog_control','HOG1_SCHPO'),(25520,'HOG1_MAPK','hog_control','HOG1_ASPNC');
-- @output A11/pfam_trait_proteins.parquet
-- one row per (genome, protein, set) for representative genomes
COPY (SELECT u.genome, u.protein_id, s.set_name FROM (SELECT split_part(protein_id,'_',1) AS genome, protein_id, unnest(string_split(pfam_arch, ';')) AS acc FROM '{A}/A3/protein_pfam_arch.parquet' WHERE split_part(protein_id,'_',1) IN (SELECT genome FROM rep)) u JOIN sets s ON s.pfam_acc = u.acc GROUP BY 1,2,3) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A11/pfam_trait_counts.parquet
COPY (SELECT genome, set_name, count(*) AS n_proteins FROM '{A}/A11/pfam_trait_proteins.parquet' GROUP BY 1,2 ORDER BY 1,2) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A11/nlr_arch_counts.parquet
COPY (SELECT genome,
   count(*) FILTER (WHERE het) AS n_het,
   count(*) FILTER (WHERE nacht OR nbarc) AS n_nacht_nbarc,
   count(*) FILTER (WHERE (nacht OR nbarc) AND (ank OR tpr OR wd40)) AS n_nlr_arch,
   count(*) FILTER (WHERE het AND (nacht OR nbarc)) AS n_het_nlr,
   count(*) FILTER (WHERE het AND NOT (nacht OR nbarc)) AS n_het_only,
   count(*) FILTER (WHERE sesa) AS n_sesa, count(*) FILTER (WHERE hets) AS n_hets,
   count(*) FILTER (WHERE te) AS n_te_proteins
 FROM (SELECT genome, protein_id, bool_or(set_name='het') het, bool_or(set_name='nacht') nacht, bool_or(set_name='nbarc') nbarc, bool_or(set_name='ank') ank, bool_or(set_name='tpr') tpr, bool_or(set_name='wd40') wd40, bool_or(set_name='sesa') sesa, bool_or(set_name='hets') hets, bool_or(set_name='te') te
       FROM '{A}/A11/pfam_trait_proteins.parquet' GROUP BY 1,2)
 GROUP BY 1 ORDER BY 1) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A11/c3_gene_presence.parquet
-- per genome x gene: proteins whose best SwissProt hit is in the gene set, with relaxed gates (pident >= 25 %, qcov and scov >= 50 %) and strict gates
COPY (SELECT c.genome, g.gene, g.module,
   count(*) FILTER (WHERE c.pident_x10 >= 250 AND c.qcov_x10 >= 500 AND c.scov_x10 >= 500) AS n_relaxed,
   count(*) FILTER (WHERE c.pident_x10 >= 400 AND c.qcov_x10 >= 800 AND c.scov_x10 >= 800) AS n_strict,
   max(c.pident_x10) AS max_pident_x10, max(c.qcov_x10) AS max_qcov_x10, max(c.scov_x10) AS max_scov_x10, count(*) AS n_any
 FROM '{A}/A3/protein_core.parquet' c JOIN c3genes g USING (sp_id)
 WHERE c.genome IN (SELECT genome FROM rep) AND c.q_model
 GROUP BY 1,2,3 ORDER BY 1,2) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A11/idp_family_genome.parquet
-- P2 within-family disorder: per (representative genome, SwissProt family sp_id) median IDP fraction of complete-model proteins with strict best-hit gates
COPY (SELECT c.genome, c.sp_id, count(*) AS n_prot, median(i.IDP_fraction) AS med_idp, median(c.length) AS med_len
 FROM '{A}/A3/protein_core.parquet' c JOIN '{T}/idp_summary.parquet' i ON i.protein_id = c.gene_id
 WHERE c.genome IN (SELECT genome FROM rep) AND c.q_model AND c.sp_id IS NOT NULL AND c.pident_x10 >= 400 AND c.qcov_x10 >= 800 AND c.scov_x10 >= 800
 GROUP BY 1,2 ORDER BY 1,2) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A11/idp_genome_model.parquet
-- per representative genome: median IDP fraction over complete models, all and intracellular-only (no SwissProt hit filter; secreted set removed via A5 secreted_proteins)
COPY (SELECT c.genome, count(*) AS n_model, median(i.IDP_fraction) AS med_idp_model,
   median(i.IDP_fraction) FILTER (WHERE s.protein_id IS NULL) AS med_idp_nonsecreted,
   avg(CASE WHEN i.IDP_fraction >= 0.5 THEN 1 ELSE 0 END) AS frac_idp_ge50_model
 FROM '{A}/A3/protein_core.parquet' c JOIN '{T}/idp_summary.parquet' i ON i.protein_id = c.gene_id
 LEFT JOIN (SELECT DISTINCT protein_id FROM '{A}/A5/secreted_proteins.parquet') s ON s.protein_id = c.protein_id
 WHERE c.genome IN (SELECT genome FROM rep) AND c.q_model
 GROUP BY 1 ORDER BY 1) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
-- @output A11/c3_anchor_presence.parquet
-- C3, anchor detector: per genome x gene, proteins whose best hit INSIDE one reference proteome (A9: YEAST, NEUCR, SCHPO, EMENI) is the gene's entry.
-- Gates use reference coverage (scov) and identity only: no query-coverage gate, because HEX-1 orthologs carry N-terminal extensions (mean qcov 58 % on the 176-aa HEX1_NEUCR).
COPY (WITH g AS (SELECT * FROM (VALUES
 (297270,'HEX1','woronin'),(556249,'HEX1','woronin'),(575006,'WSC','woronin'),
 (471618,'PRM1','fusion'),(409790,'PRM1','fusion'),(273122,'PRM1','fusion'),(556001,'PRM1','fusion'),(258977,'FIG1','fusion'),(575028,'HAM7','fusion'),
 (240526,'FUS3_MAPK','fusion_mapk'),(238647,'FUS3_MAPK','fusion_mapk'),(249837,'FUS3_MAPK','fusion_mapk'),(247828,'STE50','fusion_mapk'),(246276,'STE11','fusion_mapk'),
 (216590,'STE7','fusion_mapk'),(327839,'STE7','fusion_mapk'),(237703,'STE12','fusion_mapk'),(257630,'STE12','fusion_mapk'),
 (304112,'SLT2_MAPK','cwi_control'),(524143,'SLT2_MAPK','cwi_control'),
 (254040,'HOG1_MAPK','hog_control'),(530151,'HOG1_MAPK','hog_control'),(555829,'HOG1_MAPK','hog_control'),(316046,'HOG1_MAPK','hog_control')) AS t(anchor_sp_id, gene, module))
 SELECT a.genome, g.gene, g.module, a.org_code,
   count(*) FILTER (WHERE a.pident_x10 >= 250 AND a.scov_x10 >= 500) AS n_relaxed,
   count(*) FILTER (WHERE a.pident_x10 >= 400 AND a.scov_x10 >= 800) AS n_strict,
   max(a.pident_x10) AS max_pident_x10, max(a.qcov_x10) AS max_qcov_x10, max(a.scov_x10) AS max_scov_x10, count(*) AS n_any
 FROM '{A}/A9/protein_model_anchor.parquet' a JOIN g USING (anchor_sp_id)
 WHERE a.genome IN (SELECT genome FROM rep)
 GROUP BY 1,2,3,4 ORDER BY 1,2,4) TO '{OUT}' (FORMAT PARQUET, COMPRESSION ZSTD);
