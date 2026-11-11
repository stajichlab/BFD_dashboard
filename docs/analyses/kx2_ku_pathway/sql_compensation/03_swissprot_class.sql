-- Purpose: assign each protein to a repair-gene class from its best SwissProt hit NAME (name regex), used to subdivide Pfam families.
-- Inputs: bfd.swissprot_best (protein_name of best hit; neg_log10_evalue)
-- Output: bfd_scratch.kx2_kc_sbclass (protein_id, LOCUSTAG, hit name, evalue, coverage, cls)
-- Thresholds: neg_log10_evalue >= 10 (e-value <= 1e-10). Names verified against swissprot_subject (README gene table).
-- Produces: class labels (LIG1/LIG3/LIG4, FEN1, RAD51, DMC1, RAD55/57, XRCC2/3, RAD52/59, RAD54/RDH54, MUS81, MMS4/EME1, XPF/RAD1, RAD10/ERCC1, TERT, PARP1, POLQ-name, KU-name, BRCA2).
-- S. pombe-style names are included (rhp54=Rad54, rhp55=Rad55, rhp57=Rad57, rhm52=Rad52 homolog): the best SwissProt hit of fungal Rad54 is 'DNA repair protein rhp54'.
-- Caveat: name-based subtyping is unreliable for divergent families (see README); family counts from Pfam are the primary evidence.
CREATE OR REPLACE TABLE kx2_kc_sbclass AS
SELECT protein_id, split_part(protein_id,'_',1) AS LOCUSTAG, protein_name, neg_log10_evalue, qcov_x10, scov_x10, q_len, s_len,
 CASE
  WHEN regexp_matches(n,'ku70|ku80|x-ray repair cross-complementing protein [56]|atp-dependent dna helicase 2 subunit|xrcc[56]') THEN 'KU'
  WHEN regexp_matches(n,'dna polymerase theta') THEN 'POLQ'
  WHEN regexp_matches(n,'poly \[adp-ribose\] polymerase 1') THEN 'PARP1'
  WHEN regexp_matches(n,'xrcc1') THEN 'XRCC1'
  WHEN regexp_matches(n,'^dna ligase 1') THEN 'LIG1'
  WHEN regexp_matches(n,'^dna ligase 3') THEN 'LIG3'
  WHEN regexp_matches(n,'^dna ligase 4') THEN 'LIG4'
  WHEN regexp_matches(n,'^flap endonuclease 1|^probable flap endonuclease 1') THEN 'FEN1'
  WHEN regexp_matches(n,'^exonuclease 1') THEN 'EXO1'
  WHEN regexp_matches(n,'flap endonuclease gen') THEN 'GEN1'
  WHEN regexp_matches(n,'dmc1') THEN 'DMC1'
  WHEN regexp_matches(n,'rad55|rhp55') THEN 'RAD55'
  WHEN regexp_matches(n,'xrcc2') THEN 'XRCC2'
  WHEN regexp_matches(n,'rad57|rhp57') THEN 'RAD57'
  WHEN regexp_matches(n,'xrcc3') THEN 'XRCC3'
  WHEN regexp_matches(n,'rad51 homolog [234]|rad51[bcd]$') THEN 'RAD51_PARALOG_BCD'
  WHEN regexp_matches(n,'dna repair protein rad51( homolog (1|a))?$') THEN 'RAD51'
  WHEN regexp_matches(n,'rad59') THEN 'RAD59'
  WHEN regexp_matches(n,'(recombination|repair) protein (rad52|rhm52)') THEN 'RAD52'
  WHEN regexp_matches(n,'rdh54') THEN 'RDH54'
  WHEN regexp_matches(n,'rad54b') THEN 'RAD54B'
  WHEN regexp_matches(n,'rad54|rhp54') THEN 'RAD54'
  WHEN regexp_matches(n,'breast cancer type 2') THEN 'BRCA2'
  WHEN regexp_matches(n,'mus-?81') THEN 'MUS81'
  WHEN regexp_matches(n,'mms4|eme1') THEN 'MMS4_EME1'
  WHEN regexp_matches(n,'^dna repair protein rad1$|dna repair endonuclease xpf|ercc4|^dna repair protein rad16$') THEN 'XPF_RAD1'
  WHEN regexp_matches(n,'rad10|ercc1') THEN 'RAD10_ERCC1'
  WHEN regexp_matches(n,'telomerase reverse transcriptase') THEN 'TERT'
 END AS cls
FROM (SELECT *, lower(protein_name) AS n FROM bfd.swissprot_best WHERE neg_log10_evalue>=10 AND regexp_matches(lower(protein_name),'ku70|ku80|x-ray repair cross|dna helicase 2 subunit|xrcc|dna polymerase theta|poly \[adp-ribose\] polymerase 1|dna ligase [134]|flap endonuclease|^exonuclease 1|rad5[14579]|rhp5[457]|rhm52|rad52|rdh54|dmc1|breast cancer type 2|mus-?81|mms4|eme1|^dna repair protein rad1$|xpf|ercc[14]|rad16|rad10|telomerase reverse')) x
WHERE cls IS NOT NULL;
SELECT cls, count(*), count(distinct LOCUSTAG) FROM kx2_kc_sbclass GROUP BY 1 ORDER BY 1;
