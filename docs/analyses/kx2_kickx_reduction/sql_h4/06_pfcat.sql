-- 06_pfcat.sql
-- Purpose: map each Pfam to a broad functional category from GO terms of SwissProt entries carrying that Pfam
--          (keyword argmax; crude, NOT a GO slim or ontology walk).
-- Inputs : bfd.swissprot_subject (pfam, go_terms)
-- Output : bfd_scratch.kx2_h4_pfcat (pfam_acc, n_sp = number of SwissProt entries, cat, k = matches in best cat, tot)
-- Method : for each Pfam, concatenate lower-cased go_terms over all SwissProt entries listing it; count regex matches
--          per category; assign the category with the most matches (ties: row_number order, arbitrary); 'unclassified'
--          when no category matched. Only Pfams present in SwissProt subjects are mapped (others become 'unclassified'
--          in downstream joins).
-- GO keyword list (regex per category) is in the cats CTE below.
-- Report numbers: functional-category enrichment of losses (step 07).

create or replace table kx2_h4_pfcat as
with e as (select unnest(string_split(pfam,'; ')) p, go_terms from bfd.swissprot_subject where pfam is not null),
e2 as (select regexp_extract(p,'^(PF[0-9]+)',1) pfam_acc, lower(coalesce(go_terms,'')) go from e),
s as (select pfam_acc, count(*) n_sp, string_agg(go,' ') g from e2 where pfam_acc<>'' group by 1),
cats as (select * from (values
 ('DNA/chromatin/repair','dna repair|dna recombination|dna replication|telomere|chromatin|histone|dna damage'),
 ('RNA/translation','spliceos|mrna processing|rna splicing|trna|rrna|ribosom|translation|rna helicase|rna processing'),
 ('transcription','regulation of transcription|dna-binding transcription|transcription factor|transcription by rna polymerase'),
 ('mitochondrial/respiration','mitochondri|respiratory|electron transport|atp synth|oxidative phosph|nadh dehydrogenase'),
 ('signalling','protein kinase|phosphorylation|signal transduction|gtpase|g protein|receptor'),
 ('transport','transmembrane transport|transporter|ion transport|channel'),
 ('carbohydrate/cell wall','carbohydrate|polysaccharide|cell wall|chitin|glucan|cellulose|glycosidase|glycosyl'),
 ('lipid','fatty acid|lipid|sterol|ergosterol|phospholipid'),
 ('amino acid/N metabolism','amino acid|nitrogen|methionine|arginine|lysine|tryptophan|cysteine'),
 ('nucleotide metabolism','nucleotide|purine|pyrimidine|nucleoside'),
 ('cofactor/vitamin','vitamin|cofactor|coenzyme|heme|thiamine|riboflavin|folate|biotin|molybdenum'),
 ('oxidoreductase/secondary/catabolic','secondary metabolite|antibiotic|oxidoreductase|monooxygenase|oxidation-reduction|catabolic'),
 ('proteostasis/protein turnover','proteolysis|peptidase|ubiquitin|protein folding|chaperone|proteasome|protein import|protein targeting'),
 ('trafficking/cytoskeleton/cell cycle','vesicle|endocytosis|golgi|cytoskeleton|actin|microtubule|autophagy|cell cycle|mitosis|cytokinesis')) c(cat,rx)),
sc as (select s.pfam_acc, s.n_sp, c.cat, len(regexp_extract_all(s.g, c.rx)) k from s cross join cats c),
best as (select pfam_acc, n_sp, cat, k, row_number() over (partition by pfam_acc order by k desc) rn, sum(k) over (partition by pfam_acc) tot from sc)
select pfam_acc, n_sp, case when tot=0 then 'unclassified' else cat end cat, k, tot from best where rn=1;
