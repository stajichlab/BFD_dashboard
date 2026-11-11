-- 04_class_counts.sql
-- Purpose: per-genome functional-class protein counts (CAZy classes, TF families, secretome, Pfam-defined classes).
-- Inputs : bfd.cazy_overview, bfd.tf_inventory, bfd.signalp_signal_peptide, bfd_scratch.kw_pfam_matrix, kx2_h4_meta
-- Outputs: bfd_scratch.kx2_h4_cazy_fam (genome, fam, n), kx2_h4_tf (genome, tf_family, n), kx2_h4_sp (genome, n_sp),
--          kx2_h4_cls (genome, cls, n; zeros filled for all genomes)
-- Thresholds:
--   CAZy: dbCAN toolcount>=2, families parsed with regexp '[A-Z]{2,3}[0-9]+', distinct proteins per family
--         (a multi-domain protein counts once per family). CAZy_total excludes CBM.
--   TF: distinct proteins per tf_family; 'Fungal-specific transcription factor domain' rows excluded from TF totals.
--   Secretome: signalp probability >= 0.5, distinct proteins.
--   Pfam-defined classes (sum of n_proteins; may double count proteins with two listed Pfams):
--     Kinase PF00069,PF07714; MFS PF07690,PF00083; ABC PF00005; P450 PF00067; PKS_KS PF00109; NRPS_C PF00668;
--     AMP_ligase PF00501; SDR PF00106; AKR PF00248; MTase11 PF08241; Pept_S PF00082,PF00089,PF00450,PF00326;
--     Pept_C PF00112; Pept_A PF00026; Pept_M PF01435,PF00246,PF01433; Zn2Cys6_dom PF00172.
--   MEROPS not used as a class (bfd.merops has no family/class column). Pept_C dropped later (median 0).
-- Report numbers: all scaling classes (CAZy_*, TF_*, Secretome, Kinase, MFS, ABC, P450, SDR, AKR, Pept_*, ...).

create or replace table kx2_h4_cazy_fam as
select species_prefix genome, fam, count(distinct protein_id) n from (
 select species_prefix, protein_id, unnest(regexp_extract_all(cazyme_fam,'[A-Z]{2,3}[0-9]+')) fam
 from bfd.cazy_overview where toolcount>=2 and species_prefix in (select genome from kx2_h4_meta)) group by 1,2;

create or replace table kx2_h4_tf as
select species_prefix genome, tf_family, count(distinct protein_id) n from bfd.tf_inventory
where species_prefix in (select genome from kx2_h4_meta) group by 1,2;

create or replace table kx2_h4_sp as
select species_prefix genome, count(distinct protein_id) n_sp from bfd.signalp_signal_peptide
where species_prefix in (select genome from kx2_h4_meta) and probability>=0.5 group by 1;

create or replace table kx2_h4_cls as
with cls_def as (select * from (values
 ('PF00069','Kinase'),('PF07714','Kinase'),('PF07690','MFS'),('PF00083','MFS'),('PF00005','ABC'),('PF00067','P450'),
 ('PF00109','PKS_KS'),('PF00668','NRPS_C'),('PF00501','AMP_ligase'),('PF00106','SDR'),('PF00248','AKR'),('PF08241','MTase11'),
 ('PF00082','Pept_S'),('PF00089','Pept_S'),('PF00450','Pept_S'),('PF00326','Pept_S'),('PF00112','Pept_C'),('PF00026','Pept_A'),
 ('PF01435','Pept_M'),('PF00246','Pept_M'),('PF01433','Pept_M'),('PF00172','Zn2Cys6_dom')) t(pfam_acc,cls)),
pf as (select m.genome, d.cls, sum(m.n_proteins) n from kw_pfam_matrix m join cls_def d using(pfam_acc)
       where m.genome in (select genome from kx2_h4_meta) group by 1,2),
cz as (select genome, case when regexp_extract(fam,'^[A-Z]+')='CBM' then 'CAZy_CBM' else 'CAZy_'||regexp_extract(fam,'^[A-Z]+') end cls, sum(n) n
       from kx2_h4_cazy_fam group by 1,2),
cz_all as (select genome, 'CAZy_total' cls, sum(n) n from kx2_h4_cazy_fam where regexp_extract(fam,'^[A-Z]+')<>'CBM' group by 1),
tf as (select genome, case when tf_family like 'Zn(II)2Cys6%' then 'TF_Zn2Cys6' when tf_family like 'C2H2 zinc finger%' then 'TF_C2H2'
        when tf_family like 'Myb-like%' then 'TF_Myb' when tf_family like 'Basic leucine%' then 'TF_bZIP'
        when tf_family='Basic helix-loop-helix' then 'TF_bHLH' when tf_family='GATA zinc finger' then 'TF_GATA'
        when tf_family like 'Homeo%' then 'TF_Homeo' when tf_family='High mobility group box' then 'TF_HMG' else 'TF_other' end cls, sum(n) n
       from kx2_h4_tf where tf_family not like 'Fungal-specific%' group by 1,2),
tf_all as (select genome, 'TF_total' cls, sum(n) n from kx2_h4_tf where tf_family not like 'Fungal-specific%' group by 1),
sp as (select genome, 'Secretome' cls, n_sp n from kx2_h4_sp),
u as (select * from pf union all select * from cz union all select * from cz_all union all select * from tf union all select * from tf_all union all select * from sp),
cl as (select distinct cls from u)
select m.genome, c.cls, coalesce(u.n,0) n from kx2_h4_meta m cross join cl c left join u on u.genome=m.genome and u.cls=c.cls;
