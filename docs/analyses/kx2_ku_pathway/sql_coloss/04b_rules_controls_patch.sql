-- 04b_rules_controls_patch.sql
-- Purpose : revise the control rules AFTER first look at the matrix (documented, because it is tuning on the data):
--           * Rad51 control: the name rule failed in ~40% of Coemansia although PF08423 (Rad51 family) is in 110/111
--             Kickxellales-core genomes (best hit is often Dmc1/Rad57-like) -> use family-level Pfam rule, label 'family'.
--           * Pol delta control: PF14260 (zf-C4pol) is missed in ~10% of genomes -> use POL3_N + DNA_pol_B + exo1.
--           * Pol alpha (PF12254 N-terminal domain, 91/111), Rad54 (PF08658 Rad54_N: 0/4 Harpellales, 32/111 Kickx core)
--             and Rad52 show family/model-level sensitivity differences, so they are demoted from strict controls to
--             'weak controls' (is_ctrl=false, category CTRL2). Rad54_N absence in Harpellales is itself a possible
--             real change of the HR machinery (Q2 territory) and is NOT used to estimate false absence.
-- Inputs  : bfd_scratch.kx2_ku_rules (04)
-- Output  : bfd_scratch.kx2_ku_rules (updated in place); re-run 05 and 06 afterwards.
-- Reports : strict controls = Rad51 family, Pol delta, Pol epsilon, Mcm2.
UPDATE bfd_scratch.kx2_ku_rules SET gene='CTRL Rad51 family', mode='pfam', name_re=NULL, name_excl_re=NULL,
  note='Rad51/Dmc1/Rad57 family (PF08423+PF14520), family level' WHERE gene='CTRL Rad51';
UPDATE bfd_scratch.kx2_ku_rules SET mode='pfam', req_all=['PF24055','PF00136','PF03104'], name_re=NULL,
  note='POL3_N + DNA_pol_B + exo1 (delta catalytic subunit)' WHERE gene='CTRL Pol delta (Pol3)';
UPDATE bfd_scratch.kx2_ku_rules SET gene='WEAK Pol alpha (Pol1)', category='CTRL2', is_ctrl=false WHERE gene='CTRL Pol alpha (Pol1)';
UPDATE bfd_scratch.kx2_ku_rules SET gene='WEAK Rad54', category='CTRL2', is_ctrl=false WHERE gene='CTRL Rad54';
UPDATE bfd_scratch.kx2_ku_rules SET gene='WEAK Rad52', category='CTRL2', is_ctrl=false WHERE gene='CTRL Rad52';
