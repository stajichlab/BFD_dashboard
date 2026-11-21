#!/usr/bin/env python3
"""POST HOC diagnosis of why control C4 failed (43% of effector-family proteins pass secreted_strict; 70% required).

Usage: 03_c4_diagnosis.py OUT_DIR [REPO_ROOT]
Not pre-registered. Read-only on BFD.duckdb. Same proteins as C4: genus-representative genomes (v_genome_labels.is_genus_rep,
annotated), Pfam domain_i_evalue <= 1e-5, families CFEM PF05730, LysM PF01476, CAP PF00188, Chitin_bind_1 PF00187.
Each protein gets ONE reason, in this order:
  no_signalp      no SignalP row (no signal peptide called)
  sp_below_0.9    SignalP probability < 0.9
  tm_gt1          SignalP >= 0.9 but more than 1 TMHMM helix
  gpi             SignalP >= 0.9, <= 1 helix, GPI-anchor predicted
  pass_strict     all of the above clear (= secreted_strict, the C4 filter)
Also writes the DeepLoc primary label for each family and a looser pass rate (secreted_noTM: GPI allowed).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE.parents[2]
sys.path.insert(0, str(REPO / "lib"))
FAM = {"PF05730": "CFEM", "PF01476": "LysM", "PF00188": "CAP", "PF00187": "Chitin_bind_1"}
LABELS = ["cytoplasm", "nucleus", "extracellular", "cell_membrane", "mitochondrion", "plastid", "endoplasmic_reticulum",
          "lysosome_vacuole", "golgi", "peroxisome"]


def main():
    import bfd_connect
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    con = bfd_connect.connect(scratch=False, threads=8, memory="48GB")
    con.execute("""create temp table reps as select LOCUSTAG as species_prefix from v_genome_labels
                   where is_genus_rep and annotation_status = 'annotated'""")
    acc = ",".join(f"'{a}'" for a in FAM)
    con.execute(f"""create temp table eff as
        select distinct x.species_prefix, x.protein_id, split_part(x.pfam_acc, '.', 1) as acc
        from pfam x join reps using (species_prefix)
        where x.domain_i_evalue <= 1e-5 and split_part(x.pfam_acc, '.', 1) in ({acc})""")
    cols = ", ".join(f"p_{c}" for c in LABELS)
    case = "case " + " ".join(f"when greatest({cols}) = p_{c} then '{c}'" for c in LABELS) + " end"
    con.execute(f"""create temp table dl as select d.species_prefix, d.protein_id, {case} as loc
                    from deeploc d join (select distinct species_prefix, protein_id from eff) e using (species_prefix, protein_id)""")
    con.execute("""create temp table prot as
        select e.species_prefix, e.protein_id, e.acc,
               s.probability as sp_prob, coalesce(t.PredHel, 0) as tm, (g.protein_id is not null) as gpi, dl.loc,
               g2.length
        from eff e
        left join signalp s using (species_prefix, protein_id)
        left join tmhmm t using (species_prefix, protein_id)
        left join (select distinct species_prefix, protein_id from predgpi where feature = 'GPI-anchor') g using (species_prefix, protein_id)
        left join dl using (species_prefix, protein_id)
        left join gene_proteins g2 using (species_prefix, protein_id)""")
    con.execute("""create temp table r as select *,
        case when sp_prob is null then 'no_signalp' when sp_prob < 0.9 then 'sp_below_0.9' when tm > 1 then 'tm_gt1'
             when gpi then 'gpi' else 'pass_strict' end as reason from prot""")
    fam_case = "case acc " + " ".join(f"when '{a}' then '{n}'" for a, n in FAM.items()) + " end"
    reason = con.execute(f"""select {fam_case} as family, reason, count(*) as n from r group by all order by 1, 2""").df()
    tot = reason.groupby("family").n.sum().rename("family_total")
    reason = reason.join(tot, on="family")
    reason["share"] = reason.n / reason.family_total
    allr = con.execute("select 'ALL' as family, reason, count(*) as n from r group by all order by 2").df()
    allr["family_total"] = allr.n.sum()
    allr["share"] = allr.n / allr.family_total
    import pandas as pd
    reason = pd.concat([reason, allr])
    reason.to_csv(out / "c4_reason_by_family.tsv", sep="\t", index=False)
    loose = con.execute(f"""select {fam_case} as family, count(*) as n,
        -- coalesce: a protein with no SignalP row is a fail, not a NULL that avg() skips
        avg(coalesce(sp_prob >= 0.9 and tm <= 1 and not gpi, false)::int) as strict,
        avg(coalesce(sp_prob >= 0.9 and tm <= 1, false)::int) as noTM_gpi_allowed,
        avg(coalesce(sp_prob >= 0.9, false)::int) as signalp_ge_0_9,
        avg((sp_prob is not null)::int) as any_signalp,
        avg(coalesce(loc in ('extracellular','cell_membrane'), false)::int) as deeploc_extra_or_membrane,
        avg(coalesce(loc = 'extracellular', false)::int) as deeploc_extracellular,
        median(length) as median_length
        from r group by all order by 1""").df()
    loose.to_csv(out / "c4_pass_rates_by_family.tsv", sep="\t", index=False)
    dlt = con.execute(f"select {fam_case} as family, coalesce(loc,'none') as loc, count(*) as n from r group by all order by 1, 3 desc").df()
    dlt.to_csv(out / "c4_deeploc_by_family.tsv", sep="\t", index=False)
    print(reason.to_string(index=False))
    print(loose.to_string(index=False))


if __name__ == "__main__":
    main()
