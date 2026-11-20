#!/usr/bin/env python3
"""Build the compact tables for the DeepLoc/IDP analysis from the local BFD.duckdb (HPCC job, read-only on the database).

Usage: 01_build_cube.py OUT_DIR [REPO_ROOT]
Unit: one representative genome per genus (v_genome_labels.is_genus_rep, annotated).
Writes to OUT_DIR:
  reps.parquet        species_prefix, PHYLUM, CLASS, GENUS, gc_pct
  cube.parquet        species_prefix x length bin x primary DeepLoc label x flags -> n, n_row, sum_idr
                      flags: sec (strict secretome), tf, has_pfam, unann (no Pfam, no CAZy, no MEROPS), mito_carrier (PF00153),
                      effector_fam (CFEM PF05730, LysM PF01476, CAP PF00188, Chitin_bind_1 PF00187)
  ssp_genome.parquet  species_prefix -> n_ssp (effector-like candidates, rule in README), n_proteins
  ssp_candidates.parquet  species_prefix, protein_id, length (the candidates themselves)
Domain gate for Pfam, CAZy and MEROPS: e-value <= 1e-5. IDR measure: idp_summary.IDP_fraction (no row = 0 in sum_idr).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE.parents[2]
sys.path.insert(0, str(REPO / "lib"))

LABELS = ["cytoplasm", "nucleus", "extracellular", "cell_membrane", "mitochondrion", "plastid", "endoplasmic_reticulum",
          "lysosome_vacuole", "golgi", "peroxisome"]
EFFECTOR_FAMILIES = ("PF05730", "PF01476", "PF00188", "PF00187")


def argmax_case():
    cols = ", ".join(f"p_{c}" for c in LABELS)
    return "case " + " ".join(f"when greatest({cols}) = p_{c} then '{c}'" for c in LABELS) + " end"


def main():
    import bfd_connect
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    con = bfd_connect.connect(scratch=False, threads=16, memory="80GB")
    con.execute("""create temp table reps as
        select LOCUSTAG as species_prefix, PHYLUM, CLASS, GENUS, gc_pct from v_genome_labels
        where is_genus_rep and annotation_status = 'annotated'""")
    n = con.execute("select count(*) from reps").fetchone()[0]
    print("genus representatives:", n, flush=True)
    con.execute("""create temp table p as
        select g.species_prefix, g.protein_id, g.length from gene_proteins g join reps using (species_prefix)""")
    print("proteins:", con.execute("select count(*) from p").fetchone()[0], flush=True)
    con.execute(f"""create temp table loc as
        select d.species_prefix, d.protein_id, {argmax_case()} as loc
        from deeploc d join reps r using (species_prefix)""")
    con.execute("""create temp table idr as
        select i.species_prefix, i.protein_id, i.IDP_fraction as f from idp_summary i join reps using (species_prefix)""")
    con.execute("""create temp table sec as
        select distinct s.species_prefix, s.protein_id from v_secretome_proteins s join reps using (species_prefix)
        where s.secreted_strict""")
    con.execute("""create temp table tf as
        select distinct t.species_prefix, t.protein_id from tf_inventory t join reps using (species_prefix)""")
    con.execute("""create temp table pf as
        select distinct x.species_prefix, x.protein_id, split_part(x.pfam_acc, '.', 1) as acc from pfam x join reps using (species_prefix)
        where x.domain_i_evalue <= 1e-5""")
    con.execute("create temp table haspf as select distinct species_prefix, protein_id from pf")
    con.execute("""create temp table mito as select distinct species_prefix, protein_id from pf where acc = 'PF00153'""")
    fam = ",".join(f"'{a}'" for a in EFFECTOR_FAMILIES)
    con.execute(f"create temp table eff as select distinct species_prefix, protein_id from pf where acc in ({fam})")
    con.execute("""create temp table cz as select distinct c.species_prefix, c.protein_id from cazy c join reps using (species_prefix)
        where c.evalue <= 1e-5""")
    con.execute("""create temp table mp as select distinct m.species_prefix, m.protein_id from merops m join reps using (species_prefix)
        where m.evalue <= 1e-5""")
    print("temp tables built", flush=True)
    con.execute("""create temp table prot as
        select p.species_prefix, p.protein_id, p.length,
               case when p.length < 100 then '<100' when p.length < 300 then '100-299' when p.length < 1000 then '300-999' else '>=1000' end as lbin,
               coalesce(l.loc, 'none') as loc,
               (s.protein_id is not null) as sec, (t.protein_id is not null) as tf,
               (h.protein_id is not null) as has_pfam,
               (h.protein_id is null and c.protein_id is null and m.protein_id is null) as unann,
               (mi.protein_id is not null) as mito_carrier, (e.protein_id is not null) as effector_fam,
               (i.protein_id is not null) as has_row, coalesce(i.f, 0.0) as idr
        from p left join loc l using (species_prefix, protein_id) left join idr i using (species_prefix, protein_id)
        left join sec s using (species_prefix, protein_id) left join tf t using (species_prefix, protein_id)
        left join haspf h using (species_prefix, protein_id) left join mito mi using (species_prefix, protein_id)
        left join eff e using (species_prefix, protein_id) left join cz c using (species_prefix, protein_id)
        left join mp m using (species_prefix, protein_id)""")
    con.execute(f"""copy (select * from reps) to '{out}/reps.parquet' (format parquet, compression zstd)""")
    con.execute(f"""copy (select species_prefix, lbin, loc, sec, tf, has_pfam, unann, mito_carrier, effector_fam,
                            count(*) as n, sum(has_row::int) as n_row, sum(idr) as sum_idr
                     from prot group by all) to '{out}/cube.parquet' (format parquet, compression zstd)""")
    con.execute(f"""copy (select species_prefix, count(*) filter (where sec and length <= 300 and unann) as n_ssp, count(*) as n_proteins
                     from prot group by 1) to '{out}/ssp_genome.parquet' (format parquet, compression zstd)""")
    con.execute(f"""copy (select species_prefix, protein_id, length from prot where sec and length <= 300 and unann)
                    to '{out}/ssp_candidates.parquet' (format parquet, compression zstd)""")
    for f in ("reps", "cube", "ssp_genome", "ssp_candidates"):
        print(f, con.execute(f"select count(*) from '{out}/{f}.parquet'").fetchone()[0], "rows", flush=True)


if __name__ == "__main__":
    main()
