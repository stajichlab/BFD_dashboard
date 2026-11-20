#!/usr/bin/env python3
"""Build the tables for T4 and T6 from the local BFD.duckdb (HPCC job, read-only).

Usage: 01_build_tables.py OUT_DIR [REPO_ROOT]
Unit: one representative genome per genus (is_genus_rep, annotated). Writes OUT_DIR/reps.parquet (species_prefix, GENUS, CLASS, ORDER_,
PHYLUM, n_proteins), OUT_DIR/copies.parquet (species_prefix, acc, n = proteins carrying the family's domain; domain_i_evalue <= 1e-5,
version stripped) for every family present in at least MIN_PREV of the genus representatives plus the module families, and
OUT_DIR/names.parquet (acc, pfam_id, prevalence).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE.parents[3]
sys.path.insert(0, str(REPO / "lib"))

MIN_PREV = 0.90
MODULES = ["PF00982", "PF02358", "PF01204", "PF03632", "PF05693", "PF00343", "PF06202", "PF01704", "PF02878", "PF01210",
           "PF00082", "PF00704"]


def main():
    import bfd_connect
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    con = bfd_connect.connect(scratch=False, threads=8, memory="48GB")
    con.execute("""create temp table reps as
        select v.LOCUSTAG as species_prefix, v.GENUS, v.CLASS, v.ORDER_, v.PHYLUM
        from v_genome_labels v where v.is_genus_rep and v.annotation_status = 'annotated'""")
    con.execute("""create temp table np as
        select g.species_prefix, count(*) as n_proteins from gene_proteins g join reps using (species_prefix) group by 1""")
    con.execute("""create temp table cp as
        select p.species_prefix, split_part(p.pfam_acc, '.', 1) as acc, any_value(p.pfam_id) as pfam_id, count(distinct p.protein_id) as n
        from pfam p join reps r using (species_prefix) where p.domain_i_evalue <= 1e-5 group by 1, 2""")
    n_reps = con.execute("select count(*) from reps join np using (species_prefix)").fetchone()[0]
    mods = ",".join(f"'{m}'" for m in MODULES)
    con.execute(f"""create temp table fam as
        select acc, any_value(pfam_id) as pfam_id, count(*)::double / {n_reps} as prevalence from cp group by 1
        having count(*)::double / {n_reps} >= {MIN_PREV} or acc in ({mods})""")
    con.execute(f"copy (select r.*, np.n_proteins from reps r join np using (species_prefix)) to '{out}/reps.parquet' (format parquet)")
    con.execute(f"copy (select species_prefix, acc, n from cp where acc in (select acc from fam)) to '{out}/copies.parquet' (format parquet)")
    con.execute(f"copy (select * from fam) to '{out}/names.parquet' (format parquet)")
    for f in ("reps", "copies", "names"):
        print(f, con.execute(f"select count(*) from '{out}/{f}.parquet'").fetchone()[0], "rows", flush=True)
    print("near-universal families (prevalence >= 0.90):",
          con.execute(f"select count(*) from fam where prevalence >= {MIN_PREV}").fetchone()[0])


if __name__ == "__main__":
    main()
