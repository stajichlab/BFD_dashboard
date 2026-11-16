#!/usr/bin/env python3
"""Aggregate the local BFD.duckdb into small per-genome tables for the Quarto summary.

Run on the HPCC (see workflows/summary/build_summary.sh). Reads only; writes parquet files to OUT.
Unit: one representative genome per species (v_genome_labels.is_representative, annotated).

Outputs (one parquet each):
  genome_core      one row per representative genome: taxonomy, size, GC, BUSCO, gene count,
                   protein length, Pfam coverage, secretome / CAZy / MEROPS / TF counts
  pfam_genome      species_prefix x pfam_acc: proteins with the domain, domain copies
                   (domain_i_evalue <= 1e-5, version stripped; families in the top PREV_TOP by
                   prevalence among representatives, plus the curated panel in panel.py)
  pfam_names       pfam_acc -> pfam_id (short name)
  arch_genome      per genome: distribution of domains per protein (0,1,2,3,4,5+)
  deeploc_genome   species_prefix x primary DeepLoc label (argmax of p_*), protein counts
  wolfpsort_genome species_prefix x WoLF PSORT localization, protein counts
  idp_genome       per genome: proteins in idp_summary, mean IDP fraction, share with > 0.3
  cazy_genome      species_prefix x cazy_class x cazy_family, proteins
  tf_genome        species_prefix x tf_family, proteins
  merops_genome    species_prefix x merops class (first letter of merops_id), proteins
"""
import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bfd_connect  # noqa: E402
from panel import PANEL_ACCESSIONS  # noqa: E402

PREV_TOP = 1500
DEEPLOC_COLS = ["cytoplasm", "nucleus", "extracellular", "cell_membrane", "mitochondrion", "plastid",
                "endoplasmic_reticulum", "lysosome_vacuole", "golgi", "peroxisome"]


def copy(con, name, sql, out):
    path = Path(out) / f"{name}.parquet"
    con.execute(f"COPY ({sql}) TO '{path}' (FORMAT PARQUET, COMPRESSION ZSTD)")
    n = con.execute(f"select count(*) from '{path}'").fetchone()[0]
    print(f"{name}: {n} rows", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--memory", default="40GB")
    ap.add_argument("--db", default=None)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tmp = os.environ.get("SCRATCH")
    con = bfd_connect.connect(a.db, scratch=False, threads=a.threads, memory=a.memory, tmp=tmp)
    # the db is attached read-only as 'bfd'; temp tables live in memory
    con.execute("""
        create temp table reps as
        select v.LOCUSTAG as species_prefix, v.PHYLUM, sp.SUBPHYLUM, v.CLASS, v.ORDER_ as ORDER, v.FAMILY, v.GENUS, v.SPECIES,
               v.total_length_bp, v.gc_pct, v.busco_complete, v.N50_bp, v.contig_count
        from v_genome_labels v
        left join (select LOCUSTAG, any_value(SUBPHYLUM) as SUBPHYLUM from species group by 1) sp on sp.LOCUSTAG = v.LOCUSTAG
        where v.is_representative
    """)
    print("representatives:", con.execute("select count(*) from reps").fetchone()[0], flush=True)

    copy(con, "genome_core", """
        select r.*, g.n_genes, g.mean_len, f.n_signalp, f.n_signalp_09, f.n_secreted_noTM, f.n_secreted_strict,
               f.n_cazy, f.n_cazy_secreted, f.n_merops, f.n_merops_secreted, f.n_tf_proteins
        from reps r
        left join (select species_prefix, count(*) as n_genes, avg(length) as mean_len
                   from gene_proteins group by 1) g using (species_prefix)
        left join v_genome_function_counts f on f.LOCUSTAG = r.species_prefix
    """, out)

    # Pfam: gate and strip the version as in bin/atlas/sql/J2_pfam.sql
    con.execute("""
        create temp table pf as
        select p.species_prefix, split_part(p.pfam_acc, '.', 1) as pfam_acc, p.protein_id, p.pfam_id
        from pfam p join reps r using (species_prefix)
        where p.domain_i_evalue <= 1e-5
    """)
    con.execute("""
        create temp table prev as
        select pfam_acc, count(distinct species_prefix) as n_genomes from pf group by 1
    """)
    panel = ",".join(f"'{x}'" for x in PANEL_ACCESSIONS)
    con.execute(f"""
        create temp table keep as
        select pfam_acc from (select pfam_acc from prev order by n_genomes desc limit {PREV_TOP})
        union select pfam_acc from prev where pfam_acc in ({panel})
    """)
    copy(con, "pfam_genome", """
        select species_prefix, pfam_acc, count(distinct protein_id) as n_proteins, count(*) as n_domains
        from pf where pfam_acc in (select pfam_acc from keep) group by 1, 2
    """, out)
    copy(con, "pfam_names", "select pfam_acc, any_value(pfam_id) as pfam_id, any_value(n_genomes) as n_genomes "
         "from (select pf.pfam_acc, pf.pfam_id, prev.n_genomes from pf join prev using (pfam_acc)) group by 1", out)
    copy(con, "arch_genome", """
        with pp as (select species_prefix, protein_id, count(*) as nd from pf group by 1, 2)
        select species_prefix, least(nd, 5) as n_domains, count(*) as n_proteins from pp group by 1, 2
    """, out)

    cols = ", ".join(f"p_{c}" for c in DEEPLOC_COLS)
    cases = " ".join(f"when greatest({cols}) = p_{c} then '{c}'" for c in DEEPLOC_COLS)
    copy(con, "deeploc_genome", f"""
        select d.species_prefix, case {cases} end as label, count(*) as n_proteins
        from deeploc d join reps r using (species_prefix) group by 1, 2
    """, out)
    copy(con, "wolfpsort_genome", """
        select w.species_prefix, w.localization, count(*) as n_proteins
        from wolfpsort w join reps r using (species_prefix) group by 1, 2
    """, out)
    copy(con, "idp_genome", """
        select i.species_prefix, count(*) as n_proteins_idp_summary, avg(IDP_fraction) as mean_idp_fraction,
               avg((IDP_fraction > 0.3)::int) as frac_idp_gt_0_3
        from idp_summary i join reps r using (species_prefix) group by 1
    """, out)
    copy(con, "cazy_genome", """
        select c.species_prefix, c.cazy_class, c.cazy_family, sum(c.n_proteins) as n_proteins
        from v_cazy_family_counts c join reps r using (species_prefix) group by 1, 2, 3
    """, out)
    copy(con, "tf_genome", """
        select t.species_prefix, t.tf_family, sum(t.n_proteins) as n_proteins
        from v_tf_family_counts t join reps r using (species_prefix) group by 1, 2
    """, out)
    copy(con, "merops_genome", """
        select m.species_prefix, left(m.merops_id, 1) as merops_class, count(distinct m.protein_id) as n_proteins
        from v_merops_best m join reps r using (species_prefix) group by 1, 2
    """, out)
    print("done", flush=True)


if __name__ == "__main__":
    main()
