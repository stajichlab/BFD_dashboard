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
  genome_all       one row per genome in v_genome_labels (all, not only representatives): taxonomy,
                   assembly, BUSCO, and gene-structure summary
  gene_structure_genome  per-genome medians of gene, CDS, protein, exon, intron, intergenic lengths
  codon_genome     species_prefix x 61 sense-codon frequencies (stops dropped, renormalised to sum 1)
  aa_genome        species_prefix x 20 amino-acid frequencies (* and X dropped, renormalised to sum 1)
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


CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
          if a + b + c not in ("TAA", "TAG", "TGA")]
AMINO = list("ACDEFGHIKLMNPQRSTVWY")


def build_structure_tables(con, out):
    """Write genome_all, gene_structure_genome, codon_genome, aa_genome.

    Needs: v_genome_labels, gene_info, gene_proteins, gene_exons, gene_introns,
    gene_intergenic_distances, codon_freq, aa_freq. gene_proteins, gene_exons and gene_introns
    carry species_prefix (a view column on the real db). Each transcript is one gene model.
    """
    con.execute("""
        create temp table gs_prot as
        select species_prefix, count(distinct gene_id) as n_genes, median(length) as median_protein_len_aa,
               avg(length) as mean_len
        from gene_proteins group by 1
    """)
    con.execute("""
        create temp table gs_genelen as
        select LOCUSTAG as species_prefix, median("end" - start + 1) as median_gene_len_bp
        from gene_info where gene_id in (select gene_id from gene_proteins) group by 1
    """)
    con.execute("""
        create temp table gs_tx as
        select species_prefix, transcript_id, count(*) as n_exons, sum("end" - start + 1) as cds_len
        from gene_exons group by 1, 2
    """)
    con.execute("""
        create temp table gs_exon as
        select species_prefix, median(cds_len) as median_cds_len_bp, sum(cds_len) as total_exon_bp,
               median(n_exons) as median_exons_per_gene, avg((n_exons > 1)::int) as intron_frac_genes_with_intron
        from gs_tx group by 1
    """)
    con.execute("""
        create temp table gs_intron as
        select species_prefix, median("end" - start + 1) as median_intron_len_bp
        from gene_introns group by 1
    """)
    con.execute("""
        create temp table gs_inter as
        select species_prefix, median(distance) as median_intergenic_bp
        from gene_intergenic_distances group by 1
    """)
    con.execute("""
        create temp table gs as
        select p.species_prefix, p.n_genes, p.mean_len, g.median_gene_len_bp, e.median_cds_len_bp,
               p.median_protein_len_aa,
               p.n_genes / (l.total_length_bp / 1e6) as genes_per_mb,
               e.total_exon_bp / l.total_length_bp as coding_fraction,
               e.median_exons_per_gene, e.intron_frac_genes_with_intron,
               i.median_intron_len_bp, d.median_intergenic_bp
        from gs_prot p
        join v_genome_labels l on l.LOCUSTAG = p.species_prefix
        left join gs_genelen g using (species_prefix)
        left join gs_exon e using (species_prefix)
        left join gs_intron i using (species_prefix)
        left join gs_inter d using (species_prefix)
    """)
    copy(con, "gene_structure_genome", """
        select species_prefix, n_genes, median_gene_len_bp, median_cds_len_bp, median_protein_len_aa,
               genes_per_mb, coding_fraction, median_exons_per_gene, intron_frac_genes_with_intron,
               median_intron_len_bp, median_intergenic_bp
        from gs order by species_prefix
    """, out)
    copy(con, "genome_all", """
        select l.LOCUSTAG as species_prefix, l.PHYLUM, l."CLASS", l.ORDER_ as "ORDER", l.FAMILY, l.GENUS, l.SPECIES,
               l.is_representative, l.total_length_bp, l.gc_pct,
               l.busco_complete, l.busco_single, l.busco_duplicated, l.busco_fragmented, l.busco_missing,
               l.N50_bp, l.contig_count,
               g.n_genes, g.mean_len, g.genes_per_mb, g.coding_fraction, g.median_intron_len_bp, g.median_intergenic_bp
        from v_genome_labels l left join gs g on g.species_prefix = l.LOCUSTAG
        order by l.LOCUSTAG
    """, out)
    # composition: drop stops / ambiguous symbols and renormalise (source frequencies are rounded to 4 decimals)
    cc = ", ".join(f"max(frequency) filter (where codon = '{c}') / any_value(tot) as codon_{c}" for c in CODONS)
    copy(con, "codon_genome", f"""
        select species_prefix, {cc} from (
            select *, sum(frequency) over (partition by species_prefix) as tot
            from codon_freq where codon in ({",".join(f"'{c}'" for c in CODONS)})
        ) group by species_prefix order by species_prefix
    """, out)
    ac = ", ".join(f"max(frequency) filter (where amino_acid = '{a}') / any_value(tot) as aa_{a}" for a in AMINO)
    copy(con, "aa_genome", f"""
        select species_prefix, {ac} from (
            select *, sum(frequency) over (partition by species_prefix) as tot
            from aa_freq where amino_acid in ({",".join(f"'{a}'" for a in AMINO)})
        ) group by species_prefix order by species_prefix
    """, out)


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
    build_structure_tables(con, out)
    print("done", flush=True)


if __name__ == "__main__":
    main()
