"""Fixture test for the per-genome structure and composition tables (Task 2)."""
import itertools
import sys
from pathlib import Path

import duckdb
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_summary_tables as b  # noqa: E402

BASES = "ACGT"
STOPS = {"TAA", "TAG", "TGA"}
AAS = "ACDEFGHIKLMNPQRSTVWY"


@pytest.fixture()
def con():
    c = duckdb.connect()
    # three genomes: G1 and G2 annotated, G2 has no introns, G3 is not representative and has no genes
    c.execute("""create table v_genome_labels as select * from (values
      ('G1','Asco','Sord','Sordariales','Fam1','Gen1','Gen1 a', true,  2000000::bigint, 50.0, 95.0, 90.0, 5.0, 3.0, 2.0, 100000::bigint, 10::bigint),
      ('G2','Asco','Sord','Sordariales','Fam1','Gen1','Gen1 b', false, 1000000::bigint, 48.0, 90.0, 85.0, 5.0, 4.0, 6.0,  50000::bigint, 20::bigint),
      ('G3','Basi','Agar','Agaricales','Fam2','Gen2','Gen2 a', true,  4000000::bigint, 40.0, 80.0, 75.0, 5.0, 8.0, 12.0, 70000::bigint, 30::bigint)
    ) t(LOCUSTAG, PHYLUM, "CLASS", ORDER_, FAMILY, GENUS, SPECIES, is_representative, total_length_bp, gc_pct,
        busco_complete, busco_single, busco_duplicated, busco_fragmented, busco_missing, N50_bp, contig_count)""")
    c.execute("create table gene_info (gene_id varchar, LOCUSTAG varchar, start bigint, \"end\" bigint)")
    c.execute("create table gene_proteins (gene_id varchar, protein_id varchar, transcript_id varchar, length bigint, species_prefix varchar)")
    c.execute("create table gene_exons (exon_id varchar, transcript_id varchar, \"order\" bigint, start bigint, \"end\" bigint, species_prefix varchar)")
    c.execute("create table gene_introns (intron_id varchar, transcript_id varchar, start bigint, \"end\" bigint, species_prefix varchar)")
    c.execute("create table gene_intergenic_distances (species_prefix varchar, left_gene varchar, right_gene varchar, distance bigint)")
    # G1: 4 genes. g1 has 3 exons, 2 introns; g2 has 2 exons, 1 intron; g3, g4 single exon
    genes = {"G1": [(1, 1000, 3, 300), (1, 5000, 2, 200), (1, 9000, 1, 100), (1, 12000, 1, 100)],
             "G2": [(1, 500, 1, 90), (1, 3000, 1, 90)]}
    for sp, lst in genes.items():
        for i, (_, start, nex, exlen) in enumerate(lst, 1):
            gid = f"{sp}_{i:06d}"
            tid = gid + "-T1"
            c.execute("insert into gene_info values (?,?,?,?)", [gid, sp, start, start + exlen * 3])
            c.execute("insert into gene_proteins values (?,?,?,?,?)", [gid, tid, tid, exlen // 3, sp])
            per = exlen // nex
            for k in range(nex):
                c.execute("insert into gene_exons values (?,?,?,?,?,?)",
                          [f"{tid}.exon{k+1}", tid, k, start + k * 1000, start + k * 1000 + per - 1, sp])
            for k in range(nex - 1):
                c.execute("insert into gene_introns values (?,?,?,?,?)",
                          [f"{tid}.intron{k}", tid, start, start + 49 + k * 10 - 0, sp])  # lengths 50, 60
    c.execute("insert into gene_intergenic_distances values ('G1','a','b',100),('G1','b','c',300),('G1','c','d',200),('G2','a','b',50)")
    c.execute("create table codon_freq (species_prefix varchar, codon varchar, frequency double)")
    c.execute("create table aa_freq (species_prefix varchar, amino_acid varchar, frequency double)")
    codons = [''.join(p) for p in itertools.product(BASES, repeat=3)]
    for sp in ("G1", "G2"):
        for i, cd in enumerate(codons):
            c.execute("insert into codon_freq values (?,?,?)", [sp, cd, round((i + 1) / 2080, 4)])
        for i, a in enumerate(AAS + "*X"):
            c.execute("insert into aa_freq values (?,?,?)", [sp, a, round((i + 1) / 253, 4)])
    return c


@pytest.fixture()
def built(con, tmp_path):
    b.build_structure_tables(con, tmp_path)
    return tmp_path


def q(con, path, sql):
    return con.execute(sql.replace("T", f"'{path}'")).fetchall()


def test_genome_all_has_every_genome(con, built):
    n = con.execute("select count(*) from v_genome_labels").fetchone()[0]
    rows = q(con, built / "genome_all.parquet", "select count(*) from T")[0][0]
    assert rows == n == 3
    cols = [r[0] for r in con.execute(f"describe select * from '{built}/genome_all.parquet'").fetchall()]
    for want in ["species_prefix", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES", "is_representative",
                 "total_length_bp", "gc_pct", "busco_complete", "busco_single", "busco_duplicated", "busco_fragmented",
                 "busco_missing", "N50_bp", "contig_count", "n_genes", "mean_len", "genes_per_mb", "coding_fraction",
                 "median_intron_len_bp", "median_intergenic_bp"]:
        assert want in cols
    # non-representative genome kept
    assert q(con, built / "genome_all.parquet", "select is_representative from T where species_prefix='G2'")[0][0] is False


def test_structure_one_row_per_genome_with_genes(con, built):
    p = built / "gene_structure_genome.parquet"
    assert sorted(r[0] for r in q(con, p, "select species_prefix from T")) == ["G1", "G2"]
    n, gpm = q(con, p, "select n_genes, genes_per_mb from T where species_prefix='G1'")[0]
    assert n == 4
    assert gpm == pytest.approx(4 / (2000000 / 1e6))
    cols = [r[0] for r in con.execute(f"describe select * from '{p}'").fetchall()]
    assert cols == ["species_prefix", "n_genes", "median_gene_len_bp", "median_cds_len_bp", "median_protein_len_aa",
                    "genes_per_mb", "coding_fraction", "median_exons_per_gene", "intron_frac_genes_with_intron",
                    "median_intron_len_bp", "median_intergenic_bp"]


def test_structure_values(con, built):
    p = built / "gene_structure_genome.parquet"
    r = q(con, p, """select median_exons_per_gene, intron_frac_genes_with_intron, median_intron_len_bp,
                     median_intergenic_bp, median_cds_len_bp, coding_fraction from T where species_prefix='G1'""")[0]
    assert r[0] == pytest.approx(1.5)       # exons per gene 3,2,1,1 -> median 1.5
    assert r[1] == pytest.approx(0.5)       # 2 of 4 genes have an intron
    assert r[2] == pytest.approx(50.0)      # intron lengths 50,50,51 -> median 50
    assert r[3] == pytest.approx(200.0)
    assert r[4] == pytest.approx(150.0)     # cds 300,200,100,100
    assert r[5] == pytest.approx(700 / 2000000)


def test_no_introns_gives_null(con, built):
    r = q(con, built / "gene_structure_genome.parquet",
          "select median_intron_len_bp, intron_frac_genes_with_intron from T where species_prefix='G2'")[0]
    assert r[0] is None
    assert r[1] == 0.0
    r = q(con, built / "genome_all.parquet", "select median_intron_len_bp from T where species_prefix='G2'")[0]
    assert r[0] is None


def test_genome_without_genes_is_in_genome_all_with_null(con, built):
    r = q(con, built / "genome_all.parquet", "select n_genes, genes_per_mb from T where species_prefix='G3'")[0]
    assert r == (None, None)


def test_codon_rows_sum_to_one(con, built):
    p = built / "codon_genome.parquet"
    cols = [r[0] for r in con.execute(f"describe select * from '{p}'").fetchall()]
    assert len(cols) == 62 and cols[0] == "species_prefix"
    assert "codon_AAA" in cols and "codon_TAA" not in cols
    rows = con.execute(f"select * from '{p}'").fetchall()
    assert len(rows) == 2
    for row in rows:
        assert sum(row[1:]) == pytest.approx(1.0, abs=1e-6)


def test_aa_rows_sum_to_one(con, built):
    p = built / "aa_genome.parquet"
    cols = [r[0] for r in con.execute(f"describe select * from '{p}'").fetchall()]
    assert len(cols) == 21 and "aa_A" in cols and "aa_X" not in cols
    for row in con.execute(f"select * from '{p}'").fetchall():
        assert sum(row[1:]) == pytest.approx(1.0, abs=1e-6)
