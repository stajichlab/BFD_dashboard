"""Tests for the JGI-to-BFD map and the DAP-seq target sets (jgi_map.py, targets.py). Toy files only.

Run on the HPCC:  python -m pytest -q -W error::UserWarning bin/tf/test_tf_targets.py
The MMseqs integration test needs `module load MMseqs2/17-b804f`; it skips when mmseqs is not on PATH.
"""
import gzip
import random
import shutil
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "summary"))
import jgi_map as jm  # noqa: E402
import targets as tg  # noqa: E402

NAMES_TSV = Path(__file__).resolve().parents[2] / "config" / "paper_tf_names.tsv"


# ---------------------------------------------------------------- readers
def test_read_jgi_proteins_parses_headers(tmp_path):
    p = tmp_path / "p.aa.fasta.gz"
    with gzip.open(p, "wt") as f:
        f.write(">jgi|Rhimi59_2|101|fgenesh1_pg.1_#_1\nMKV\nLLA\n>jgi|Rhimi59_2|102|name two\nMMM\n>jgi|Rhimi59_2|103|x\nAAA\n")
    d = jm.read_jgi_proteins(p)
    assert d == {"101": ("fgenesh1_pg.1_#_1", "MKVLLA"), "102": ("name two", "MMM"), "103": ("x", "AAA")}


def test_read_catalog_gff3(tmp_path):
    p = tmp_path / "c.gff3.gz"
    with gzip.open(p, "wt") as f:
        f.write("##gff-version 3\n"
                "scaffold_1\tprediction\tgene\t10\t90\t0\t-\t.\tID=gene_1;Name=jgi.p|Rhimi59_2|555;product_name=Carbonic anhydrase;proteinId=555;transcriptId=9\n"
                "scaffold_1\tprediction\tmRNA\t10\t90\t.\t-\t.\tID=mRNA_1;Parent=gene_1;proteinId=555;track=FilteredModels1;transcriptId=9\n"
                "scaffold_1\tprediction\tCDS\t10\t90\t.\t-\t0\tID=CDS_1;Parent=mRNA_1\n"
                "scaffold_2\tprediction\tgene\t5\t50\t0\t+\t.\tID=gene_2;Name=jgi.p|Rhimi59_2|556;product_name=Hypothetical;proteinId=556\n"
                "scaffold_2\tprediction\tmRNA\t5\t50\t.\t+\t.\tID=mRNA_2;Parent=gene_2;proteinId=556;track=FilteredModels2\n")
    df = jm.read_jgi_catalog(p)
    assert list(df.columns) == ["protein_id", "scaffold", "start", "end", "strand", "product_name", "track"]
    assert df["protein_id"].tolist() == ["555", "556"]
    assert df.loc[0, "product_name"] == "Carbonic anhydrase" and df.loc[0, "track"] == "FilteredModels1"
    assert df.loc[1, "start"] == 5 and df.loc[1, "strand"] == "+"


def _xlsx(path, sheets):
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for name, ids in sheets.items():
        ws = wb.create_sheet(name)
        ws.append(["Target (JGI Gene ID)"])
        for i in ids:
            ws.append([i])
    wb.save(path)


def test_read_targets_sheet_names_and_dedup(tmp_path):
    x = tmp_path / "t.xlsx"
    _xlsx(x, {"HSF1": [1, 2, 2, 3], "GATA Znf10": [3, 4]})
    df = jm.read_targets(x)
    assert list(df.columns) == ["tf_name", "jgi_protein_id"]
    assert df[df.tf_name == "HSF1"]["jgi_protein_id"].tolist() == ["1", "2", "3"]
    assert df[df.tf_name == "GATA Znf10"]["jgi_protein_id"].tolist() == ["3", "4"]


def test_read_targets_joins_names_table(tmp_path):
    x = tmp_path / "t.xlsx"
    _xlsx(x, {"HSF1": [1], "GATA Znf10": [2], "CP2": [3]})
    names = pd.DataFrame({"jgi_id": ["11", "22", "33"], "name": ["HSF1", "GATA_Znf10", "CP2_TF1"],
                          "tf_family": ["HSF", "GATA", "CP2"], "success_frip5_any": [1, 1, 0]})
    df = jm.read_targets(x, names=names)
    d = df.set_index("tf_name")
    assert d.loc["GATA Znf10", "tf_jgi_id"] == "22" and d.loc["GATA Znf10", "tf_family"] == "GATA"
    assert d.loc["CP2", "tf_jgi_id"] == "33"


def test_read_targets_unmatched_sheet_fails(tmp_path):
    x = tmp_path / "t.xlsx"
    _xlsx(x, {"NoSuchTF": [1]})
    names = pd.DataFrame({"jgi_id": ["11"], "name": ["HSF1"], "tf_family": ["HSF"], "success_frip5_any": [1]})
    with pytest.raises(ValueError, match="NoSuchTF"):
        jm.read_targets(x, names=names)


def test_read_tested_tfs_splits_fragment_suffix(tmp_path):
    x = tmp_path / "m6.xlsx"
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["DNA sequences were reverse translated"])
    ws.append([])
    ws.append(["Protein ID", "Protein sequence", "Codon refactored DNA sequence"])
    ws.append(["Rhimi59_2_1926449a", "MKV", "ATG"])
    ws.append(["Rhimi59_2_1926449b", "LLL", "CTG"])
    ws.append(["Rhimi59_2_1762580", "AAA", "GCC"])
    wb.save(x)
    df = jm.read_tested_tfs(x)
    assert list(df.columns) == ["tested_id", "jgi_protein_id", "fragment"]
    assert df["jgi_protein_id"].tolist() == ["1926449", "1926449", "1762580"]
    assert df["fragment"].tolist() == ["a", "b", ""]


def test_read_tf_names_real_file():
    df = jm.read_tf_names(NAMES_TSV)
    assert len(df) == 88 and int(df["success_frip5_any"].sum()) == 58
    assert df["jgi_id"].map(type).eq(str).all() and df["jgi_id"].iloc[0] == "1762580"


# ---------------------------------------------------------------- mapping rate and target sets
def test_mapping_rate_counts_filtered_models_separately():
    m = pd.DataFrame({"jgi_protein_id": ["1", "2", "3", "4"], "bfd_protein_id": ["B1", "", "B3", None]})
    r = jm.mapping_rate(m, ["1", "2", "3", "9"], best_ids={"1", "2"})
    assert r["n_ids"] == 4 and r["n_mapped"] == 2 and r["rate"] == pytest.approx(0.5)
    assert r["n_ids_best"] == 2 and r["n_mapped_best"] == 1 and r["rate_best"] == pytest.approx(0.5)


def test_target_sets_keep_unmapped_and_count_them():
    t = pd.DataFrame({"tf_name": ["A", "A", "B"], "jgi_protein_id": ["1", "2", "1"]})
    m = pd.DataFrame({"jgi_protein_id": ["1"], "bfd_protein_id": ["B1"]})
    ts = tg.target_sets(t, m)
    assert list(ts.columns) == ["tf_name", "jgi_protein_id", "bfd_protein_id"]
    assert len(ts) == 3
    assert ts.loc[ts.jgi_protein_id == "2", "bfd_protein_id"].iloc[0] == ""
    c = tg.set_sizes(ts)
    assert c.set_index("tf_name").loc["A", "n_targets"] == 2 and c.set_index("tf_name").loc["A", "n_mapped"] == 1
    assert c.set_index("tf_name").loc["A", "n_unmapped"] == 1


# ---------------------------------------------------------------- conservation
TAX = pd.DataFrame({"species_prefix": ["G1", "G2", "G3", "G4"],
                    "PHYLUM": ["P1", "P1", "P2", "P2"], "ORDER": ["O1", "O1", "O2", "O3"],
                    "GENUS": ["g1", "g2", "g3", "g4"]})


def test_count_taxa_thresholds_and_zero_hits():
    hits = pd.DataFrame({"query": ["q1", "q1", "q1", "q1", "q2", "q2"],
                         "target": ["G1_1-T1", "G2_5-T1", "G3_2-T1", "G4_2-T1", "G1_1-T1", "G4_9-T1"],
                         "evalue": [1e-50, 1e-30, 1e-6, 1e-4, 1e-10, 1e-20],
                         "qcov": [0.9, 0.9, 0.5, 0.9, 0.4, 0.8]})
    out = tg.count_taxa(hits, TAX, ["q1", "q2", "q3"], min_cov=0.5, evalue=1e-5)
    d = out.set_index("bfd_protein_id")
    assert d.loc["q1", ["n_phyla", "n_orders", "n_genera"]].tolist() == [2, 2, 3]   # G4 fails e-value
    assert d.loc["q2", ["n_phyla", "n_orders", "n_genera"]].tolist() == [1, 1, 1]   # G1 fails coverage
    assert d.loc["q3", ["n_phyla", "n_orders", "n_genera"]].tolist() == [0, 0, 0]   # no hit kept
    assert out.attrs["min_cov"] == 0.5 and out.attrs["evalue"] == 1e-5


def test_conservation_counts_phyla_and_orders(tmp_path):
    hits = pd.DataFrame({"query": ["q1", "q1"], "target": ["G1_1-T1", "G3_2-T1"], "evalue": [1e-9, 1e-9], "qcov": [1.0, 1.0]})
    out = tg.conservation(["q1"], tmp_path / "db", tmp_path / "cons.tsv", tmp_path, taxonomy=TAX, hits=hits)
    assert out.loc[0, "n_phyla"] == 2 and out.loc[0, "n_orders"] == 2
    text = (tmp_path / "cons.tsv").read_text()
    assert text.startswith("# min_cov=0.5 evalue=1e-05")


# ---------------------------------------------------------------- comparison to background
def _planted():
    rng = np.random.default_rng(1)
    bg_ids = [f"B{i}" for i in range(500)]
    cons = pd.DataFrame({"bfd_protein_id": bg_ids, "n_phyla": 1, "n_genera": 5,
                         "n_orders": rng.integers(0, 20, 500)})
    ts = []
    hi = bg_ids[:60]
    cons.loc[cons.bfd_protein_id.isin(hi), "n_orders"] = rng.integers(40, 60, 60)
    for t in hi:
        ts.append(("TFhi", f"j{t}", t))
    same = [f"B{i}" for i in range(100, 160)]
    for t in same:
        ts.append(("TFmid", f"j{t}", t))
    for t in bg_ids[200:219]:
        ts.append(("TFsmall", f"j{t}", t))
    for k in range(5):
        ts.append(("TFsmall", f"u{k}", ""))
    return cons, pd.DataFrame(ts, columns=["tf_name", "jgi_protein_id", "bfd_protein_id"]), bg_ids


def test_compare_to_background_known_answer():
    cons, ts, bg = _planted()
    r = tg.compare_to_background(cons, ts, bg).set_index("tf_name")
    assert r.loc["TFhi", "p"] < 0.01 and r.loc["TFhi", "effect"] > 0.5 and r.loc["TFhi", "p_adj"] < 0.01
    assert r.loc["TFhi", "median_n_orders"] > r.loc["TFhi", "background_median_n_orders"]
    assert r.loc["TFsmall", "n_targets_bfd"] == 19
    assert np.isnan(r.loc["TFsmall", "p"]) and "fewer than 20" in r.loc["TFsmall", "reason"]
    assert r.loc["TFmid", "reason"] == ""


# ---------------------------------------------------------------- profile
def test_profile_targets_top_families_and_swissprot():
    ts = pd.DataFrame({"tf_name": ["A", "A", "A", "B"], "jgi_protein_id": list("1234"),
                       "bfd_protein_id": ["B1", "B2", "B3", "B4"]})
    pf = pd.DataFrame({"protein_id": ["B1", "B2", "B2", "B3", "B4"],
                       "pfam_id": ["HLH", "HLH", "Myc", "Zn", "HLH"],
                       "pfam_acc": ["PF00010.31", "PF00010.31", "PF99999.1", "PF00096.28", "PF00010.31"]})
    sp = pd.DataFrame({"protein_id": ["B1", "B2", "B4"], "accession": ["P1", "P1", "P2"],
                       "entry_name": ["X_YEAST", "X_YEAST", "Y_HUMAN"], "protein_name": ["Xp", "Xp", "Yp"]})
    out = tg.profile_targets(ts, pf, sp, top=2)
    a = out[(out.tf_name == "A") & (out.source == "pfam")].set_index("term")
    assert a.loc["PF00010", "n_proteins"] == 2 and a.loc["PF00010", "rank"] == 1 and a.loc["PF00010", "n_targets_bfd"] == 3
    assert len(out[(out.tf_name == "A") & (out.source == "pfam")]) == 2
    s = out[(out.tf_name == "A") & (out.source == "swissprot")]
    assert s["term"].tolist() == ["P1"] and s["n_proteins"].tolist() == [2]


# ---------------------------------------------------------------- MMseqs integration
def _rand_prot(rng, n=120):
    return "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=n))


def _write_fa(path, recs):
    path.write_text("".join(f">{k}\n{v}\n" for k, v in recs.items()))


@pytest.mark.skipif(shutil.which("mmseqs") is None, reason="mmseqs is not on PATH (module load MMseqs2/17-b804f)")
def test_map_jgi_to_bfd_rbh_and_one_way_hit(tmp_path):
    rng = random.Random(3)
    base = {f"p{i}": _rand_prot(rng) for i in range(20)}
    jgi = dict(base)
    bfd = {f"B{i}": s for i, s in zip(range(20), base.values())}
    # planted one-way hit: jgi pX is a close copy of p0, but B0's best partner is p0 (the identical one)
    jgi["pX"] = base["p0"][:100] + _rand_prot(rng, 20)
    # a jgi protein with no partner at all
    jgi["pLone"] = _rand_prot(rng, 150)
    jf, bf = tmp_path / "jgi.fa", tmp_path / "bfd.fa"
    _write_fa(jf, jgi)
    _write_fa(bf, bfd)
    scratch = tmp_path / "scr"
    scratch.mkdir()
    out = tmp_path / "map.tsv"
    m = jm.map_jgi_to_bfd(jf, bf, out, scratch, min_pident=90.0, min_cov=0.8).set_index("jgi_protein_id")
    assert m.loc["p3", "bfd_protein_id"] == "B3" and bool(m.loc["p3", "rbh"])
    assert m.loc["pLone", "bfd_protein_id"] == ""
    assert m.loc["pX", "bfd_protein_id"] == ""      # one-way hit does not count
    assert out.read_text().startswith("# min_pident=90.0 min_cov=0.8")


def test_page_tf_targets_draws_both_figures(tmp_path):
    td = tmp_path / "tables_tf"
    td.mkdir()
    rows = [(t, "Fam" + t[-1], f"j{i}", f"B{i}", 1, o, 3) for t in ("TFa", "TFb") for i, o in enumerate(range(5, 30))]
    pd.DataFrame(rows, columns=["tf_name", "tf_family", "jgi_protein_id", "bfd_protein_id", "n_phyla", "n_orders", "n_genera"]).to_csv(
        td / "T21_target_conservation.tsv", sep="\t", index=False)
    pd.DataFrame({"tf_name": ["TFa", "TFb"], "n_targets_bfd": [25, 25], "median_n_orders": [17, 17], "background_median_n_orders": [10, 10],
                  "effect": [.1, .1], "p": [.01, .01], "p_adj": [.01, .01], "reason": ["", ""]}).to_csv(
        td / "T21_target_background_tests.tsv", sep="\t", index=False)
    pd.DataFrame({"tf_name": ["TFa", "TFb"], "n_targets": [40, 30], "n_mapped": [25, 25], "n_unmapped": [15, 5]}).to_csv(
        td / "T21_tf_set_sizes.tsv", sep="\t", index=False)
    made = tg.page_tf_targets(tmp_path, td)
    assert [m.kind for m in made] == ["figure", "figure"]
    assert (tmp_path / "figures/tf/T21_target_conservation.png").exists()
    assert (tmp_path / "figures/tf/T21_target_set_sizes.png").exists()
