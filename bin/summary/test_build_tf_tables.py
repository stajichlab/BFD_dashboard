"""Tests for the per-genome TF family tables (TF survey, Task 1). Tiny parquet fixtures."""
import math
import sys
from pathlib import Path

import duckdb
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_tf_tables as t  # noqa: E402

MAP_HEADER = "pfam_acc\tfamily\tpaper_family\tnote\tpfam_id\n"


def write_map(tmp_path, body):
    p = tmp_path / "map.tsv"
    p.write_text("# comment\n\n" + MAP_HEADER + body)
    return p


def fmap_df(rows):
    return pd.DataFrame(rows, columns=["pfam_acc", "family", "paper_family", "note"])


def to_parquet(con, df, path):
    def typ(c):
        return "double" if c.startswith("domain_") else "bigint" if c == "length" else "varchar"
    sel = ", ".join(f'cast("{c}" as {typ(c)}) as "{c}"' for c in df.columns)
    con.register("tmp_df", df)
    con.execute(f"copy (select {sel} from tmp_df) to '{path}' (format parquet)")
    con.unregister("tmp_df")


@pytest.fixture()
def con():
    return duckdb.connect()


FMAP = fmap_df([
    ["PF00172", "Zn", "Zn cluster", ""],       # existing: from tf_inventory
    ["PF03106", "WRKY", "DNA-binding WRKY", ""],  # added: from pfam
    ["", "", "STAT", "unmapped"],
])


def make_tables(con, tmp_path, inv, pfam, genes, genomes):
    tdir = tmp_path / "tables"
    sdir = tmp_path / "summary"
    tdir.mkdir()
    sdir.mkdir()
    to_parquet(con, pd.DataFrame(inv, columns=["species_prefix", "protein_id", "pfam_acc", "pfam_name", "tf_family",
                                               "domain_i_evalue", "domain_score"]), tdir / "tf_inventory.parquet")
    to_parquet(con, pd.DataFrame(pfam, columns=["protein_id", "pfam_id", "pfam_acc", "domain_i_evalue", "domain_score"]),
               tdir / "pfam.parquet")
    to_parquet(con, pd.DataFrame(genes, columns=["gene_id", "protein_id", "transcript_id", "length"]),
               tdir / "gene_proteins.parquet")
    to_parquet(con, pd.DataFrame({"species_prefix": genomes}), sdir / "genome_all.parquet")
    return tdir, sdir


def run(con, tmp_path, inv, pfam, genes, genomes, fmap=FMAP):
    tdir, sdir = make_tables(con, tmp_path, inv, pfam, genes, genomes)
    out = tmp_path / "out"
    out.mkdir()
    res = t.build_tf_tables(con, tdir, sdir, fmap, out)
    fam = pd.read_parquet(out / "tf_family_genome.parquet")
    summ = pd.read_parquet(out / "tf_summary_genome.parquet").set_index("species_prefix")
    return res, fam, summ


def genes_for(prefix, n):
    return [(f"{prefix}_{i:06d}", f"{prefix}_{i:06d}-T1", f"{prefix}_{i:06d}-T1", 100) for i in range(1, n + 1)]


def test_family_map_rejects_duplicate_accessions(tmp_path):
    p = write_map(tmp_path, "PF00172\tZn\tZn cluster\t\tZn_clus\nPF00172\tZn2\tZn cluster\t\tZn_clus\n")
    with pytest.raises(ValueError, match="PF00172"):
        t.load_family_map(p)


def test_family_map_ignores_comments_and_keeps_unmapped_rows(tmp_path):
    p = write_map(tmp_path, "PF00172\tZn\tZn cluster\t\tZn_clus\n\t\tSTAT\tunmapped\t\n")
    m = t.load_family_map(p)
    assert list(m.columns) == ["pfam_acc", "family", "paper_family", "note", "pfam_id"]
    assert len(m) == 2 and m.loc[1, "pfam_acc"] == ""


def test_validate_accessions_lists_all_missing():
    fm = fmap_df([["PF00001", "a", "", ""], ["PF00002", "b", "", ""], ["PF00003", "c", "", ""], ["", "", "STAT", "unmapped"]])
    names = pd.DataFrame({"pfam_acc": ["PF00003.5"], "pfam_id": ["x"], "n_genomes": [1]})
    with pytest.raises(ValueError) as e:
        t.validate_accessions(fm, names)
    assert "PF00001" in str(e.value) and "PF00002" in str(e.value) and "PF00003" not in str(e.value)
    t.validate_accessions(fm.iloc[2:], names)  # no accession rows missing: no error


def test_tables_cover_every_genome_including_no_tf(con, tmp_path):
    inv = [("G1", "G1_000001-T1", "PF00172", "Zn_clus", "Zn", 1e-10, 50.0)]
    genes = genes_for("G1", 4) + genes_for("G2", 2)  # G3 has no proteome at all
    res, fam, summ = run(con, tmp_path, inv, [], genes, ["G1", "G2", "G3"])
    assert set(summ.index) == {"G1", "G2", "G3"}
    assert summ.loc["G1", "n_tf_proteins"] == 1 and summ.loc["G1", "tf_fraction"] == 0.25
    assert summ.loc["G2", "n_tf_proteins"] == 0 and summ.loc["G2", "tf_fraction"] == 0.0
    assert summ.loc["G3", "n_tf_proteins"] == 0 and summ.loc["G3", "n_proteins_total"] == 0
    assert math.isnan(summ.loc["G3", "tf_fraction"])
    assert set(fam["species_prefix"]) == {"G1"}  # no zero rows in the family table
    assert res["n_genomes"] == 3 and res["n_zero_tf"] == 2


def test_counts_proteins_not_domains(con, tmp_path):
    inv = [("G1", "G1_000001-T1", "PF00172", "Zn_clus", "Zn", 1e-10, 50.0),
           ("G1", "G1_000001-T1", "PF00172", "Zn_clus", "Zn", 1e-9, 40.0),
           ("G1", "G1_000002-T1", "PF00172", "Zn_clus", "Zn", 1e-9, 40.0)]
    _, fam, summ = run(con, tmp_path, inv, [], genes_for("G1", 5), ["G1"])
    assert fam.set_index("family").loc["Zn", "n_proteins"] == 2
    assert summ.loc["G1", "n_tf_proteins"] == 2


def test_species_prefix_from_protein_id(con, tmp_path):
    pfam = [("F0000001_000001-T1", "WRKY", "PF03106.20", 1e-5, 30.0),
            ("FF38F669A_000001-T1", "WRKY", "PF03106.20", 1e-5, 30.0),
            ("FF38F669A_000002-T1", "WRKY", "PF03106.20", 1e-5, 30.0)]
    genes = genes_for("F0000001", 2) + genes_for("FF38F669A", 3)
    _, fam, summ = run(con, tmp_path, [], pfam, genes, ["F0000001", "FF38F669A"])
    f = fam.set_index(["species_prefix", "family"])["n_proteins"]
    assert f[("F0000001", "WRKY")] == 1 and f[("FF38F669A", "WRKY")] == 2
    assert summ.loc["FF38F669A", "n_proteins_total"] == 3


def test_shannon_values(con, tmp_path):
    inv = [("G1", "G1_000001-T1", "PF00172", "Zn_clus", "Zn", 1e-10, 50.0),
           ("G2", "G2_000001-T1", "PF00172", "Zn_clus", "Zn", 1e-10, 50.0)]
    pfam = [("G2_000002-T1", "WRKY", "PF03106.20", 1e-5, 30.0)]
    genes = genes_for("G1", 3) + genes_for("G2", 3) + genes_for("G3", 3)
    _, fam, summ = run(con, tmp_path, inv, pfam, genes, ["G1", "G2", "G3"])
    assert summ.loc["G1", "shannon"] == 0.0 and summ.loc["G1", "n_families"] == 1
    assert summ.loc["G2", "shannon"] == pytest.approx(math.log(2)) and summ.loc["G2", "n_families"] == 2
    assert math.isnan(summ.loc["G3", "shannon"]) and summ.loc["G3", "n_families"] == 0


def test_same_gate_as_pipeline(con, tmp_path):
    """The pipeline has no e-value gate after hmmsearch --cut_ga, so no extra gate is applied.

    A weak but reported hit (the weakest in tf_inventory: i-evalue 0.087, score 15.8) is kept.
    A hit to an accession that is not in the map is excluded.
    """
    pfam = [("G1_000001-T1", "WRKY", "PF03106.20", 0.087, 15.8),
            ("G1_000002-T1", "other", "PF99999.1", 1e-30, 300.0)]
    _, fam, summ = run(con, tmp_path, [], pfam, genes_for("G1", 3), ["G1"])
    assert fam.set_index("family").loc["WRKY", "n_proteins"] == 1
    assert summ.loc["G1", "n_tf_proteins"] == 1


def test_protein_in_both_sources_counts_once_and_unlisted_genome_dropped(con, tmp_path):
    inv = [("G1", "G1_000001-T1", "PF00172", "Zn_clus", "Zn", 1e-10, 50.0),
           ("GX", "GX_000001-T1", "PF00172", "Zn_clus", "Zn", 1e-10, 50.0)]  # GX not in genome_all
    pfam = [("G1_000001-T1", "WRKY", "PF03106.20", 1e-5, 30.0)]
    res, fam, summ = run(con, tmp_path, inv, pfam, genes_for("G1", 2), ["G1"])
    assert summ.loc["G1", "n_tf_proteins"] == 1 and summ.loc["G1", "n_families"] == 2
    assert "GX" not in set(fam["species_prefix"])
    assert res["tf_genomes_not_in_genome_all"] == 1


def test_found_accessions_reported(con, tmp_path):
    inv = [("G1", "G1_000001-T1", "PF00172", "Zn_clus", "Zn", 1e-10, 50.0)]
    res, _, _ = run(con, tmp_path, inv, [], genes_for("G1", 2), ["G1"])
    assert res["accessions_with_hits"] == ["PF00172"]
    assert res["accessions_without_hits"] == ["PF03106"]


DAT = """#=GF ID   Zn_clus
#=GF AC   PF00172.24
#=GF DE   Fungal Zn(2)-Cys(6) binuclear cluster domain
//
#=GF ID   WRKY
#=GF AC   PF03106.20
//
#=GF ID   Ovate
#=GF AC   PF04844.18
//
"""


def pf_map(rows):
    return pd.DataFrame(rows, columns=["pfam_acc", "family", "paper_family", "note", "pfam_id"])


def test_check_pfam_ids_match_wrong_and_missing(tmp_path):
    dat = tmp_path / "Pfam-A.hmm.dat"
    dat.write_text(DAT)
    ok = pf_map([["PF00172", "Zn", "", "", "Zn_clus"], ["PF03106", "WRKY", "", "", "WRKY"], ["", "", "STAT", "unmapped", ""]])
    assert t.check_pfam_ids(ok, dat) == []
    bad = pf_map([["PF04844", "APSES", "", "", "DNA_binding_1"], ["PF99999", "X", "", "", "X"]])
    msgs = t.check_pfam_ids(bad, dat)
    assert len(msgs) == 2
    assert "PF04844" in msgs[0] and "Ovate" in msgs[0] and "DNA_binding_1" in msgs[0]
    assert "PF99999" in msgs[1]


REAL_DAT = Path("/srv/projects/db/pfam/current/Pfam-A.hmm.dat")


@pytest.mark.skipif(not REAL_DAT.exists(), reason="Pfam-A.hmm.dat not available on this host")
def test_committed_map_matches_pfam_release():
    fm = t.load_family_map(t.DEFAULT_MAP)
    assert t.check_pfam_ids(fm, REAL_DAT) == []


def test_excluded_accession_in_tf_inventory_is_not_counted(con, tmp_path):
    """PF08880 (QLQ) is absent from the map. Its tf_inventory rows must not count or leak a family label."""
    inv = [("G1", "G1_000001-T1", "PF08880", "Rrn7", "RNA Pol I transcription factor Rrn7-type", 1e-10, 50.0),
           ("G1", "G1_000002-T1", "PF08880", "Rrn7", "RNA Pol I transcription factor Rrn7-type", 1e-10, 50.0),
           ("G1", "G1_000002-T1", "PF00172", "Zn_clus", "Zn", 1e-10, 50.0)]
    _, fam, summ = run(con, tmp_path, inv, [], genes_for("G1", 4), ["G1"])
    assert list(fam["family"]) == ["Zn"] and list(fam["n_proteins"]) == [1]
    assert summ.loc["G1", "n_tf_proteins"] == 1 and summ.loc["G1", "n_families"] == 1


def test_committed_map_excludes_pf08880():
    fm = t.load_family_map(t.DEFAULT_MAP)
    assert "PF08880" not in set(fm["pfam_acc"])
    row = fm[fm["note"].str.contains("PF08880")]
    assert len(row) == 1 and row.iloc[0]["pfam_acc"] == "" and row.iloc[0]["pfam_id"] == ""
    assert row.iloc[0]["note"].startswith("excluded")


def cli_setup(con, tmp_path):
    inv = [("G1", "G1_000001-T1", "PF00172", "Zn_clus", "Zn", 1e-10, 50.0)]
    tdir, sdir = make_tables(con, tmp_path, inv, [], genes_for("G1", 2), ["G1"])
    pd.DataFrame({"pfam_acc": ["PF00172.24"], "pfam_id": ["Zn_clus"], "n_genomes": [1]}).to_parquet(sdir / "pfam_names.parquet")
    mp = write_map(tmp_path, "PF00172\tZn\tZn cluster\t\tZn_clus\n")
    return tdir, sdir, mp


def run_cli(tdir, sdir, mp, tmp_path, dat):
    t.main([str(tdir), str(sdir), str(tmp_path / "cli_out"), "--map", str(mp), "--pfam-dat", str(dat), "--tmp", str(tmp_path)])


def test_cli_explicit_missing_pfam_dat_fails(con, tmp_path, capsys):
    tdir, sdir, mp = cli_setup(con, tmp_path)
    with pytest.raises(SystemExit) as e:
        run_cli(tdir, sdir, mp, tmp_path, tmp_path / "typo.dat")
    assert e.value.code not in (0, None)
    assert "typo.dat" in str(e.value.code) or "typo.dat" in capsys.readouterr().err


def test_cli_default_dat_absent_notes_and_continues(con, tmp_path, capsys, monkeypatch):
    tdir, sdir, mp = cli_setup(con, tmp_path)
    missing = str(tmp_path / "absent_default.dat")
    monkeypatch.setattr(t, "DEFAULT_DAT", missing)
    t.main([str(tdir), str(sdir), str(tmp_path / "cli_out"), "--map", str(mp), "--tmp", str(tmp_path)])
    assert f"pfam_id check: skipped (dat file not found: {missing})" in capsys.readouterr().err
    assert (tmp_path / "cli_out" / "tf_summary_genome.parquet").exists()


def test_cli_none_disables_and_ok_reports(con, tmp_path, capsys):
    tdir, sdir, mp = cli_setup(con, tmp_path)
    run_cli(tdir, sdir, mp, tmp_path, "none")
    assert "pfam_id check: disabled" in capsys.readouterr().err
    dat = tmp_path / "Pfam-A.hmm.dat"
    dat.write_text(DAT)
    run_cli(tdir, sdir, mp, tmp_path, dat)
    assert "pfam_id check: ok (1 accessions)" in capsys.readouterr().err
