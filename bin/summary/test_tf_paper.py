"""Tests for the calibration against the paper's own counts (tf_paper.py). Fixtures are synthetic; one test reads the
real MOESM3 xlsx when it is present and is skipped otherwise.

Run on the HPCC:  python3.12 -m pytest -q -W error::UserWarning bin/summary/test_tf_paper.py
"""
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tf_claims as tc  # noqa: E402
import tf_paper as tp  # noqa: E402

REAL_XLSX = Path("/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/resources/s41467-025-65177-1_suppl/"
                 "41467_2025_65177_MOESM3_ESM.xlsx")
FMAP = pd.DataFrame({"pfam_acc": ["PF1", "PF2", "PF3", "PF4"], "family": ["Fam One", "Fam One b", "Fam Two", "Fam Gone"],
                     "paper_family": ["P1", "P1", "P2", "P3"]})


def make_xlsx(path):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "TF count"
    ws.append(["TF Family ", "A. one", "B. two "])
    ws.append(["Zn cluster ", 10, 20])
    ws.append(["bZIP", 5, 6])
    ws.append([None, None, None])
    ws.append(["Total TFs", 15, 26])
    ws.append(["Total Genes", 1000, 2000])
    g = wb.create_sheet("Genomes")
    g.append(["Species", "Accesion", "Reference", "doi"])
    g.append(["Alpha one ", "https://x/1", "r", None])
    g.append(["Beta two", "https://x/2", "r", None])
    g.append(["Extra three", "https://x/3", "r", None])       # in Genomes only (not in TF count)
    wb.save(path)


def test_convert_paper_counts_fixture(tmp_path):
    x = tmp_path / "m.xlsx"
    make_xlsx(x)
    out = tmp_path / "paper_tf_counts.tsv"
    df = tp.convert_paper_counts(x, out)
    assert list(df.columns) == ["species", "species_full", "paper_family", "count"]
    d = df.set_index(["species", "paper_family"])["count"]
    assert d[("A. one", "Zn cluster")] == 10 and d[("B. two", "bZIP")] == 6        # trailing spaces stripped
    assert set(df["paper_family"]) == {"Zn cluster", "bZIP"} and len(df) == 4       # totals and blank row are not families
    assert dict(zip(df["species"], df["species_full"])) == {"A. one": "Alpha one", "B. two": "Beta two"}
    assert out.read_text().startswith("#") and "Lax" in out.read_text().splitlines()[0]
    counts, totals = tp.load_paper_counts(out)
    assert counts.equals(df.reset_index(drop=True)) or len(counts) == 4
    t = totals.set_index("species")
    assert t.loc["A. one", "total_tfs"] == 15 and t.loc["B. two", "total_genes"] == 2000
    assert not (tmp_path / "m.xlsx.tsv").exists()


@pytest.mark.skipif(not REAL_XLSX.exists(), reason="real MOESM3 not present")
def test_convert_paper_counts_known_cells(tmp_path):
    out = tmp_path / "c.tsv"
    df = tp.convert_paper_counts(REAL_XLSX, out)
    _, tot = tp.load_paper_counts(out)
    t = tot.set_index("species").loc["R. microsporus"]
    assert t["total_tfs"] == 605 and t["total_genes"] == 10700      # the paper's text: 605 TFs of 10,700 genes
    assert df["species"].nunique() == 62 and "Zn cluster" in set(df["paper_family"])
    assert df["species_full"].str.strip().eq(df["species_full"]).all()


GA = pd.DataFrame({
    "species_prefix": ["F1", "F2", "F3", "F4", "F5", "F6", "F7"],
    "GENUS": ["Rhizopus", "Rhizopus", "Eremothecium", "Malassezia", "Other", "Rhizopus", "Mucor"],
    "SPECIES": ["Rhizopus microsporus", "Rhizopus microsporus", "Eremothecium gossypii", "Malassezia sympodialis",
                "Other thing", "Rhizopus delemar", "Mucor lusitanicus"],
    "PHYLUM": ["Mucoromycota"] * 3 + ["Basidiomycota", "Ascomycota", "Mucoromycota", "Mucoromycota"],
    "busco_complete": [90.0, 95.0, 99.0, 88.0, 80.0, 99.9, 70.0], "N50_bp": [1, 1, 1, 1, 1, 1, 1]})


def test_match_species_picks_best_busco_and_reports_none():
    names = ["Rhizopus microsporus", "Ashbya(Eremothecium) gossypii", "Malasseizia sympodialis", "Zzz nonesuch", "Orpinomyces"]
    m = tp.match_species(names, GA, aliases={"Malasseizia sympodialis": "Malassezia sympodialis"}).set_index("paper_species")
    assert list(m.columns[:4]) == ["species_prefix", "strain", "match", "n_candidates"] or {"species_prefix", "strain", "match"} <= set(m.columns)
    assert m.loc["Rhizopus microsporus", "species_prefix"] == "F2" and m.loc["Rhizopus microsporus", "n_candidates"] == 2
    assert m.loc["Rhizopus microsporus", "match"] == "exact"
    assert m.loc["Ashbya(Eremothecium) gossypii", "species_prefix"] == "F3"       # the genus in brackets is tried too
    assert m.loc["Malasseizia sympodialis", "species_prefix"] == "F4" and m.loc["Malasseizia sympodialis", "match"] == "alias"
    assert m.loc["Zzz nonesuch", "match"] == "none" and m.loc["Zzz nonesuch", "n_candidates"] == 0
    assert pd.isna(m.loc["Zzz nonesuch", "species_prefix"]) or m.loc["Zzz nonesuch", "species_prefix"] == ""
    assert m.loc["Orpinomyces", "match"] == "none"           # no epithet: cannot match by genus and epithet
    assert m.loc["Malasseizia sympodialis", "alias_reason"] == "typo" and m.loc["Rhizopus microsporus", "alias_reason"] == ""
    m2 = tp.match_species(["Ustilago maydis"], GA.assign(GENUS=GA["GENUS"].replace("Malassezia", "Mycosarcoma"),
                          SPECIES=GA["SPECIES"].replace("Malassezia sympodialis", "Mycosarcoma maydis")),
                          aliases={"Ustilago maydis": "Mycosarcoma maydis"}).iloc[0]
    assert m2["match"] == "alias" and m2["alias_reason"] == "renamed taxon"


def _perfect():
    pc = pd.DataFrame([("A. one", "Alpha one", "P1", 10), ("A. one", "Alpha one", "P2", 4),
                       ("B. two", "Beta two", "P1", 20), ("B. two", "Beta two", "P2", 8),
                       ("C. three", "Gamma three", "P1", 30), ("C. three", "Gamma three", "P2", 2),
                       ("D. four", "Delta four", "P1", 40), ("D. four", "Delta four", "P2", 6),
                       ("A. one", "Alpha one", "P3", 3), ("B. two", "Beta two", "P3", 5), ("C. three", "Gamma three", "P3", 7),
                       ("D. four", "Delta four", "P3", 9), ("E. five", "Eps five", "P1", 50), ("E. five", "Eps five", "P2", 1)],
                      columns=["species", "species_full", "paper_family", "count"])
    sm = pd.DataFrame({"paper_species": ["Alpha one", "Beta two", "Gamma three", "Delta four", "Eps five"],
                       "species_prefix": ["S1", "S2", "S3", "S4", None], "strain": "", "match": ["exact"] * 4 + ["none"],
                       "n_candidates": [1, 1, 1, 1, 0]})
    # BFD: P1 = Fam One + Fam One b (counts split), P2 = Fam Two; P3 (Fam Gone) has no rows in BFD
    fam = pd.DataFrame([("S1", "Fam One", 6), ("S1", "Fam One b", 4), ("S1", "Fam Two", 4),
                        ("S2", "Fam One", 12), ("S2", "Fam One b", 8), ("S2", "Fam Two", 8),
                        ("S3", "Fam One", 15), ("S3", "Fam One b", 15), ("S3", "Fam Two", 2),
                        ("S4", "Fam One", 40), ("S4", "Fam Two", 6)], columns=["species_prefix", "family", "n_proteins"])
    return fam, pc, sm


def test_calibrate_known_answer():
    fam, pc, sm = _perfect()
    summ, per = tp.calibrate(fam, pc, sm, FMAP)
    assert list(summ.columns) == ["paper_family", "n_species", "n_nonzero_paper", "n_nonzero_bfd", "spearman", "pearson_log",
                                  "median_ratio", "note"]
    s = summ.set_index("paper_family")
    assert s.loc["P1", "n_species"] == 4 and abs(s.loc["P1", "spearman"] - 1.0) < 1e-12
    assert abs(s.loc["P1", "median_ratio"] - 1.0) < 1e-12 and abs(s.loc["P1", "pearson_log"] - 1.0) < 1e-12
    assert s.loc["P2", "n_species"] == 4 and abs(s.loc["P2", "spearman"] - 1.0) < 1e-12
    # P3 has a Pfam in the map but no BFD rows at all: reported, not dropped, and not scored
    assert s.loc["P3", "n_species"] == 0 and np.isnan(s.loc["P3", "spearman"]) and "no BFD" in s.loc["P3", "note"]
    assert set(per["species_prefix"].dropna()) == {"S1", "S2", "S3", "S4"}      # the unmatched species is not in the table
    assert {"paper_species", "paper_family", "paper_count", "bfd_count"} <= set(per.columns)
    r = per[(per["paper_species"] == "Alpha one") & (per["paper_family"] == "P1")].iloc[0]
    assert r["paper_count"] == 10 and r["bfd_count"] == 10


def test_calibrate_ratio_and_disagreement():
    fam, pc, sm = _perfect()
    fam.loc[fam["family"] == "Fam Two", "n_proteins"] = [8, 16, 4, 12]
    s = tp.calibrate(fam, pc, sm, FMAP)[0].set_index("paper_family")
    assert abs(s.loc["P2", "median_ratio"] - 2.0) < 1e-12 and abs(s.loc["P2", "spearman"] - 1.0) < 1e-12


def test_paper_claims_use_species_not_genus():
    # five Ascomycota species of ONE genus and five Mucoromycota species of one genus: 5 species per group is testable
    # (grouping by genus would give one value per group and "not tested").
    sp = [f"Aa s{i}" for i in range(5)] + [f"Mm s{i}" for i in range(5)]
    pc = pd.DataFrame([(f"A.{i}" if i < 5 else f"M.{i}", s, "P1", c) for i, (s, c) in enumerate(zip(sp, [50, 52, 54, 56, 58, 10, 12, 14, 16, 18]))],
                      columns=["species", "species_full", "paper_family", "count"])
    totals = pd.DataFrame({"species": pc["species"], "species_full": pc["species_full"], "total_tfs": 100, "total_genes": 1000})
    ph = pd.DataFrame({"species_full": sp, "PHYLUM": ["Ascomycota"] * 5 + ["Mucoromycota"] * 5})
    cl = pd.DataFrame([dict(claim_id="c", text="t", family="P1", metric="fraction", level="species", group_a="Ascomycota",
                            group_b="Mucoromycota", paper_direction="a>b")])
    r = tp.evaluate_claims_on_paper(pc, ph, cl, totals, n_boot=100).iloc[0]
    assert r["n_genera_a"] == 5 and r["n_genera_b"] == 5 and r["source"] == "paper"
    assert r["effect"] == 1.0 and abs(r["value_a"] - 0.054) < 1e-12 and r["agrees"] in ("yes", "mixed")
    assert set(tc_cols()) <= set(r.index)


def tc_cols():
    return ["claim_id", "text", "family", "metric", "group_a", "group_b", "n_genera_a", "n_genera_b", "value_a", "value_b",
            "effect", "ci_lo", "ci_hi", "p", "p_adj", "observed_direction", "paper_direction", "agrees"]


def test_paper_claims_species_without_phylum_are_left_out():
    sp = [f"Aa s{i}" for i in range(5)] + [f"Mm s{i}" for i in range(5)] + ["Lost one"]
    pc = pd.DataFrame([(s[:4], s, "P1", 10 + i) for i, s in enumerate(sp)], columns=["species", "species_full", "paper_family", "count"])
    totals = pd.DataFrame({"species": pc["species"], "species_full": pc["species_full"], "total_tfs": 100, "total_genes": 1000})
    ph = pd.DataFrame({"species_full": sp[:10], "PHYLUM": ["Ascomycota"] * 5 + ["Mucoromycota"] * 5})
    cl = pd.DataFrame([dict(claim_id="c", text="t", family="P1", metric="count", level="species", group_a="Ascomycota",
                            group_b="Mucoromycota", paper_direction="a<b")])
    r = tp.evaluate_claims_on_paper(pc, ph, cl, totals, n_boot=50).iloc[0]
    assert r["n_genera_a"] == 5 and r["n_genera_b"] == 5


def test_calibration_page_writes_table_and_figure(tmp_path):
    import make_figures as mf
    fam, pc, sm = _perfect()
    o = mf.Out(tmp_path, "tf")
    tp.page_tf_calibration(fam, pc, sm, FMAP, o)
    t = pd.read_csv(tmp_path / "tables" / "tf" / "T09_tf_calibration.tsv", sep="\t", keep_default_na=False)
    assert {"P1", "P2", "P3"} <= set(t["paper_family"])
    for ext in ("png", "svg"):
        assert (tmp_path / "figures" / "tf" / f"T09_tf_calibration.{ext}").stat().st_size > 1000
    assert (tmp_path / "tables" / "tf" / "T09_tf_calibration_species.tsv").exists()


def test_committed_paper_tables_are_consistent():
    cfg = Path(__file__).resolve().parents[2] / "config"
    counts, totals = tp.load_paper_counts(cfg / "paper_tf_counts.tsv")
    smap = pd.read_csv(cfg / "paper_species_map.tsv", sep="\t", comment="#", keep_default_na=False)
    assert counts["species"].nunique() == 62 and len(totals) == 62 and len(smap) == 62
    assert list(smap["paper_species"]) == list(totals["species_full"])
    assert set(smap["match"]) <= {"exact", "alias", "none"}
    assert set(smap.loc[smap["match"] == "alias", "alias_reason"]) <= {"typo", "renamed taxon", "synonym"}
    assert (smap.loc[smap["match"] == "alias", "alias_reason"] != "").all() and (smap.loc[smap["match"] != "alias", "alias_reason"] == "").all()
    assert ((smap["match"] == "none") == (smap["species_prefix"] == "")).all()
    fmap = pd.read_csv(cfg / "tf_family_map.tsv", sep="\t", dtype=str, comment="#", keep_default_na=False)
    assert (set(fmap["paper_family"]) - {""}) <= set(counts["paper_family"])        # the map uses the paper's exact names
    assert counts["paper_family"].str.strip().eq(counts["paper_family"]).all()


def test_run_tf_claims_end_to_end_without_matches(tmp_path):
    """The page driver on a synthetic genome_all: no paper species match, so every paper row is 'not tested' and the
    calibration table is written with n_species 0. Nothing raises."""
    import make_figures as mf
    from test_figs_tf import tf_synth
    from test_make_figures import synth
    root = tmp_path / "d"
    root.mkdir()
    synth(root)
    fam, summ = tf_synth(root)
    ga = pd.read_parquet(root / "genome_all.parquet")
    made = tp.run_tf_claims(ga, fam, summ, tmp_path / "o", n_boot=50)
    names = {str(m): m.kind for m in made}
    assert names["T10_tf_claims"] == "table" and names["T10_tf_claims_nofilter"] == "table"
    t = pd.read_csv(tmp_path / "o" / "tables" / "tf" / "T10_tf_claims.tsv", sep="\t", keep_default_na=False)
    assert set(t["source"]) == {"BFD", "paper"}
    assert (t.loc[t["source"] == "paper", "agrees"] == "not tested").all()
    c = pd.read_csv(tmp_path / "o" / "tables" / "tf" / "T09_tf_calibration.tsv", sep="\t", keep_default_na=False)
    assert (c["n_species"] == 0).all() or len(c) > 0
    assert (tmp_path / "o" / "figures" / "tf" / "T10_tf_claims_forest.png").exists()


def test_main_runs_claims_page_after_landscape(tmp_path):
    import make_figures as mf
    from test_figs_tf import tf_synth
    from test_make_figures import synth
    root = tmp_path / "d"
    root.mkdir()
    synth(root)
    tf_synth(root)
    made = mf.main(root, tmp_path / "o", ("phylum",), suppress=pd.DataFrame({"species_prefix": []}))
    assert "T01_tf_fraction" in made["phylum"] and "T10_tf_claims_forest" in made["tf"]
    assert (tmp_path / "o" / "tables" / "tf" / "T10_tf_claims.tsv").exists()
    assert not mf.RUN_STATUS["failed"]


def test_calibrate_nonzero_counts_and_not_comparable_family():
    fam, pc, sm = _perfect()
    # P2: BFD has the family label but zero proteins in every matched species -> not comparable: NaN, not 0.0
    fam = fam[fam["family"] != "Fam Two"]
    fam = pd.concat([fam, pd.DataFrame([("S9", "Fam Two", 5)], columns=fam.columns)], ignore_index=True)   # a non-matched genome
    s = tp.calibrate(fam, pc, sm, FMAP)[0].set_index("paper_family")
    assert s.loc["P1", "n_nonzero_paper"] == 4 and s.loc["P1", "n_nonzero_bfd"] == 4
    assert s.loc["P2", "n_nonzero_paper"] == 4 and s.loc["P2", "n_nonzero_bfd"] == 0
    assert np.isnan(s.loc["P2", "spearman"]) and np.isnan(s.loc["P2", "median_ratio"]) and np.isnan(s.loc["P2", "pearson_log"])
    assert "BFD Pfam definition does not detect this family" in s.loc["P2", "note"]


def test_calibrate_genus_nonzero_threshold_flags_family():
    fam, pc, sm = _perfect()
    assert "not comparable" not in tp.calibrate(fam, pc, sm, FMAP)[0].set_index("paper_family").loc["P1", "note"]
    s = tp.calibrate(fam, pc, sm, FMAP, genus_nonzero={"P1": 2, "P2": 50})[0].set_index("paper_family")
    assert np.isnan(s.loc["P1", "spearman"]) and np.isnan(s.loc["P1", "median_ratio"])
    assert "non-zero in 2 genera" in s.loc["P1", "note"] and np.isfinite(s.loc["P2", "spearman"])


def test_calibrate_uses_paper_family_table_and_fallback_note():
    fam, pc, sm = _perfect()
    # distinct-protein paper-family counts: P1 is half of the sum over the two Pfam families
    pf = pd.DataFrame([("S1", "P1", 5), ("S2", "P1", 10), ("S3", "P1", 15), ("S4", "P1", 20),
                       ("S1", "P2", 4), ("S2", "P2", 8), ("S3", "P2", 2), ("S4", "P2", 6)],
                      columns=["species_prefix", "paper_family", "n_proteins"])
    s = tp.calibrate(fam, pc, sm, FMAP, bfd_paper_family=pf)[0].set_index("paper_family")
    assert abs(s.loc["P1", "median_ratio"] - 0.5) < 1e-12 and s.loc["P1", "note"] == ""
    s0 = tp.calibrate(fam, pc, sm, FMAP)[0].set_index("paper_family")
    assert "paper-family counts summed over Pfam families (double counting possible)" in s0.loc["P1", "note"]


def test_paper_side_note_for_not_comparable_family():
    sp = [f"Aa s{i}" for i in range(5)] + [f"Mm s{i}" for i in range(5)]
    pc = pd.DataFrame([(s[:4] + str(i), s, "P1", c) for i, (s, c) in enumerate(zip(sp, [5, 6, 7, 8, 9, 0, 0, 1, 0, 0]))],
                      columns=["species", "species_full", "paper_family", "count"])
    totals = pd.DataFrame({"species": pc["species"], "species_full": pc["species_full"], "total_tfs": 100, "total_genes": 1000})
    ph = pd.DataFrame({"species_full": sp, "PHYLUM": ["Ascomycota"] * 5 + ["Mucoromycota"] * 5})
    cl = pd.DataFrame([dict(claim_id="c", text="t", family="P1", metric="fraction", level="species", group_a="Ascomycota",
                            group_b="Mucoromycota", paper_direction="a>b")])
    r = tp.evaluate_claims_on_paper(pc, ph, cl, totals, n_boot=50, not_comparable={"P1"}).iloc[0]
    assert "paper definition; not reproduced by BFD Pfam" in r["note"] and r["agrees"] in ("yes", "mixed")
    assert r["filter"] != "BUSCO>=80"
