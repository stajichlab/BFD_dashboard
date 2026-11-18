"""Tests for sql/bfd_traits.sql and the v_genome_labels join, on small synthetic tables."""
import csv
import re
from pathlib import Path

import duckdb

ROOT = Path(__file__).resolve().parents[1]
TRAITS_SQL = ROOT / "sql" / "bfd_traits.sql"
EXTRA_SQL = ROOT / "sql" / "bfd_extra_views.sql"

FT_COLS = ["jrk_template", "GENUS", "primary_lifestyle", "Secondary_lifestyle",
           "Plant_pathogenic_capacity_template", "Decay_substrate_template",
           "Animal_biotrophic_capacity_template", "Growth_form_template"]


def _write(path, header, rows, delim):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=delim)
        w.writerow(header)
        w.writerows(rows)


def ft(jrk, genus, life=None, growth=None):
    return [jrk, genus, life, None, None, None, None, growth]


def load(tmp_path, species, ft_rows=(), ovr_rows=(), cls_rows=()):
    """species: list of (LOCUSTAG, GENUS, CLASS, SUBPHYLUM). Returns a connection with the traits SQL applied."""
    _write(tmp_path / "fungaltraits_1.2.utf8.csv", FT_COLS, ft_rows, ",")
    _write(tmp_path / "growth_form_overrides.tsv", ["GENUS", "growth_form", "source", "n_genomes", "CLASS"],
           ovr_rows, "\t")
    _write(tmp_path / "class_subphylum.tsv", ["CLASS", "SUBPHYLUM", "source", "note"], cls_rows, "\t")
    sql = TRAITS_SQL.read_text()
    first = sql.splitlines()[0]
    assert first.startswith("-- REF_DIR="), "line 1 of bfd_traits.sql must be '-- REF_DIR=<path>'"
    default = first[len("-- REF_DIR="):]
    con = duckdb.connect()
    con.execute("CREATE TABLE species (LOCUSTAG VARCHAR, GENUS VARCHAR, CLASS VARCHAR, SUBPHYLUM VARCHAR)")
    con.executemany("INSERT INTO species VALUES (?,?,?,?)", species)
    con.execute(sql.replace(default, str(tmp_path)))
    return con


def one(con, locustag):
    cur = con.execute("SELECT * FROM v_genome_traits WHERE LOCUSTAG = ?", [locustag])
    cols = [d[0] for d in cur.description]
    rows = cur.fetchall()
    assert len(rows) == 1, f"{locustag}: {len(rows)} rows"
    return dict(zip(cols, rows[0]))


def test_duplicate_genus_keeps_first_by_jrk_template(tmp_path):
    con = load(tmp_path, [("A1", "Mucor", None, None), ("A2", "Mucor", None, None)],
               ft_rows=[ft("jrk2", "Mucor", "wrong"), ft("jrk1", "Mucor", "right")])
    assert con.execute("SELECT count(*) FROM v_genome_traits").fetchone()[0] == 2
    assert one(con, "A1")["primary_lifestyle"] == "right"


def test_override_wins_and_sets_source(tmp_path):
    con = load(tmp_path, [("A1", "Candida", None, None), ("A2", "Mucor", None, None)],
               ft_rows=[ft("j1", "Candida", "x", "filamentous_mycelium"), ft("j2", "Mucor", "y", "filamentous_mycelium")],
               ovr_rows=[("Candida", "yeast", "class_rule", "9", "Saccharomycetes"),
                         ("Mucor", "", "class_rule", "1", "Mucoromycetes")])
    a1, a2 = one(con, "A1"), one(con, "A2")
    assert (a1["growth_form"], a1["growth_form_source"]) == ("yeast", "override")
    assert (a2["growth_form"], a2["growth_form_source"]) == ("filamentous_mycelium", "fungaltraits")


def test_unmatched_genus_and_null_genus(tmp_path):
    con = load(tmp_path, [("A1", "Nogenus", None, None), ("A2", None, None, None)], ft_rows=[ft("j1", "Mucor", "y")])
    for lt in ("A1", "A2"):
        r = one(con, lt)
        assert r["in_fungaltraits"] is False
        assert r["primary_lifestyle"] is None and r["growth_form"] is None and r["growth_form_source"] is None


def test_subphylum_from_class_table(tmp_path):
    cls = [("Eurotiomycetes", "Pezizomycotina", "taxonomy_py", "")]
    con = load(tmp_path, [("A1", "G", "Eurotiomycetes", None), ("A2", "G", "Unlisted", None), ("A3", "G", None, None)],
               cls_rows=cls)
    r = one(con, "A1")
    assert (r["SUBPHYLUM"], r["subphylum_source"]) == ("Pezizomycotina", "class_table")
    for lt in ("A2", "A3"):
        r = one(con, lt)
        assert r["SUBPHYLUM"] is None and r["subphylum_source"] is None


def test_species_subphylum_wins(tmp_path):
    cls = [("Eurotiomycetes", "Pezizomycotina", "taxonomy_py", "")]
    con = load(tmp_path, [("A1", "G", "Eurotiomycetes", "Other-sub")], cls_rows=cls)
    r = one(con, "A1")
    assert (r["SUBPHYLUM"], r["subphylum_source"]) == ("Other-sub", "species")


def test_yeast_origin_matches_atlas_rules(tmp_path):
    spec = [("A1", "Gy", "Saccharomycetes", None), ("A2", "Gd", "Eurotiomycetes", None),
            ("A3", "Gf", "Sordariomycetes", None), ("A4", "Gn", "Sordariomycetes", None)]
    con = load(tmp_path, spec, ft_rows=[ft("j1", "Gy", None, "yeast"), ft("j2", "Gd", None, "dimorphic_yeast"),
                                        ft("j3", "Gf", None, "filamentous_mycelium"), ft("j4", "Gn", None, "other_form")])
    got = {lt: one(con, lt)["yeast_origin"] for lt in ("A1", "A2", "A3", "A4")}
    assert got == {"A1": "Saccharomycotina", "A2": "Pezizomycotina-dimorphic", "A3": "filamentous", "A4": "other"}
