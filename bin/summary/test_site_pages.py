"""Tests for the Quarto site: navbar targets, the gallery helper and the required statements in the pages."""
import html as htmllib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.parse import quote

import yaml

SITE = Path(__file__).resolve().parents[2] / "docs" / "quarto_summary"
sys.path.insert(0, str(SITE))
import _helpers  # noqa: E402

INDEX = ("name\tlevel\ttopic\tkind\tscope\tn_genomes\tn_groups\tpooled_share\teffect_size\tscore\tsignal\tsignal_kind\n"
         "S01_gene_len\tphylum\tstructure\tfigure\trepresentatives\t240\t5\t0.0\t0.2\t0.2\t0.2\tepsilon2\n"
         "A01_busco_vs_genome_size\tall\tscatter\tfigure\tall_genomes\t360\t\t\t\t\t0.31\tspearman_r\n"
         "O01_codon_pcoa_phylum\tordination\tordination\tfigure\tgenus_reps\t89\t5\t\t\t\t0.12\tpermanova_R2\n")


def _read(name):
    return (SITE / name).read_text()


def test_navbar_hrefs_exist():
    cfg = yaml.safe_load(_read("_quarto.yml"))
    hrefs = [i["href"] for i in cfg["website"]["navbar"]["left"]]
    assert len(hrefs) >= 12
    for h in hrefs:
        assert (SITE / h).is_file(), h
    for page in ("07_gene_structure", "08_composition", "09_comparisons", "10_ordination", "11_gallery"):
        assert f"{page}.qmd" in hrefs


def test_gallery_has_every_row_and_filters(tmp_path):
    p = tmp_path / "figures_index.tsv"
    p.write_text(INDEX)
    h = _helpers.gallery_html(index_path=p)
    for name in ("S01_gene_len", "A01_busco_vs_genome_size", "O01_codon_pcoa_phylum"):
        assert f"figures/{'phylum' if name[0] == 'S' else 'all' if name[0] == 'A' else 'ordination'}/{name}.png" in h
        assert f"{name}.svg" in h
    for col in ("topic", "level", "scope", "signal_kind"):
        assert f'data-col="{col}"' in h
    assert "<script>" in h and "http" not in h
    assert h.index("S01_gene_len") < h.index("A01_busco") < h.index("O01_codon")  # index order kept
    assert "A01_busco" not in _helpers.gallery_html(topic="structure", index_path=p)
    top2 = _helpers.gallery_html(top=2, index_path=p)
    assert "O01_codon" not in top2 and top2.count("<tr data-") == 2


def test_gallery_missing_index(tmp_path):
    h = _helpers.gallery_html(index_path=tmp_path / "none.tsv")
    assert "has not been built" in h and "<table" not in h


def test_required_statements():
    comp = _read("09_comparisons.qmd")
    assert "do not control for phylogeny" in comp and "no multiple-testing correction" in comp
    ordn = _read("10_ordination.qmd")
    assert "one genome per genus" in ordn and "best BUSCO, then N50, then lowest id" in ordn
    assert "complete codon, amino-acid and gene-structure data" in ordn
    assert "999 permutations" in ordn and "by phylum and by order" in ordn
    for page in ("07_gene_structure.qmd", "08_composition.qmd"):
        t = _read(page)
        assert "protein-coding transcripts" in t, page
    assert "exon counts may include UTR exons" in _read("07_gene_structure.qmd")
    gal = _read("11_gallery.qmd")
    assert "effect_size" in gal and "pooled_share" in gal and "not comparable" in gal


LV_INDEX = (INDEX +
            "G1\tphylum\tstructure\tfigure\trepresentatives\t10\t3\t0.0\t0.1\t0.1\t0.1\tepsilon2\n"
            "G2\tsubphylum\tstructure\tfigure\trepresentatives\t10\t3\t0.0\t0.1\t0.1\t0.1\tepsilon2\n"
            "G3\torder\tcomposition\tfigure\trepresentatives\t10\t3\t0.0\t0.1\t0.1\t0.1\tepsilon2\n"
            "G4\tordination\tordination\tfigure\tgenus_reps\t10\t3\t\t\t\t0.1\tpermanova_R2\n")


LV_ONLY = INDEX.splitlines()[0] + "\n" + LV_INDEX.split("\n", 4)[4]


def test_gallery_escapes_paths(tmp_path):
    p = tmp_path / "i.tsv"
    bad = 'x"><script>alert(1)</script>'
    p.write_text(INDEX.splitlines()[0] + "\n" + f"{bad}\tphylum\tstructure\tfigure\trepresentatives\t1\t1\t0\t0.1\t0.1\t0.1\tepsilon2\n")
    h = _helpers.gallery_html(index_path=p)
    body = h.split("<tbody>")[1].split("</tbody>")[0]
    assert "<script>" not in body and 'x">' not in body
    assert f'href="figures/phylum/{quote(bad, safe="")}.png"' in htmllib.unescape(body) or quote(bad, safe="") in body
    assert quote(bad, safe="") in body


def test_gallery_structure_and_a11y(tmp_path):
    p = tmp_path / "i.tsv"
    p.write_text(LV_ONLY)
    h = _helpers.gallery_html(index_path=p)
    assert re.search(r'<div class="gal-scroll" style="overflow-x:auto">\s*<table', h)
    assert "aria-label=" in h.split("<table")[1].split(">")[0] or "<caption" in h
    assert 'aria-live="polite"' in h
    assert 'width="90" height="53"' in h and 'decoding="async"' in h and 'tabindex="-1"' in h
    assert 'aria-label="PNG of G1"' in h and 'aria-label="SVG of G1"' in h
    assert 'alt="G1 (structure, phylum)"' in h
    # closed-set columns are selects with exactly the unique values; name is free text
    for col, vals in (("level", ["phylum", "subphylum", "order", "ordination"]), ("topic", ["structure", "composition", "ordination"]),
                      ("scope", ["representatives", "genus_reps"]), ("signal_kind", ["epsilon2", "permanova_R2"])):
        m = re.search(rf'<select class="gal-filter" data-col="{col}" data-mode="exact"[^>]*>(.*?)</select>', h, re.S)
        assert m, col
        opts = re.findall(r'<option value="([^"]*)"', m.group(1))
        assert opts[0] == "" and sorted(opts[1:]) == sorted(vals), (col, opts)
    assert 'data-col="name"' in h and 'data-mode="contains"' in h


def _run_filter(h, settings):
    """Execute the gallery script in node against a minimal DOM stub; returns the visible row names."""
    script = h.split("<script>")[1].split("</script>")[0]
    rows = [dict(re.findall(r'data-(\w+)="([^"]*)"', m)) for m in re.findall(r"<tr (data-[^>]*)>", h)]
    harness = """
const rows = %s.map(d => ({dataset: d, style: {}}));
const ctrls = %s.map(c => ({dataset: {col: c[0], mode: c[1]}, value: c[2], handlers: [], addEventListener(e, f) { this.handlers.push(f); }}));
const count = {textContent: ''};
const root = {querySelectorAll(s) { return s.indexOf('tbody') >= 0 ? rows : ctrls; }, querySelector() { return count; }};
const document = {currentScript: {parentNode: root}};
%s
console.log(JSON.stringify(rows.filter(r => r.style.display !== 'none').map(r => r.dataset.name)));
""" % (json.dumps(rows), json.dumps(settings), script)
    out = subprocess.run(["node", "-e", harness], capture_output=True, text=True, check=True).stdout
    return json.loads(out.strip().splitlines()[-1])


import pytest  # noqa: E402


@pytest.mark.skipif(shutil.which("node") is None, reason="node not available")
def test_gallery_js_exact_filters(tmp_path):
    p = tmp_path / "i.tsv"
    p.write_text(LV_INDEX)
    h = _helpers.gallery_html(index_path=p)
    mk = lambda col, mode, v: [col, mode, v]  # noqa: E731
    assert _run_filter(h, [mk("level", "exact", "phylum")]) == ["S01_gene_len", "G1"]
    assert _run_filter(h, [mk("level", "exact", "order")]) == ["G3"]
    assert _run_filter(h, [mk("topic", "exact", "ordination"), mk("level", "exact", "")]) == ["O01_codon_pcoa_phylum", "G4"]
    assert _run_filter(h, [mk("name", "contains", "busco")]) == ["A01_busco_vs_genome_size"]
    assert _run_filter(h, [mk("name", "contains", "GENE")]) == ["S01_gene_len"]
    assert len(_run_filter(h, [mk("level", "exact", "")])) == 7


def test_facts_tolerates_nan(tmp_path, monkeypatch):
    d = tmp_path / "figures" / "phylum"
    d.mkdir(parents=True)
    (d / "z.json").write_text('{"n_genomes": NaN, "n_groups": "x", "pooled_share": NaN}')
    monkeypatch.setattr(_helpers, "HERE", tmp_path)
    assert _helpers._facts("phylum", "z") == ""
    (d / "y.json").write_text('{"n_genomes": 5, "n_groups": NaN, "pooled_share": 0.5}')
    assert "n = 5 genomes" in _helpers._facts("phylum", "y")


def test_index_and_notes_caveats():
    for page in ("index.qmd", "06_data_notes.qmd"):
        t = _read(page)
        assert "no statistical tests" not in t and "do not test hypotheses" not in t, page
    idx = _read("index.qmd")
    assert "do not control for phylogeny" in idx and "multiple testing" in idx and "genus representatives" in idx


def test_explore_page_wired():
    cfg = yaml.safe_load(_read("_quarto.yml"))
    nav = {i["href"]: i["text"] for i in cfg["website"]["navbar"]["left"]}
    assert nav.get("12_explore.qmd") == "Explore"
    assert "data/genomes.json" in cfg["project"]["resources"] and "CNAME" in cfg["project"]["resources"]
    page = _read("12_explore.qmd")
    assert "{ojs}" in page and 'FileAttachment("data/genomes.json")' in page
    assert "<script" not in page
    assert "phase B" in page and "Each point is one genome" in page and "controlled for phylogeny" in page


def test_committed_genomes_json():
    p = SITE / "data" / "genomes.json"
    j = json.loads(p.read_text())
    assert len(j["rows"]) == j["meta"]["n_genomes"]
    assert p.stat().st_size < 5_000_000
    assert "NaN" not in p.read_text()
    assert j["meta"]["n_suppressed"] == 3
    assert "species" in j["fields"] and "strain" in j["fields"] and "id" in j["fields"]
    i = j["fields"].index("species")
    assert all(r[i] for r in j["rows"])


def test_explore_shows_species_not_id():
    page = _read("12_explore.qmd")
    assert 'columns: ["species", "strain"' in page
    assert "${d.id}" not in page
    assert "d.species" in page and "d.strain" in page


def test_explore_log10_and_defaults():
    page = _read("12_explore.qmd")
    assert "Math.log10" in page and "(log10)" in page
    assert "untransformed" not in page
    assert "logDefault" in page and "viewof logx = Inputs.toggle({label: \"Log x\"" in page and ": logDefault.has(xvar)" in page
    assert "genome_mb" in page.split("logDefault =")[1].split("```")[0]
    assert "pooled small groups" in page and "aria-label" in page
    notes = _read("06_data_notes.qmd")
    assert "python3.12 bin/summary/export_clade_json.py" in notes and "docs/quarto_summary/data/genomes.json" in notes


# ---- final-review text fixes ----
def test_composition_captions_match_code():
    t = _read("08_composition.qmd")
    assert "Amino-acid frequency against genome GC." not in t
    assert "GC3" in t.split('fig("C05_aa_vs_gc"')[1].split("\n")[0]
    assert "standard genetic code" in t and "CTG" in t


def test_ordination_page_shows_negative_eigen_columns():
    t = _read("10_ordination.qmd")
    for c in ("n_negative_eigen", "neg_eigen_over_abs_total", "neg_eigen_over_pos_total", "n_permutations_run"):
        assert t.count(c) >= 3, c


def test_index_scope_wording():
    t = _read("index.qmd")
    assert "Every figure uses **one representative genome per species**" not in t
    assert "Each figure is shown at four taxonomic levels" not in t
    assert "all genomes" in t and "one genome per genus" in t and "largest" in t


def test_explore_scope_wording():
    t = _read("12_explore.qmd")
    assert "The static pages use representatives only" not in t
    assert "level pages use one representative per species" in t and "Comparisons page use all genomes" in t


def test_data_notes_counts_and_build():
    t = _read("06_data_notes.qmd")
    assert "meta" in t and "n_genomes" in t and "genome_all" in t and "gene-structure" in t
    assert "23,680" not in t
    assert "sbatch workflows/summary_tables.sh" in t and "sbatch workflows/summary_figures.sh" in t
    assert "workflows/release_site_assets.sh --upload" in t and "export_clade_json.py" in t
    assert "figures_release.json" in t and "not in git" in t
    assert "build_summary.sh` is superseded" in t


def test_gene_structure_definitions_match_sql():
    t = _read("07_gene_structure.qmd")
    assert "come from protein-coding transcripts only" not in t
    assert "per transcript" in t and "isoform" in t and "intergenic" in t and "not filtered" in t


def test_gallery_note_wording():
    t = _read("11_gallery.qmd")
    assert "n_genomes is the number of genomes in the figure" not in t
    assert "group table" in t and "O02" in t and "by-phylum" in t


def test_ordination_page_adjusted_r2_and_clipping():
    t = _read("10_ordination.qmd")
    assert t.count('"R2", "R2_adj"') >= 3
    assert "R2_adj" in t and "inflated by the number of groups" in t
    assert "clipped at" in t and "n_clipped_values" in t and "n_genera_clipped" in t
