"""Tests for the Quarto site: navbar targets, the gallery helper and the required statements in the pages."""
import sys
from pathlib import Path

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
    assert "O01_codon" not in _helpers.gallery_html(top=2, index_path=p)


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
