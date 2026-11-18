"""Tests for the ortholog tracking modules (seeds, ortho_search, ortho_tree, ortho_report) and the make_figures hook.
No network. The MMseqs2/MAFFT/FastTree integration test skips when the tools are not on PATH.

Run:  module load MMseqs2/17-b804f ; module load mafft/7.505 ; module load fasttree/2.1.11   (each on its own line)
      python -m pytest -q -W error::UserWarning bin/tf/test_tf_ortho.py
"""
import json
import os
import random
import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "bin" / "summary"))
import ortho_report as orp  # noqa: E402
import ortho_search as ors  # noqa: E402
import ortho_tree as ort  # noqa: E402
import seeds as sd  # noqa: E402

SEED_CFG = HERE.parents[1] / "config" / "tf_seeds.tsv"
PFAM_NAMES = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/summary_tables/pfam_names.parquet"
PAS = "PF00989|PF08447|PF13426"
PAS_ALL = ["PF00989", "PF08447", "PF13426", "PF08446", "PF14598", "PF08448", "PF08348", "PF07310", "PF12860", "PF13188"]


# ---------------------------------------------------------------- seeds
def _uni_seeds():
    return pd.DataFrame({"target": ["Zap1-like", "bZip5"], "seed_name": ["ZAP1", "CAD1"], "source": ["uniprot"] * 2,
                         "key": ["ZAP1", "CAD1"], "taxon_id": ["559292", "559292"]})


def test_seed_fetch_uses_cache_and_reports_missing(tmp_path):
    calls = []

    def fetch(url):
        calls.append(url)
        if "CAD1" in url:
            raise OSError("network down")
        return ">sp|P47988|ZAP1_YEAST Zinc-responsive\nMKLV\nAAGG\n"
    df = sd.fetch_seeds(_uni_seeds(), tmp_path, fetch=fetch)
    assert list(df.columns) == ["target", "seed_name", "accession", "sequence", "status"]
    z = df[df.seed_name == "ZAP1"].iloc[0]
    assert z.status == "ok" and z.accession == "P47988" and z.sequence == "MKLVAAGG"
    c = df[df.seed_name == "CAD1"].iloc[0]
    assert c.status.startswith("missing") and "network down" in c.status and c.sequence == ""
    n = len(calls)
    df2 = sd.fetch_seeds(_uni_seeds(), tmp_path, fetch=fetch)
    assert len(calls) == n + 1          # ZAP1 came from the cache; only the failing key is tried again
    assert df2[df2.seed_name == "ZAP1"].iloc[0].sequence == "MKLVAAGG"


def test_seed_fetch_empty_answer_is_missing(tmp_path):
    df = sd.fetch_seeds(_uni_seeds().head(1), tmp_path, fetch=lambda url: "")
    assert df.iloc[0].status.startswith("missing") and "not found" in df.iloc[0].status


def test_seed_config_accessions_exist_in_pfam_names():
    cfg = sd.read_seed_config(SEED_CFG)
    have = set(pd.read_parquet(PFAM_NAMES)["pfam_acc"])
    for _, r in cfg.iterrows():
        for a in sd.rule_accessions(r["target_pfam"], r["required_pfam"]):
            assert a in have, a
    for a in ors.PAS_ACCS:
        assert a in have, a
    assert set(cfg["target"]) == {"WC-2", "WC-1", "Zap1-like", "bZip5"}
    assert set(cfg.loc[cfg.source == "paper", "key"]) == {"1799646", "1921537", "1909793"}
    wc1 = cfg[cfg.target == "WC-1"]
    assert "wc-1" in set(wc1.key) and (wc1.required_pfam == "PAS_regions>=2").all()
    assert (cfg[cfg.target == "WC-2"].required_pfam == "PAS_regions==1").all()
    bfd = cfg[cfg.source == "bfd"]
    assert len(bfd) and not bfd.key.str.contains(r"GATA\+PAS$").any()      # the new region rule, not 'GATA+PAS'


def test_bfd_seed_selector_uses_pas_regions():
    pf = pd.DataFrame({
        "protein_id": ["G_1-T1"] * 2 + ["G_2-T1"] * 4 + ["G_3-T1"] * 2 + ["G_4-T1"],
        "pfam_acc": ["PF00320.25", "PF00989.30"] + ["PF00320.25", "PF00989.30", "PF13426.12", "PF13426.12"]
                    + ["PF00320.25", "PF08447.14"] + ["PF00989.30"],
        "ali_from": [300, 10, 700, 10, 20, 400, 300, 10, 10],
        "ali_to": [350, 120, 750, 120, 130, 520, 350, 120, 120]})
    # G_1: GATA + 1 PAS. G_2: GATA + PAS/PAS_9 overlapping (one region) + a second PAS_9 (two regions). G_3: 1 region.
    assert sd.bfd_selector_proteins("G", "pfam:GATA+PAS_regions==1", pf) == ["G_1-T1", "G_3-T1"]
    assert sd.bfd_selector_proteins("G", "pfam:GATA+PAS_regions>=2", pf) == ["G_2-T1"]


def test_rule_parser_one_place():
    assert sd.parse_rule("PF00096>=4") == ("count", {"PF00096"}, ">=", 4)
    assert sd.parse_rule("PF00989|PF08447") == ("any", {"PF00989", "PF08447"})
    assert sd.parse_rule("PAS_regions==1") == ("pas", "==", 1)
    assert sd.parse_rule("PAS_regions>=2") == ("pas", ">=", 2)
    assert sd.rule_accessions("PF00320", "PAS_regions==1") == ["PF00320"]
    with pytest.raises(ValueError):
        sd.parse_rule("PAS_regions<3")


def test_resolve_paper_seeds_names_and_architecture_warning(tmp_path):
    names = tmp_path / "n.tsv"
    names.write_text("# c\njgi_id\tname\n1799646\tZnf1\n1921537\tGATA_Znf10\n1909793\tbZip5\n")
    paper = {"Rhimi59_2_1799646": "MKKK", "Rhimi59_2_1921537": "MLLL", "Rhimi59_2_1909793": "MPPP", "Rhimi59_2_5": "MQQQ"}
    bfd = tmp_path / "bfd.faa"
    bfd.write_text(">B_1-T1\nMKKK\n>B_2-T1\nMLLL\n>B_3-T1\nMPPP\n")

    def search(q, t, out, tmp, pident, cov):
        return pd.DataFrame({"query": ["Rhimi59_2_1799646", "Rhimi59_2_1921537", "Rhimi59_2_1909793"],
                             "target": ["B_1-T1", "B_2-T1", "B_3-T1"], "pident": [99.0, 98.0, 100.0],
                             "qcov": [1.0, 1.0, 1.0], "tcov": [1.0, 1.0, 1.0], "evalue": [0.0] * 3, "bits": [900.0] * 3})
    arch = {"B_1-T1": "PF00096>PF00096>PF00096>PF00096", "B_2-T1": "PF00320", "B_3-T1": "PF00170"}   # WC-2C lacks PAS
    df = sd.resolve_paper_seeds(paper, bfd, names_tsv=names, arch=arch, tmp=tmp_path, search=search)
    r = df.set_index("paper_id")
    assert r.loc["Rhimi59_2_1799646", "name"] == "Znf1" and r.loc["Rhimi59_2_1799646", "warning"] == ""
    assert r.loc["Rhimi59_2_1921537", "name"] == "GATA_Znf10"
    assert "PAS" in r.loc["Rhimi59_2_1921537", "warning"]
    assert r.loc["Rhimi59_2_1909793", "warning"] == ""
    assert r.loc["Rhimi59_2_5", "bfd_protein"] == "" and r.loc["Rhimi59_2_5", "name"] == ""


# ---------------------------------------------------------------- search and architecture
def test_architecture_ok_rules():
    assert ors.architecture_ok("PF00320>PF00989", [PAS])
    assert ors.architecture_ok("PF00320>PF13426>PF00989", [PAS])
    assert not ors.architecture_ok("PF00320", [PAS])
    assert ors.architecture_ok("PF00320", [])
    assert ors.architecture_ok(">".join(["PF00096"] * 4), ["PF00096>=4"])
    assert not ors.architecture_ok(">".join(["PF00096"] * 3), ["PF00096>=4"])
    assert not ors.architecture_ok("", ["PF00096>=4"])
    assert ors.architecture_ok("PF00320>PF00989", ["PAS_regions==1"], n_pas=1)
    assert not ors.architecture_ok("PF00320>PF00989", ["PAS_regions==1"], n_pas=3)
    assert ors.architecture_ok("PF00320", ["PAS_regions>=2"], n_pas=2)
    assert not ors.architecture_ok("PF00320", ["PAS_regions>=2"])      # no region count given: fails


def _dom(rows):
    return pd.DataFrame(rows, columns=["pfam_acc", "ali_from", "ali_to"])


def test_pas_regions_merges_overlapping_family_hits():
    one = _dom([("PF00989.30", 10, 100), ("PF08447.14", 20, 110), ("PF13426.12", 95, 200), ("PF00320.25", 400, 450)])
    assert ors.pas_regions(one, ors.PAS_ACCS) == [(10, 200)]
    three = _dom([("PF00989.30", 10, 100), ("PF08447.14", 20, 110), ("PF13426.12", 300, 400), ("PF14598.5", 385, 410),
                  ("PF13426.12", 600, 700), ("PF00320.25", 800, 850)])
    assert ors.pas_regions(three, ors.PAS_ACCS) == [(10, 110), (300, 410), (600, 700)]
    touching = _dom([("PF00989.30", 10, 100), ("PF08447.14", 95, 180)])      # overlap 6 aa < 10: two regions
    assert len(ors.pas_regions(touching, ors.PAS_ACCS)) == 2
    assert ors.pas_regions(_dom([("PF00320.25", 1, 50)]), ors.PAS_ACCS) == []


REPS = pd.DataFrame({"species_prefix": ["AAA", "BBB", "CCC"], "PHYLUM": ["Ascomycota", "Ascomycota", "Mucoromycota"],
                     "ORDER": ["Eurotiales", "Saccharomycetales", "Mucorales"], "GENUS": ["g1", "g2", "g3"]})


def _pf():
    rows = [("AAA_1-T1", "PF00320.25", 10, 60), ("AAA_1-T1", "PF00989.30", 100, 200), ("AAA_2-T1", "PF00320.25", 5, 55),
            ("BBB_1-T1", "PF00096.30", 1, 30), ("CCC_1-T1", "PF00320.25", 1, 50), ("CCC_1-T1", "PF08447.14", 90, 190),
            ("CCC_1-T1", "PF13426.12", 400, 500)]
    return pd.DataFrame(rows, columns=["protein_id", "pfam_acc", "ali_from", "ali_to"]).assign(pfam_id="x")


def test_select_candidates_includes_genera_without_hits_as_absent():
    seeds = pd.DataFrame({"target": ["WC-2"], "target_pfam": ["PF00320"], "required_pfam": [PAS]})
    c = ors.select_candidates(REPS, _pf(), "WC-2", seeds)
    prot = c[c.protein_id.notna()]
    assert sorted(prot.protein_id) == ["AAA_1-T1", "AAA_2-T1", "CCC_1-T1"]
    assert prot.set_index("protein_id").loc["AAA_1-T1", "architecture"] == "PF00320>PF00989"
    ps = prot.set_index("protein_id")
    assert ps.loc["AAA_1-T1", "n_pas_regions"] == 1 and ps.loc["AAA_2-T1", "n_pas_regions"] == 0
    assert ps.loc["CCC_1-T1", "n_pas_regions"] == 2
    assert ps.loc["AAA_1-T1", "gata_rel_pos"] == pytest.approx(10 / 200)
    bbb = c[c.species_prefix == "BBB"]
    assert len(bbb) == 1 and not bool(bbb.iloc[0]["present"]) and pd.isna(bbb.iloc[0]["protein_id"])
    assert c[c.species_prefix == "AAA"]["present"].all()


def test_search_skipped_when_mmseqs_missing(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(shutil, "which", lambda name: None)
    f = tmp_path / "x.faa"
    f.write_text(">a\nMKV\n")
    df = ors.search(f, f, tmp_path / "o.tsv", tmp_path)
    assert df.empty and "mmseqs" in df.attrs["skipped"]
    assert "mmseqs" in capsys.readouterr().out
    assert df.attrs["min_pident"] == 30.0 and df.attrs["min_cov"] == 0.5


def test_thresholds_recorded_in_output(tmp_path):
    df = pd.DataFrame({"query": ["s"], "target": ["t"], "pident": [50.0]})
    p = tmp_path / "h.tsv"
    ors.write_commented(df, p, {"min_pident": 30.0, "min_cov": 0.5})
    first = p.read_text().splitlines()[0]
    assert first.startswith("# ") and "min_pident=30.0" in first and "min_cov=0.5" in first
    back = pd.read_csv(p, sep="\t", comment="#")
    assert list(back.columns) == ["query", "target", "pident"]


def test_pass_hits_applies_architecture_after_search():
    cand = pd.DataFrame({"species_prefix": ["AAA", "AAA", "CCC"], "protein_id": ["AAA_1-T1", "AAA_2-T1", "CCC_1-T1"],
                         "architecture": ["PF00320>PF00989", "PF00320", "PF00320>PF08447"], "present": True,
                         "n_pas_regions": [1, 0, 1]})
    hits = pd.DataFrame({"query": ["s", "s", "s"], "target": cand.protein_id, "pident": [60.0, 55.0, 40.0],
                         "qcov": 0.9, "tcov": 0.9, "evalue": 1e-9, "bits": [300.0, 250.0, 120.0]})
    out = ors.pass_hits(hits, cand, [PAS])
    assert sorted(out.protein_id) == ["AAA_1-T1", "CCC_1-T1"]      # AAA_2 hit the seed but has no PAS domain
    assert out.set_index("protein_id").loc["AAA_1-T1", "seed"] == "s"


def test_assign_best_seed_reports_disagreement_instead_of_dropping():
    cand = pd.DataFrame({"species_prefix": ["A", "A", "B", "B"], "protein_id": ["A_1-T1", "A_2-T1", "B_1-T1", "B_2-T1"],
                         "architecture": "PF00320>PF00989", "present": True, "n_pas_regions": [1, 3, 1, 3]})
    hits = pd.DataFrame({"query": ["wc2seed", "wc2seed", "wc1seed", "wc2seed", "wc1seed"],
                         "target": ["A_1-T1", "A_2-T1", "A_2-T1", "B_1-T1", "B_1-T1"],
                         "pident": 50.0, "qcov": 0.9, "tcov": 0.9, "evalue": 1e-9, "bits": [300.0, 200.0, 250.0, 100.0, 400.0]})
    rules = {"WC-2": ["PAS_regions==1"], "WC-1": ["PAS_regions>=2"]}
    st = {"wc2seed": "WC-2", "wc1seed": "WC-1"}
    a = ors.assign_best_seed(hits, cand, st, rules).set_index("protein_id")
    assert a.loc["A_1-T1", "seed_target"] == "WC-2" and a.loc["A_1-T1", "rule_ok"]
    assert a.loc["A_2-T1", "seed_target"] == "WC-1" and a.loc["A_2-T1", "rule_ok"]          # best seed is WC-1 (250 > 200)
    assert a.loc["B_1-T1", "seed_target"] == "WC-1" and not a.loc["B_1-T1", "rule_ok"]      # WC-1 seed but one PAS region
    assert a.loc["B_1-T1", "arch_type"] == "WC-2" and a.loc["B_1-T1", "disagree"]
    assert "B_2-T1" not in a.index                                                          # no hit: not assigned
    sel = ors.pass_for_target(a.reset_index(), "WC-1")
    assert sorted(sel.protein_id) == ["A_2-T1"]
    assert ors.count_disagreements(a.reset_index(), "WC-1") == 1


def test_stamp_changes_with_config_and_invalidates_outputs(tmp_path):
    root = tmp_path / "root"
    (root / "config").mkdir(parents=True)
    cfg = root / "config" / "tf_seeds.tsv"
    cfg.write_text("a\n")
    stage = tmp_path / "stage"
    stage.mkdir()
    assert ors.invalidate_if_changed(stage, root, code_files=[]) is True        # first run: no stamp yet
    for n in ("seeds.tsv", "pass_wc2.tsv", "hits_wc.tsv", "tree_wc2.nwk", "counts.tsv", "pfam_candidates.parquet"):
        (stage / n).write_text("x")
    assert ors.invalidate_if_changed(stage, root, code_files=[]) is False       # same config: outputs stay
    assert (stage / "seeds.tsv").exists()
    cfg.write_text("b\n")
    assert ors.invalidate_if_changed(stage, root, code_files=[]) is True
    assert not (stage / "seeds.tsv").exists() and not (stage / "pass_wc2.tsv").exists() and not (stage / "tree_wc2.nwk").exists()
    assert (stage / "pfam_candidates.parquet").exists()                          # the Pfam cache is not an output of the config
    assert ors.invalidate_if_changed(stage, root, code_files=[]) is False


def test_stamp_includes_code_and_thresholds(tmp_path):
    f = tmp_path / "code.py"
    f.write_text("x=1")
    h1 = ors.config_hash(tmp_path / "none.tsv", [f])
    f.write_text("x=2")
    assert ors.config_hash(tmp_path / "none.tsv", [f]) != h1
    assert ors.config_hash(tmp_path / "none.tsv", [f], min_cov=0.6) != ors.config_hash(tmp_path / "none.tsv", [f])


def test_link_reuse_replaces_dangling_symlink(tmp_path):
    real = tmp_path / "real.tsv"
    real.write_text("r")
    link = tmp_path / "link.tsv"
    link.symlink_to(tmp_path / "gone.tsv")
    assert link.is_symlink() and not link.exists()
    ors.link_reuse(link, real)
    assert link.read_text() == "r"
    ors.link_reuse(link, real)          # second call: no error


# ---------------------------------------------------------------- report
def test_report_prevalence_and_copies_known_answer():
    hits = pd.DataFrame({"species_prefix": ["AAA", "AAA", "AAA", "CCC"],
                         "protein_id": ["AAA_1-T1", "AAA_2-T1", "AAA_3-T1", "CCC_1-T1"],
                         "architecture": ["PF00320>PF00989"] * 3 + ["PF00320>PF08447"]})
    s = orp.summarise(hits, REPS, target="WC-2")
    assert list(s.columns) == ["target", "lineage", "parent", "n_genera", "n_with_ortholog", "prevalence",
                               "median_copies", "architecture_counts"]
    r = s.set_index("lineage")
    assert r.loc["Ascomycota", "n_genera"] == 2 and r.loc["Ascomycota", "n_with_ortholog"] == 1
    assert r.loc["Ascomycota", "prevalence"] == 0.5 and r.loc["Ascomycota", "median_copies"] == 3
    assert r.loc["Mucoromycota", "prevalence"] == 1.0 and r.loc["Mucoromycota", "median_copies"] == 1
    assert r.loc["Saccharomycetales", "n_with_ortholog"] == 0 and r.loc["Saccharomycetales", "parent"] == "Ascomycota"
    assert r.loc["Eurotiales", "median_copies"] == 3
    assert "PF00320>PF00989:3" in r.loc["Eurotiales", "architecture_counts"]
    assert (s["target"] == "WC-2").all()


def test_summarise_empty_hits_all_absent():
    s = orp.summarise(pd.DataFrame(columns=["species_prefix", "protein_id", "architecture"]), REPS, target="bZip5")
    assert (s["n_with_ortholog"] == 0).all() and (s["prevalence"] == 0).all() and len(s) == 5


# ---------------------------------------------------------------- tree
def test_build_tree_skips_with_clear_message_when_mafft_missing(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(shutil, "which", lambda name: None)
    f = tmp_path / "a.faa"
    f.write_text(">a\nMKV\n>b\nMKI\n")
    assert ort.build_tree(f, tmp_path / "t", tmp_path) is None
    assert "mafft" in capsys.readouterr().out


def test_trim_columns_drops_gappy_columns():
    aln = {"a": "AC-GT", "b": "AC-GT", "c": "A--GA", "d": "A--G-"}
    out = ort.trim_gappy(aln, max_gap=0.5)
    assert out == {"a": "ACGT", "b": "ACGT", "c": "A-GA", "d": "A-G-"}   # the all-gap column goes; 50% gaps stays


def test_sample_tips_one_per_order_plus_seeds_and_anchors():
    rows = [{"tip": f"p{i}", "phylum": "Ascomycota", "order": f"O{i % 40}", "anchor": i == 7, "is_seed": False, "bits": float(i)}
            for i in range(300)]
    rows += [{"tip": "S1", "phylum": "Ascomycota", "order": "", "anchor": False, "is_seed": True, "bits": float("inf")}]
    m = ort.sample_tips(pd.DataFrame(rows), max_tips=150)
    fig = m[m.in_figure]
    assert len(m) == 301 and len(fig) == 42            # 40 order representatives + the seed + the anchor tip p7
    assert "S1" in set(fig.tip) and "p7" in set(fig.tip)                      # seeds and anchors always drawn
    non_fixed = fig[~fig.is_seed & (fig.tip != "p7")]
    assert non_fixed["order"].is_unique                                         # one tip per order
    assert fig[fig["order"] == "O3"].tip.tolist() == ["p283"]                   # the highest bit score of the order
    small = ort.sample_tips(pd.DataFrame(rows), max_tips=20)
    assert small.in_figure.sum() == 20


def _toy_fasta(path, n=12, seed=1):
    rnd = random.Random(seed)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    base = "".join(rnd.choice(aa) for _ in range(120))
    with open(path, "w") as f:
        for i in range(n):
            s = list(base)
            for j in rnd.sample(range(120), 12 + i):
                s[j] = rnd.choice(aa)
            f.write(f">t{i}\n{''.join(s)}\n")


@pytest.mark.skipif(not all(shutil.which(x) for x in ("mmseqs", "mafft", "FastTree")),
                    reason="mmseqs, mafft or FastTree is not on PATH (module load MMseqs2/17-b804f mafft/7.505 fasttree/2.1.11)")
def test_integration_search_align_tree_12_tips(tmp_path):
    fa = tmp_path / "toy.faa"
    _toy_fasta(fa)
    hits = ors.search(fa, fa, tmp_path / "hits.tsv", tmp_path, min_pident=30.0, min_cov=0.5)
    assert set(hits["query"]) == {f"t{i}" for i in range(12)}
    nwk = ort.build_tree(fa, tmp_path / "toy", tmp_path)
    from Bio import Phylo
    tree = Phylo.read(str(nwk), "newick")
    assert len(tree.get_terminals()) == 12
    meta = pd.DataFrame({"tip": [f"t{i}" for i in range(12)], "phylum": ["Ascomycota"] * 6 + ["Mucoromycota"] * 6,
                         "anchor": [False] * 11 + [True]})
    out = tmp_path / "tree.png"
    ort.draw_tree(nwk, meta, {"Ascomycota": "#0072B2", "Mucoromycota": "#D55E00"}, out)
    assert out.stat().st_size > 1000


# ---------------------------------------------------------------- make_figures hook
def _prebuilt(root, fig="T21_x", level="tf"):
    (root / "figures" / level).mkdir(parents=True)
    (root / "figures" / level / f"{fig}.png").write_bytes(b"png")
    (root / "figures" / level / f"{fig}.svg").write_text("<svg/>")
    (root / "figures" / level / f"{fig}.json").write_text(json.dumps(
        {"name": fig, "level": level, "kind": "figure", "topic": "tf", "scope": "genus_reps", "n_genomes": None,
         "n_groups": 3, "pooled_share": 0.0, "signal": 4.0, "signal_kind": "n_tfs_tested"}))
    (root / "tables" / level).mkdir(parents=True)
    (root / "tables" / level / "T20_ortho_wc2.tsv").write_text("a\tb\n1\t2\n")
    (root / "tables_tf").mkdir()
    (root / "tables_tf" / "T21_target_conservation.tsv").write_text("c\n1\n")


def test_make_figures_hook_copies_prebuilt_and_indexes(tmp_path):
    import make_figures as mf
    from figure_index import build_index
    pre = tmp_path / "analysis" / "tf_targets"
    _prebuilt(pre)
    out = tmp_path / "out"
    made = mf.import_tf_extra(out, [pre])
    assert (out / "figures" / "tf" / "T21_x.png").exists() and (out / "figures" / "tf" / "T21_x.json").exists()
    assert (out / "tables" / "tf" / "T20_ortho_wc2.tsv").exists()
    assert (out / "tables" / "tf" / "T21_target_conservation.tsv").exists()
    assert [str(m) for m in made] == ["T21_x"]
    idx = build_index(out)
    row = idx[idx["name"] == "T21_x"].iloc[0]
    assert row["level"] == "tf" and row["topic"] == "tf" and row["scope"] == "genus_reps"


def test_make_figures_hook_missing_dir_prints_one_note(tmp_path, capsys):
    import make_figures as mf
    made = mf.import_tf_extra(tmp_path / "out", [tmp_path / "nope1", tmp_path / "nope2"])
    assert made == []
    assert capsys.readouterr().out.count("not found") == 1


def test_make_figures_extra_dirs_resolution(tmp_path, monkeypatch):
    import make_figures as mf
    data = tmp_path / "analysis" / "summary_tables"
    monkeypatch.delenv("BFD_TF_EXTRA_DIRS", raising=False)
    assert [d.name for d in mf.tf_extra_dirs(data)] == ["tf_targets", "tf_orthologs"]
    monkeypatch.setenv("BFD_TF_EXTRA_DIRS", f"{tmp_path}/a:{tmp_path}/b")
    assert [d.name for d in mf.tf_extra_dirs(data)] == ["a", "b"]
    assert [d.name for d in mf.tf_extra_dirs(data, "none")] == []
    assert [d.name for d in mf.tf_extra_dirs(data, f"{tmp_path}/c")] == ["c"]


# ---------------------------------------------------------------- config and workflow
def test_seed_config_columns_and_workflow_script():
    cfg = sd.read_seed_config(SEED_CFG)
    assert list(cfg.columns) == ["target", "seed_name", "source", "key", "taxon_id", "target_pfam", "required_pfam", "note"]
    assert set(cfg["source"]) == {"paper", "uniprot", "bfd"}
    sh = HERE.parents[1] / "workflows" / "tf_orthologs.sh"
    txt = sh.read_text()
    assert "BASH_SOURCE" not in txt and "${SCRATCH:?}" in txt and "-J tf_ortho" in txt and "--mem 32G" in txt
    for m in ("MMseqs2/17-b804f", "mafft/7.505", "fasttree/2.1.11"):
        assert f"module load {m}\n" in txt
    assert os.access(sh, os.X_OK)


def test_hook_skips_symlinks_and_warns_on_differing_clash(tmp_path, capsys):
    import make_figures as mf
    pre = tmp_path / "pre"
    _prebuilt(pre)
    outside = tmp_path / "secret.txt"
    outside.write_text("secret")
    (pre / "tables" / "tf" / "link.tsv").symlink_to(outside)
    out = tmp_path / "out"
    (out / "tables" / "tf").mkdir(parents=True)
    (out / "tables" / "tf" / "T20_ortho_wc2.tsv").write_text("old bytes\n")       # clash with different content
    (out / "tables" / "tf" / "extra_keep.tsv").write_text("keep")
    mf.import_tf_extra(out, [pre])
    msg = capsys.readouterr().out
    assert not (out / "tables" / "tf" / "link.tsv").exists()                       # symlink not copied
    assert "T20_ortho_wc2.tsv" in msg and "differs" in msg                          # clash warning names the file
    assert (out / "tables" / "tf" / "T20_ortho_wc2.tsv").read_text() == "a\tb\n1\t2\n"   # overwritten with the new bytes
    assert (out / "tables" / "tf" / "extra_keep.tsv").read_text() == "keep"        # nothing removed
    mf.import_tf_extra(out, [pre])                                                  # identical bytes: no new warning
    assert "differs" not in capsys.readouterr().out


def test_gap_check_overlap_with_gene_models():
    import ortho_gap_check as gc
    genes = pd.DataFrame({"scaffold": ["s1", "s1", "s2"], "start": [100, 1000, 50], "end": [500, 1500, 90],
                          "gene_id": ["g1", "g2", "g3"], "protein_id": ["g1-T1", "g2-T1", "g3-T1"]})
    hits = pd.DataFrame({"target": ["s1", "s1", "s2", "s3"], "tstart": [450, 700, 90, 10], "tend": [300, 800, 120, 40]})
    out = gc.overlap_genes(hits, genes)
    assert out["gene_id"].tolist() == ["g1", "", "g3", ""]            # in model g1; between models; 1 aa overlap with g3; other scaffold
    assert out["protein_id"].tolist() == ["g1-T1", "", "g3-T1", ""]


def test_cov_mode_per_target_zap1_strict_others_seed_coverage():
    assert ors.cov_mode_for("Zap1-like") == 0              # both-sequence coverage: C2H2 repeats alone must not pass
    assert ors.cov_mode_for("WC-2") == 2 and ors.cov_mode_for("WC-1") == 2 and ors.cov_mode_for("bZip5") == 2
    assert ors.LOOSE_MODE == 2                              # sensitivity run only
    assert ors.variants("zap1") == [("", 0), ("_loose", 2)]
    assert ors.variants("wc") == [("", 2)] and ors.variants("bzip5") == [("", 2)]


def test_loose_zap1_table_has_note_and_no_figure_slug(tmp_path):
    hits = pd.DataFrame({"species_prefix": ["AAA"], "protein_id": ["AAA_1-T1"], "architecture": ["PF00096>PF00096>PF00096>PF00096"]})
    p = tmp_path / "T20_ortho_zap1_loose.tsv"
    orp.write_loose_table(hits, REPS, p, cov_mode=2)
    first = p.read_text().splitlines()[0]
    assert first.startswith("# ") and "cov_mode=2" in first and "C2H2 repeat" in first and "sensitivity" in first
    t = pd.read_csv(p, sep="\t", comment="#")
    assert list(t.columns) == orp.SUMMARY_COLS and t.loc[t.lineage == "Ascomycota", "n_with_ortholog"].iloc[0] == 1
