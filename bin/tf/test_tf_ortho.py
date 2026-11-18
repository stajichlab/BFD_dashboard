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
    assert set(cfg["target"]) == {"WC-2", "Zap1-like", "bZip5"}
    assert set(cfg.loc[cfg.source == "paper", "key"]) == {"1799646", "1921537", "1909793"}


def test_bfd_seed_selector_gata_plus_pas():
    pf = pd.DataFrame({"protein_id": ["G_1-T1", "G_1-T1", "G_2-T1", "G_3-T1", "G_3-T1"],
                       "pfam_acc": ["PF00320.25", "PF08447.14", "PF00320.25", "PF00989.30", "PF00320.25"]})
    assert sd.bfd_selector_proteins("G", "pfam:GATA+PAS", pf) == ["G_1-T1", "G_3-T1"]


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


REPS = pd.DataFrame({"species_prefix": ["AAA", "BBB", "CCC"], "PHYLUM": ["Ascomycota", "Ascomycota", "Mucoromycota"],
                     "ORDER": ["Eurotiales", "Saccharomycetales", "Mucorales"], "GENUS": ["g1", "g2", "g3"]})


def _pf():
    rows = [("AAA_1-T1", "PF00320.25", 10), ("AAA_1-T1", "PF00989.30", 100), ("AAA_2-T1", "PF00320.25", 5),
            ("BBB_1-T1", "PF00096.30", 1), ("CCC_1-T1", "PF00320.25", 1), ("CCC_1-T1", "PF08447.14", 90)]
    return pd.DataFrame(rows, columns=["protein_id", "pfam_acc", "ali_from"]).assign(pfam_id="x")


def test_select_candidates_includes_genera_without_hits_as_absent():
    seeds = pd.DataFrame({"target": ["WC-2"], "target_pfam": ["PF00320"], "required_pfam": [PAS]})
    c = ors.select_candidates(REPS, _pf(), "WC-2", seeds)
    prot = c[c.protein_id.notna()]
    assert sorted(prot.protein_id) == ["AAA_1-T1", "AAA_2-T1", "CCC_1-T1"]
    assert prot.set_index("protein_id").loc["AAA_1-T1", "architecture"] == "PF00320>PF00989"
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
                         "architecture": ["PF00320>PF00989", "PF00320", "PF00320>PF08447"], "present": True})
    hits = pd.DataFrame({"query": ["s", "s", "s"], "target": cand.protein_id, "pident": [60.0, 55.0, 40.0],
                         "qcov": 0.9, "tcov": 0.9, "evalue": 1e-9, "bits": [300.0, 250.0, 120.0]})
    out = ors.pass_hits(hits, cand, [PAS])
    assert sorted(out.protein_id) == ["AAA_1-T1", "CCC_1-T1"]      # AAA_2 hit the seed but has no PAS domain
    assert out.set_index("protein_id").loc["AAA_1-T1", "seed"] == "s"


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
