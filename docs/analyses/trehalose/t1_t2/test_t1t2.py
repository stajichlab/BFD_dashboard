"""Fixture tests for the T1+T2 summary rules (run on the HPCC: python3.12 -m pytest -q test_t1t2.py)."""
import importlib.util
from pathlib import Path

import pandas as pd

spec = importlib.util.spec_from_file_location("summ", Path(__file__).with_name("04_summarize.py"))
summ = importlib.util.module_from_spec(spec)
spec.loader.exec_module(summ)


def hit(q, scaf, e, bits, qs, qe, qlen):
    return {"qseqid": q, "sseqid": scaf, "evalue": e, "bitscore": bits, "qstart": qs, "qend": qe, "qlen": qlen}


def test_merged_span_counts_split_gene_once():
    assert summ.merged_span([(1, 100), (50, 150), (300, 400)]) == 150 + 101


def test_call_thresholds():
    assert summ.call(1e-30, 0.95) == "full-length"
    assert summ.call(1e-30, 0.50) == "strong"
    assert summ.call(1e-30, 0.20) == "weak"          # e-value good but coverage below 30%
    assert summ.call(1e-6, 0.90) == "weak"
    assert summ.call(1e-2, 0.90) == "none"
    assert summ.call(None, 0) == "none"


def test_best_hit_merges_hsps_on_best_scaffold():
    rows = [hit("TPS|A", "s1", 1e-40, 200, 1, 200, 500), hit("TPS|A", "s1", 1e-30, 150, 250, 450, 500),
            hit("TPS|A", "s2", 1e-5, 30, 1, 40, 500)]
    q, e, cov = summ.best_hit(rows)
    assert q == "TPS|A" and e == 1e-40 and abs(cov - 401 / 500) < 1e-9   # 200 + 201 residues (inclusive intervals)


def test_genome_outcomes():
    ok = {"TPS": True, "TPP": True}
    gap = {"TPS": "full-length", "TPP": "none", "Rad51": "strong"}
    assert summ.genome_outcome(gap, 1_000_000, ok)[0] == "annotation_gap"
    loss = {"TPS": "none", "TPP": "none", "Rad51": "strong"}
    assert summ.genome_outcome(loss, 1_000_000, ok) == ("candidate_loss", True)
    assert summ.genome_outcome(loss, 10_000, ok) == ("ambiguous", False)               # fragmented assembly
    assert summ.genome_outcome({**loss, "Rad51": "none"}, 1_000_000, ok)[0] == "ambiguous"   # search did not work
    assert summ.genome_outcome({**loss, "TPS": "weak"}, 1_000_000, ok)[0] == "ambiguous"
    assert summ.genome_outcome(loss, 1_000_000, {"TPS": False, "TPP": True})[0] == "ambiguous"  # assay failed for TPS


def test_clade_verdicts():
    sup = [("candidate_loss", True)] * 3 + [("ambiguous", True)]
    assert summ.clade_verdict(sup).startswith("supported")
    assert summ.clade_verdict(sup + [("annotation_gap", True)]).startswith("mixed")       # a gap blocks "supported"
    assert summ.clade_verdict([("annotation_gap", True)] * 2 + [("candidate_loss", True)]).startswith("killed")
    assert summ.clade_verdict([("candidate_loss", False)]).startswith("not evaluable")


def test_end_to_end(tmp_path):
    g = [("C1", "Cladonia", "control_pos"), ("C2", "Umbilicaria", "control_pos"), ("L1", "Physcia", "lecan_primary"),
         ("S1", "Scedosporium", "scedo")]
    pd.DataFrame([{"locustag": a, "genus": b, "grp": c, "min_priority": 1, "asmid": "x", "species": b + " sp", "n50_bp": 500000}
                  for a, b, c in g]).to_csv(tmp_path / "genomes.tsv", sep="\t", index=False)
    (tmp_path / "hits").mkdir()

    def write(loc, lines):
        (tmp_path / "hits" / f"{loc}.tsv").write_text("".join(lines))
    row = lambda q, e, b, qs, qe, ql: f"{q}\ts1\t60\t100\t{e}\t{b}\t{qs}\t{qe}\t{ql}\t1\t300\t+1\n"  # noqa: E731
    for loc in ("C1", "C2"):                       # positive controls: TPS, TPP and Rad51 all strong
        write(loc, [row("TPS|P1", 1e-80, 400, 1, 500, 515), row("TPP_bfd|F1_1", 1e-60, 300, 1, 230, 240),
                    row("Rad51|P25454", 1e-90, 500, 1, 380, 400)])
    write("L1", [row("Rad51|P25454", 1e-90, 500, 1, 380, 400)])           # TPS, TPP none -> candidate loss
    write("S1", [row("Rad51|P25454", 1e-90, 500, 1, 380, 400), row("TPS|P1", 1e-50, 300, 1, 450, 515)])  # TPS full length
    summ.main(tmp_path, str(tmp_path / "out"))
    res = pd.read_csv(tmp_path / "out_genomes.tsv", sep="\t").set_index("locustag")
    assert res.loc["L1", "outcome"] == "candidate_loss" and res.loc["S1", "outcome"] == "annotation_gap"
    assert "T1 overall: NOT supported" in (tmp_path / "out_verdicts.txt").read_text()


def test_t2_rule():
    spec2 = importlib.util.spec_from_file_location("life", Path(__file__).with_name("05_lifestyle.py"))
    life = importlib.util.module_from_spec(spec2)
    spec2.loader.exec_module(life)
    ft = pd.DataFrame({"GENUS": ["Peltigera", "Scopulariopsis"], "primary_lifestyle": ["lichenized", "soil_saprotroph"],
                       "Secondary_lifestyle": [None, None]})
    cand = pd.DataFrame({"genus": ["Peltigera", "Scopulariopsis", "Nolabel"]})
    lab = life.label_genomes(cand, ft)
    assert list(lab.primary_lifestyle) == ["lichenized", "soil_saprotroph", "no label"]
    assert list(lab.free_living_saprotroph) == [False, True, False]
    assert life.t2_verdict(lab).startswith("T2 KILLED")
    assert life.t2_verdict(lab.iloc[[0, 2]]).startswith("T2 not killed")
