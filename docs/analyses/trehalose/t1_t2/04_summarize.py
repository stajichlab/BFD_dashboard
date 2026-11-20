#!/usr/bin/env python3
"""Summarize the T1+T2 tblastn hits by the rules fixed in README.md.

Usage: 04_summarize.py BENCH_DIR OUT_PREFIX      (BENCH_DIR holds hits/<locustag>.tsv, genomes.tsv, queries.faa)
Writes OUT_PREFIX_calls.tsv (per genome and gene), OUT_PREFIX_genomes.tsv (per genome outcome),
OUT_PREFIX_candidate_losses.tsv and OUT_PREFIX_verdicts.txt.

Calls (README): strong = best e-value <= 1e-10 and merged query coverage >= 30% on one scaffold; full-length = strong and
coverage >= 70%; weak = e-value <= 1e-5 but not strong; none otherwise.
Genes: TPS = keys TPS + TPS_bfd; TPP = key TPP_bfd; Rad51 = key Rad51. TPS2like is recorded but never used in a call.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

GENES = {"TPS": ("TPS", "TPS_bfd"), "TPP": ("TPP_bfd",), "Rad51": ("Rad51",)}
STRONG_E, WEAK_E, STRONG_COV, FULL_COV = 1e-10, 1e-5, 0.30, 0.70
N50_MIN, ASSAY_PASS_FRACTION, CLADE_FRACTION = 50_000, 0.8, 0.5
CLADE_STRATA = ("lecan_primary", "scedo")
TESTED = ("lecan_primary", "lecan_fragmented", "scedo", "other_loss")


def merged_span(intervals):
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def best_hit(rows):
    """rows: list of dicts (qseqid, sseqid, evalue, bitscore, qstart, qend, qlen) for one query set in one genome.
    Returns (query, evalue, coverage) of the query with the highest coverage x total bit score on its best scaffold."""
    by_q = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by_q[r["qseqid"]][r["sseqid"]].append(r)
    best = None
    for q, scaf in by_q.items():
        sid, hs = max(scaf.items(), key=lambda kv: sum(h["bitscore"] for h in kv[1]))
        qlen = hs[0]["qlen"]
        cov = merged_span([(min(h["qstart"], h["qend"]), max(h["qstart"], h["qend"])) for h in hs]) / qlen
        score = cov * sum(h["bitscore"] for h in hs)
        cand = (score, q, min(h["evalue"] for h in hs), cov)
        if best is None or cand > best:
            best = cand
    return None if best is None else (best[1], best[2], best[3])


def call(evalue, cov):
    if evalue is None:
        return "none"
    if evalue <= STRONG_E and cov >= STRONG_COV:
        return "full-length" if cov >= FULL_COV else "strong"
    return "weak" if evalue <= WEAK_E else "none"


def genome_outcome(calls, n50, assay_ok):
    """calls: {'TPS': call, 'TPP': call, 'Rad51': call} for a genome that Pfam called TPS- and TPP-absent."""
    informative = calls["Rad51"] in ("strong", "full-length") and n50 >= N50_MIN
    if "full-length" in (calls["TPS"], calls["TPP"]):
        return "annotation_gap", informative
    if informative and assay_ok["TPS"] and assay_ok["TPP"] and calls["TPS"] == "none" and calls["TPP"] == "none":
        return "candidate_loss", informative
    return "ambiguous", informative


def clade_verdict(outcomes):
    """outcomes: list of (outcome, informative). Pre-registered rule: denominator = informative genomes."""
    inf = [o for o, i in outcomes if i]
    if not inf:
        return "not evaluable (no informative genome)"
    loss, gap = inf.count("candidate_loss"), inf.count("annotation_gap")
    if loss >= CLADE_FRACTION * len(inf) and gap == 0:
        return f"supported ({loss} of {len(inf)} informative genomes are candidate losses, none an annotation gap)"
    if gap >= CLADE_FRACTION * len(inf):
        return f"killed ({gap} of {len(inf)} informative genomes show an annotation gap)"
    return f"mixed ({loss} candidate loss, {gap} annotation gap, {len(inf) - loss - gap} ambiguous of {len(inf)} informative)"


def read_hits(path):
    out = []
    if not Path(path).exists():
        return out
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        out.append({"qseqid": f[0], "sseqid": f[1], "evalue": float(f[4]), "bitscore": float(f[5]),
                    "qstart": int(f[6]), "qend": int(f[7]), "qlen": int(f[8])})
    return out


def main(bench, prefix):
    bench = Path(bench)
    genomes = list(csv.DictReader(open(bench / "genomes.tsv"), delimiter="\t"))
    rows = []
    calls = {}
    for g in genomes:
        hits = read_hits(bench / "hits" / f"{g['locustag']}.tsv")
        calls[g["locustag"]] = {}
        for gene, keys in GENES.items():
            sel = [h for h in hits if h["qseqid"].split("|")[0] in keys]
            b = best_hit(sel)
            c = call(b[1], b[2]) if b else "none"
            calls[g["locustag"]][gene] = c
            rows.append({"locustag": g["locustag"], "genus": g["genus"], "stratum": g["grp"], "gene": gene,
                         "call": c, "best_query": b[0] if b else "", "best_evalue": b[1] if b else "",
                         "query_coverage": round(b[2], 3) if b else ""})
    with open(f"{prefix}_calls.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    ctrl = [g for g in genomes if g["grp"] == "control_pos"]
    assay_ok = {}
    lines = []
    for gene in ("TPS", "TPP"):
        n_strong = sum(calls[g["locustag"]][gene] in ("strong", "full-length") for g in ctrl)
        assay_ok[gene] = bool(ctrl) and n_strong >= ASSAY_PASS_FRACTION * len(ctrl)
        lines.append(f"assay {gene}: strong in {n_strong} of {len(ctrl)} positive-control genomes -> "
                     f"{'passes' if assay_ok[gene] else 'FAILS: negatives for this gene are uninformative'}")
    rad = sum(calls[g["locustag"]]["Rad51"] in ("strong", "full-length") for g in genomes)
    lines.append(f"Rad51 strong in {rad} of {len(genomes)} genomes")

    per = []
    by_stratum = defaultdict(list)
    for g in genomes:
        if g["grp"] not in TESTED:
            continue
        out, inf = genome_outcome(calls[g["locustag"]], int(g["n50_bp"]), assay_ok)
        per.append({"locustag": g["locustag"], "genus": g["genus"], "species": g["species"], "stratum": g["grp"],
                    "n50_bp": g["n50_bp"], "TPS": calls[g["locustag"]]["TPS"], "TPP": calls[g["locustag"]]["TPP"],
                    "Rad51": calls[g["locustag"]]["Rad51"], "informative": inf, "outcome": out})
        by_stratum[g["grp"]].append((out, inf))
    with open(f"{prefix}_genomes.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per[0]), delimiter="\t")
        w.writeheader()
        w.writerows(per)
    cand = [p for p in per if p["outcome"] == "candidate_loss"]
    with open(f"{prefix}_candidate_losses.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per[0]), delimiter="\t")
        w.writeheader()
        w.writerows(cand)

    lines.append("")
    verdicts = {}
    for s in TESTED:
        if s in by_stratum:
            n = len(by_stratum[s])
            tally = {k: sum(o == k for o, _ in by_stratum[s]) for k in ("candidate_loss", "annotation_gap", "ambiguous")}
            lines.append(f"{s}: {n} genomes, {sum(i for _, i in by_stratum[s])} informative; {tally}")
            if s in CLADE_STRATA:
                verdicts[s] = clade_verdict(by_stratum[s])
                lines.append(f"  clade verdict {s}: {verdicts[s]}")
    both = all(verdicts.get(s, "").startswith("supported") for s in CLADE_STRATA)
    lines.append("")
    lines.append(f"T1 overall: {'supported' if both else 'NOT supported'} (supported only if both clades are supported; "
                 "a candidate loss is not proof)")
    lines.append(f"T2: evaluated only on the {len(cand)} candidate-loss genomes; lifestyle labels are joined in a separate step")
    open(f"{prefix}_verdicts.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
