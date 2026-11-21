#!/usr/bin/env python3
"""Extract the domain-trimmed protein sequences of the carbon-storage families for one genome per genus (HPCC job).

Usage: 01_extract_domains.py OUT_DIR [REPO_ROOT]
Reads the local BFD.duckdb (read-only) and the proteome FASTA files (read-only). Writes OUT_DIR/fam/<FAM>.faa (header
>protein_id|locustag, the best-scoring domain envelope +-MARGIN aa) and OUT_DIR/tips.tsv
(family, protein_id, locustag, genus, class, phylum, trimmed_length, copy_index).
A protein shorter than MIN_REL of the family's median trimmed length is dropped.
"""
import csv
import gzip
import statistics
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE.parents[3]
sys.path.insert(0, str(REPO / "lib"))

FAMILIES = {"TPS": "PF00982", "TPP": "PF02358", "NTH": "PF01204", "UGP": "PF01704", "GSY": "PF05693", "GPH": "PF00343",
            "GDB": "PF06202", "GPD": "PF01210"}
MARGIN, MIN_REL = 10, 0.60
MAP = "/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/locustag_to_pep.tsv"


def trim_window(env_from, env_to, length, margin=MARGIN):
    """1-based inclusive window of the domain envelope +- margin, clipped to the protein."""
    return max(1, int(env_from) - margin), min(int(length), int(env_to) + margin)


def keep_lengths(lengths, min_rel=MIN_REL):
    """Indices of sequences at least min_rel of the median length."""
    med = statistics.median(lengths)
    return [i for i, n in enumerate(lengths) if n >= min_rel * med]


def read_wanted(path, ids):
    op = gzip.open if str(path).endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name in ids:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name in ids:
        seqs[name] = "".join(buf)
    return seqs


def main():
    import bfd_connect
    out = Path(sys.argv[1])
    (out / "fam").mkdir(parents=True, exist_ok=True)
    con = bfd_connect.connect(scratch=False, threads=8, memory="40GB")
    accs = ",".join(f"'{a}'" for a in FAMILIES.values())
    hits = con.execute(f"""
        with reps as (select LOCUSTAG as species_prefix, GENUS, CLASS, PHYLUM from v_genome_labels
                      where is_genus_rep and annotation_status = 'annotated')
        select p.species_prefix, p.protein_id, split_part(p.pfam_acc, '.', 1) as acc, p.env_from, p.env_to, r.GENUS, r.CLASS, r.PHYLUM
        from pfam p join reps r using (species_prefix)
        where p.domain_i_evalue <= 1e-5 and split_part(p.pfam_acc, '.', 1) in ({accs})
        qualify row_number() over (partition by p.protein_id, split_part(p.pfam_acc, '.', 1) order by p.domain_score desc) = 1
    """).df()
    print("domain hits:", len(hits), flush=True)
    acc2fam = {a: f for f, a in FAMILIES.items()}
    hits["fam"] = hits.acc.map(acc2fam)
    pep = {r["LOCUSTAG"]: r["pep_path"] for r in csv.DictReader(open(MAP), delimiter="\t")}
    recs = defaultdict(list)
    for loc, sub in hits.groupby("species_prefix"):
        seqs = read_wanted(pep[loc], set(sub.protein_id))
        for r in sub.itertuples():
            if r.protein_id not in seqs:
                raise SystemExit(f"{loc}: {r.protein_id} not found in {pep[loc]}")
            s = seqs[r.protein_id].rstrip("*")
            a, b = trim_window(r.env_from, r.env_to, len(s))
            recs[r.fam].append((r.protein_id, loc, r.GENUS, r.CLASS, r.PHYLUM, s[a - 1:b]))
    with open(out / "tips.tsv", "w") as tf:
        w = csv.writer(tf, delimiter="\t")
        w.writerow(["family", "protein_id", "locustag", "genus", "class", "phylum", "trimmed_length", "copy_index"])
        for fam in FAMILIES:
            rows = recs[fam]
            keep = keep_lengths([len(x[5]) for x in rows])
            dropped = len(rows) - len(keep)
            copies = defaultdict(int)
            with open(out / "fam" / f"{fam}.faa", "w") as fa:
                for i in keep:
                    pid, loc, genus, cls, ph, seq = rows[i]
                    copies[loc] += 1
                    fa.write(f">{pid}|{loc}\n{seq}\n")
                    w.writerow([fam, pid, loc, genus, cls, ph, len(seq), copies[loc]])
            print(f"{fam}: {len(keep)} sequences written, {dropped} dropped as short, median length "
                  f"{statistics.median([len(x[5]) for x in rows]):.0f}", flush=True)


if __name__ == "__main__":
    main()
