#!/usr/bin/env python3
"""Select the genomes and queries for the T1+T2 tblastn check, by the rules in README.md (no hit information is used).

Usage (HPCC job, local BFD.duckdb):  01_select_genomes.py OUT_DIR [REPO_ROOT]
Writes genomes.tsv, queries.tsv and selection_log.txt to OUT_DIR.

genomes.tsv columns 1,2,3,4,5 are locustag, genus, stratum, priority, asmid (the layout run_tblastn.sh of the Ku check
reads); further columns carry the species, class, order, BUSCO, N50, contig count and the Pfam TPS/TPP counts.
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE.parents[3]
sys.path.insert(0, str(REPO / "lib"))

N50_MIN = 50_000
N_LECAN_PRIMARY, N_LECAN_FRAG, N_SCEDO, N_OTHER = 10, 4, 5, 4
MAX_PER_GENUS = 2
UNIPROT_TPS = ["TPS1A_ASPFU", "TPS1_BOTFB", "TPS1_EMENI", "TPS1_CANAL", "TPS1_YEAST", "TPS1_SCHPO"]
UNIPROT_TPP = ["TPS2_YEAST", "TPS2_SCHPO", "TPS2_ZYGRO"]


def pick(df, n, per_genus=MAX_PER_GENUS):
    """Best BUSCO then N50 first, at most per_genus genomes of a genus, first n rows."""
    df = df.sort_values(["busco_complete", "N50_bp"], ascending=False)
    return df.groupby("GENUS", sort=False).head(per_genus).head(n)


def select(cand):
    """cand: annotated genomes with columns GENUS, SPECIES, CLASS, ORDER_, busco_complete, N50_bp, tps, tpp, credible."""
    out = []
    lec = cand[(cand.CLASS == "Lecanoromycetes") & cand.credible & (cand.tps == 0) & (cand.tpp == 0)]
    out.append(("lecan_primary", pick(lec[lec.N50_bp >= N50_MIN], N_LECAN_PRIMARY)))
    out.append(("lecan_fragmented", pick(lec[lec.N50_bp < N50_MIN], N_LECAN_FRAG)))
    sc = cand[cand.SPECIES.str.match(r"^(Scedosporium|Lomentospora) [a-z]+$") & cand.credible]
    sc = sc.sort_values(["busco_complete", "N50_bp"], ascending=False).groupby("SPECIES", sort=False).head(1)
    out.append(("scedo", sc.head(N_SCEDO)))
    used = set(pd.concat([d for _, d in out]).LOCUSTAG)
    oth = cand[cand.credible & ~cand.LOCUSTAG.isin(used) & (cand.CLASS != "Lecanoromycetes")
               & ~cand.SPECIES.str.match(r"^(Scedosporium|Lomentospora)")]
    out.append(("other_loss", pick(oth, N_OTHER, per_genus=1)))
    ctrl = cand[(cand.tps >= 1) & (cand.tpp >= 1) & (cand.N50_bp >= N50_MIN) & ~cand.credible]
    parts = [pick(ctrl[ctrl.CLASS == "Lecanoromycetes"], 5, per_genus=1)]
    sord = ctrl[ctrl.CLASS == "Sordariomycetes"]
    micro = sord[sord.ORDER_ == "Microascales"]
    parts.append(pick(micro if len(micro) >= 3 else sord, 3, per_genus=1))
    parts.append(pick(ctrl[ctrl.CLASS == "Eurotiomycetes"], 2, per_genus=1))
    out.append(("control_pos", pd.concat(parts)))
    return out


def main():
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    import bfd_connect
    con = bfd_connect.connect(scratch=False, threads=4, memory="24GB")
    cred = pd.read_csv(REPO / "docs/analyses/trehalose/data_j13/th_credible_loss_genomes.tsv", sep="\t")
    con.register("cred", cred[["SPECIES"]].drop_duplicates())
    con.execute("""
        create temp table c as
        select v.LOCUSTAG, v.ASMID, v.GENUS, v.SPECIES, v.CLASS, v.ORDER_, v.busco_complete, v.N50_bp, v.contig_count,
               (v.SPECIES in (select SPECIES from cred)) as credible
        from v_genome_labels v where v.annotation_status = 'annotated'
          and (v.CLASS in ('Lecanoromycetes', 'Sordariomycetes', 'Eurotiomycetes', 'Candelariomycetes', 'Coniocybomycetes'))""")
    con.execute("""
        create temp table pc as
        select species_prefix as LOCUSTAG,
               count(distinct case when split_part(pfam_acc, '.', 1) = 'PF00982' then protein_id end) as tps,
               count(distinct case when split_part(pfam_acc, '.', 1) = 'PF02358' then protein_id end) as tpp
        from pfam where domain_i_evalue <= 1e-5 and split_part(pfam_acc, '.', 1) in ('PF00982', 'PF02358')
          and species_prefix in (select LOCUSTAG from c) group by 1""")
    cand = con.execute("select c.*, coalesce(pc.tps, 0) as tps, coalesce(pc.tpp, 0) as tpp from c left join pc using (LOCUSTAG)").df()
    strata = select(cand)
    rows, log = [], []
    for name, d in strata:
        log.append(f"{name}: {len(d)} genomes, {d.GENUS.nunique()} genera, N50 median {d.N50_bp.median():,.0f}")
        for _, r in d.iterrows():
            rows.append({"locustag": r.LOCUSTAG, "genus": r.GENUS, "grp": name, "min_priority": 1, "asmid": r.ASMID,
                         "species": r.SPECIES, "class": r.CLASS, "order": r.ORDER_, "busco": r.busco_complete,
                         "n50_bp": int(r.N50_bp), "contigs": int(r.contig_count), "tps_pfam": int(r.tps), "tpp_pfam": int(r.tpp)})
    g = pd.DataFrame(rows)
    assert g.locustag.is_unique, "a genome fell in two strata"
    g.to_csv(out_dir / "genomes.tsv", sep="\t", index=False)

    # queries: UniProt entries are looked up by entry name in bfd.swissprot_subject (accession and exact length);
    # BFD proteins are the best-scoring PF00982 / PF02358 protein of the first control genomes
    names = ",".join(f"'{n}'" for n in UNIPROT_TPS + UNIPROT_TPP)
    sp = con.execute(f"select entry_name, accession, s_len from swissprot_subject where entry_name in ({names})").df().set_index("entry_name")
    missing = [n for n in UNIPROT_TPS + UNIPROT_TPP if n not in sp.index]
    assert not missing, f"not in swissprot_subject: {missing}"
    q = [("TPS", "uniprot", sp.loc[n, "accession"], int(sp.loc[n, "s_len"])) for n in UNIPROT_TPS]
    q += [("TPP", "uniprot", sp.loc[n, "accession"], int(sp.loc[n, "s_len"])) for n in UNIPROT_TPP]
    rad = pd.read_csv(REPO / "docs/analyses/kx2_ku_pathway/tblastn/queries.tsv", sep="\t")
    rad = rad[rad.key == "Rad51"].iloc[0]
    q.append(("Rad51", "uniprot", rad.accession, int(rad.expected_len)))
    ctrl = g[g.grp == "control_pos"]
    bfd_genomes = list(ctrl[ctrl["class"] == "Lecanoromycetes"].locustag[:2]) + \
        list(ctrl[ctrl["class"] == "Sordariomycetes"].locustag[:2]) + list(ctrl[ctrl["class"] == "Eurotiomycetes"].locustag[:1])
    for loc in bfd_genomes:
        for key, acc in (("TPS_bfd", "PF00982"), ("TPP_bfd", "PF02358")):
            p = con.execute(f"""select protein_id from pfam where species_prefix = '{loc}' and split_part(pfam_acc, '.', 1) = '{acc}'
                                and domain_i_evalue <= 1e-5 order by domain_score desc limit 1""").fetchone()
            assert p, f"no {acc} protein in {loc}"
            q.append((key, "bfd", p[0], "NA"))
    pd.DataFrame(q, columns=["key", "source", "accession", "expected_len"]).to_csv(out_dir / "queries.tsv", sep="\t", index=False)
    log.append(f"queries: {len(q)} ({sum(x[1] == 'uniprot' for x in q)} UniProt, {sum(x[1] == 'bfd' for x in q)} BFD)")
    (out_dir / "selection_log.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log))
    print(g[["locustag", "genus", "grp", "species", "busco", "n50_bp", "tps_pfam", "tpp_pfam"]].to_string())


if __name__ == "__main__":
    main()
