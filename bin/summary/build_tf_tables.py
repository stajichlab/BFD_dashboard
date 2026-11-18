#!/usr/bin/env python3
"""Per-genome transcription factor (TF) tables for the TF survey.

Outputs (ZSTD parquet, built for EVERY genome in genome_all; the suppress list is applied by the readers):
  tf_family_genome   species_prefix, family, n_proteins   (one row per genome and family with n >= 1)
  tf_summary_genome  species_prefix, n_proteins_total, n_tf_proteins, tf_fraction, n_families, shannon

Inputs: <tables_dir>/{tf_inventory,pfam,gene_proteins}.parquet, <summary_dir>/{genome_all,pfam_names}.parquet,
and the family map config/tf_family_map.tsv (pfam_acc, family, paper_family, note, pfam_id).

Counting rule. A protein counts once per family if it has at least one domain of the family.
Pfams that are in tf_inventory are counted from tf_inventory. The added Pfams are counted from pfam.parquet,
filtered by accession (version suffix removed). Both tables come from the same hmmsearch domtblout files,
so a protein is never in both for one accession.

Domain gate. nextflow/bin/merge_tf_domains.py (and sql/merge/tf.sql) keep every domain hit whose Pfam
accession is in the list. They set no e-value or score cut-off. The only gate is upstream: RUN_PFAM runs
`hmmsearch --cut_ga --noali` (Pfam gathering thresholds), and pfam.parquet holds the same hits. So the
same gate applies to the added Pfams: none beyond --cut_ga. In the data, tf_inventory has domain_i_evalue
up to 0.087 and domain_score from 15.8.

Genome key. gene_proteins has no species_prefix. The genome is split_part(protein_id, '_', 1)
(for example FBF4BF29_000001-T1 gives FBF4BF29). Never use left(protein_id, 8): some LOCUSTAGs have 9 characters.

n_proteins_total is the number of distinct proteins in gene_proteins (0 when the genome has none).
tf_fraction = n_tf_proteins / n_proteins_total; NaN when the proteome is empty.
n_tf_proteins counts distinct proteins with a TF domain of any mapped family.
shannon is the natural-log Shannon index of the family distribution p_f = n_f / sum(n_f) over the families
of the genome. A protein with domains of two families (for example Zn_clus and Fungal_trans) adds to both
n_f, so sum(n_f) can exceed n_tf_proteins. The index is 0 for one family and NaN for none.
"""
import argparse
import os
import sys
from pathlib import Path

import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MAP = ROOT / "config" / "tf_family_map.tsv"
COLS = ["pfam_acc", "family", "paper_family", "note", "pfam_id"]
DEFAULT_DAT = "/srv/projects/db/pfam/current/Pfam-A.hmm.dat"
ACC_RE = "^(PF[0-9]+)"


def load_family_map(path) -> pd.DataFrame:
    """Read the family map. '#' lines and blank lines are ignored. A duplicate accession is an error.

    Rows with an empty pfam_acc are kept (they list unmapped families) but are not counted.
    """
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"{path}: empty family map")
    header = lines[0].split("\t")
    if header != COLS:
        raise ValueError(f"{path}: header must be {COLS}, got {header}")
    rows = []
    for i, ln in enumerate(lines[1:], 2):
        f = ln.split("\t")
        if len(f) > len(COLS):
            raise ValueError(f"{path}: data row {i} has {len(f)} columns, expected {len(COLS)}")
        rows.append((f + [""] * len(COLS))[:len(COLS)])
    df = pd.DataFrame(rows, columns=COLS).apply(lambda s: s.str.strip())
    acc = df.loc[df["pfam_acc"] != "", "pfam_acc"]
    dup = sorted(acc[acc.duplicated()].unique())
    if dup:
        raise ValueError(f"{path}: duplicate pfam_acc {dup}")
    bad = df[(df["pfam_acc"] != "") & (df["family"] == "")]
    if len(bad):
        raise ValueError(f"{path}: empty family for {sorted(bad['pfam_acc'])}")
    return df.reset_index(drop=True)


def validate_accessions(fmap: pd.DataFrame, pfam_names: pd.DataFrame) -> None:
    """Raise ValueError listing every mapped accession that is not in pfam_names (versions stripped)."""
    known = set(pfam_names["pfam_acc"].astype(str).str.extract(ACC_RE)[0].dropna())
    want = [a for a in fmap["pfam_acc"] if a]
    missing = sorted(a for a in want if a not in known)
    if missing:
        raise ValueError(f"{len(missing)} accession(s) in the family map are missing from pfam_names: {missing}")


def check_pfam_ids(fmap: pd.DataFrame, dat_path) -> list:
    """Compare pfam_id of every mapped accession with `#=GF ID` in a Pfam-A.hmm.dat file.

    Returns one message per accession that is missing from the file or has a different ID. Empty list = all match.
    """
    ids, cur = {}, None
    with open(dat_path) as fh:
        for ln in fh:
            if ln.startswith("#=GF ID"):
                cur = ln.split(None, 2)[2].strip()
            elif ln.startswith("#=GF AC") and cur is not None:
                ids[ln.split()[2].split(".")[0]] = cur
                cur = None
    msgs = []
    for acc, pid in zip(fmap["pfam_acc"], fmap["pfam_id"]):
        if not acc:
            continue
        if acc not in ids:
            msgs.append(f"{acc}: not found in {dat_path}")
        elif ids[acc] != pid:
            msgs.append(f"{acc}: map pfam_id '{pid}' but Pfam ID is '{ids[acc]}'")
    return msgs


def build_tf_tables(con, tables_dir, summary_dir, fmap: pd.DataFrame, out) -> dict:
    """Write tf_family_genome.parquet and tf_summary_genome.parquet to out. Return a dict of counts."""
    tables_dir, summary_dir, out = Path(tables_dir), Path(summary_dir), Path(out)
    m = fmap.loc[fmap["pfam_acc"] != "", ["pfam_acc", "family"]].drop_duplicates()
    con.register("tf_map", m)
    inv = f"'{tables_dir / 'tf_inventory.parquet'}'"
    pfam = f"'{tables_dir / 'pfam.parquet'}'"
    gp = f"'{tables_dir / 'gene_proteins.parquet'}'"
    ga = f"'{summary_dir / 'genome_all.parquet'}'"
    # existing Pfams: tf_inventory (its species_prefix is split_part(protein_id, '_', 1) too)
    con.execute(f"""
        create or replace temp table tf_inv_acc as select distinct pfam_acc from {inv}""")
    con.execute(f"""
        create or replace temp table tf_hits as
        select distinct species_prefix, protein_id, family from (
            select i.species_prefix, i.protein_id, m.family
            from {inv} i join tf_map m on m.pfam_acc = i.pfam_acc
            union all
            select split_part(p.protein_id, '_', 1) as species_prefix, p.protein_id, m.family
            from (select protein_id, regexp_extract(pfam_acc, '{ACC_RE}', 1) as acc from {pfam}
                  where regexp_extract(pfam_acc, '{ACC_RE}', 1) in
                        (select pfam_acc from tf_map where pfam_acc not in (select pfam_acc from tf_inv_acc))) p
            join tf_map m on m.pfam_acc = p.acc
        )""")
    con.execute(f"""
        create or replace temp table tf_proteome as
        select split_part(protein_id, '_', 1) as species_prefix, count(distinct protein_id) as n_proteins_total
        from {gp} group by 1""")
    con.execute(f"create or replace temp table tf_genomes as select distinct species_prefix from {ga}")
    fam_path, sum_path = out / "tf_family_genome.parquet", out / "tf_summary_genome.parquet"
    con.execute(f"""
        copy (select h.species_prefix, h.family, count(*)::bigint as n_proteins
              from tf_hits h join tf_genomes g using (species_prefix)
              group by 1, 2 order by 1, 2) to '{fam_path}' (format parquet, compression zstd)""")
    con.execute(f"""
        copy (
          with fam as (select species_prefix, family, count(*)::double as n from tf_hits
                       where species_prefix in (select species_prefix from tf_genomes) group by 1, 2),
          tot as (select species_prefix, count(*)::bigint as n_families, sum(n) as tot from fam group by 1),
          sh as (select f.species_prefix, -sum((f.n / t.tot) * ln(f.n / t.tot)) as shannon
                 from fam f join tot t using (species_prefix) group by 1),
          tfp as (select species_prefix, count(distinct protein_id)::bigint as n_tf_proteins from tf_hits group by 1)
          select g.species_prefix,
                 coalesce(p.n_proteins_total, 0)::bigint as n_proteins_total,
                 coalesce(t.n_tf_proteins, 0)::bigint as n_tf_proteins,
                 case when coalesce(p.n_proteins_total, 0) > 0
                      then coalesce(t.n_tf_proteins, 0)::double / p.n_proteins_total else 'NaN'::double end as tf_fraction,
                 coalesce(tot.n_families, 0)::bigint as n_families,
                 case when tot.n_families is null then 'NaN'::double else sh.shannon + 0.0 end as shannon
          from tf_genomes g
          left join tf_proteome p using (species_prefix)
          left join tfp t using (species_prefix)
          left join tot using (species_prefix)
          left join sh using (species_prefix)
          order by 1) to '{sum_path}' (format parquet, compression zstd)""")
    n_genomes = con.execute("select count(*) from tf_genomes").fetchone()[0]
    n_zero = con.execute(f"select count(*) from '{sum_path}' where n_tf_proteins = 0").fetchone()[0]
    n_empty = con.execute(f"select count(*) from '{sum_path}' where n_proteins_total = 0").fetchone()[0]
    outside = con.execute(
        "select count(distinct species_prefix) from tf_hits where species_prefix not in (select species_prefix from tf_genomes)"
    ).fetchone()[0]
    hit_acc = _accessions_with_hits(con, tables_dir, m)
    want = sorted(m["pfam_acc"])
    con.unregister("tf_map")
    return {
        "n_genomes": int(n_genomes), "n_zero_tf": int(n_zero), "n_empty_proteome": int(n_empty),
        "tf_genomes_not_in_genome_all": int(outside),
        "n_family_rows": int(con.execute(f"select count(*) from '{fam_path}'").fetchone()[0]),
        "accessions_with_hits": [a for a in want if a in hit_acc],
        "accessions_without_hits": [a for a in want if a not in hit_acc],
    }


def _accessions_with_hits(con, tables_dir, m):
    inv = f"'{tables_dir / 'tf_inventory.parquet'}'"
    pfam = f"'{tables_dir / 'pfam.parquet'}'"
    a = con.execute(f"select distinct pfam_acc from {inv} where pfam_acc in (select pfam_acc from tf_map)").fetchall()
    b = con.execute(f"""select distinct regexp_extract(pfam_acc, '{ACC_RE}', 1) from {pfam}
                        where regexp_extract(pfam_acc, '{ACC_RE}', 1) in
                              (select pfam_acc from tf_map where pfam_acc not in (select pfam_acc from tf_inv_acc))""").fetchall()
    return {r[0] for r in a} | {r[0] for r in b}


def check_dat(fmap, pfam_dat):
    """Run the pfam_id guard as the CLI does. Exit with a message on a mismatch or on a missing explicit path."""
    if str(pfam_dat).lower() == "none":
        print("pfam_id check: disabled", file=sys.stderr)
    elif not Path(pfam_dat).exists():
        if str(pfam_dat) == str(DEFAULT_DAT):
            print(f"pfam_id check: skipped (dat file not found: {pfam_dat})", file=sys.stderr)
        else:
            sys.exit(f"--pfam-dat file not found: {pfam_dat} (use --pfam-dat none to disable the check)")
    else:
        msgs = check_pfam_ids(fmap, pfam_dat)
        if msgs:
            sys.exit("family map does not match the Pfam release:\n" + "\n".join(msgs))
        print(f"pfam_id check: ok ({int((fmap['pfam_acc'] != '').sum())} accessions)", file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("tables_dir", help="dir with tf_inventory, pfam, gene_proteins parquet")
    ap.add_argument("summary_dir", help="dir with genome_all.parquet and pfam_names.parquet")
    ap.add_argument("out_dir")
    ap.add_argument("--map", default=str(DEFAULT_MAP))
    ap.add_argument("--pfam-dat", default=str(DEFAULT_DAT),
                    help="Pfam-A.hmm.dat to check pfam_id against ('none' disables; the default path is skipped with a note if absent; an explicit missing path is an error)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--memory", default="24GB")
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH"), help="DuckDB temp dir (default $SCRATCH)")
    a = ap.parse_args(argv)
    fmap = load_family_map(a.map)
    names = pd.read_parquet(Path(a.summary_dir) / "pfam_names.parquet")
    validate_accessions(fmap, names)
    check_dat(fmap, a.pfam_dat)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(f"SET threads={a.threads}")
    con.execute(f"SET memory_limit='{a.memory}'")
    if a.tmp:
        con.execute(f"SET temp_directory='{a.tmp}'")
    res = build_tf_tables(con, a.tables_dir, a.summary_dir, fmap, out)
    for k, v in res.items():
        print(f"{k}: {v}", flush=True)
    unm = fmap[fmap["pfam_acc"] == ""]
    print(f"unmapped/not-counted map rows: {len(unm)}", flush=True)
    if res["accessions_without_hits"]:
        print(f"WARNING: no hits for {res['accessions_without_hits']}", file=sys.stderr)


if __name__ == "__main__":
    main()
