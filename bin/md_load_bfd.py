#!/usr/bin/env python3
"""Load BFD parquet tables into a MotherDuck database.

Each parquet file (or directory of parquet parts) in --tables becomes one table.
Table names: file stem with '.' -> '_' (cazy.overview -> cazy_overview).

Restartable: a table that already exists in MotherDuck with the same row count
as its source is skipped, so a failed or timed-out job can simply be resubmitted.

Table comments come from --tables/_manifest/<name>.json ("description" key)
when present, else from a short default.

Token: read from $motherduck_token, else from --token-file
(default ~/.config/motherduck/token, must be chmod 600). Never pass it on the
command line.

Usage:
    bin/md_load_bfd.py --tables /bigdata/.../Fungi_BFD_runs/tables --db bfd --dry-run
    bin/md_load_bfd.py --tables /bigdata/.../Fungi_BFD_runs/tables --db bfd
    bin/md_load_bfd.py ... --only pfam species        # subset
    bin/md_load_bfd.py ... --include-swissprot         # also load swissprot_hits (~37 GB)
"""
import argparse
import fnmatch
import json
import os
import stat
import sys
import time
from pathlib import Path

import duckdb

# Excluded unless --include-swissprot / explicitly listed with --only.
DEFAULT_EXCLUDE = ["swissprot_hits", "*.bak*", "*.csv*", "*.md", "_*"]

# Genome keys: ASMID (assembly accession) and LOCUSTAG (= species_prefix = the
# prefix of every gene/transcript/protein id, e.g. 'F2D1C788_000724-T1' -> 'F2D1C788').
_G = "Genome key species_prefix = species.LOCUSTAG."
_ID = "Genome = split_part(<id>,'_',1) = species.LOCUSTAG."
DEFAULT_COMMENTS = {
    # per-genome
    "species": "One row per genome. Taxonomy PHYLUM..SPECIES (SUBPHYLUM is mostly empty; use CLASS). Keys: ASMID; LOCUSTAG = protein_id prefix = species_prefix.",
    "samples": "Per-genome sample metadata: input species name, strain, BioProject, NCBI taxid, BUSCO lineage, translation table, taxonomy. Keys ASMID, LOCUSTAG.",
    "asm_stats": "Assembly statistics per genome (ASMID): contig count, total length, N50/L50, N90/L90, GC, gaps, repeat-masked bases, T2T scaffolds, telomere counts, annotation status.",
    "busco_genome": "BUSCO completeness per genome (ASMID): complete/single/duplicated/fragmented/missing %. All scored with fungi_odb12; budding-yeast classes (Saccharomycetes, Pichiomycetes) score ~75-80% from real gene loss, not poor assembly.",
    "telomere_summary": "Telomere repeat summary per genome (ASMID, LOCUSTAG): scaffolds with telomeres at one/both ends, tract counts and lengths, 5'/3' and terminal/internal counts, dominant repeat monomer.",
    "telomere_tracts": "One row per telomeric repeat tract: scaffold, end (5'/3'), strand, monomer, repeat count, length, coordinates, whether terminal, distance to scaffold end, sequence. Keys ASMID, LOCUSTAG.",
    "aa_freq": "Proteome amino-acid composition: one row per genome x amino acid, frequency = fraction of all residues. " + _G,
    "codon_freq": "Codon usage: one row per genome x codon, frequency = fraction of all codons in CDS. " + _G,
    "wgd_ks_summary": "Paralog Ks summary per genome: pair/family/node counts, hybrid flag, identical-copy groups. " + _G,
    # gene models
    "gene_info": "One row per gene: coordinates (chrom, start, end, strand) and gene_type. LOCUSTAG joins species.",
    "gene_transcripts": "One row per transcript: gene_id, coordinates, partial flag, start/stop codon presence, CDS length, md5. " + _ID,
    "gene_proteins": "One row per protein: gene_id, protein_id, transcript_id, length (aa), md5 of sequence. " + _ID,
    "gene_exons": "One row per exon: transcript_id, exon order, coordinates, GC content. " + _ID,
    "gene_CDS": "One row per CDS segment: transcript_id, order, coordinates, GC content. " + _ID,
    "gene_introns": "One row per intron: transcript_id, intron number, coordinates, 5'/3' splice dinucleotides, GC, codon phase. " + _ID,
    "gene_intergenic_distances": "Distance (bp) between each pair of adjacent genes (left_gene, right_gene). " + _G,
    "gene_trnas": "tRNA genes: gene_id, amino acid, anticodon. " + _ID,
    "infernal_rfam": "Infernal cmscan hits to Rfam non-coding RNA models: model name/accession, contig, coordinates, strand, e-value, score. " + _G,
    # protein domains and families
    "pfam": "Pfam hmmsearch domain hits, one row per domain: pfam_id (name), pfam_acc, e-values, scores, HMM/alignment/envelope coordinates. Genome = split_part(protein_id,'_',1) = species.LOCUSTAG.",
    "tf_inventory": "Transcription-factor inventory: proteins with TF-defining Pfam domains, assigned tf_family, domain e-value and score. " + _G,
    "cazy_overview": "dbCAN CAZyme calls, one row per protein: CAZy family, subfamily, EC, predicted substrate, number of dbCAN tools agreeing (toolcount). " + _G,
    "cazy_cazymes_hmm": "dbCAN HMM hits to CAZy families: HMM_id, profile/protein lengths, e-value, coordinates, coverage. " + _G,
    "merops": "MEROPS peptidase BLAST hits: merops_id, % identity, alignment stats, e-value, bitscore. " + _G,
    # localization and protein properties
    "signalp_signal_peptide": "SignalP predicted signal peptides (secretion): peptide start/end, probability. Only proteins with a predicted signal peptide are listed. " + _G,
    "targetP": "TargetP targeting predictions: only proteins predicted SP (secretory signal peptide) or mTP (mitochondrial transit peptide) are listed; absent = no targeting peptide. Probability, cleavage site, motif. " + _G,
    "wolfpsort": "WoLF PSORT localization predictions: localization class and score per protein. " + _G,
    "predgpi": "PredGPI GPI-anchor predictions (GFF-style rows: feature, coordinates, score, attributes). " + _G,
    "tmhmm": "TMHMM transmembrane helix predictions per protein: length, expected TM residues (ExpAA), PredHel count, topology string. " + _G,
    "idp": "Intrinsically disordered regions: one row per region with start, end, length, mean disorder score. " + _G,
    "idp_summary": "Per-protein disorder summary: disordered residues, IDP_fraction, protein length. " + _G,
    # paralogs / WGD
    "wgd_ks": "Paralog pair Ks estimates (wgd): family, genes, dN, dS, dN/dS, alignment stats. dS is the synonymous divergence used for WGD peaks. " + _G,
    "wgd_ks_nodes": "wgd gene-family tree nodes: pair counts, pass/sample flags, mean/min/max dS, mean identity. " + _G,
    # SwissProt
    "swissprot_best": "Best SwissProt DIAMOND hit per protein, with subject annotation: accession, protein name, % identity/coverage x10 (divide by 10), bitscore, -log10 e-value. " + _G,
    "swissprot_hits": "All ranked SwissProt DIAMOND hits per protein (rk = rank). pident/qcov/scov stored x10. Join sp_id to swissprot_subject. " + _G,
    "swissprot_subject": "SwissProt subject proteins referenced by hits (sp_id): accession, name, gene, organism, GO terms, EC, InterPro, Pfam.",
    "swissprot_subject_summary": "Per SwissProt subject (sp_id): number of hits, number of genomes hit, best hits, max bitscore and identity (x10).",
    "swissprot_annot": "SwissProt annotation by accession: protein and gene names, organism, GO terms, EC numbers, InterPro, Pfam, length.",
}


def table_name(p: Path) -> str:
    stem = p.name.removesuffix(".parquet")
    return stem.replace(".", "_").replace("-", "_")


def discover(tables: Path, only, include_swissprot):
    exclude = [e for e in DEFAULT_EXCLUDE if not (include_swissprot and e == "swissprot_hits")]
    out = []
    for p in sorted(tables.iterdir()):
        is_file = p.is_file() and p.name.endswith(".parquet")
        is_dir = p.is_dir() and any(p.glob("*.parquet"))
        if not (is_file or is_dir):
            continue
        name = table_name(p)
        if only:
            if name not in only and p.name.removesuffix(".parquet") not in only:
                continue
        elif any(fnmatch.fnmatch(p.name, pat) for pat in exclude):
            continue
        src = str(p) if is_file else str(p / "*.parquet")
        size = p.stat().st_size if is_file else sum(f.stat().st_size for f in p.glob("*.parquet"))
        out.append((name, src, size, p))
    return out


def load_token(token_file: Path):
    if os.environ.get("motherduck_token"):
        return
    if not token_file.exists():
        sys.exit(f"No $motherduck_token and no token file at {token_file}")
    mode = token_file.stat().st_mode
    if mode & (stat.S_IRWXG | stat.S_IRWXO):
        sys.exit(f"{token_file} is readable by group/others; run: chmod 600 {token_file}")
    os.environ["motherduck_token"] = token_file.read_text().strip()


def comment_for(tables: Path, name: str, stem: str):
    """Default description plus build provenance from _manifest/<stem>.json if present."""
    text = DEFAULT_COMMENTS.get(name, "")
    m = tables / "_manifest" / f"{stem}.json"
    if m.exists():
        try:
            d = json.loads(m.read_text())
            text = d.get("description") or text
            prov = ", ".join(f"{k}={d[k]}" for k in ("built_at", "built_by", "schema_version") if k in d)
            if prov:
                text = f"{text} [{prov}]".strip()
        except json.JSONDecodeError:
            pass
    return text or None


def set_comment(c, tables: Path, name: str, p: Path):
    txt = comment_for(tables, name, p.name.removesuffix(".parquet"))
    if txt:
        c.sql(f"COMMENT ON TABLE \"{name}\" IS '" + txt.replace("'", "''") + "'")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", type=Path, required=True)
    ap.add_argument("--db", default="bfd")
    ap.add_argument("--only", nargs="*", default=None, help="table names to load (overrides excludes)")
    ap.add_argument("--include-swissprot", action="store_true")
    ap.add_argument("--token-file", type=Path, default=Path("~/.config/motherduck/token").expanduser())
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--memory", default="24GB")
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH", "/tmp") + "/md_duckdb_tmp")
    ap.add_argument("--dry-run", action="store_true", help="list what would load and the local row counts")
    a = ap.parse_args()

    todo = discover(a.tables, a.only, a.include_swissprot)
    total = sum(s for _, _, s, _ in todo)
    print(f"{len(todo)} tables, {total/1e9:.1f} GB parquet -> md:{a.db}", flush=True)

    local = duckdb.connect()
    if a.dry_run:
        for name, src, size, _ in todo:
            n = local.sql(f"SELECT count(*) FROM read_parquet('{src}')").fetchone()[0]
            print(f"  {name:36s} {size/1e9:7.2f} GB {n:>14,d} rows")
        return

    load_token(a.token_file)
    c = duckdb.connect("md:")
    c.sql(f"SET threads={a.threads}; SET memory_limit='{a.memory}'; SET temp_directory='{a.tmp}'")
    c.sql(f"CREATE DATABASE IF NOT EXISTS {a.db}")
    c.sql(f"USE {a.db}")
    existing = {r[0] for r in c.sql(
        f"SELECT table_name FROM duckdb_tables() WHERE database_name='{a.db}'").fetchall()}

    failed = []
    for name, src, size, p in todo:
        want = local.sql(f"SELECT count(*) FROM read_parquet('{src}')").fetchone()[0]
        if name in existing:
            have = c.sql(f"SELECT count(*) FROM {a.db}.main.\"{name}\"").fetchone()[0]
            if have == want:
                set_comment(c, a.tables, name, p)  # backfill/refresh comment on reruns
                print(f"skip  {name:36s} {want:>14,d} rows (already loaded)", flush=True)
                continue
            print(f"redo  {name:36s} have {have:,d} != {want:,d}", flush=True)
        t0 = time.time()
        try:
            c.sql(f"CREATE OR REPLACE TABLE \"{name}\" AS SELECT * FROM read_parquet('{src}')")
            have = c.sql(f"SELECT count(*) FROM \"{name}\"").fetchone()[0]
            if have != want:
                raise RuntimeError(f"row count {have:,d} != source {want:,d}")
            set_comment(c, a.tables, name, p)
            dt = time.time() - t0
            print(f"load  {name:36s} {want:>14,d} rows {size/1e9:7.2f} GB {dt:7.0f}s "
                  f"({size/1e6/max(dt,1):.0f} MB/s)", flush=True)
        except Exception as e:  # keep going; resubmit retries only what's missing
            failed.append(name)
            print(f"FAIL  {name}: {str(e)[:300]}", flush=True)

    if failed:
        sys.exit(f"{len(failed)} table(s) failed: {' '.join(failed)} -- resubmit to retry")
    print("done", flush=True)


if __name__ == "__main__":
    main()
