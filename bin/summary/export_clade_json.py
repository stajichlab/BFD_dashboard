#!/usr/bin/env python3
"""Export the per-genome table as plain JSON for the Explore page (Observable JS).

    python3 export_clade_json.py TABLES_DIR OUT_JSON [--suppress PATH|none]

Genomes in config/summary_suppress.tsv are left out (see suppress.py). --suppress none keeps all genomes.

Regenerate the committed file (from the repository root):

    python3.12 bin/summary/export_clade_json.py <tables_dir> docs/quarto_summary/data/genomes.json

<tables_dir> is Fungi_BFD_runs/analysis/summary_tables, written by workflows/summary_tables.sh
(build_summary_tables.py).

TABLES_DIR holds genome_all.parquet (one row per genome) and gene_structure_genome.parquet.
Output: {"meta", "levels", "columns", "colors", "fields", "rows"}. "rows" is a list of arrays, one per genome,
in the order of "fields" (this keeps the file under 5 MB; the page turns each row into an object). Floats are rounded to 3 decimals,
NaN becomes null, missing taxonomy becomes "Unknown". Colours come from make_figures.assign_groups.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_figures import OTHER, assign_groups, phylum_order  # noqa: E402
from suppress import drop_suppressed, load_suppress  # noqa: E402

LEVELS = ("phylum", "subphylum", "class", "order")
# field -> (label, unit)
COLUMNS = {
    "genome_mb": ("Genome size", "Mb"),
    "n_genes": ("Gene count", "genes"),
    "gc": ("GC content", "%"),
    "genes_per_mb": ("Gene density", "genes/Mb"),
    "median_intron_len_bp": ("Median intron length", "bp"),
    "median_intergenic_bp": ("Median intergenic length", "bp"),
    "coding_fraction": ("Coding fraction", "fraction"),
    "busco_complete": ("BUSCO complete", "%"),
    "n50_bp": ("Assembly N50", "bp"),
    "median_protein_len_aa": ("Median protein length", "aa"),
}


def _num(s):
    return pd.to_numeric(s, errors="coerce").round(3).astype(object).where(lambda x: x.notna(), None)


REQUIRED = ["species_prefix", "PHYLUM", "CLASS", "ORDER", "GENUS", "is_representative", "total_length_bp",
            "gc_pct", "busco_complete", "N50_bp", "n_genes", "genes_per_mb", "coding_fraction",
            "median_intron_len_bp", "median_intergenic_bp"]


def export_clade_json(core: pd.DataFrame, structure: pd.DataFrame, out: Path,
                      generated_from: str = "summary_tables/genome_all.parquet", suppress=None) -> Path:
    """Write the genome records and the colour map; return the path.
    suppress: None = the repo list, a path, or a DataFrame (an empty one disables). Suppressed genomes are dropped
    from core and structure, and listed in meta."""
    out = Path(out)
    sup = suppress if isinstance(suppress, pd.DataFrame) else load_suppress(suppress)
    if "species_prefix" in core.columns:
        dropped = sorted(set(core["species_prefix"].astype(str)) & set(sup["species_prefix"].astype(str)))
    else:
        dropped = []
    core, _ = drop_suppressed(core, sup)
    structure, _ = drop_suppressed(structure, sup)
    miss = [c for c in REQUIRED if c not in core.columns]
    if "species_prefix" not in structure.columns or "median_protein_len_aa" not in structure.columns:
        miss += [c for c in ("species_prefix", "median_protein_len_aa") if c not in structure.columns]
    if miss:
        raise ValueError(f"export_clade_json: missing required columns: {', '.join(miss)}")
    d = core.reset_index(drop=True).copy()
    colors, group_of = {}, {}
    cols_by_level = {}
    porder = phylum_order(d)    # the same phylum colours as the figures
    for lv in LEVELS:
        g, _levels, cmap = assign_groups(d, lv, phylum_order=porder)
        # the plotted value per genome: the label itself when it has a colour, else pooled
        raw = g["subphylum_or_phylum"] if lv == "subphylum" else g["PHYLUM"] if lv == "phylum" else \
            d["CLASS" if lv == "class" else "ORDER"].fillna("Unknown")
        cols_by_level[lv] = raw.astype(str)
        colors[lv] = {v: cmap.get(v, cmap[OTHER]) if grp != OTHER else cmap[OTHER]
                      for v, grp in zip(raw.astype(str), g["grp"])}
    prot = structure.drop_duplicates("species_prefix").set_index("species_prefix")["median_protein_len_aa"]
    n = len(d)
    rec = pd.DataFrame({
        "id": d["species_prefix"].astype(str),
        "phylum": cols_by_level["phylum"], "subphylum": cols_by_level["subphylum"],
        "class": cols_by_level["class"], "order": cols_by_level["order"],
        "genus": d["GENUS"].fillna("Unknown").astype(str),
        "genome_mb": _num(d["total_length_bp"] / 1e6),
        "n_genes": _num(d["n_genes"]), "gc": _num(d["gc_pct"]), "genes_per_mb": _num(d["genes_per_mb"]),
        "median_intron_len_bp": _num(d["median_intron_len_bp"]),
        "median_intergenic_bp": _num(d["median_intergenic_bp"]),
        "coding_fraction": _num(d["coding_fraction"]),
        "busco_complete": _num(d["busco_complete"]), "n50_bp": _num(d["N50_bp"]),
        "median_protein_len_aa": _num(d["species_prefix"].map(prot)),
        "is_representative": d["is_representative"].map(lambda v: None if pd.isna(v) else bool(v)),
    })
    fields = list(rec.columns)

    def clean(v):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return None
        if isinstance(v, (np.floating, np.integer)):
            v = v.item()
        if isinstance(v, float) and v.is_integer():
            return int(v)          # whole numbers (gene counts, N50) are written without ".0"
        return v
    rows = [[clean(v) for v in r] for r in rec.itertuples(index=False, name=None)]
    obj = {
        "meta": {"n_genomes": n, "generated_from": generated_from, "n_suppressed": len(dropped), "suppressed": dropped,
                 "unit": "one record per genome",
                 "rows_layout": "array of arrays; column order = fields",
                 "floats": "rounded to 3 decimals; null = missing"},
        "levels": list(LEVELS),
        "columns": {k: {"label": l, "unit": u} for k, (l, u) in COLUMNS.items()},
        "colors": colors,
        "fields": fields,
        "rows": rows,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(obj, separators=(",", ":"), allow_nan=False))
    return out


def main(argv):
    argv = list(argv)
    suppress = None
    if "--suppress" in argv:
        i = argv.index("--suppress")
        val = argv[i + 1] if i + 1 < len(argv) else sys.exit(__doc__)
        suppress = pd.DataFrame({"species_prefix": []}) if val == "none" else val
        del argv[i:i + 2]
    if len(argv) != 3:
        sys.exit(__doc__)
    t = Path(argv[1])
    core = pd.read_parquet(t / "genome_all.parquet")
    st = pd.read_parquet(t / "gene_structure_genome.parquet", columns=["species_prefix", "median_protein_len_aa"])
    p = export_clade_json(core, st, Path(argv[2]), suppress=suppress)
    n = json.loads(p.read_text())["meta"]
    print(f"{p}: {p.stat().st_size} bytes, {n['n_genomes']} records, {n['n_suppressed']} suppressed")


if __name__ == "__main__":
    main(sys.argv)
