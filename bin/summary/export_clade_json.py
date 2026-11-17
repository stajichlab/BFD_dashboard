#!/usr/bin/env python3
"""Export the per-genome table as plain JSON for the Explore page (Observable JS).

    python3 export_clade_json.py TABLES_DIR OUT_JSON

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
from make_figures import OTHER, assign_groups  # noqa: E402

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


def export_clade_json(core: pd.DataFrame, structure: pd.DataFrame, out: Path) -> Path:
    """Write the genome records and the colour map; return the path."""
    out = Path(out)
    d = core.reset_index(drop=True).copy()
    colors, group_of = {}, {}
    cols_by_level = {}
    for lv in LEVELS:
        g, _levels, cmap = assign_groups(d, lv)
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
        "meta": {"n_genomes": n, "unit": "one record per genome", "floats": "rounded to 3 decimals; null = missing"},
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
    if len(argv) != 3:
        sys.exit(__doc__)
    t = Path(argv[1])
    core = pd.read_parquet(t / "genome_all.parquet")
    st = pd.read_parquet(t / "gene_structure_genome.parquet", columns=["species_prefix", "median_protein_len_aa"])
    p = export_clade_json(core, st, Path(argv[2]))
    print(f"{p}: {p.stat().st_size} bytes, {len(core)} records")


if __name__ == "__main__":
    main(sys.argv)
