"""Suppress list for the summary site.

Genomes listed in config/summary_suppress.tsv are left out of every summary figure, table, test and the
Explore page JSON. The list is applied when the parquet tables are read, so removing a line and re-running
make_figures.py / export_clade_json.py restores the genome. No table rebuild is needed.

File format: tab-separated, header species_prefix, species, reason, added, until.
Lines that start with # are comments. Blank lines are ignored.
"""
import io
from pathlib import Path

import pandas as pd

COLUMNS = ["species_prefix", "species", "reason", "added", "until"]
DEFAULT_PATH = Path(__file__).resolve().parents[2] / "config" / "summary_suppress.tsv"


def load_suppress(path=None) -> pd.DataFrame:
    """Read the suppress list. A missing file gives an empty list."""
    p = Path(path) if path is not None else DEFAULT_PATH
    if not p.exists():
        return pd.DataFrame(columns=COLUMNS)
    lines = [ln for ln in p.read_text().splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return pd.DataFrame(columns=COLUMNS)
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", dtype=str, keep_default_na=False)
    if "species_prefix" not in df.columns:
        raise ValueError(f"{p}: header must contain the column species_prefix")
    for c in COLUMNS:
        if c not in df.columns:
            df[c] = ""
    df = df[COLUMNS].apply(lambda s: s.str.strip())
    bad = df.index[df["species_prefix"] == ""]
    if len(bad):
        raise ValueError(f"{p}: empty species_prefix in data row(s) {[int(i) + 1 for i in bad]}")
    return df.reset_index(drop=True)


def drop_suppressed(df: pd.DataFrame, suppress: pd.DataFrame, key: str = "species_prefix"):
    """Return (df without suppressed rows, number of rows dropped). No-op when df has no key column."""
    if key not in df.columns or suppress is None or len(suppress) == 0:
        return df, 0
    mask = df[key].astype(str).isin(set(suppress["species_prefix"].astype(str)))
    n = int(mask.sum())
    return (df[~mask] if n else df), n
