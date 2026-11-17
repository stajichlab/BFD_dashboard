"""Figure index for the gallery: one row per figures/<level>/*.png, ranked by effect size.

    build_index(out_root) -> DataFrame, also written to OUT_ROOT/figures_index.tsv

Each PNG may have a sidecar JSON (written by make_figures.Out.save) with topic, stat, kind, scope, n_genomes,
n_groups and pooled_share. A PNG without a sidecar still gets a row, with those fields empty.
effect_size is the epsilon^2 of tables/<level>/T_anova.tsv for the figure's stat. score = effect_size * (1 - pooled_share):
a figure that mostly shows the pooled 'Other' group scores lower. Scored rows sort first, by score descending;
unscored rows follow.

signal / signal_kind give a second, coarser ranking for figures that have no ANOVA effect size:
    epsilon2      the T_anova effect_size itself (figures that have one)
    spearman_r    |Spearman r| of an A or R scatter figure
    permanova_R2  PERMANOVA R^2 by phylum for an ordination matrix
Signals of different kinds are NOT comparable (a correlation, a rank-based effect size and a variance fraction
measure different things): use them to order figures within one kind. Sort order: score descending (scored rows
first), then signal descending, then level, then name.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

COLUMNS = ["name", "level", "topic", "kind", "scope", "n_genomes", "n_groups", "pooled_share", "effect_size", "score", "signal", "signal_kind"]


def _sidecar(png):
    try:
        d = json.loads(png.with_suffix(".json").read_text())
        return d if isinstance(d, dict) else {}
    except (OSError, ValueError):
        return {}


def _effects(root, level):
    """stat -> finite effect size from tables/<level>/T_anova.tsv ({} when absent or unreadable)."""
    f = Path(root) / "tables" / level / "T_anova.tsv"
    try:
        t = pd.read_csv(f, sep="\t")
        s = pd.to_numeric(t["effect_size"], errors="coerce")
        return {k: float(v) for k, v in zip(t["stat"], s) if np.isfinite(v)}
    except (OSError, KeyError, ValueError, pd.errors.EmptyDataError):
        return {}


def _num(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return np.nan
    return v if np.isfinite(v) else np.nan


def _signal(side, effect):
    """(signal, signal_kind): the effect size when there is one, else the sidecar's own signal (empty if none)."""
    if np.isfinite(effect):
        return {"signal": effect, "signal_kind": "epsilon2"}
    sig = _num(side.get("signal"))
    if np.isfinite(sig) and side.get("signal_kind"):
        return {"signal": sig, "signal_kind": str(side["signal_kind"])}
    return {"signal": np.nan, "signal_kind": None}


def build_index(out_root):
    root = Path(out_root)
    rows = []
    for png in sorted((root / "figures").glob("*/*.png")):
        level = png.parent.name
        side = _sidecar(png)
        eff = _effects(root, level).get(side.get("stat")) if side.get("stat") else None
        ps = _num(side.get("pooled_share"))
        if np.isfinite(ps):
            ps = min(max(ps, 0.0), 1.0)
        e = np.nan if eff is None else eff
        rows.append({"name": png.stem, "level": level, "topic": side.get("topic"), "kind": side.get("kind") or "figure",
                     "scope": side.get("scope"), "n_genomes": _num(side.get("n_genomes")),
                     "n_groups": _num(side.get("n_groups")), "pooled_share": ps, "effect_size": e,
                     "score": e * (1 - (ps if np.isfinite(ps) else 0.0)) if np.isfinite(e) else np.nan,
                     **_signal(side, e)})
    df = pd.DataFrame(rows, columns=COLUMNS)
    df = df.sort_values(["score", "signal", "level", "name"], ascending=[False, False, True, True], na_position="last", kind="stable").reset_index(drop=True)
    for c in ("n_genomes", "n_groups"):
        df[c] = df[c].astype("Int64")
    df.to_csv(root / "figures_index.tsv", sep="\t", index=False)
    return df


if __name__ == "__main__":
    print(build_index(Path(sys.argv[1])).to_string())
