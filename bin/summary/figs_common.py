"""Shared helpers for the figs_* page modules: safe figure wrapper and per-genome join."""
import sys
import traceback

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def safe(o, name, topic, fn):
    """Run fn() -> Figure (or None to skip) and save it as `name`. A failure logs and continues; the name is
    recorded as skipped in o.made."""
    try:
        fig = fn()
        if fig is None:
            o.skip(name, topic)
            return False
        o.save(fig, name, topic=topic)
    except Exception as e:  # one broken figure (or a failed save) must not stop the page
        traceback.print_exc(file=sys.stderr)
        print(f"  skipped {name}: {type(e).__name__}: {e}", flush=True)
        plt.close("all")
        o.skip(name, topic)
        return False
    return True


def join_groups(core, table):
    """Per-genome table joined to the group column of core (inner join on species_prefix)."""
    return core[["species_prefix", "grp", "gc_pct"]].merge(table, on="species_prefix", how="inner")


def positive(d, col):
    v = d[col].dropna()
    return len(v) > 0 and bool((v > 0).all())
