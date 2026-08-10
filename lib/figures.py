"""Palette/style constants shared by every BFD_dashboard HTML view and any later
PDF export, so the whole output reads as one visual system.

Self-contained by design (no import of another project's lib/): the palette is the
CVD-safe Okabe-Ito-derived set used across this lab's comparative-genomics reports.
An optional env var BFD_FIGURES_MODULE can point at an external figures.py to reuse
a canonical copy, but the bundled constants are the default so this package works
standalone.
"""
import os
import sys

# Canonical palette (matches the lab's shared report figures.py).
BLUE = "#2a78d6"
GREEN = "#008300"
ORANGE = "#e69f00"
PURPLE = "#cc79a7"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
INK = "#0b0b0b"
MUTED = "#898781"
ABSENT = "#ecebe4"
CATEGORICAL = [BLUE, ORANGE, PURPLE, GREEN]

# Extra series colors needed for up to ~7 GENUS/group categories (still CVD-safe,
# Okabe-Ito derived), appended after the 4 already validated above.
CATEGORICAL_EXTENDED = CATEGORICAL + ["#56b4e9", "#009e73", "#d55e00"]

# Optional external override: point BFD_FIGURES_MODULE at a directory containing a
# figures.py with the same names, e.g. to reuse NovInvenio's canonical constants.
_external = os.environ.get("BFD_FIGURES_MODULE")
if _external:
    try:
        if _external not in sys.path:
            sys.path.insert(0, _external)
        import figures as _ext  # noqa: E402

        for _name in ("BLUE", "GREEN", "ORANGE", "PURPLE", "GRID", "AXIS", "INK", "MUTED", "ABSENT", "CATEGORICAL"):
            globals()[_name] = getattr(_ext, _name)
        CATEGORICAL_EXTENDED = CATEGORICAL + ["#56b4e9", "#009e73", "#d55e00"]
    except Exception:
        pass  # bundled copy is fine

PALETTE_JSON = {
    "blue": BLUE,
    "green": GREEN,
    "orange": ORANGE,
    "purple": PURPLE,
    "grid": GRID,
    "axis": AXIS,
    "ink": INK,
    "muted": MUTED,
    "absent": ABSENT,
    "categorical": CATEGORICAL_EXTENDED,
}
