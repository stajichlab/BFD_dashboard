"""Helpers shared by the .qmd pages. Figures and tables are written by bin/summary/make_figures.py into
figures/<level>/ and tables/<level>/ for level in phylum, subphylum, class; each is shown as three tabs."""
from pathlib import Path

import pandas as pd
from IPython.display import Markdown, display

HERE = Path(__file__).resolve().parent
PFAM = "https://www.ebi.ac.uk/interpro/entry/pfam/{}/"
LEVELS = [("phylum", "Phylum"), ("subphylum", "Subphylum"), ("class", "Class")]


def _tabs(blocks):
    out = ["::: {.panel-tabset}"]
    for title, body in blocks:
        out += [f"## {title}", "", body, ""]
    out.append(":::")
    display(Markdown("\n".join(out)))


def fig(name, caption, width="95%", levels=None):
    """Embed figures/<level>/<name>.png as one tab per level; the SVG sits next to each PNG."""
    blocks = []
    for lv, title in LEVELS:
        if levels and lv not in levels:
            continue
        png = HERE / "figures" / lv / f"{name}.png"
        if png.exists():
            blocks.append((title, f"![{caption}](figures/{lv}/{name}.png){{width={width}}}\n\n*{caption}* "
                                  f"([SVG](figures/{lv}/{name}.svg))"))
        else:
            blocks.append((title, f"*Figure `{name}` has not been built for this level.*"))
    _tabs(blocks)


def table(name, link_pfam=None, digits=3, max_rows=None, levels=None):
    """Show tables/<level>/<name>.tsv as one tab per level; link_pfam = column holding Pfam accessions."""
    blocks = []
    for lv, title in LEVELS:
        if levels and lv not in levels:
            continue
        path = HERE / "tables" / lv / f"{name}.tsv"
        if not path.exists():
            blocks.append((title, f"*Table `{name}` has not been built for this level.*"))
            continue
        df = pd.read_csv(path, sep="\t")
        if max_rows:
            df = df.head(max_rows)
        if link_pfam and link_pfam in df.columns:
            df[link_pfam] = df[link_pfam].map(lambda a: f"[{a}]({PFAM.format(str(a).split('.')[0])})")
        blocks.append((title, df.round(digits).to_markdown(index=False)))
    _tabs(blocks)
