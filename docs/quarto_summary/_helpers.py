"""Small helpers shared by the .qmd pages (tables are TSV files written by bin/summary/make_figures.py)."""
from pathlib import Path

import pandas as pd
from IPython.display import Markdown, display

HERE = Path(__file__).resolve().parent
PFAM = "https://www.ebi.ac.uk/interpro/entry/pfam/{}/"


def table(name, link_pfam=None, digits=3, max_rows=None):
    """Show tables/<name>.tsv as a Markdown table; link_pfam = column holding Pfam accessions."""
    path = HERE / "tables" / f"{name}.tsv"
    if not path.exists():
        display(Markdown(f"*Table `{name}` has not been built yet (run bin/summary/make_figures.py).*"))
        return
    df = pd.read_csv(path, sep="\t")
    if max_rows:
        df = df.head(max_rows)
    if link_pfam and link_pfam in df.columns:
        df[link_pfam] = df[link_pfam].map(lambda a: f"[{a}]({PFAM.format(str(a).split('.')[0])})")
    display(Markdown(df.round(digits).to_markdown(index=False)))


def fig(name, caption, width="95%"):
    """Embed figures/<name>.png (the SVG sits next to it) with a caption; a note if it is missing."""
    png = HERE / "figures" / f"{name}.png"
    if not png.exists():
        display(Markdown(f"*Figure `{name}` has not been built yet.*"))
        return
    display(Markdown(f"![{caption}](figures/{name}.png){{width={width}}}\n\n*{caption}* ([SVG](figures/{name}.svg))"))
