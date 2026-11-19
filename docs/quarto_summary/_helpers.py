"""Helpers shared by the .qmd pages. Figures and tables are written by bin/summary/make_figures.py into
figures/<level>/ and tables/<level>/ for level in phylum, subphylum, class; each is shown as three tabs."""
import html
import json
from pathlib import Path
from urllib.parse import quote

import pandas as pd
try:
    from IPython.display import HTML, Markdown, display
except ImportError:  # unit tests run without IPython; the helpers then only return strings
    def Markdown(x): return x
    def HTML(x): return x
    def display(x): print(x)

HERE = Path(__file__).resolve().parent
PFAM = "https://www.ebi.ac.uk/interpro/entry/pfam/{}/"
LEVELS = [("phylum", "Phylum"), ("subphylum", "Subphylum"), ("class", "Class"), ("order", "Order")]


def _tabs(blocks):
    out = ["::: {.panel-tabset}"]
    for title, body in blocks:
        out += [f"## {title}", "", body, ""]
    out.append(":::")
    display(Markdown("\n".join(out)))


def _count(v):
    """Integer text for a sidecar value, or None when it is missing, NaN or not numeric."""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f or f in (float("inf"), float("-inf")) else f"{int(f):,}"


def _facts(folder, name):
    """One factual sentence from the sidecar figures/<folder>/<name>.json (n, groups, pooled share), or ''."""
    path = HERE / "figures" / folder / f"{name}.json"
    try:
        side = json.loads(path.read_text())
    except (OSError, ValueError):
        return ""
    parts = []
    n, g = _count(side.get("n_genomes")), _count(side.get("n_groups"))
    if n:
        parts.append(f"n = {n} genomes")
    if g:
        parts.append(f"{g} groups")
    ps = side.get("pooled_share")
    if isinstance(ps, (int, float)) and ps == ps:
        parts.append(f"{ps * 100:.0f}% of genomes are in the pooled \"Other\" group")
    return " (" + "; ".join(parts) + ")" if parts else ""


def fig(name, caption, width="95%", levels=None):
    """Embed figures/<level>/<name>.png as one tab per level; the SVG sits next to each PNG.
    The sidecar JSON, when present, adds n, number of groups and the pooled share to the caption."""
    blocks = []
    for lv, title in LEVELS:
        if levels and lv not in levels:
            continue
        png = HERE / "figures" / lv / f"{name}.png"
        if png.exists():
            cap = caption + _facts(lv, name)
            blocks.append((title, f"![{cap}](figures/{lv}/{name}.png){{width={width}}}\n\n*{cap}* "
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


def single(name, caption, folder, width="85%"):
    """Embed one figure that has no level tabs: figures/<folder>/<name>.png (folder = all or ordination)."""
    png = HERE / "figures" / folder / f"{name}.png"
    if not png.exists():
        display(Markdown(f"*Figure `{name}` has not been built (figures/{folder}/).*"))
        return
    cap = caption + _facts(folder, name)
    display(Markdown(f"![{cap}](figures/{folder}/{name}.png){{width={width}}}\n\n*{cap}* "
                     f"([SVG](figures/{folder}/{name}.svg))"))


def _read_commented(path):
    """(DataFrame, comment lines): leading lines that start with '#' are returned as text (without '#') and skipped."""
    comments = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            comments.append(line[1:].strip())
    return pd.read_csv(path, sep="\t", skiprows=len(comments)), comments


def read_comments(name, folder):
    """Leading '#' lines of tables/<folder>/<name>.tsv (thresholds written by the workflows), or [] when absent."""
    path = HERE / "tables" / folder / f"{name}.tsv"
    return _read_commented(path)[1] if path.exists() else []


def flat_table(name, folder, digits=3, max_rows=None, columns=None, where=None):
    """Show tables/<folder>/<name>.tsv (folder = ordination or tf); says so when the table is absent.
    Leading '#' lines are shown as a note above the table. where = optional function DataFrame -> DataFrame (row filter);
    an empty result prints a short note. The table sits in a scrollable div."""
    path = HERE / "tables" / folder / f"{name}.tsv"
    if not path.exists():
        display(Markdown(f"*Table `{name}` has not been built (tables/{folder}/).*"))
        return
    df, comments = _read_commented(path)
    if where is not None:
        df = where(df)
    if not len(df):
        display(Markdown(f"*Table `{name}` has no rows that match this view.*"))
        return
    if columns:
        df = df[[c for c in columns if c in df.columns]]
    if max_rows:
        df = df.head(max_rows)
    note = "".join(f"*{html.escape(c)}*\n\n" for c in comments)
    display(Markdown(note + '<div style="overflow-x:auto">\n\n' + df.round(digits).to_markdown(index=False) + "\n\n</div>"))


def read_flat(name, folder):
    """DataFrame of tables/<folder>/<name>.tsv (leading '#' lines skipped), or None when absent."""
    path = HERE / "tables" / folder / f"{name}.tsv"
    return _read_commented(path)[0] if path.exists() else None


def read_family_map(path=None):
    """data/tf_family_map.tsv (a copy of config/tf_family_map.tsv; lines that start with # are comments) as a DataFrame of
    strings, or None when the file cannot be read."""
    try:
        return pd.read_csv(Path(path) if path else HERE / "data" / "tf_family_map.tsv", sep="\t", comment="#", dtype=str,
                           keep_default_na=False)
    except (OSError, ValueError):
        return None


def say(text):
    """Show one Markdown paragraph."""
    display(Markdown(text))


INDEX_COLS = ["name", "level", "topic", "kind", "scope", "n_genomes", "n_groups", "pooled_share", "effect_size",
              "score", "signal", "signal_kind"]
_NUM = {"n_genomes": "{:,.0f}", "n_groups": "{:,.0f}", "pooled_share": "{:.2f}", "effect_size": "{:.3f}",
        "score": "{:.3f}", "signal": "{:.3f}"}
_FOLDER = {"all_genomes": "all", "genus_reps": "ordination"}

_JS = """<script>
(function(){
  var root=document.currentScript.parentNode, rows=root.querySelectorAll('tbody tr'), ctl=root.querySelectorAll('.gal-filter');
  function apply(){
    var shown=0;
    rows.forEach(function(r){
      var ok=true;
      ctl.forEach(function(c){
        var v=c.value, d=r.dataset[c.dataset.col]||'';
        if(!v){return;}
        if(c.dataset.mode==='exact'){ if(d!==v){ok=false;} }
        else if(d.toLowerCase().indexOf(v.toLowerCase())<0){ok=false;}
      });
      r.style.display=ok?'':'none'; if(ok){shown++;}
    });
    root.querySelector('.gal-count').textContent=shown+' of '+rows.length+' figures shown';
  }
  ctl.forEach(function(c){c.addEventListener('input',apply); c.addEventListener('change',apply);}); apply();
})();
</script>"""
_SELECTS = ["topic", "level", "scope", "signal_kind"]
_CSS = ("<style>.gallery table{font-size:.8rem}.gallery img{width:90px;height:auto}"
        ".gallery .gal-filter{margin:0 .5rem .5rem 0;max-width:100%}"
        ".gallery .gal-scroll{overflow-x:auto;max-width:100%}"
        "@media (max-width:600px){.gallery .lowval{display:none}}</style>")
_LOWVAL = {"kind", "n_groups", "pooled_share", "scope"}


def _fmt(col, v):
    if v is None or (isinstance(v, float) and v != v) or v == "":
        return ""
    if col in _NUM:
        try:
            return _NUM[col].format(float(v))
        except (TypeError, ValueError):
            pass
    return str(v)


def _fig_dir(row):
    scope = str(row.get("scope", ""))
    return _FOLDER.get(scope, str(row.get("level", "")))


def _url(row, ext):
    """Relative figure URL with every path segment percent-encoded (then HTML-escaped by the caller)."""
    return f"figures/{quote(_fig_dir(row), safe='')}/{quote(str(row['name']), safe='')}.{ext}"


def gallery_html(topic=None, level=None, top=None, index_path=None):
    """HTML for the figure gallery built from figures_index.tsv (index order kept; top = first N rows after the
    topic/level filter). Returns a short note instead when the index has not been built. Closed-set columns
    get exact-match <select> filters; the figure name gets a substring search box."""
    path = Path(index_path) if index_path else HERE / "figures_index.tsv"
    if not path.exists():
        return ("<p><em>figures_index.tsv has not been built yet, so the gallery is empty. "
                "It is written by bin/summary/make_figures.py.</em></p>")
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    if topic:
        df = df[df["topic"] == topic]
    if level:
        df = df[df["level"] == level]
    if top:
        df = df.head(int(top))
    e = html.escape
    cols = [c for c in INDEX_COLS if c in df.columns]
    head = ["thumbnail"] + cols + ["files"]
    sel = [c for c in _SELECTS if c in df.columns]
    out = ['<div class="gallery">' + _CSS]
    ctl = []
    for c in sel:
        vals = sorted({str(v) for v in df[c] if str(v) != ""})
        ctl.append(f'<label>{c} <select class="gal-filter" data-col="{c}" data-mode="exact">'
                   '<option value="">all</option>' + "".join(f'<option value="{e(v)}">{e(v)}</option>' for v in vals)
                   + "</select></label>")
    ctl.append('<label>name <input class="gal-filter" data-col="name" data-mode="contains" type="text" '
               'placeholder="contains..."></label>')
    out.append("<p>" + " ".join(ctl) + '<span class="gal-count" aria-live="polite"></span></p>')
    out.append('<div class="gal-scroll" style="overflow-x:auto">'
               '<table aria-label="Figure gallery"><thead><tr>'
               + "".join(f'<th class="{"lowval" if h in _LOWVAL else ""}">{e(h)}</th>' for h in head) + "</tr></thead><tbody>")
    for _, r in df.iterrows():
        png, svg = e(_url(r, "png"), quote=True), e(_url(r, "svg"), quote=True)
        nm = e(str(r["name"]), quote=True)
        data = " ".join(f'data-{c}="{e(str(r[c]), quote=True)}"' for c in sel + ["name"])
        alt = e(f"{r['name']} ({r.get('topic', '')}, {r.get('level', '')})", quote=True)
        cells = [f'<td><a href="{png}" tabindex="-1"><img loading="lazy" decoding="async" width="90" height="53" '
                 f'src="{png}" alt="{alt}"></a></td>']
        cells += [f'<td class="{"lowval" if c in _LOWVAL else ""}">{e(_fmt(c, r[c]))}</td>' for c in cols]
        cells.append(f'<td><a href="{png}" aria-label="PNG of {nm}">PNG</a> <a href="{svg}" aria-label="SVG of {nm}">SVG</a></td>')
        out.append(f"<tr {data}>" + "".join(cells) + "</tr>")
    out.append("</tbody></table></div>" + _JS + "</div>")
    return "\n".join(out)


def gallery(topic=None, level=None, top=None):
    """Display the filterable figure gallery (see gallery_html)."""
    display(HTML(gallery_html(topic, level, top)))
