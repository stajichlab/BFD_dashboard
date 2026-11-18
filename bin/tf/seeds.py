"""Seeds for ortholog tracking (WC-2, Zap1-like, bZip5).

Config: config/tf_seeds.tsv, columns target, seed_name, source, key, taxon_id, target_pfam, required_pfam, note.
    source paper    key = JGI protein id of R. microsporus ATCC 11559 (Rhimi59_2); sequence from the JGI protein file
    source uniprot  key = gene symbol, taxon_id = organism; fetched from UniProt REST, cached by key
    source bfd      key = <species_prefix>|pfam:<DOMAINS joined by +>, e.g. FBF4BF29|pfam:GATA+PAS

Rule syntax in required_pfam: accessions joined by | are alternatives (any one satisfies); an entry
`PF00096>=4` needs at least 4 domains of that accession. See ortho_search.architecture_ok.
"""
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import jgi_map as jm  # noqa: E402

CONFIG_COLS = ["target", "seed_name", "source", "key", "taxon_id", "target_pfam", "required_pfam", "note"]
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{key}+AND+organism_id:{taxon}&format=fasta&size=1"
PAPER_PREFIX = "Rhimi59_2_"
# Short names of the domain groups usable in a bfd selector.
SELECTOR_GROUPS = {"GATA": ["PF00320"], "PAS": ["PF00989", "PF08447", "PF13426"], "C2H2": ["PF00096"], "BZIP": ["PF00170"]}
# Sanity rules for the named TFs of the paper: name -> (rule list, text for the warning).
NAMED_RULES = {
    "Znf1": (["PF00096"], "no C2H2 (PF00096) domain"),
    "GATA_Znf10": (["PF00320", "PF00989|PF08447|PF13426"], "needs a GATA domain (PF00320) and a PAS domain (PF00989, PF08447 or PF13426)"),
    "bZip5": (["PF00170"], "no bZIP_1 (PF00170) domain"),
}
TAXON_PHYLUM = {"559292": "Ascomycota", "367110": "Ascomycota"}


def read_seed_config(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment="#")
    missing = [c for c in CONFIG_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    return df[CONFIG_COLS]


def strip_version(acc):
    return re.sub(r"\.\d+$", "", str(acc))


def rule_accessions(target_pfam, required_pfam):
    """Every Pfam accession named in a rule row (target Pfam and required entries), without '>=N' counts."""
    out = [a for a in re.split(r"[|;,]", str(target_pfam)) if a]
    for entry in str(required_pfam).split(";"):
        for a in entry.split("|"):
            a = a.split(">=")[0].strip()
            if a:
                out.append(a)
    return out


def rules_for_target(cfg, target):
    """(target_pfam, [required entries]) of a target from the config (first row of the target)."""
    r = cfg[cfg["target"] == target].iloc[0]
    req = [x for x in str(r["required_pfam"]).split(";") if x] if r["required_pfam"] else []
    return r["target_pfam"], req


# ---------------------------------------------------------------- UniProt
def _default_fetch(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": "tf-survey-seed-fetch"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode()


def _parse_fasta_first(text):
    lines = text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        return None
    head = lines[0][1:]
    parts = head.split("|")
    acc = parts[1] if len(parts) >= 3 else head.split()[0]
    return acc, "".join(x.strip() for x in lines[1:] if not x.startswith(">"))


def fetch_seeds(seeds: pd.DataFrame, cache_dir: Path, fetch=None) -> pd.DataFrame:
    """UniProt seeds of `seeds` (rows with source == 'uniprot'): columns target, seed_name, accession, sequence, status.
    status is 'ok' or 'missing: <reason>'. No exception for a network or not-found error. A good answer is cached as
    <cache_dir>/<taxon>_<key>.fasta; a cached key is not fetched again. Failures are not cached."""
    fetch = fetch or _default_fetch
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, r in seeds[seeds["source"] == "uniprot"].iterrows():
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", f"{r['taxon_id']}_{r['key']}")
        cf = cache_dir / f"{safe}.fasta"
        status, acc, seq = "ok", "", ""
        try:
            if cf.exists() and cf.stat().st_size > 0:
                text = cf.read_text()
            else:
                url = UNIPROT_URL.format(key=urllib.parse.quote(r["key"]), taxon=r["taxon_id"])
                text = fetch(url)
                if _parse_fasta_first(text or "") is None:
                    raise LookupError("not found in UniProt")
                cf.write_text(text)
            acc, seq = _parse_fasta_first(text)
        except Exception as e:     # network, timeout, not found: never stop the run
            status = f"missing: {type(e).__name__ if not str(e) else str(e)}"
            print(f"  seed {r['seed_name']} ({r['key']}): {status}", flush=True)
        rows.append({"target": r["target"], "seed_name": r["seed_name"], "accession": acc, "sequence": seq, "status": status})
    return pd.DataFrame(rows, columns=["target", "seed_name", "accession", "sequence", "status"])


# ---------------------------------------------------------------- paper seeds
def paper_sequences(seeds: pd.DataFrame, jgi_fasta: Path) -> pd.DataFrame:
    """Full-length sequences of the paper seeds from the JGI protein file. Columns as fetch_seeds."""
    prot = jm.read_jgi_proteins(jgi_fasta)
    rows = []
    for _, r in seeds[seeds["source"] == "paper"].iterrows():
        hit = prot.get(str(r["key"]))
        rows.append({"target": r["target"], "seed_name": r["seed_name"], "accession": PAPER_PREFIX + str(r["key"]),
                     "sequence": hit[1] if hit else "",
                     "status": "ok" if hit else f"missing: JGI protein {r['key']} not in the JGI protein file"})
    return pd.DataFrame(rows, columns=["target", "seed_name", "accession", "sequence", "status"])


def check_named_architecture(name, arch):
    """'' when the BFD protein architecture `arch` satisfies the sanity rule of the named TF, else a warning text."""
    from ortho_search import architecture_ok
    if name not in NAMED_RULES:
        return ""
    rule, text = NAMED_RULES[name]
    return "" if architecture_ok(arch, rule) else f"architecture check failed: {text}; matched BFD protein has '{arch}'"


def resolve_paper_seeds(paper_seqs, bfd_proteome_fasta, names_tsv, arch, tmp, min_pident=95.0, min_cov=0.9, search=None):
    """Match the paper sequences (dict paper id -> sequence) to the proteome of FBF4BF29 with MMseqs2.
    Returns one row per paper id: paper_id, name, bfd_protein, pident, qcov, tcov, architecture, warning.
    `arch` maps a BFD protein id to its Pfam architecture string. `search(query_fa, target_fa, out_tsv, tmp, pident, cov)`
    returns the hit table (default: MMseqs2 through jgi_map.easy_search). The name comes from Supplementary Table 1
    (config/paper_tf_names.tsv). A named TF without the expected Pfam architecture in BFD gets a warning."""
    tmp = Path(tmp)
    tmp.mkdir(parents=True, exist_ok=True)
    names = jm.read_tf_names(names_tsv)
    name_of = dict(zip(names["jgi_id"].astype(str), names["name"]))
    q = tmp / "paper_seeds.faa"
    q.write_text("".join(f">{k}\n{v}\n" for k, v in paper_seqs.items()))
    search = search or (lambda a, b, o, t, p, c: jm.easy_search(a, b, o, t, p, c))
    hits = search(q, Path(bfd_proteome_fasta), tmp / "paper_vs_bfd.tsv", tmp / "mm", min_pident, min_cov)
    best = {}
    if len(hits):
        for k, g in hits.groupby("query"):
            best[k] = g.sort_values("bits", ascending=False).iloc[0]
    rows = []
    for pid in paper_seqs:
        jid = pid[len(PAPER_PREFIX):] if pid.startswith(PAPER_PREFIX) else pid
        name = name_of.get(jid, "")
        b = best.get(pid)
        row = {"paper_id": pid, "name": name, "bfd_protein": "", "pident": float("nan"), "qcov": float("nan"),
               "tcov": float("nan"), "architecture": "", "warning": ""}
        if b is not None:
            a = arch.get(b["target"], "")
            row.update(bfd_protein=b["target"], pident=float(b["pident"]), qcov=float(b["qcov"]), tcov=float(b["tcov"]),
                       architecture=a, warning=check_named_architecture(name, a))
        elif name:
            row["warning"] = f"{name}: no BFD match at identity {min_pident} and coverage {min_cov}"
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- BFD seeds
def bfd_selector_proteins(species_prefix, selector, pfam_hits: pd.DataFrame):
    """Proteins of one genome with every domain group of a selector 'pfam:GATA+PAS' (a group matches when the protein
    has any accession of it). pfam_hits has columns protein_id, pfam_acc (versions allowed)."""
    tokens = selector.split(":", 1)[1].split("+")
    groups = []
    for t in tokens:
        if t.upper() not in SELECTOR_GROUPS:
            raise ValueError(f"unknown domain group {t!r} in selector {selector!r}")
        groups.append(set(SELECTOR_GROUPS[t.upper()]))
    ph = pfam_hits[pfam_hits["protein_id"].astype(str).str.split("_").str[0] == species_prefix]
    have = ph.assign(acc=ph["pfam_acc"].map(strip_version)).groupby("protein_id")["acc"].agg(set)
    return sorted(p for p, accs in have.items() if all(accs & g for g in groups))


def read_fasta_dict(path):
    seqs, name = {}, None
    with jm._open_text(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}
