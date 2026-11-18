"""Seeds for ortholog tracking (WC-2, Zap1-like, bZip5).

Config: config/tf_seeds.tsv, columns target, seed_name, source, key, taxon_id, target_pfam, required_pfam, note.
    source paper    key = JGI protein id of R. microsporus ATCC 11559 (Rhimi59_2); sequence from the JGI protein file
    source uniprot  key = gene symbol, taxon_id = organism; fetched from UniProt REST, cached by key
    source bfd      key = <species_prefix>|pfam:<DOMAINS joined by +>, e.g. FBF4BF29|pfam:GATA+PAS

Rule syntax in required_pfam (entries separated by ';', all must hold; parsed only in parse_rule):
    A|B          at least one domain of any listed accession
    A>=N, A==N   at least / exactly N domains of A (A may be A|B)
    PAS_regions==N, PAS_regions>=N   number of merged PAS regions (all PAS-family hits merged by coordinates)
See ortho_search.architecture_ok.
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
# PAS-family Pfams (Pfam has no LOV family): PAS, PAS_3, PAS_9, PAS_2, PAS_11, PAS_4, PAS_6, PAS_5, PAS_7, PAS_8.
PAS_ACCS = ["PF00989", "PF08447", "PF13426", "PF08446", "PF14598", "PF08448", "PF08348", "PF07310", "PF12860", "PF13188"]
SELECTOR_GROUPS = {"GATA": ["PF00320"], "PAS": PAS_ACCS, "C2H2": ["PF00096"], "BZIP": ["PF00170"]}
# Sanity rules for the named TFs of the paper: name -> (rule list, text for the warning).
NAMED_RULES = {
    "Znf1": (["PF00096"], "no C2H2 (PF00096) domain"),
    "GATA_Znf10": (["PF00320", "|".join(PAS_ACCS)], "needs a GATA domain (PF00320) and a PAS-family domain"),
    "bZip5": (["PF00170"], "no bZIP_1 (PF00170) domain"),
}
TAXON_PHYLUM = {"559292": "Ascomycota", "367110": "Ascomycota", "426428": "Ascomycota", "746128": "Ascomycota",
                "235443": "Basidiomycota"}


def read_seed_config(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment="#")
    missing = [c for c in CONFIG_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    return df[CONFIG_COLS]


def strip_version(acc):
    return re.sub(r"\.\d+$", "", str(acc))


def parse_rule(entry):
    """Parse one rule entry. Returns ('any', {acc}), ('count', {acc}, op, n) or ('pas', op, n); op is '==' or '>='.
    ValueError for anything else. This is the only place that knows the syntax."""
    entry = entry.strip()
    m = re.fullmatch(r"PAS_regions(==|>=)(\d+)", entry)
    if m:
        return ("pas", m.group(1), int(m.group(2)))
    m = re.fullmatch(r"([A-Za-z0-9_|]+)(==|>=)(\d+)", entry)
    if m:
        return ("count", set(m.group(1).split("|")), m.group(2), int(m.group(3)))
    if re.fullmatch(r"[A-Za-z0-9_]+(\|[A-Za-z0-9_]+)*", entry):
        return ("any", set(entry.split("|")))
    raise ValueError(f"bad architecture rule entry: {entry!r}")


def rule_accessions(target_pfam, required_pfam):
    """Every Pfam accession named in a rule row (target Pfam and required entries); PAS_regions entries name none."""
    out = [a for a in re.split(r"[|;,]", str(target_pfam)) if a]
    for entry in str(required_pfam).split(";"):
        if entry.strip():
            r = parse_rule(entry)
            if r[0] in ("any", "count"):
                out += sorted(r[1])
    return out


def pas_regions(domain_rows, pas_accs=None, min_overlap=10):
    """Merge the PAS-family domain hits of ONE protein into regions by alignment coordinates. Two hits join one region
    when they overlap by at least `min_overlap` residues (transitive). domain_rows has pfam_acc, ali_from, ali_to.
    Returns the sorted list of (start, end) regions."""
    pas = set(pas_accs or PAS_ACCS)
    d = domain_rows[domain_rows["pfam_acc"].map(strip_version).isin(pas)].sort_values(["ali_from", "ali_to"])
    regions = []
    for a, b in zip(d["ali_from"], d["ali_to"]):
        if regions and regions[-1][1] - a + 1 >= min_overlap:
            regions[-1] = (regions[-1][0], max(regions[-1][1], int(b)))
        else:
            regions.append((int(a), int(b)))
    return regions


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
    """Proteins of one genome that satisfy a selector 'pfam:GATA+PAS_regions==1'. Tokens joined by '+': a domain group
    name (GATA, PAS, C2H2, BZIP; the protein needs any accession of it) or a PAS_regions rule (merged PAS regions,
    see pas_regions). pfam_hits has protein_id, pfam_acc (versions allowed) and, for PAS_regions, ali_from, ali_to."""
    tokens = selector.split(":", 1)[1].split("+")
    groups, region_rules = [], []
    for t in tokens:
        if t.startswith("PAS_regions"):
            region_rules.append(parse_rule(t))
        elif t.upper() in SELECTOR_GROUPS:
            groups.append(set(SELECTOR_GROUPS[t.upper()]))
        else:
            raise ValueError(f"unknown domain group {t!r} in selector {selector!r}")
    ph = pfam_hits[pfam_hits["protein_id"].astype(str).str.split("_").str[0] == species_prefix]
    out = []
    for pid, g in ph.groupby("protein_id"):
        accs = set(g["pfam_acc"].map(strip_version))
        if not all(accs & grp for grp in groups):
            continue
        ok = True
        for _, op, n in region_rules:
            k = len(pas_regions(g))
            ok = ok and (k == n if op == "==" else k >= n)
        if ok:
            out.append(pid)
    return sorted(out)


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
