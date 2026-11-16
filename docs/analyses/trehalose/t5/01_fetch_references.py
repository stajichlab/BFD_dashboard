#!/usr/bin/env python3
"""T5 step 1: fetch GH65 reference and context proteins from UniProt and verify them.

Usage (HPCC login node, light; only public queries are sent):  01_fetch_references.py <BENCH>
Writes $BENCH/refs/refs.tsv (one row per protein: accession, group, role, EC, organism, taxid, domain, kingdom, reviewed, length, name),
$BENCH/refs/refs.faa (headers  >GROUP|ACCESSION|Organism_name  ), and $BENCH/refs/fetch_log.txt (query strings and counts).

Selection (fixed before any tree was built; the random seed is part of the rule):
  references (the gate): every REVIEWED UniProtKB entry that has a CAZy GH65 cross-reference OR the Pfam PF03632 cross-reference
      (query (xref:cazy-GH65 OR xref:pfam-PF03632) AND reviewed:true), length 500-1500 aa, no fragments. Their clade label comes from the
      UniProt EC number and, for EC 3.2.1.28, from kingdom Fungi (group_of()).
  context taxa (not part of the gate, never used to define a clade): UNREVIEWED entries with the Pfam PF03632 cross-reference, non-fragment,
      length 500-1500 aa, not Fungi (fungal proteins come from BFD). Bacteria are sampled per phylum (up to PER_PHYLUM genera each, so that
      one well-sampled phylum does not dominate), Archaea (up to 20 genera) and non-fungal Eukaryota in three pools (Metazoa 15, plants 10, all other 25 genera, so animals do not dominate). One protein per genus
      (first in a shuffled order seeded with 20261113), and no genus that already has a reference. UniProt returns at most MAX_PAGES x 500
      entries per query in its own order, so the pool is not a random sample of all GH65 sequences. Their function is automatic annotation
      and is NOT used.
Every sequence is verified: the sequence returned in the TSV must have the length UniProt reports for the entry. A mismatch stops the script.
"""
import csv
import os
import random
import re
import sys
import urllib.parse
import urllib.request

SEED = 20261113
MAX_PAGES = 4
PER_PHYLUM = 14
ARCHAEA_CAP = 20
EUK_POOLS = {'euk_metazoa': (33208, 15), 'euk_plants': (33090, 10), 'euk_other': (None, 25)}   # taxid (None = all other non-fungal), cap
BACT_PHYLA = {'Firmicutes': 1239, 'Bacteroidetes': 976, 'Actinobacteria': 201174, 'Proteobacteria': 1224, 'Chloroflexi': 200795,
              'Spirochaetes': 203691, 'Planctomycetes': 203682, 'Verrucomicrobia': 74201, 'Acidobacteria': 57723, 'Thermotogae': 200918}
FIELDS = 'accession,protein_name,organism_name,organism_id,length,ec,reviewed,lineage,sequence,fragment'
BASE = 'https://rest.uniprot.org/uniprotkb/search'
EC_NAMES = {'3.2.1.216': 'KOJIBIOSE_HYDROLASE', '2.4.1.64': 'TREHALOSE_PHOSPHORYLASE', '2.4.1.216': 'TREHALOSE6P_PHOSPHORYLASE',
            '2.4.1.230': 'KOJIBIOSE_PHOSPHORYLASE', '2.4.1.8': 'MALTOSE_PHOSPHORYLASE', '2.4.1.279': 'NIGEROSE_PHOSPHORYLASE',
            '2.4.1.282': 'GLUCOSYLRHAMNOSE_PHOSPHORYLASE', '2.4.1.332': 'GLUCOSYLGLYCEROL_PHOSPHORYLASE', '3.2.1.28': 'TREHALASE_NONFUNGAL'}
LEN = 'length:[500 TO 1500] AND fragment:false'


def fetch(query, log, max_pages=None):
    """Rows of a UniProt TSV search, following the paging Link header (at most max_pages pages of 500 if given)."""
    url = BASE + '?' + urllib.parse.urlencode({'query': query, 'fields': FIELDS, 'format': 'tsv', 'size': 500})
    rows, page = [], 0
    while url and (max_pages is None or page < max_pages):
        req = urllib.request.Request(url, headers={'User-Agent': 'bfd-dashboard-t5'})
        with urllib.request.urlopen(req, timeout=120) as r:
            text = r.read().decode()
            link = r.headers.get('Link', '')
        lines = text.splitlines()
        header = lines[0].split('\t')
        rows += [dict(zip(header, ln.split('\t'))) for ln in lines[1:]]
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = m.group(1) if m else None
        page += 1
    log.append(f'{query}\t{len(rows)}')
    return rows


def lineage_of(row):
    """(domain, kingdom) from the UniProt lineage column ('Taxonomic lineage')."""
    lin = row.get('Taxonomic lineage') or row.get('Lineage') or ''
    parts = [re.sub(r'\s*\(.*\)$', '', p.strip()) for p in lin.split(',')]
    dom = next((d for d in ('Bacteria', 'Archaea', 'Eukaryota') if d in parts), 'Other')
    king = next((k for k in ('Fungi', 'Metazoa', 'Viridiplantae') if k in parts), '')
    return dom, king


def group_of(ec_field, kingdom):
    """Clade label of a reviewed reference from its UniProt EC number(s) and kingdom."""
    ecs = sorted(set(filter(None, ec_field.replace(' ', '').split(';'))))
    if '3.2.1.28' in ecs and kingdom == 'Fungi':
        return 'TREHALASE_FUNGI'
    if '3.2.1.107' in ecs:
        return 'ATHL1_METAZOA' if kingdom == 'Metazoa' else 'GH65_GLUCOSIDASE_3.2.1.107'
    for e in ecs:
        if e in EC_NAMES:
            return EC_NAMES[e]
    ecs = [e for e in ecs if not e.endswith('-')]
    return 'GH65_OTHER' if not ecs else 'GH65_OTHER_' + ecs[0]


def genus_of(org):
    m = re.match(r'([A-Z][a-z]+)', org)
    return m.group(1) if m else org


def collapse_near_duplicates(rows, min_id=0.98):
    """Drop a reviewed row that is the same species, same length and >= min_id identical (ungapped) to a row with a lower accession
    (for example two strains of one yeast), so that near-identical references cannot make the leave-one-out test trivial."""
    kept = []
    for r in sorted(rows, key=lambda x: x['Entry']):
        sp = ' '.join(re.sub(r'\s*\(.*', '', r['Organism']).split()[:2])
        dup = False
        for k in kept:
            if k['_sp'] == sp and len(k['Sequence']) == len(r['Sequence']):
                same = sum(a == b for a, b in zip(k['Sequence'], r['Sequence'])) / len(r['Sequence'])
                dup = dup or same >= min_id
        if not dup:
            r = dict(r)
            r['_sp'] = sp
            kept.append(r)
    return kept


def select(reviewed, pools, seed=SEED):
    """Pure selection logic (tested). reviewed: rows (all become references). pools: {name: (rows, cap)} of unreviewed rows; one per genus
    (seeded shuffle), at most cap per pool, never a genus that already has a reference or was taken in another pool, never Fungi."""
    reviewed = collapse_near_duplicates(reviewed)
    out = [(r, *lineage_of(r), 'reference') for r in reviewed]
    taken = {genus_of(r['Organism']) for r in reviewed}
    for name in sorted(pools):
        rows, cap = pools[name]
        pool = list(rows)
        random.Random(f'{seed}-{name}').shuffle(pool)
        n = 0
        for r in pool:
            dom, king = lineage_of(r)
            g = genus_of(r['Organism'])
            if king == 'Fungi' or g in taken or n >= cap:
                continue
            taken.add(g)
            n += 1
            out.append((r, dom, king, 'context'))
    return out


def main(bench):
    os.makedirs(f'{bench}/refs', exist_ok=True)
    log = []
    reviewed = fetch(f'(xref:cazy-GH65 OR xref:pfam-PF03632) AND reviewed:true AND {LEN}', log)
    base = f'xref:pfam-PF03632 AND reviewed:false AND {LEN}'
    pools = {}
    for name, tid in BACT_PHYLA.items():
        pools['bact_' + name] = (fetch(f'{base} AND taxonomy_id:{tid}', log, MAX_PAGES), PER_PHYLUM)
    pools['archaea'] = (fetch(f'{base} AND taxonomy_id:2157', log, MAX_PAGES), ARCHAEA_CAP)
    for name, (tid, cap) in EUK_POOLS.items():
        q = f'{base} AND taxonomy_id:{tid}' if tid else f'{base} AND taxonomy_id:2759 NOT taxonomy_id:4751 NOT taxonomy_id:33208 NOT taxonomy_id:33090'
        pools[name] = (fetch(q, log, MAX_PAGES), cap)
    chosen = select(reviewed, pools)
    rows_out = []
    with open(f'{bench}/refs/refs.faa', 'w') as fa:
        for r, dom, king, role in chosen:
            seq = r['Sequence']
            if len(seq) != int(r['Length']):
                raise SystemExit(f"LENGTH MISMATCH {r['Entry']}: sequence {len(seq)} vs UniProt length {r['Length']}")
            if role == 'reference':
                grp = group_of(r['EC number'], king)
            else:
                grp = 'CONTEXT_' + (king.upper() if king else dom.upper())
            org = re.sub(r'[^A-Za-z0-9]+', '_', r['Organism']).strip('_')[:40]
            fa.write(f">{grp}|{r['Entry']}|{org}\n")
            fa.writelines(seq[i:i + 80] + '\n' for i in range(0, len(seq), 80))
            rows_out.append([r['Entry'], grp, role, r['EC number'], r['Organism'], r.get('Organism (ID)', ''), dom, king,
                             r['Reviewed'], r['Length'], r['Protein names']])
    with open(f'{bench}/refs/refs.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['accession', 'group', 'role', 'EC', 'organism', 'taxid', 'domain', 'kingdom', 'reviewed', 'length', 'name'])
        w.writerows(rows_out)
    nref = sum(1 for x in rows_out if x[2] == 'reference')
    log.append(f'selected {len(rows_out)}: {nref} references, {len(rows_out) - nref} context; seed {SEED}')
    open(f'{bench}/refs/fetch_log.txt', 'w').write('\n'.join(log) + '\n')
    print('\n'.join(log))


if __name__ == '__main__':
    main(sys.argv[1])
