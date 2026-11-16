#!/usr/bin/env python3
"""Step 2: annotated candidate proteins and tree inputs for the four paralog-prone meiotic families (MutS, MutL, kleisin, HORMA).

Usage: 02_build_inputs.py <BENCH> <backbone_taxa.tsv> <coemansia_strains.txt> [tables_dir]
Run on the HPCC (python3.12 with duckdb). Reads the BFD Pfam table and the proteomes read-only.

Genomes: the 63 backbone genera (backbone_taxa.tsv, column LOCUSTAG) plus the Coemansia strains listed in coemansia_strains.txt
(one LOCUSTAG per line), the same 73-genome design as the RecA tree in ../../family_trees/.
Candidate rule (fixed before any result was seen): a protein is a candidate for a family if it has the family domain at domain_i_evalue <= 1e-5
and length >= MIN_LEN:
  muts    : PF00488 MutS_V         (the MutS-family ATPase domain shared by Msh1-6)           MIN_LEN 500
  mutl    : PF01119 DNA_mis_repair (the MutL N-terminal domain shared by Mlh1, Pms1, Mlh2, Mlh3) MIN_LEN 400
  kleisin : PF04824 Rad21_Rec8 or PF04825 Rad21_Rec8_N (Rad21 and Rec8)                        MIN_LEN 250
  horma   : PF02301 HORMA          (Hop1, Mad2, Rev7)                                          MIN_LEN 150
Outputs in $BENCH/trees/: <fam>_input.faa (labels REF_n / ANN_n, safe for tree software) and <fam>_members.tsv
(label, kind ref|ann, locustag, id, clade_or_best_ref, length, info). References come from $BENCH/refs/<fam>.faa (step 1).
"""
import csv, os, sys, gzip
import duckdb

bench, bb_tsv, strains_txt = sys.argv[1], sys.argv[2], sys.argv[3]
tables = sys.argv[4] if len(sys.argv) > 4 else '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables'
MAP = '/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/locustag_to_pep.tsv'   # read-only
FAM = {'muts': ['PF00488'], 'mutl': ['PF01119'], 'kleisin': ['PF04824', 'PF04825'], 'horma': ['PF02301']}
MIN_LEN = {'muts': 500, 'mutl': 400, 'kleisin': 250, 'horma': 150}
EVAL = 1e-5

loci = [r['LOCUSTAG'] for r in csv.DictReader(open(bb_tsv), delimiter='\t')]
loci += [l.strip() for l in open(strains_txt) if l.strip()]
loci = sorted(set(loci))
print('genomes:', len(loci), flush=True)
rng = ' OR '.join(f"(protein_id >= '{l}_' AND protein_id < '{l}`')" for l in loci)
accs = sorted({a for v in FAM.values() for a in v})
acc_in = ','.join(f"'{a}'" for a in accs)
con = duckdb.connect()
con.execute("SET threads=8; SET memory_limit='24GB'; SET preserve_insertion_order=false")
hits = con.execute(f"""
SELECT protein_id, split_part(protein_id,'_',1) AS loc, split_part(pfam_acc,'.',1) AS acc, pfam_id, domain_i_evalue
FROM '{tables}/pfam.parquet'
WHERE ({rng}) AND domain_i_evalue <= {EVAL} AND split_part(pfam_acc,'.',1) IN ({acc_in})""").fetchall()
print('domain hits:', len(hits), flush=True)

cand = {f: {} for f in FAM}                       # family -> protein_id -> (locustag, set of 'pfam_id:acc')
for pid, loc, acc, pfid, ev in hits:
    for f, a in FAM.items():
        if acc in a:
            cand[f].setdefault(pid, (loc, set()))[1].add(f'{pfid}:{acc}')

pep = {}
for line in open(MAP):
    p = line.rstrip('\n').split('\t')
    if p[0] in loci:
        pep[p[0]] = p[2]
need = {}                                          # locustag -> protein ids needed
for f in cand:
    for pid, (loc, _) in cand[f].items():
        need.setdefault(loc, set()).add(pid)
seqs = {}
for loc, ids in need.items():
    name, buf = None, []
    with gzip.open(pep[loc], 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name in ids:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if name in ids:
            seqs[name] = ''.join(buf)

os.makedirs(f'{bench}/trees', exist_ok=True)
for f in FAM:
    refs = []
    name, buf = None, []
    for line in open(f'{bench}/refs/{f}.faa'):
        if line.startswith('>'):
            if name:
                refs.append((name, ''.join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    refs.append((name, ''.join(buf)))
    n_ann, rows, out = 0, [], []
    for i, (hdr, s) in enumerate(refs, 1):
        clade, acc, gene, org = hdr.split('|')
        out.append(f'>REF_{i}\n{s}\n')
        rows.append((f'REF_{i}', 'ref', '', acc, clade, len(s), f'{gene} {org}'))
    for pid in sorted(cand[f]):
        s = seqs.get(pid, '')
        if len(s) < MIN_LEN[f]:
            continue
        n_ann += 1
        loc, doms = cand[f][pid]
        out.append(f'>ANN_{n_ann}\n{s}\n')
        rows.append((f'ANN_{n_ann}', 'ann', loc, pid, '', len(s), ';'.join(sorted(doms))))
    open(f'{bench}/trees/{f}_input.faa', 'w').write(''.join(out))
    with open(f'{bench}/trees/{f}_members.tsv', 'w') as fh:
        fh.write('label\tkind\tlocustag\tid\tclade_or_best_ref\tlength\tinfo\n')
        for r in rows:
            fh.write('\t'.join(str(x) for x in r) + '\n')
    print(f'{f}: {len(refs)} references + {n_ann} annotated candidates ({len(cand[f])} before the length filter)', flush=True)
