"""SwissProt best-hit inputs for 08_swissprot_contam.py, read from the raw BFD table swissprot_best (run after 01_build_set.py).

Writes
  data/feat/sprot_org_long.parquet  per genome x SwissProt organism code (comparison set only): n (best hits),
                                    pid (mean % identity) and pid_cov (mean % identity of the hits whose alignment covers
                                    >= 80 % of both query and subject; reproduces the 2026-11-10 laptop values exactly)
  data/ao_swissprot_best.parquet    all best hits of the Antarctolithica MAGs (per-contig contamination screen)
  data/orgcode_map0.csv             org_code -> organism, genus (first word of the organism name), fungal
                                    (genus is a BFD genus in genome_labels, or listed in fungal_genus_extra.txt)
Environment: BFD_DATA / AO_TABLES as in 00_prep_data.py; NCPU threads for DuckDB. Run from AO_BASE."""
import os, duckdb, pandas as pd
from pathlib import Path

COV_MIN_X10 = 800   # 80 % query and subject coverage for pid_cov


def tables_dir(env=os.environ):
    bfd = Path(env.get('BFD_DATA', '/bigdata/stajichlab/shared/projects/BFD'))
    return Path(env.get('AO_TABLES', bfd / 'Fungi_BFD_runs' / 'tables'))


def org_long_sql(src, cov=COV_MIN_X10):
    """Per genome x org_code aggregation of a swissprot_best-like relation `src` (needs species_prefix, org_code, pident_x10,
    qcov_x10, scov_x10). Columns are cast first: the x10 columns are USMALLINT and overflow in products."""
    return f"""SELECT species_prefix AS g, org_code, count(*) AS n, avg(p) / 10 AS pid,
                      avg(CASE WHEN qc >= {cov} AND sc >= {cov} THEN p END) / 10 AS pid_cov
               FROM (SELECT species_prefix, org_code, pident_x10::DOUBLE p, qcov_x10::DOUBLE qc, scov_x10::DOUBLE sc FROM {src})
               GROUP BY ALL"""


def orgcode_map(subject, bfd_genera, extra):
    m = subject[['org_code', 'organism']].drop_duplicates('org_code').copy()
    m['genus'] = m.organism.str.split().str[0]
    m['fungal'] = m.genus.isin(set(bfd_genera) | set(extra))
    return m.sort_values('org_code').reset_index(drop=True)


if __name__ == '__main__':
    T = tables_dir(); d = duckdb.connect(); d.sql(f"SET threads={int(os.environ.get('NCPU', '4'))}")
    U = pd.read_parquet('data/ao_set.parquet'); d.register('U', U[['LOCUSTAG']])
    sb = f"(SELECT s.* FROM read_parquet('{T}/swissprot_best.parquet') s JOIN U ON s.species_prefix = U.LOCUSTAG)"
    os.makedirs('data/feat', exist_ok=True)
    L = d.sql(org_long_sql(sb)).df(); L.to_parquet('data/feat/sprot_org_long.parquet', index=False)
    ao = tuple(U.LOCUSTAG[U.role == 'AO'])
    A = d.sql(f"SELECT species_prefix, protein_id, org_code, pident_x10, qcov_x10, scov_x10, bitscore FROM read_parquet('{T}/swissprot_best.parquet') WHERE species_prefix IN {ao}").df()
    A.to_parquet('data/ao_swissprot_best.parquet', index=False)
    S = d.sql(f"SELECT DISTINCT org_code, organism FROM read_parquet('{T}/swissprot_subject.parquet')").df()
    labels = pd.read_parquet('data/labels/genome_labels.parquet')
    extra = [l.strip() for l in open(Path(__file__).parent / 'fungal_genus_extra.txt') if l.strip() and not l.startswith('#')]
    M = orgcode_map(S, labels.GENUS.dropna(), extra); M.to_csv('data/orgcode_map0.csv', index=False)
    print(f'sprot_org_long: {len(L):,} rows, {L.g.nunique()} genomes; AO best hits: {len(A):,} ({len(ao)} MAGs); org codes: {len(M):,} ({int(M.fungal.sum())} fungal)')
