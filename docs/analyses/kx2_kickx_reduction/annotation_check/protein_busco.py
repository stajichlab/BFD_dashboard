import os, csv, glob, re
B = os.path.realpath(os.path.expanduser('~/projects/BFD_dashboard/bench/kickx_tree'))
G = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/results/genome_stats/BUSCO_protein'
locs = []
for t in ('backbone', 'coemansia'):
    for r in csv.DictReader(open(f'{B}/taxa/{t}_taxa.tsv'), delimiter='\t'):
        if r['LOCUSTAG'] not in locs:
            locs.append(r['LOCUSTAG'])
pat = re.compile(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)')
print('locustag\tp_C\tp_S\tp_D\tp_F\tp_M')
missing = []
for loc in locs:
    f = glob.glob(f'{G}/*/{loc}.BUSCO_summary.fungi_odb12.txt')
    if not f:
        missing.append(loc)
        continue
    m = pat.search(open(f[0]).read())
    print('\t'.join([loc] + list(m.groups()[:5])) if m else f'{loc}\tPARSE_FAIL')
import sys
print('# missing protein BUSCO summary:', missing, file=sys.stderr)
print('# loci:', len(locs), file=sys.stderr)
