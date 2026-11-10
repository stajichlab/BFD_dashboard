"""Evaluate placement of Antarctolithica in the main and robustness trees."""
import glob, pandas as pd
from Bio import Phylo
from mdtab import md
T = pd.read_csv('phylo/taxa.tsv', sep='\t').set_index('label')
taph = {'Taphrinomycetes','Schizosaccharomycetes','Pneumocystomycetes','Saitoellomycetes','Neolectomycetes','Novakomycetes'}; sacc = {'Saccharomycetes','Pichiomycetes','Dipodascomycetes','Lipomycetes','Trigonopsidomycetes','Alloascoideomycetes','Sporopachydermiomycetes'}
def grp(n):
    r = T.loc[n]
    if 'Antarctolithica' in n: return 'AO'
    if r.PHYLUM == 'Mucoromycota': return 'Mucoro'
    if r.PHYLUM != 'Ascomycota': return 'Basidio'
    if pd.isna(r.CLASS): return 'unplaced'
    return 'Taphrino' if r.CLASS in taph else 'Saccharo' if r.CLASS in sacc else 'Pezizo'
rows = []
for f, desc in [('phylo/concat.fasttree.nwk', 'main: 300 markers, 159 taxa, LG+G')] + [(f'phylo/jk/jk{i}.nwk', f'jackknife {i+1}: 100 random markers, LG+G') for i in range(4)] + [('phylo/jk/wag.nwk', '100 random markers, WAG+G'), ('phylo/jk/reduced.nwk', '150 random markers, long-branch taxa removed, LG+G')]:
    try: tr = Phylo.read(f, 'newick')
    except Exception as e: rows.append(dict(tree=desc, note='not finished')); continue
    og = [t for t in tr.get_terminals() if grp(t.name) == 'Mucoro'] or [t for t in tr.get_terminals() if grp(t.name) == 'Basidio' and 'Saccharomycotina_sp' not in t.name]
    tr.root_with_outgroup(tr.common_ancestor(og)); ao = [t for t in tr.get_terminals() if grp(t.name) == 'AO']; m = tr.common_ancestor(ao); path = tr.get_path(m); par = path[-2] if len(path) > 1 else tr.root
    sis = [t.name for t in par.get_terminals() if grp(t.name) != 'AO']; tot = pd.Series([grp(t.name) for t in tr.get_terminals()]).value_counts(); sc = pd.Series([grp(n) for n in sis]).value_counts()
    full = {g: f"{sc.get(g, 0)}/{tot.get(g, 0)}" for g in ['Taphrino', 'Saccharo', 'Pezizo']}
    if sc.get('Saccharo', 0) == tot.get('Saccharo', 0) and sc.get('Pezizo', 0) == tot.get('Pezizo', 0) and sc.get('Taphrino', 0) == 0: place = 'sister to Saccharomycotina + Pezizomycotina'
    elif sc.get('Taphrino', 0) > 0 and sc.get('Pezizo', 0) == 0 and sc.get('Saccharo', 0) == 0: place = 'sister to / within Taphrinomycotina'
    elif sc.get('Pezizo', 0) == tot.get('Pezizo', 0) and sc.get('Saccharo', 0) == 0 and sc.get('Taphrino', 0) == 0: place = 'sister to Pezizomycotina'
    elif sc.get('Saccharo', 0) == tot.get('Saccharo', 0) and sc.get('Pezizo', 0) == 0 and sc.get('Taphrino', 0) == 0: place = 'sister to Saccharomycotina'
    else: place = 'other'
    tm = tr.common_ancestor([t for t in tr.get_terminals() if grp(t.name) == 'Taphrino']); taph_mono = all(grp(t.name) in ('Taphrino',) for t in tm.get_terminals())
    rows.append(dict(tree=desc, AO_monophyletic=len(m.get_terminals()) == len(ao), placement=place, support_of_AO_plus_sister=par.confidence, sister_Taphrino=full['Taphrino'], sister_Saccharo=full['Saccharo'], sister_Pezizo=full['Pezizo'], Taphrinomycotina_monophyletic=taph_mono))
R = pd.DataFrame(rows); R.to_csv('results/phylo_robustness.csv', index=False); print(R.to_string()); open('report/jk.md', 'w').write(md(R, 3))
