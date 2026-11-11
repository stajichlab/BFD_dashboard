"""Root the FastTree tree on Mucoromycota, report the placement of Antarctolithica, draw the tree."""
import sys, re, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from Bio import Phylo
from common import COL, FIG
nwk = sys.argv[1] if len(sys.argv) > 1 else 'phylo/concat.fasttree.nwk'; out = sys.argv[2] if len(sys.argv) > 2 else f'{FIG}/fable_C_ao_fig8_phylogeny.png'
T = pd.read_csv('phylo/taxa.tsv', sep='\t').set_index('label'); tr = Phylo.read(nwk, 'newick')
og = [t for t in tr.get_terminals() if T.PHYLUM.get(t.name) == 'Mucoromycota']
tr.root_with_outgroup(tr.common_ancestor(og) if len(og) > 1 else og[0]); tr.ladderize()
taph = {'Taphrinomycetes','Schizosaccharomycetes','Pneumocystomycetes','Saitoellomycetes','Neolectomycetes','Novakomycetes'}; sacc = {'Saccharomycetes','Pichiomycetes','Dipodascomycetes','Lipomycetes','Trigonopsidomycetes','Alloascoideomycetes','Sporopachydermiomycetes'}
def grp(n):
    r = T.loc[n]
    if 'Antarctolithica' in n: return 'Antarctolithica'
    if r.PHYLUM != 'Ascomycota': return 'Basidiomycota (outgroup)'
    if pd.isna(r.CLASS): return 'unplaced Ascomycota'
    if r.CLASS in taph: return 'Taphrinomycotina'
    if r.CLASS in sacc: return 'Saccharomycotina'
    return r.CLASS if r.CLASS in COL else 'other Pezizomycotina'
ao = [t for t in tr.get_terminals() if 'Antarctolithica' in t.name]; mrca = tr.common_ancestor(ao)
print('AO monophyletic:', len(mrca.get_terminals()) == len(ao), 'support', mrca.confidence)
path = tr.get_path(mrca); 
for depth, node in enumerate(reversed(path[:-1])):
    tips = [t.name for t in node.get_terminals() if 'Antarctolithica' not in t.name]
    cl = pd.Series([f"{T.CLASS.get(t)}" for t in tips]).value_counts().to_dict()
    print(f'ancestor -{depth+1}: support={node.confidence} n_other_tips={len(tips)} classes={cl}')
    if depth == 0: print('   sister tips:', [f"{T.SPECIES[t]} [{T.CLASS[t]}]" for t in tips][:12])
    if depth >= 5: break
# patristic distances from AO representative to all tips
rep = [t for t in ao if 'FF1FC3B3' in t.name][0]
D = pd.Series({t.name: tr.distance(rep, t) for t in tr.get_terminals() if 'Antarctolithica' not in t.name}).sort_values()
print('closest tips by patristic distance:'); print('\n'.join(f"  {d:.3f}  {T.SPECIES[n]} [{T.CLASS[n]}]" for n, d in D.head(12).items()))
root_tip = pd.Series({t.name: tr.distance(tr.root, t) for t in tr.get_terminals()}); gs = root_tip.groupby([grp(n) for n in root_tip.index]).median(); print('median root-to-tip by group:'); print(gs.round(3).to_string())
pd.DataFrame({'patristic_from_AO': D, 'CLASS': T.CLASS.reindex(D.index), 'SPECIES': T.SPECIES.reindex(D.index)}).to_csv('results/phylo_patristic.csv')
# draw
tips = tr.get_terminals(); y = {t: i for i, t in enumerate(tips)}; x = tr.depths()
def ypos(c):
    if c.is_terminal(): return y[c]
    v = np.mean([ypos(k) for k in c.clades]); y[c] = v; return v
ypos(tr.root); fig, a = plt.subplots(figsize=(13, 0.155 * len(tips) + 2))
def draw(c):
    for k in c.clades:
        a.plot([x[c], x[c]], [y[c], y[k]], color='#55554f', lw=.7); a.plot([x[c], x[k]], [y[k], y[k]], color='#55554f', lw=.7); draw(k)
    if not c.is_terminal() and c.confidence is not None and c.confidence < 0.9: a.text(x[c], y[c], f'{c.confidence:.2f}', fontsize=5, color='#e34948', ha='right', va='bottom')
draw(tr.root)
for t in tips:
    g = grp(t.name); r = T.loc[t.name]; lab = f"{r.SPECIES} {'' if pd.isna(r.STRAIN) or g != 'Antarctolithica' else r.STRAIN}  [{r.CLASS if not pd.isna(r.CLASS) else 'unplaced'}]"
    a.scatter([x[t]], [y[t]], s=22 if g != 'Antarctolithica' else 90, c=COL[g], marker='o' if g != 'Antarctolithica' else '*', zorder=3, edgecolors='white', linewidths=.4)
    a.text(x[t] + .01, y[t], lab, fontsize=6.3, va='center', fontweight='bold' if g == 'Antarctolithica' else 'normal', color='#111111')
a.set_ylim(len(tips), -1); a.set_yticks([]); a.spines[['left', 'top', 'right']].set_visible(False); a.set_xlabel('substitutions per site (LG+G, FastTree); red numbers = SH-like support < 0.90')
import matplotlib.lines as ml
a.legend(handles=[ml.Line2D([], [], marker='o' if g != 'Antarctolithica' else '*', ls='', color=c, label=g, markersize=7) for g, c in COL.items() if g in set(grp(t.name) for t in tips)], frameon=False, fontsize=8, loc='lower left')
a.set_xlim(0, max(x.values()) * 1.45); plt.tight_layout(); plt.savefig(out, dpi=130)
