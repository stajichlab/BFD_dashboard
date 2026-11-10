"""Random-forest classification of CLASS / subphylum from each feature block and all blocks.
CV: 5-fold StratifiedGroupKFold, groups = genus (species of one genus never split across folds).
Novelty calibration: leave-one-class-out (LOCO). Incompleteness: binomial thinning of count blocks."""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import balanced_accuracy_score, f1_score
from common import *
RF = dict(n_estimators=500, max_features='sqrt', min_samples_leaf=1, class_weight='balanced_subsample', n_jobs=12, random_state=0)
nprot = 10 ** feat('struct').log10_n_genes
B = blocks()
raw = {}
pf = feat('pfam'); pf = pf.loc[:, (pf > 0).mean() >= 0.01]; raw['pfam_ab'] = pf
for n in ['cazy', 'merops']:
    x = feat(n); raw[n] = x.loc[:, (x > 0).mean() >= 0.01]
raw['tf'] = feat('tf')
def tf_tr(x): return clr(x, 0.5)
def rel(x, npx): return np.sqrt(x.div(npx, axis=0) * 1000)
M = {k: v[0] for k, v in B.items() if k not in ('pfam_pa',)}
M['pfam_ab'] = rel(raw['pfam_ab'], nprot); M['cazy'] = rel(raw['cazy'], nprot); M['merops'] = rel(raw['merops'], nprot)
M['pfam_pa'] = (raw['pfam_ab'] > 0).astype(float)
allb = ['struct', 'cell', 'codon', 'aa', 'trna', 'tf', 'cazy', 'merops', 'pfam_ab']
M['all'] = pd.concat([M[b].add_prefix(b + ':') for b in allb], axis=1)
M['all_nonpfam'] = pd.concat([M[b].add_prefix(b + ':') for b in allb[:-1]], axis=1)
tr = U[U.role == 'train'].copy()
cn = tr.CLASS.value_counts(); keep = cn[cn >= 8].index
tr = tr[tr.CLASS.isin(keep)]
print('training species:', len(tr), 'classes:', len(keep)); print(cn[keep].to_string())
targets = {'class': tr.CLASS, 'subphylum': tr.subphylum}
perf, preds, imps, thin = [], [], [], []
rng = np.random.default_rng(0)
for tname, y in targets.items():
    for b in ['struct', 'struct_all', 'cell', 'codon', 'aa', 'trna', 'tf', 'cazy', 'merops', 'pfam_pa', 'pfam_ab', 'all_nonpfam', 'all']:
        X = M[b]; Xt = X.loc[tr.index].values; yv = y.values
        cv = StratifiedGroupKFold(5, shuffle=True, random_state=0)
        oof = np.empty(len(yv), dtype=object); oofp = np.zeros(len(yv)); pimp = np.zeros(X.shape[1]); thin_pred = {0.93: oof.copy(), 0.75: oof.copy()}
        for f, (a, te) in enumerate(cv.split(Xt, yv, tr.genus_key.values)):
            m = RandomForestClassifier(**RF).fit(Xt[a], yv[a]); P = m.predict_proba(Xt[te])
            oof[te] = m.classes_[P.argmax(1)]; oofp[te] = P.max(1)
            if b in raw and tname == 'class':     # MAG-like incompleteness: thin held-out genomes' protein counts
                for p in thin_pred:
                    r = raw[b].loc[tr.index[te]]; rt = pd.DataFrame(rng.binomial(r.values.astype(int), p), index=r.index, columns=r.columns)
                    xt = rel(rt, nprot.loc[r.index] * p) if b != 'tf' else tf_tr(rt)
                    thin_pred[p][te] = m.classes_[m.predict_proba(xt.values).argmax(1)]
        row = dict(target=tname, block=b, n_features=X.shape[1], balanced_acc=balanced_accuracy_score(yv, oof), macro_F1=f1_score(yv, oof, average='macro'), acc=(yv == oof).mean(), median_maxprob_correct=np.median(oofp[yv == oof]))
        if b in raw and tname == 'class':
            for p in thin_pred: row[f'balanced_acc_thin{p}'] = balanced_accuracy_score(yv, thin_pred[p])
        perf.append(row); print(row, flush=True)
        pd.DataFrame({'LOCUSTAG': tr.index, 'true': yv, 'pred': oof, 'maxprob': oofp}).to_csv(f'results/rf_oof_{tname}_{b}.csv', index=False)
        m = RandomForestClassifier(**RF).fit(Xt, yv)
        # predict AO MAGs + unplaced Ascomycota
        q = U.index[U.role.isin(['AO', 'unplaced'])]; P = pd.DataFrame(m.predict_proba(X.loc[q].values), index=q, columns=m.classes_)
        P.insert(0, 'block', b); P.insert(0, 'target', tname); preds.append(P.reset_index())
        im = pd.DataFrame({'target': tname, 'block': b, 'feature': X.columns, 'gini': m.feature_importances_, 'perm_bal_acc': np.nan})
        if tname == 'class':   # permutation importance (3 repeats) for the top-30 Gini features on a genus-grouped 25% hold-out
            a, te = next(StratifiedGroupKFold(4, shuffle=True, random_state=1).split(Xt, yv, tr.genus_key.values)); m2 = RandomForestClassifier(**RF).fit(Xt[a], yv[a])
            top = np.argsort(m.feature_importances_)[::-1][:30]; base = balanced_accuracy_score(yv[te], m2.predict(Xt[te])); pv = np.full(X.shape[1], np.nan)
            for j in top:
                d = []
                for r_ in range(3):
                    Xp = Xt[te].copy(); Xp[:, j] = rng.permutation(Xp[:, j]); d.append(base - balanced_accuracy_score(yv[te], m2.predict(Xp)))
                pv[j] = np.mean(d)
            im['perm_bal_acc'] = pv
        imps.append(im)
pd.DataFrame(perf).to_csv('results/rf_performance.csv', index=False)
pd.concat(preds).to_csv('results/rf_predictions.csv', index=False)
pd.concat(imps).to_csv('results/rf_importance.csv', index=False)
# ---- leave-one-class-out novelty calibration (class target)
loco = []
for b in ['all', 'pfam_ab', 'all_nonpfam', 'codon', 'struct']:
    X = M[b]
    for c in keep:
        a = tr.index[tr.CLASS != c]; te = tr.index[tr.CLASS == c]
        m = RandomForestClassifier(**{**RF, 'n_estimators': 300}).fit(X.loc[a].values, tr.CLASS[a].values); P = m.predict_proba(X.loc[te].values)
        top = pd.Series(m.classes_[P.argmax(1)]).value_counts()
        Pa = m.predict_proba(X.loc[HQ].values)
        loco.append(dict(block=b, left_out=c, n=len(te), median_maxprob=np.median(P.max(1)), q90_maxprob=np.quantile(P.max(1), .9), top_assigned=top.index[0], top_assigned_frac=top.iloc[0] / len(te),
                         AO_top=m.classes_[Pa.mean(0).argmax()], AO_maxprob=Pa.mean(0).max()))
    print(b, 'LOCO done', flush=True)
pd.DataFrame(loco).to_csv('results/rf_loco.csv', index=False)
