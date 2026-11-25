import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
os.chdir(WORK_DIR)
import pandas as pd, numpy as np, json
meta=pd.read_csv('meta_raw.csv'); E={g:pd.read_pickle(f'{g}_expr.pkl') for g in ['GSE63060','GSE63061','GSE140829']}
log=[]; keepall=[]
for g,ex in E.items():
    m=meta[meta.cohort==g].set_index('sample').loc[ex.columns]
    n0=len(m); ok=m.dx.notna()&m.age.notna()&m.sex.notna()
    log.append((g,'raw',n0)); log.append((g,'excluded: diagnosis not CN/MCI/AD or missing',int((~ok).sum())))
    y=[x for x in ['RPS4Y1','DDX3Y','KDM5D','EIF1AY','UTY'] if x in ex.index]
    sc=ex.loc['XIST']-ex.loc[y].mean()
    thr=(sc.max()+sc.min())/2
    # 2-means 1D
    c=np.array([sc.min(),sc.max()])
    for _ in range(50):
        lab=np.abs(sc.values[:,None]-c[None,:]).argmin(1); c=np.array([sc.values[lab==k].mean() for k in (0,1)])
    pred=np.where(lab==1,'F','M'); mis=(pred!=m.sex.values)&ok.values
    log.append((g,'excluded: sex mismatch (XIST vs Y genes)',int(mis.sum())))
    ok=ok&~mis
    X=ex.loc[:,ok.values]; v=X.var(axis=1).sort_values(ascending=False).index[:5000]
    Z=X.loc[v].T; Z=(Z-Z.mean())/Z.std()
    U,S,Vt=np.linalg.svd(Z.values-Z.values.mean(axis=0),full_matrices=False); pc=U[:,:2]*S[:2]
    pcz=(pc-pc.mean(axis=0))/pc.std(axis=0); out=(np.abs(pcz)>3).any(1)
    log.append((g,'excluded: PCA outlier (>3 SD on PC1/PC2)',int(out.sum())))
    keep=X.columns[~out]; keepall+=list(keep)
    print(g, 'sex-score separation centers', c.round(2))
# duplicates across AddNeuroMed batches
a,b=E['GSE63060'],E['GSE63061']; gs=a.index.intersection(b.index)
A=a.loc[gs,[s for s in a.columns if s in keepall]]; B=b.loc[gs,[s for s in b.columns if s in keepall]]
v=pd.concat([A,B],axis=1).var(axis=1).sort_values(ascending=False).index[:5000]
Az=A.loc[v]; Bz=B.loc[v]
Az=(Az-Az.mean())/Az.std(); Bz=(Bz-Bz.mean())/Bz.std()
r=(Az.T.values@Bz.values)/(len(v)-1)
dups=np.argwhere(r>0.99); print('cross-batch pairs r>0.99:',len(dups),' max r',r.max().round(3))
# within-cohort duplicates
for g in E:
    X=E[g][[s for s in E[g].columns if s in keepall]]; vv=X.var(axis=1).sort_values(ascending=False).index[:5000]
    Xz=(X.loc[vv]-X.loc[vv].mean())/X.loc[vv].std(); R=(Xz.T.values@Xz.values)/(len(vv)-1); np.fill_diagonal(R,0)
    print(g,'within-cohort max r',R.max().round(3),'pairs>0.99',int((R>0.99).sum()/2))
drop=set(B.columns[dups[:,1]]) if len(dups) else set()
log.append(('GSE63061','excluded: duplicate of GSE63060 sample',len(drop)))
keepall=[s for s in keepall if s not in drop]
json.dump(keepall,open('keep_samples.json','w'))
qc=pd.DataFrame(log,columns=['cohort','step','n']); qc.to_csv('qc_log.csv',index=False); print(qc.to_string())
mk=meta[meta['sample'].isin(keepall)]
print(mk.groupby(['cohort','dx']).size().unstack())

# within-GSE140829 near-duplicate arrays (r > 0.99): keep one of each pair
ex=pd.read_pickle('GSE140829_expr.pkl'); k=json.load(open('keep_samples.json')); X=ex[[s for s in ex.columns if s in k]]
v=X.var(axis=1).sort_values(ascending=False).index[:5000]; Z=(X.loc[v]-X.loc[v].mean())/X.loc[v].std()
R=(Z.T.values@Z.values)/(len(v)-1); np.fill_diagonal(R,0); pr=np.argwhere(np.triu(R)>0.99)
drop=[X.columns[j] for i,j in pr]; k=[s for s in k if s not in drop]; json.dump(k,open('keep_samples.json','w'))
q=pd.read_csv('qc_log.csv'); q.loc[len(q)]=['GSE140829','excluded: near-duplicate array (r>0.99)',len(drop)]; q.to_csv('qc_log.csv',index=False)
print('final n',len(k))
