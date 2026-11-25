import os
from config import WORK_DIR as W
import pandas as pd, numpy as np, json
from scipy import stats, optimize

COH=['GSE63060','GSE63061','GSE140829']
def load():
    meta=pd.read_csv(W+'meta_raw.csv'); keep=set(json.load(open(W+'keep_samples.json')))
    E={}; M={}
    for g in COH:
        ex=pd.read_pickle(W+f'{g}_expr.pkl'); cols=[c for c in ex.columns if c in keep]
        E[g]=ex[cols]; m=meta[meta.cohort==g].set_index('sample').loc[cols].copy()
        if g=='GSE140829': m['batch']=m['batch'].astype(str)
        M[g]=m
    return E,M
def rankmat(ex):
    # descending position 0..N-1
    return (-ex).rank(axis=0,method='first').values-1
def ssgsea(R,genes_idx,alpha=0.25):
    N=R.shape[0]; P=R[genes_idx,:]; k=len(genes_idx)
    w=(N-P)**alpha; Wt=w.sum(axis=0)
    hit=(w*(N-P)).sum(axis=0)/Wt
    miss=(N*(N+1)/2-(N-P).sum(axis=0))/(N-k)
    return hit-miss
def score_sets(ex,sets,R=None):
    if R is None: R=rankmat(ex)
    idx={g:i for i,g in enumerate(ex.index)}
    return pd.DataFrame({k:ssgsea(R,[idx[x] for x in v if x in idx]) for k,v in sets.items()},index=ex.columns)
def meanz(ex,genes):
    g=[x for x in genes if x in ex.index]; X=ex.loc[g]
    return ((X.T-X.mean(axis=1))/X.std(axis=1)).mean(axis=1)
def ols(y,X):
    X=np.asarray(X,float); y=np.asarray(y,float)
    b,res,rk,sv=np.linalg.lstsq(X,y,rcond=None); n,p=X.shape; r=y-X@b; s2=r@r/(n-p)
    cov=s2*np.linalg.pinv(X.T@X); se=np.sqrt(np.diag(cov)); t=b/se; pv=2*stats.t.sf(np.abs(t),n-p)
    return b,se,pv,np.sqrt(s2),n-p
def design(m,terms,dxcontrast=None,batch=False):
    cols={'const':np.ones(len(m))}
    if dxcontrast: cols['dx']=(m.dx==dxcontrast).astype(float).values
    for t in terms:
        if t=='sex': cols['sexM']=(m.sex=='M').astype(float).values
        else: cols[t]=m[t].astype(float).values
    D=pd.DataFrame(cols,index=m.index)
    if batch and 'batch' in m and m['batch'].notna().all():
        bd=pd.get_dummies(m['batch'],drop_first=True,dtype=float); bd=bd.loc[:,bd.sum()>=1]
        D=pd.concat([D,bd],axis=1)
    return D
def reml(y,se):
    y=np.asarray(y); v=np.asarray(se)**2; k=len(y)
    if k==1: return dict(est=y[0],se=se[0],lo=y[0]-1.96*se[0],hi=y[0]+1.96*se[0],p=2*stats.norm.sf(abs(y[0]/se[0])),I2=np.nan,tau2=0,Q=np.nan,Qp=np.nan)
    def nll(t2):
        w=1/(v+t2); mu=(w*y).sum()/w.sum()
        return 0.5*(np.log(v+t2).sum()+np.log(w.sum())+(w*(y-mu)**2).sum())
    t2=optimize.minimize_scalar(nll,bounds=(0,10),method='bounded').x
    if nll(0)<=nll(t2): t2=0.0
    w=1/(v+t2); mu=(w*y).sum()/w.sum(); s=np.sqrt(1/w.sum())
    w0=1/v; mu0=(w0*y).sum()/w0.sum(); Q=(w0*(y-mu0)**2).sum(); I2=max(0,(Q-(k-1))/Q)*100 if Q>0 else 0
    return dict(est=mu,se=s,lo=mu-1.96*s,hi=mu+1.96*s,p=2*stats.norm.sf(abs(mu/s)),I2=I2,tau2=t2,Q=Q,Qp=stats.chi2.sf(Q,k-1))
