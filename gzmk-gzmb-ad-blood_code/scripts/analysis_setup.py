import sys, warnings; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings('ignore')
import json, pandas as pd, numpy as np
from scipy import stats
from core import *
try:
    E, M
except NameError:
    E,M=load()
SIG=json.load(open(W+'signatures.json'))
PROXY={'Neut':['FCGR3B','CSF3R','CXCR2','FPR1','FPR2','MME','CMTM2','PROK2','MGAM','CXCR1'],
 'Mono':['CD14','CSF1R','FCN1','VCAN','CD300E','LILRB2','MS4A7'],
 'B':['CD19','MS4A1','CD79A','CD79B','CD22','BANK1','FCRLA','PAX5'],
 'Tcell':['CD3E','CD3D','CD3G','CD2','CD5','CD6','CD28','LCK','TRAT1'],
 'CD8':['CD8B'],
 'NK':['NCR1','SH2D1B','KIR2DL1','KIR2DL3','KIR3DL2','NCAM1','KIR2DL4']}
# CD8A is in GZMB set -> use CD8B only for proxy; remove any overlap
allsig=set(SIG['GZMK_Tem'])|set(SIG['GZMB_Tem'])
for k in PROXY: PROXY[k]=[g for g in PROXY[k] if g not in allsig]
R={}
for g in COH:
    ex=E[g]; m=M[g]; Rm=rankmat(ex)
    sc=score_sets(ex,{'K':SIG['GZMK_Tem'],'B':SIG['GZMB_Tem'],'Temra':SIG['Temra']},Rm)
    cn=(m.dx=='CN').values
    z=lambda x:(x-x[cn].mean())/x[cn].std()
    m['zK']=z(sc.K); m['zB']=z(sc.B); m['zR']=z(sc.Temra)
    m['KB']=m.zK-m.zB; m['KR']=m.zK-m.zR
    mk=z(meanz(ex,SIG['GZMK_Tem'])); mb=z(meanz(ex,SIG['GZMB_Tem'])); m['KB_mz']=mk-mb
    for p,gl in PROXY.items():
        v=ex.loc[[x for x in gl if x in ex.index]].mean(axis=0); m[p]=(v-v.mean())/v.std()
    m['age10']=m.age/10
    M[g]=m; R[g]=Rm
pd.concat(M.values()).to_csv(W+'scores_all.csv')
COV=['age10','sex','Tcell','CD8','NK','Neut']; COV6=['age10','sex','Tcell','CD8','NK','Neut','Mono','B']
def smd(m,y,contrast,cov=COV,ref='CN',batch=True):
    d=m[m.dx.isin([ref,contrast])]
    D=design(d,cov,None,batch=batch and 'batch' in d.columns and d['batch'].notna().all())
    D.insert(1,'dx',(d.dx==contrast).astype(float).values)
    b,se,p,s,df=ols(d[y].values,D.values)
    return dict(est=b[1]/s,se=se[1]/s,p=p[1],raw=b[1],n_case=int((d.dx==contrast).sum()),n_ref=int((d.dx==ref).sum()))
def run(y,contrast,cohorts=COH,Mx=None,**kw):
    Mx=Mx or M; per={g:smd(Mx[g],y,contrast,**kw) for g in cohorts}
    pooled=reml([per[g]['est'] for g in cohorts],[per[g]['se'] for g in cohorts])
    return dict(per=per,pooled=pooled)
