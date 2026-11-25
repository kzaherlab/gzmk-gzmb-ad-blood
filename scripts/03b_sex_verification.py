import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
os.chdir(WORK_DIR)
"""Verify sex annotations using RPS4Y1 (decisive, bimodal) and XIST, per cohort and diagnostic group."""
import pandas as pd, numpy as np, json
from scipy import stats
meta=pd.read_csv('meta_raw.csv'); out={}
for g in ['GSE63060','GSE63061','GSE140829']:
    ex=pd.read_pickle(f'{g}_expr.pkl'); m=meta[meta.cohort==g].set_index('sample').loc[ex.columns].copy()
    r=ex.loc['RPS4Y1']; thr=(r.min()+r.max())/2
    m['callR']=np.where(r>thr,'M','F'); m['RPS4Y1']=r.values; m['XIST']=ex.loc['XIST'].values
    gap=(float(r[r<thr].max()), float(r[r>thr].min()))
    disc=pd.crosstab(m.dx,m.callR!=m.sex)
    out[g]=dict(rps4y1_gap=gap, discordant_by_dx={dx:int(disc.loc[dx].get(True,0)) for dx in disc.index}, n_by_dx=m.dx.value_counts().to_dict())
    if g=='GSE63060':
        cn=m[m.dx=='CN']; tab=pd.crosstab(cn.sex,cn.callR); OR,p=stats.fisher_exact(tab.values)
        pa=(cn.sex=='F').mean(); pi=(cn.callR=='F').mean()
        out[g].update(cn_crosstab=tab.to_dict(), cn_fisher_OR=OR, cn_fisher_p=p, cn_observed_concordance=float((cn.sex==cn.callR).mean()),
                      cn_expected_concordance_if_random=float(pa*pi+(1-pa)*(1-pi)))
        m.to_csv('GSE63060_sexcheck.csv')
        json.dump(cn[cn.callR!=cn.sex].index.tolist(),open('gse63060_discordant_cn.json','w'))
print(json.dumps(out,indent=1,default=float)); json.dump(out,open('sex_verification.json','w'),indent=1,default=float)
