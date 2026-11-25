import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
os.chdir(WORK_DIR)
import sys, warnings; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings('ignore')
import json, pandas as pd, numpy as np
from scipy import stats
from core import *
sc=pd.read_csv(W+'GSE63060_sexcheck.csv',index_col=0); disc=json.load(open(W+'gse63060_discordant_cn.json'))
def build(scenario):
    E,M=load()
    if scenario=='inclusive':
        ex=pd.read_pickle(W+'GSE63060_expr.pkl'); meta=pd.read_csv(W+'meta_raw.csv').set_index('sample')
        add=[s for s in disc]
        E['GSE63060']=ex[list(E['GSE63060'].columns)+add]
        mm=meta.loc[E['GSE63060'].columns].copy(); mm['sex']=sc.loc[mm.index,'callR'].values  # expression-inferred sex for all
        M['GSE63060']=mm
    if scenario=='no_gse63060_cn':
        keep=M['GSE63060'].dx!='CN'; M['GSE63060']=M['GSE63060'][keep]; E['GSE63060']=E['GSE63060'][M['GSE63060'].index]
    return E,M
out={}
for scen in ['primary','inclusive','no_gse63060_cn']:
    E,M=build(scen)
    NS=dict(globals()); NS['E']=E; NS['M']=M
    exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'analysis_setup.py')).read(),NS)
    cohorts=COH if scen!='no_gse63060_cn' else ['GSE63061','GSE140829']
    r={}
    for y in ['KB','zK','zB']:
        for c in ['MCI','AD']:
            rr=NS['run'](y,c,cohorts=cohorts); r[f'{y}_{c}']=dict(pooled=rr['pooled'],per={k:v for k,v in rr['per'].items()})
    # age positive control in CN
    per={}
    for c0 in cohorts:
        d=NS['M'][c0][NS['M'][c0].dx=='CN']; D=design(d,['age10','sex','Tcell','CD8','NK','Neut'],None,batch='batch' in d.columns and d['batch'].notna().all())
        b,se,p,s,df=ols(d['KB'].values,D.values); per[c0]=dict(est=b[1],se=se[1])
    r['age_KB']=reml([per[c]['est'] for c in cohorts],[per[c]['se'] for c in cohorts])
    r['n']={c0:NS['M'][c0].dx.value_counts().to_dict() for c0 in cohorts}
    out[scen]=r
    f=lambda q:f"{q['est']:.3f} [{q['lo']:.3f},{q['hi']:.3f}] p={q['p']:.2g} I2={q['I2']:.0f}"
    print('==',scen,r['n']); [print('  ',k,f(v['pooled'])) for k,v in r.items() if k not in('age_KB','n')]; print('   age KB',f(r['age_KB']))
    if scen=='inclusive': print('   GSE63060 per-cohort KB:',{c:round(r[f'KB_{c}']['per']['GSE63060']['est'],3) for c in ['MCI','AD']})
json.dump(out,open(W+'qc_scenarios.json','w'),default=float,indent=1)
