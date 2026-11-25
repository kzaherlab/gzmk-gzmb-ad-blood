import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
os.chdir(WORK_DIR)
import sys, warnings, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings('ignore')
import json, pandas as pd, numpy as np
from core import *
E,M=load(); NS=dict(globals()); NS['E']=E; NS['M']=M
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'analysis_setup.py')).read(),NS)
R=NS['R']; M=NS['M']; smd=NS['smd']; SIG=NS['SIG']
obs={c:json.load(open(W+'results.json'))['primary'][f'KB_{c}']['pooled']['est'] for c in ['MCI','AD']}
rng=np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 7)
N=int(sys.argv[1]); pre={}
for g in COH:
    ex=E[g]; mu=ex.mean(axis=1); bins=pd.qcut(mu,10,labels=False); pools={b:np.where(bins.values==b)[0] for b in range(10)}
    pre[g]=(bins,pools,{x:bins.loc[x] for x in SIG['GZMK_Tem']+SIG['GZMB_Tem']})
vals={'MCI':[],'AD':[]}; t=time.time()
for r in range(N):
    per={'MCI':[],'AD':[]}
    for g in COH:
        bins,pools,bmap=pre[g]; m=M[g].copy(); cn=(m.dx=='CN').values
        dk=[rng.choice(pools[bmap[x]]) for x in SIG['GZMK_Tem']]; db=[rng.choice(pools[bmap[x]]) for x in SIG['GZMB_Tem']]
        sk=ssgsea(R[g],dk); sb=ssgsea(R[g],db); z=lambda x:(x-x[cn].mean())/x[cn].std()
        m['KBr']=z(sk)-z(sb)
        for c in ['MCI','AD']: per[c].append(smd(m,'KBr',c))
    for c in ['MCI','AD']: vals[c].append(reml([p['est'] for p in per[c]],[p['se'] for p in per[c]])['est'])
    if r==49: print('50 iters',round(time.time()-t,1),'s',flush=True)
out={}
for c in ['MCI','AD']:
    v=np.array(vals[c]); ob=obs[c]; n=len(v)
    out[c]=dict(N=n,obs=ob,null_mean=v.mean(),null_sd=v.std(),n_le=int((v<=ob).sum()),n_abs=int((np.abs(v)>=abs(ob)).sum()),
                p_one=float(((v<=ob).sum()+1)/(n+1)),p_two=float(((np.abs(v)>=abs(ob)).sum()+1)/(n+1)),vals=v.tolist())
    print(c,{k:(round(x,4) if isinstance(x,float) else x) for k,x in out[c].items() if k!='vals'})
json.dump(out,open(W+f'null_{N}.json','w'))
