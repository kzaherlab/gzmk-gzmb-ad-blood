import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
os.chdir(WORK_DIR)
import sys, os, warnings; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings('ignore')
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'analysis_setup.py')).read())
OUT={}
# positive control: age in CN
def agebeta(y,cov=['sex','Tcell','CD8','NK','Neut']):
    per={}
    for g in COH:
        d=M[g][M[g].dx=='CN']; D=design(d,['age10']+cov,None,batch='batch' in d.columns)
        b,se,p,s,df=ols(d[y].values,D.values); per[g]=dict(est=b[1],se=se[1],p=p[1],n=len(d))
    return dict(per=per,pooled=reml([per[g]['est'] for g in COH],[per[g]['se'] for g in COH]))
OUT['age']={y:agebeta(y) for y in ['zK','zB','KB']}
OUT['age_unadj']={y:agebeta(y,cov=['sex']) for y in ['zK','zB','KB']}
# primary
OUT['primary']={f'{y}_{c}':run(y,c) for y in ['KB','zK','zB'] for c in ['AD','MCI']}
OUT['AD_vs_MCI']={y:run(y,'AD',**{}) if False else None for y in []}
OUT['AD_vs_MCI']={y:dict(per={g:smd(M[g],y,'AD',ref='MCI') for g in COH}) for y in ['KB','zK','zB']}
for y in OUT['AD_vs_MCI']:
    p=OUT['AD_vs_MCI'][y]['per']; OUT['AD_vs_MCI'][y]['pooled']=reml([p[g]['est'] for g in COH],[p[g]['se'] for g in COH])
# robustness
rob={}
for c in ['AD','MCI']:
    rob[f'six_proxies_{c}']=run('KB',c,cov=COV6)['pooled']
    rob[f'no_composition_{c}']=run('KB',c,cov=['age10','sex'])['pooled']
    rob[f'no_batch_{c}']=run('KB',c,batch=False)['pooled']
    rob[f'meanz_score_{c}']=run('KB_mz',c)['pooled']
    rob[f'Temra_instead_{c}']=run('KR',c)['pooled']
    for g in COH: rob[f'leave_out_{g}_{c}']=run('KB',c,cohorts=[x for x in COH if x!=g])['pooled']
    for s in ['F','M']:
        Ms={g:M[g][M[g].sex==s] for g in COH}
        rob[f'sex_{s}_{c}']=run('KB',c,Mx=Ms,cov=[x for x in COV if x!='sex'])['pooled']
    # inferred sex, keep all GSE63060 controls? (sex-mismatched excluded upstream) -> covered by leave-out GSE63060
OUT['robust']=rob
# fixed-effect pooled on merged data (sanity)
# random gene-set null (matched mean-expression deciles), K/B AD and MCI
rng=np.random.default_rng(20261116)
def null_dist(contrast,nrep=1000):
    obs=OUT['primary'][f'KB_{contrast}']['pooled']['est']; vals=[]
    pre={}
    for g in COH:
        ex=E[g]; mu=ex.mean(axis=1); bins=pd.qcut(mu,10,labels=False); idx={x:i for i,x in enumerate(ex.index)}
        pools={b:np.where(bins.values==b)[0] for b in range(10)}
        pre[g]=(bins,idx,pools)
    for r in range(nrep):
        per=[]
        for g in COH:
            bins,idx,pools=pre[g]; m=M[g].copy(); cn=(m.dx=='CN').values
            def draw(G):
                out=[]
                for x in G:
                    b=bins.loc[x]; out.append(rng.choice(pools[b]))
                return out
            sk=ssgsea(R[g],draw(SIG['GZMK_Tem'])); sb=ssgsea(R[g],draw(SIG['GZMB_Tem']))
            z=lambda x:(x-x[cn].mean())/x[cn].std()
            m['KBr']=z(sk)-z(sb); per.append(smd(m,'KBr',contrast))
        vals.append(reml([p['est'] for p in per],[p['se'] for p in per])['est'])
    vals=np.array(vals); return dict(obs=obs,null_mean=vals.mean(),null_sd=vals.std(),emp_p=float((np.abs(vals)>=abs(obs)).mean()),pct=float((vals<obs).mean()*100),vals=vals.tolist())
OUT['null']={c:null_dist(c) for c in ['AD','MCI']}
# leave-one-gene-out for K/B AD
logo=[]
for sname,key in [('GZMK_Tem','K'),('GZMB_Tem','B')]:
    for gene in SIG[sname]:
        per=[]
        for g in COH:
            ex=E[g]; m=M[g].copy(); cn=(m.dx=='CN').values; idx={x:i for i,x in enumerate(ex.index)}
            s=ssgsea(R[g],[idx[x] for x in SIG[sname] if x!=gene]); z=(s-s[cn].mean())/s[cn].std()
            m['KBl']= (z - m.zB) if key=='K' else (m.zK - z)
            per.append(smd(m,'KBl','AD'))
        logo.append(dict(signature=sname,gene=gene,est=reml([p['est'] for p in per],[p['se'] for p in per])['est']))
OUT['logo']=logo
# gene-level effects AD vs CN (pooled SMD per gene)
genes=[]
for sname in ['GZMK_Tem','GZMB_Tem']:
    for gene in SIG[sname]:
        for c in ['AD','MCI']:
            per=[]
            for g in COH:
                m=M[g].copy(); x=E[g].loc[gene]; cn=(m.dx=='CN').values; m['gx']=(x-x[cn].mean())/x[cn].std()
                per.append(smd(m,'gx',c))
            pr=reml([p['est'] for p in per],[p['se'] for p in per]); genes.append(dict(signature=sname,gene=gene,contrast=c,**{k:pr[k] for k in ['est','lo','hi','p']}))
gdf=pd.DataFrame(genes); gdf['fdr']=np.nan
for c in ['AD','MCI']:
    i=gdf.contrast==c; p=gdf.loc[i,'p'].values; o=np.argsort(p); q=p[o]*len(p)/np.arange(1,len(p)+1); q=np.minimum.accumulate(q[::-1])[::-1]; f=np.empty_like(q); f[o]=np.minimum(q,1); gdf.loc[i,'fdr']=f
gdf.to_csv(W+'gene_level.csv',index=False)
# immune-age gap
for g in COH:
    m=M[g]; cn=m.dx=='CN'
    D=design(m[cn],['zK','zB','sex'],None); b,_,_,_,_=ols(m.loc[cn,'age'].values,D.values)
    Dall=design(m,['zK','zB','sex'],None); m['imm_age']=Dall.values@b; m['gap_raw']=m.imm_age-m.age
    D2=design(m[cn],['age'],None); b2,_,_,_,_=ols(m.loc[cn,'gap_raw'].values,D2.values)
    m['gap']=m.gap_raw-(design(m,['age'],None).values@b2); M[g]=m
def gapyears(c):
    per={}
    for g in COH:
        d=M[g][M[g].dx.isin(['CN',c])]; D=design(d,['sex','Tcell','CD8','NK','Neut'],None,batch='batch' in d.columns); D.insert(1,'dx',(d.dx==c).astype(float).values)
        b,se,p,s,df=ols(d.gap.values,D.values); per[g]=dict(est=b[1],se=se[1],p=p[1])
    return dict(per=per,pooled=reml([per[g]['est'] for g in COH],[per[g]['se'] for g in COH]))
OUT['gap']={c:gapyears(c) for c in ['AD','MCI']}
OUT['age_model_r']={g:float(np.corrcoef(M[g].loc[M[g].dx=='CN','imm_age'],M[g].loc[M[g].dx=='CN','age'])[0,1]) for g in COH}
# APOE (GSE140829)
m=M['GSE140829'].copy(); m['e4']=m.apoe_status.astype(str).str.contains('E4').astype(float)
mm=m[m.apoe_status.notna()&(m.apoe_status.astype(str)!='nan')]
D=design(mm,['e4','age10','sex','Tcell','CD8','NK','Neut'],None,batch=True)
for c in ['MCI','AD']: D.insert(1,'dx_'+c,(mm.dx==c).astype(float).values)
b,se,p,s,df=ols(mm.KB.values,D.values); cols=list(D.columns)
OUT['apoe']=dict(n=len(mm),e4_carriers=int(mm.e4.sum()),e4_est=b[cols.index('e4')]/s,e4_se=se[cols.index('e4')]/s,e4_p=p[cols.index('e4')])
# interaction e4 x AD
D2=D.copy(); D2['e4xAD']=D2['e4']*D2['dx_AD']; b,se,p,s,df=ols(mm.KB.values,D2.values); c2=list(D2.columns)
OUT['apoe']['int_est']=b[c2.index('e4xAD')]/s; OUT['apoe']['int_p']=p[c2.index('e4xAD')]
# composition proxies by dx (for supplementary) + correlation of signatures with proxies
comp=[]
for pname in ['Neut','Mono','B','Tcell','CD8','NK']:
    for c in ['AD','MCI']:
        r=run(pname,c,cov=['age10','sex'])['pooled']; comp.append(dict(proxy=pname,contrast=c,**{k:r[k] for k in ['est','lo','hi','p']}))
OUT['composition']=comp
cors={g:M[g][['zK','zB','KB','Neut','Mono','B','Tcell','CD8','NK','age']].corr(method='spearman').round(3).to_dict() for g in COH}
OUT['cors']=cors
pd.concat(M.values()).to_csv(W+'scores_all.csv')
json.dump(OUT,open(W+'results.json','w'),default=float,indent=1)
def f(r): return f"{r['est']:.3f} [{r['lo']:.3f},{r['hi']:.3f}] p={r['p']:.2g} I2={r.get('I2',float('nan')):.0f}"
print('AGE (CN, per decade, adjusted):'); [print(' ',y,f(OUT['age'][y]['pooled']),{g:round(OUT['age'][y]['per'][g]['est'],3) for g in COH}) for y in OUT['age']]
print('AGE unadj:'); [print(' ',y,f(OUT['age_unadj'][y]['pooled'])) for y in OUT['age_unadj']]
print('PRIMARY:'); [print(' ',k,f(v['pooled']),{g:(round(v['per'][g]['est'],3),round(v['per'][g]['p'],3)) for g in COH}) for k,v in OUT['primary'].items()]
print('AD vs MCI:'); [print(' ',k,f(v['pooled'])) for k,v in OUT['AD_vs_MCI'].items()]
print('ROBUST:'); [print(' ',k,f(v)) for k,v in rob.items()]
print('NULL:',{c:(round(v['obs'],3),round(v['null_mean'],3),round(v['null_sd'],3),v['emp_p'],v['pct']) for c,v in OUT['null'].items()})
print('GAP:',{c:f(v['pooled']) for c,v in OUT['gap'].items()}, OUT['age_model_r'])
print('APOE:',OUT['apoe'])
print('COMP:'); [print(' ',d) for d in comp]
print('LOGO range:',min(x['est'] for x in logo),max(x['est'] for x in logo))
