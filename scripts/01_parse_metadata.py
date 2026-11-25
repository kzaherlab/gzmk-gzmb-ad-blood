import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
os.chdir(WORK_DIR)
import gzip, re, sys, pandas as pd, numpy as np
files={g: RAW[g] for g in ['GSE63060','GSE63061','GSE140829']}
out=[]
for gse,f in files.items():
    gsm=None; chars=[]
    with gzip.open(f,'rt') as h:
        for line in h:
            if line.startswith('!Sample_geo_accession'): gsm=[x.strip('"\n') for x in line.split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1'): chars.append([x.strip('"\n') for x in line.split('\t')[1:]])
            elif line.startswith('!series_matrix_table_begin'): break
    d={'sample':gsm}
    for row in chars:
        for i,v in enumerate(row):
            if ':' in v:
                k,val=v.split(':',1); d.setdefault(k.strip().lower(),[None]*len(gsm))[i]=val.strip()
    df=pd.DataFrame(d); df['cohort']=gse
    dxcol='status' if 'status' in df else 'diagnosis'
    df['dx']=df[dxcol].map(lambda x: {'CTL':'CN','Control':'CN','MCI':'MCI','AD':'AD'}.get(x,None))
    df['age']=pd.to_numeric(df.get('age',df.get('age_at_draw')),errors='coerce')
    sx=df.get('gender',df.get('sex')).str.lower().str[0]; df['sex']=sx.map({'f':'F','m':'M'})
    out.append(df)
    print(gse, df[dxcol].value_counts().to_dict(), 'age NA', df.age.isna().sum())
meta=pd.concat(out,ignore_index=True); meta.to_csv('meta_raw.csv',index=False)
print(meta.columns.tolist())
print(meta.groupby(['cohort','dx']).agg(n=('age','size'),age=('age','mean'),sd=('age','std'),fem=('sex',lambda s:(s=='F').mean()*100)).round(1))
if 'apoe_status' in meta: print(meta.apoe_status.value_counts())
