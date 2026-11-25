import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
os.chdir(WORK_DIR)
import gzip, pandas as pd, numpy as np, io, sys
files={'GSE63060':(RAW['GSE63060'],RAW['GPL6947'],'ID','Symbol'),
       'GSE63061':(RAW['GSE63061'],RAW['GPL10558'],'ID','Symbol'),
       'GSE140829':(RAW['GSE140829'],RAW['GPL15988'],'ID','GENE SYMBOL')}
for gse,(mf,gpl,idc,symc) in files.items():
    lines=[]; on=False
    with gzip.open(mf,'rt') as h:
        for l in h:
            if l.startswith('!series_matrix_table_begin'): on=True; continue
            if l.startswith('!series_matrix_table_end'): break
            if on: lines.append(l)
    ex=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',index_col=0)
    ex.columns=[c.strip('"') for c in ex.columns]; ex.index=[i.strip('"') for i in ex.index]
    ann=pd.read_csv(gpl,sep='\t',comment='#',dtype=str,usecols=[idc,symc],low_memory=False).dropna()
    ann=ann[ann[symc].str.strip()!='']; m=dict(zip(ann[idc],ann[symc].str.strip()))
    sym=pd.Series(ex.index.map(lambda x:m.get(x)),index=ex.index)
    ex=ex.loc[sym.notna()]; sym=sym[sym.notna()]
    ex=ex.apply(pd.to_numeric,errors='coerce')
    print(gse,'raw range',np.nanmin(ex.values).round(2),np.nanmax(ex.values).round(2),'probes mapped',len(ex))
    if np.nanmax(ex.values)>50: ex=np.log2(ex.clip(lower=1))
    mu=ex.mean(axis=1); o=mu.sort_values(ascending=False).index
    ex=ex.loc[o]; s=sym.loc[o]; keep=~s.duplicated().values
    ex=ex.loc[keep]; ex.index=s[keep].values
    ex.to_pickle(f'{gse}_expr.pkl'); print(gse, ex.shape)
