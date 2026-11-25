import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
os.chdir(WORK_DIR)
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import json, pandas as pd, numpy as np, re
from core import *
E,M=load()
T=json.load(open(W+'terekhova_markers.json'))
K0=T['CD8_conventional_T_cells|Tem GZMK+']; B0=T['CD8_conventional_T_cells|Tem GZMB+']; R0=T['CD8_conventional_T_cells|Temra']
IEG={'FOS','FOSB','JUN','JUNB','JUND','IER2','IER3','IER5','DUSP1','ZFP36','ZFP36L1','ZFP36L2','EGR1','HSPA1A','HSPA1B','NR4A1','NR4A2','ATF3','KLF6','CD69','PPP1R15A','RGS1','RGS2','TXNIP'}
nonT=set(T['All_PBMC|Myeloid cells'])|set(T['All_PBMC|B cells'])
def clean(L):
    return [g for g in L if not re.match(r'^(IG[HKL][VJC]|RP[LS]\d|MT-|HLA-D)',g) and g not in IEG and g not in nonT]
shared=set(K0)&set(B0)
K1=[g for g in clean(K0) if g not in shared]; B1=[g for g in clean(B0) if g not in shared]
R1=[g for g in clean(R0) if g not in set(K0)]
common=set.intersection(*[set(E[g].index) for g in COH])
K2=[g for g in K1 if g in common]; B2=[g for g in B1 if g in common]; R2=[g for g in R1 if g in common]
# whole-blood coherence filter (label-agnostic): leave-one-out correlation with rest of signature >=0.2 in every cohort,
# and stronger correlation with own signature than with neutrophil proxy
NEUT=['FCGR3B','CSF3R','CXCR2','FPR1','FPR2','MME','CMTM2','PROK2','MGAM','CXCR1']
def coherent(G):
    keep=[]; rows=[]
    for gene in G:
        ok=True; rr=[]
        for c in COH:
            ex=E[c]; z=((ex.loc[G].T-ex.loc[G].mean(axis=1))/ex.loc[G].std(axis=1))
            loo=z.drop(columns=gene).mean(axis=1); r=np.corrcoef(z[gene],loo)[0,1]
            nz=meanz(ex,NEUT); rn=np.corrcoef(ex.loc[gene],nz)[0,1]
            rr.append((r,rn)); ok&= (r>=0.2) and (r>rn)
        rows.append((gene,ok,*[x for p in rr for x in p])); 
        if ok: keep.append(gene)
    return keep, pd.DataFrame(rows,columns=['gene','kept']+[f'{c}_{k}' for c in COH for k in ('r_loo','r_neut')])
K3,tk=coherent(K2); B3,tb=coherent(B2); R3,tr=coherent(R2)
print('K',len(K0),'->clean',len(K1),'->platform',len(K2),'->coherent',len(K3)); print(K3)
print('B',len(B0),'->clean',len(B1),'->platform',len(B2),'->coherent',len(B3)); print(B3)
print('Temra',len(R3)); print(R3)
tk['signature']='GZMK-Tem'; tb['signature']='GZMB-Tem'; tr['signature']='Temra'
S1=pd.concat([tk,tb,tr]); S1.to_csv(W+'TableS1_signature_filtering.csv',index=False)
json.dump({'GZMK_Tem':K3,'GZMB_Tem':B3,'Temra':R3,'K_precoherence':K2,'B_precoherence':B2,'shared_removed':sorted(shared)},open(W+'signatures.json','w'),indent=1)
