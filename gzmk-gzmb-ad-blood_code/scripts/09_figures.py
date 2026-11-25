import os
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
os.chdir(WORK_DIR)
import sys,warnings; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); warnings.filterwarnings('ignore')
import json, pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from config import WORK_DIR as W; from config import FIG_DIR as F
O=json.load(open(W+'results.json')); S=pd.read_csv(W+'scores_all.csv',index_col=0); SIG=json.load(open(W+'signatures.json'))
COH=['GSE63060','GSE63061','GSE140829']; LAB={'GSE63060':'AddNeuroMed-1\n(GSE63060)','GSE63061':'AddNeuroMed-2\n(GSE63061)','GSE140829':'UCSF/UCLA\n(GSE140829)'}
LAB1={'GSE63060':'AddNeuroMed-1','GSE63061':'AddNeuroMed-2','GSE140829':'UCSF/UCLA'}
C_CN,C_MCI,C_AD='#2a78d6','#eb6834','#1baf7a'; CK,CB='#2a78d6','#eb6834'; INK='#0b0b0b'; INK2='#52514e'; GRID='#d9d8d4'
DXC={'CN':C_CN,'MCI':C_MCI,'AD':C_AD}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':INK2,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2,'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,'legend.frameon':False,'savefig.dpi':300})
def tag(ax,t): ax.text(-0.12,1.06,t,transform=ax.transAxes,fontsize=11,fontweight='bold',va='bottom',ha='left',color=INK)
def forest(ax,rows,xlabel='Adjusted SMD vs. CN (95% CI)',color=INK):
    # rows: list of (label, est, lo, hi, kind)  kind in {'cohort','pooled','head'}
    y=len(rows)
    for i,(lab,e,lo,hi,kind) in enumerate(rows):
        yy=y-i
        if kind=='head': continue
        if kind=='pooled':
            ax.plot([lo,hi],[yy,yy],color=color,lw=1.2,solid_capstyle='round')
            ax.fill([lo,e,hi,e],[yy,yy+0.28,yy,yy-0.28],color=color,zorder=3)
        else:
            ax.plot([lo,hi],[yy,yy],color=INK2,lw=1,solid_capstyle='round'); ax.plot(e,yy,'s',ms=4.5,color=INK2,zorder=3)
    ax.axvline(0,color=INK2,lw=0.6,ls=(0,(3,2)))
    ax.set_yticks([y-i for i,r in enumerate(rows) if r[4]!='head']); ax.set_yticklabels([r[0] for r in rows if r[4]!='head'],fontsize=7)
    ax.set_ylim(0.3,y+0.7); ax.set_xlabel(xlabel); ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False)
    ax.grid(axis='x',color=GRID,lw=0.4); ax.set_axisbelow(True)
def rows_for(key,src='primary'):
    r=O[src][key]; out=[]
    for g in COH:
        p=r['per'][g]; out.append((f"{LAB1[g]}  ({p['n_case']}/{p['n_ref']})",p['est'],p['est']-1.96*p['se'],p['est']+1.96*p['se'],'cohort'))
    q=r['pooled']; out.append((f"Pooled (RE), I²={q['I2']:.0f}%",q['est'],q['lo'],q['hi'],'pooled')); return out
# ---------- Figure 1: signatures ----------
T1=pd.read_csv(W+'TableS1_signature_filtering.csv')
fig=plt.figure(figsize=(7.2,5.2)); gs=fig.add_gridspec(1,3,width_ratios=[1.05,1.15,1.0],wspace=0.6)
ax=fig.add_subplot(gs[0]); steps=['Atlas markers\n(top 100)','After artefact/\nlineage filter','On all three\nplatforms','Coherent in\nwhole blood']
k=[100,len(json.load(open(W+'terekhova_markers.json'))['CD8_conventional_T_cells|Tem GZMK+'])]
cntK=[100,None,len(SIG['K_precoherence']),len(SIG['GZMK_Tem'])]; cntB=[100,None,len(SIG['B_precoherence']),len(SIG['GZMB_Tem'])]
cntK[1]=63; cntB[1]=75
yy=np.arange(4)[::-1]; h=0.36
ax.barh(yy+h/2+0.01,cntK,h,color=CK,label='GZMK-Tem'); ax.barh(yy-h/2-0.01,cntB,h,color=CB,label='GZMB-Tem')
for i in range(4): ax.text(cntK[i]+2,yy[i]+h/2,str(cntK[i]),va='center',fontsize=7,color=INK); ax.text(cntB[i]+2,yy[i]-h/2,str(cntB[i]),va='center',fontsize=7,color=INK)
ax.set_yticks(yy); ax.set_yticklabels(steps,fontsize=7); ax.set_xlabel('Genes'); ax.set_xlim(0,118); ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.12),ncol=1,fontsize=7); ax.tick_params(axis='y',length=0); tag(ax,'A')
ax=fig.add_subplot(gs[1]); kept=T1[T1.kept & T1.signature.isin(['GZMK-Tem','GZMB-Tem'])].copy()
kept['r']=kept[[f'{c}_r_loo' for c in COH]].mean(axis=1); kept['rmin']=kept[[f'{c}_r_loo' for c in COH]].min(1); kept['rmax']=kept[[f'{c}_r_loo' for c in COH]].max(1)
kept=kept.sort_values(['signature','r'],ascending=[False,True]).reset_index(drop=True)
for i,rw in kept.iterrows():
    c=CK if rw.signature=='GZMK-Tem' else CB; ax.plot([rw.rmin,rw.rmax],[i,i],color=c,lw=0.9,alpha=0.6); ax.plot(rw.r,i,'o',ms=3.2,color=c)
ax.set_yticks(range(len(kept))); ax.set_yticklabels(kept.gene,fontsize=5.6,style='italic'); ax.set_ylim(-1,len(kept))
ax.set_xlabel('Leave-one-out correlation with\nrest of signature (mean, range)'); ax.axvline(0.2,color=INK2,lw=0.6,ls=(0,(3,2))); ax.tick_params(axis='y',length=0); tag(ax,'B')
ax=fig.add_subplot(gs[2]); rowsv=['zK','zB','KB']; colsv=['Tcell','CD8','NK','Neut','Mono','B','age']
mat=np.array([[np.mean([O['cors'][g][c][r] for g in COH]) for r in rowsv] for c in colsv])
from matplotlib.colors import LinearSegmentedColormap
cm=LinearSegmentedColormap.from_list('div',['#2a78d6','#f0efec','#eb6834'])
im=ax.imshow(mat,cmap=cm,vmin=-1,vmax=1,aspect='auto')
for i in range(len(colsv)):
    for j in range(3): ax.text(j,i,f'{mat[i,j]:.2f}',ha='center',va='center',fontsize=7,color=INK)
ax.set_xticks(range(3)); ax.set_xticklabels(['GZMK-\nTem','GZMB-\nTem','K/B\nindex'],fontsize=7); ax.set_yticks(range(len(colsv))); ax.set_yticklabels(['T cells','CD8 (CD8B)','NK cells','Neutrophils','Monocytes','B cells','Age'],fontsize=7)
ax.tick_params(length=0); [sp.set_visible(False) for sp in ax.spines.values()]; ax.xaxis.tick_top()
cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=0.05,pad=0.04); cb.set_label('Spearman ρ (mean of 3 cohorts)',fontsize=6.5); cb.ax.tick_params(labelsize=6); tag(ax,'C')
fig.savefig(F+'Figure1.png',bbox_inches='tight'); fig.savefig(F+'Figure1.pdf',bbox_inches='tight'); plt.close(fig)
# ---------- Figure 2: age ----------
fig=plt.figure(figsize=(7.2,5.0)); gsa=fig.add_gridspec(2,3,height_ratios=[1,1.25],hspace=0.5,wspace=0.35)
axs=[fig.add_subplot(gsa[0,i]) for i in range(3)]+[fig.add_subplot(gsa[1,1:])]
for ax,g,t in zip(axs[:3],COH,'ABC'):
    d=S[(S.cohort==g)&(S.dx=='CN')]; ax.scatter(d.age,d.KB,s=7,color=C_CN,alpha=0.55,edgecolors='white',linewidths=0.3)
    b=np.polyfit(d.age,d.KB,1); xs=np.linspace(d.age.min(),d.age.max(),20); ax.plot(xs,np.polyval(b,xs),color=INK,lw=1.4)
    ax.set_title(f'{LAB1[g]}\nCN, n = {len(d)}',fontsize=7.5,color=INK); ax.set_xlabel('Age (years)'); 
    if g==COH[0]: ax.set_ylabel('K/B index (SD units)')
    ax.grid(color=GRID,lw=0.4); ax.set_axisbelow(True); tag(ax,t)
ax=axs[3]; rows=[]
for y,lab in [('zK','GZMK-Tem'),('zB','GZMB-Tem'),('KB','K/B index')]:
    rows.append((lab,0,0,0,'head'))
    for g in COH:
        p=O['age'][y]['per'][g]; rows.append((LAB1[g],p['est'],p['est']-1.96*p['se'],p['est']+1.96*p['se'],'cohort'))
    q=O['age'][y]['pooled']; rows.append((f'Pooled, I²={q["I2"]:.0f}%',q['est'],q['lo'],q['hi'],'pooled'))
forest(ax,rows,xlabel='β per decade of age (SD units, 95% CI)'); tag(ax,'D')
for i,r in enumerate(rows):
    if r[4]=='head': ax.text(ax.get_xlim()[0],len(rows)-i,r[0],ha='left',va='center',fontsize=7,fontweight='bold',color=INK)
fig.savefig(F+'Figure2.png',bbox_inches='tight'); fig.savefig(F+'Figure2.pdf',bbox_inches='tight'); plt.close(fig)
# ---------- Figure 3: primary ----------
fig=plt.figure(figsize=(7.2,5.6)); gs=fig.add_gridspec(2,3,height_ratios=[1,1.15],hspace=0.55,wspace=0.95)
rng=np.random.default_rng(1)
for j,g in enumerate(COH):
    ax=fig.add_subplot(gs[0,j]); d=S[S.cohort==g]
    for i,dx in enumerate(['CN','MCI','AD']):
        v=d.loc[d.dx==dx,'KB'].values
        bp=ax.boxplot(v,positions=[i],widths=0.5,showfliers=False,patch_artist=True,medianprops=dict(color=INK,lw=1.2),whiskerprops=dict(color=INK2,lw=0.7),capprops=dict(color=INK2,lw=0.7),boxprops=dict(facecolor='white',edgecolor=INK2,lw=0.7))
        ax.scatter(i+rng.uniform(-0.18,0.18,len(v)),v,s=4,color=DXC[dx],alpha=0.5,linewidths=0,zorder=3)
        ax.text(i,ax.get_ylim()[0] if False else np.nan,'')
    ax.set_xticks([0,1,2]); ax.set_xticklabels([f'CN\n{(d.dx=="CN").sum()}',f'MCI\n{(d.dx=="MCI").sum()}',f'AD\n{(d.dx=="AD").sum()}'])
    ax.set_title(LAB1[g],fontsize=8,color=INK); 
    if j==0: ax.set_ylabel('K/B index (SD units)')
    ax.grid(axis='y',color=GRID,lw=0.4); ax.set_axisbelow(True); tag(ax,'ABC'[j])
gs2=gs[1,:].subgridspec(1,3,wspace=0.12,width_ratios=[1,0.62,0.62])
for j,(y,ttl) in enumerate([('KB','K/B index'),('zK','GZMK-Tem score'),('zB','GZMB-Tem score')]):
    ax=fig.add_subplot(gs2[j]); rows=[('MCI vs. CN',0,0,0,'head')]+rows_for(f'{y}_MCI')+[('AD vs. CN',0,0,0,'head')]+rows_for(f'{y}_AD')
    forest(ax,rows,xlabel='Adjusted SMD (95% CI)',color=INK)
    if j>0: ax.set_yticklabels([])
    for i,r in enumerate(rows):
        if r[4]=='head': ax.text(ax.get_xlim()[0],len(rows)-i,r[0],ha='left',va='center',fontsize=7,fontweight='bold',color=INK)
    ax.set_title(ttl,fontsize=8,color=INK); tag(ax,'DEF'[j])
fig.savefig(F+'Figure3.png',bbox_inches='tight'); fig.savefig(F+'Figure3.pdf',bbox_inches='tight'); plt.close(fig)
# ---------- Figure 4: robustness ----------
R=O['robust']; P=O['primary']
names=[('Primary model',None),('Six composition proxies','six_proxies'),('No composition adjustment','no_composition'),('No batch term (GSE140829)','no_batch'),('Mean-z scoring','meanz_score'),('Temra in place of GZMB-Tem','Temra_instead'),('Without AddNeuroMed-1','leave_out_GSE63060'),('Without AddNeuroMed-2','leave_out_GSE63061'),('Without UCSF/UCLA','leave_out_GSE140829'),('Women only','sex_F'),('Men only','sex_M')]
fig=plt.figure(figsize=(7.2,3.9)); gs=fig.add_gridspec(2,2,width_ratios=[1.6,1],hspace=0.6,wspace=0.6)
axL=fig.add_subplot(gs[:,0]); rows=[]
for c,col in [('MCI',C_MCI),('AD',C_AD)]:
    rows.append((f'{c} vs. CN',0,0,0,'head'))
    for lab,key in names:
        q=P[f'KB_{c}']['pooled'] if key is None else R[f'{key}_{c}']
        rows.append((lab,q['est'],q['lo'],q['hi'],'pooled' if key is None else 'cohort'))
forest(axL,rows,xlabel='Pooled adjusted SMD of K/B index (95% CI)')
for i,r in enumerate(rows):
    if r[4]=='head': axL.text(axL.get_xlim()[0],len(rows)-i,r[0],ha='left',va='center',fontsize=7,fontweight='bold',color=INK)
tag(axL,'A')
NULLF=W+'null_10000.json'
NL=json.load(open(NULLF)) if os.path.exists(NULLF) else None
for j,c in enumerate(['MCI','AD']):
    ax=fig.add_subplot(gs[j,1])
    if NL: nv=np.array(NL[c]['vals']); ob=NL[c]['obs']; ptxt=f"two-sided p = {NL[c]['p_two']:.3f}"; nlab=f"{len(nv):,}"
    else: nv=np.array(O['null'][c]['vals']); ob=O['null'][c]['obs']; ptxt=f"two-sided p = {((np.abs(nv)>=abs(ob)).sum()+1)/(len(nv)+1):.3f}"; nlab=f"{len(nv):,}"
    ax.hist(nv,bins=50,color='#c3c2b7',edgecolor='white',lw=0.3); ax.axvline(ob,color=DXC[c],lw=1.8)
    ax.text(ob,ax.get_ylim()[1]*0.92,f'  observed\n  {ob:.2f}',color=INK,fontsize=6.5,va='top',ha='left')
    ax.set_title(f'{c} vs. CN ({nlab} random pairs)\n{ptxt}',fontsize=7,color=INK,loc='center')
    ax.set_xlabel('Pooled SMD'); ax.set_ylabel('Count'); ax.text(-0.42,1.2,'BC'[j],transform=ax.transAxes,fontsize=11,fontweight='bold',color=INK)
fig.savefig(F+'Figure4.png',bbox_inches='tight'); fig.savefig(F+'Figure4.pdf',bbox_inches='tight'); plt.close(fig)
# ---------- Figure 5: gene level ----------
G=pd.read_csv(W+'gene_level.csv')
fig,axs=plt.subplots(1,2,figsize=(7.2,6.2),sharey=True,gridspec_kw={'wspace':0.08})
order=G[G.contrast=='MCI'].sort_values(['signature','est'],ascending=[False,True])
genes=list(order.gene); sigs=list(order.signature)
for ax,c in zip(axs,['MCI','AD']):
    d=G[G.contrast==c].set_index('gene').loc[genes]
    for i,gn in enumerate(genes):
        col=CK if sigs[i]=='GZMK_Tem' else CB; r=d.loc[gn]
        ax.plot([r.lo,r.hi],[i,i],color=col,lw=0.8,alpha=0.7); ax.plot(r.est,i,'o',ms=3.4 if r.fdr<0.05 else 3.0,color=col if r.fdr<0.05 else 'white',markeredgecolor=col,mew=0.8,zorder=3)
    ax.axvline(0,color=INK2,lw=0.6,ls=(0,(3,2))); ax.set_title(f'{c} vs. CN',fontsize=8,color=INK); ax.set_xlabel('Pooled adjusted SMD (95% CI)')
    ax.grid(axis='x',color=GRID,lw=0.4); ax.set_axisbelow(True); ax.axhline(sum(1 for s in sigs if s=='GZMK_Tem')-0.5,color=GRID,lw=0.8)
axs[0].set_yticks(range(len(genes))); axs[0].set_yticklabels(genes,fontsize=5.2,style='italic'); axs[0].tick_params(axis='y',length=0); axs[1].tick_params(axis='y',length=0)
h=[Line2D([],[],marker='o',ls='',color=CK,label='GZMK-Tem gene'),Line2D([],[],marker='o',ls='',color=CB,label='GZMB-Tem gene'),Line2D([],[],marker='o',ls='',color=INK2,markerfacecolor='white',label='FDR ≥ 0.05 (open)')]
axs[1].legend(handles=h,loc='upper right',fontsize=6.5,bbox_to_anchor=(1.0,1.0)); tag(axs[0],'A'); tag(axs[1],'B')
fig.savefig(F+'Figure5.png',bbox_inches='tight'); fig.savefig(F+'Figure5.pdf',bbox_inches='tight'); plt.close(fig)
# ---------- Figure S1: composition proxies ; S2: sex QC ----------
fig,axs=plt.subplots(1,2,figsize=(7.2,2.8),gridspec_kw={'width_ratios':[1.2,1],'wspace':0.6})
rows=[]; 
for c in ['MCI','AD']:
    rows.append((f'{c} vs. CN',0,0,0,'head'))
    for d in O['composition']:
        if d['contrast']==c: rows.append(({'Neut':'Neutrophils','Mono':'Monocytes','B':'B cells','Tcell':'T cells','CD8':'CD8 (CD8B)','NK':'NK cells'}[d['proxy']],d['est'],d['lo'],d['hi'],'pooled'))
forest(axs[0],rows,xlabel='Pooled SMD, adjusted for age and sex')
for i,r in enumerate(rows):
    if r[4]=='head': axs[0].text(axs[0].get_xlim()[0],len(rows)-i,r[0],ha='left',va='center',fontsize=7,fontweight='bold',color=INK)
tag(axs[0],'A')
import pickle
ex=pd.read_pickle(W+'GSE63060_expr.pkl'); meta=pd.read_csv(W+'meta_raw.csv').set_index('sample'); m=meta.loc[ex.columns]
ax=axs[1]; pos=0; ticks=[]; tl=[]
for dx in ['CN','MCI','AD']:
    for sx,col in [('F',C_MCI),('M',C_CN)]:
        v=ex.loc['RPS4Y1',(m.dx==dx)&(m.sex==sx)].values; ax.scatter(pos+rng.uniform(-0.2,0.2,len(v)),v,s=5,color=col,alpha=0.7,linewidths=0); ticks.append(pos); tl.append(f'{dx}\n{sx}'); pos+=1
    pos+=0.6
ax.set_xticks(ticks); ax.set_xticklabels(tl,fontsize=6.5); ax.set_ylabel('RPS4Y1 (log2)'); ax.set_title('GSE63060: Y-gene expression by annotated sex',fontsize=7.5,color=INK); tag(ax,'B')
fig.savefig(F+'FigureS1.png',bbox_inches='tight'); fig.savefig(F+'FigureS1.pdf',bbox_inches='tight'); plt.close(fig)
print('ok')

# ---------- Figure S2: verification of GSE63060 sex annotations ----------
SC=pd.read_csv(W+'GSE63060_sexcheck.csv',index_col=0); QS=json.load(open(W+'qc_scenarios.json'))
fig=plt.figure(figsize=(7.2,5.6)); gs0=fig.add_gridspec(2,1,height_ratios=[1,0.9],hspace=0.55); gs=gs0[0].subgridspec(1,2,wspace=0.3)
for j,(grp,ttl) in enumerate([(['MCI','AD'],'MCI + AD (n = 225)'),(['CN'],'Controls (n = 104)')]):
    ax=fig.add_subplot(gs[j]); d=SC[SC.dx.isin(grp)]
    for sx,col,lab in [('F',C_MCI,'annotated female'),('M',C_CN,'annotated male')]:
        e=d[d.sex==sx]; ax.scatter(e.XIST,e.RPS4Y1,s=8,color=col,alpha=0.7,linewidths=0,label=lab)
    ax.axhline(9.5,color=INK2,lw=0.6,ls=(0,(3,2))); ax.set_xlabel('XIST (log2)')
    if j==0: ax.set_ylabel('RPS4Y1 (log2)')
    ax.set_title(f'GSE63060 {ttl}',fontsize=7.5,color=INK); tag(ax,'AB'[j])
    if j==1: ax.legend(fontsize=6.5,loc='upper left',bbox_to_anchor=(1.02,1.0))
ax=fig.add_subplot(gs0[1].subgridspec(1,3,width_ratios=[0.55,1,0.25])[1]); rows=[]
lab={'primary':'Prespecified (46 excluded)','inclusive':'All 104 kept (inferred sex)','no_gse63060_cn':'All GSE63060 controls excluded'}
for c in ['MCI','AD']:
    rows.append((f'{c} vs. CN',0,0,0,'head'))
    for s in ['primary','inclusive','no_gse63060_cn']:
        q=QS[s][f'KB_{c}']['pooled']; rows.append((lab[s],q['est'],q['lo'],q['hi'],'pooled' if s=='primary' else 'cohort'))
forest(ax,rows,xlabel='Pooled SMD of K/B index (95% CI)')
for i,r in enumerate(rows):
    if r[4]=='head': ax.text(ax.get_xlim()[0],len(rows)-i,r[0],ha='left',va='center',fontsize=7,fontweight='bold',color=INK)
tag(ax,'C')
fig.savefig(F+'FigureS2.png',bbox_inches='tight'); fig.savefig(F+'FigureS2.pdf',bbox_inches='tight'); plt.close(fig)
print('S2 ok')
