import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
os.chdir(WORK_DIR)
"""Extract ranked cluster markers (top 100 per cluster) from the Terekhova et al. 2023 supplementary marker table."""
import openpyxl, json
wb=openpyxl.load_workbook(RAW['TEREKHOVA_MARKERS'],read_only=True)
res={}
for sh in ['CD8_conventional_T_cells','All_PBMC','NK_cells']:
    rows=list(wb[sh].iter_rows(values_only=True)); hdr=rows[0]
    for j,h in enumerate(hdr):
        if h: res[sh+'|'+h]=[r[j] for r in rows[1:] if r[j]]
json.dump(res,open('terekhova_markers.json','w'))
print({k:len(v) for k,v in res.items() if k.startswith('CD8_conventional')})
