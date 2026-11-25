"""Paths and file names. Edit RAW_DIR if your downloads live elsewhere."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.environ.get('GZMK_RAW_DIR', os.path.join(ROOT, 'data', 'raw')) + os.sep
WORK_DIR = os.path.join(ROOT, 'results', 'intermediate') + os.sep
FIG_DIR = os.path.join(ROOT, 'results', 'figures') + os.sep
for d in (WORK_DIR, FIG_DIR): os.makedirs(d, exist_ok=True)
RAW = {
    'GSE63060': RAW_DIR + 'GSE63060_series_matrix.txt.gz',
    'GSE63061': RAW_DIR + 'GSE63061_series_matrix.txt.gz',
    'GSE140829': RAW_DIR + 'GSE140829_series_matrix.txt.gz',
    'GPL6947': RAW_DIR + 'GPL6947_full_table.txt',
    'GPL10558': RAW_DIR + 'GPL10558_full_table.txt',
    'GPL15988': RAW_DIR + 'GPL15988_full_table.txt',
    'TEREKHOVA_MARKERS': RAW_DIR + 'Terekhova2023_TableS5_cluster_markers.xlsx',
}
SEED = 20261116
