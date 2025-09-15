import pandas as pd
from pathlib import Path
from itertools import product

# --- Pipeline Paths --- #
# Input paths
xenium_processed_data_dir = Path(config['xenium_processed_data_dir'])
xenium_count_correction_dir = Path(config['xenium_count_correction_dir'])
std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
scrnaseq_processed_data_dir = Path(config['scrnaseq_processed_data_dir'])
palette_dir = Path(config['xenium_metadata_dir'])

# Output paths
results_dir = Path(config['results_dir'])
figures_dir = Path(config['figures_dir'])
seurat_to_h5_dir = results_dir / 'seurat_to_h5'
silhouette_dir = results_dir / 'silhouette'
scib_metrics_results_dir = results_dir / "scib_metrics_panel"

# Palettes
count_correction_palette = palette_dir / 'col_palette_correction_method.csv'
cell_type_palette = palette_dir / 'col_palette_cell_types_combo.csv'
panel_palette = palette_dir / 'col_palette_panel.csv'
sample_palette = palette_dir / 'col_palette_sample.csv'
segmentation_palette = palette_dir / 'col_palette_segmentation.csv'

# --- Params --- #
# QC params
min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5

# resolvi params
max_epochs = 100
num_samples = 30
mixture_k = 50
batch_size = 1000
use_batch = True

# UMAP params
n_comps = 50
n_neighbors = 50
min_dist = 0.5
metric = 'euclidean'

# plot params
s=0.5
alpha=0.5
points_only = True
extension = 'png'
dpi=300
annotation_mode='reference_based'
annotation_normalisation='lognorm'

# list params
signal_integrity_thresholds = [0.5,0.7]
correction_methods = ['raw','split_fully_purified',f'resolvi_panel_{use_batch=}',f'resolvi_panel_supervised_{use_batch=}'] 
correction_methods += [f'ovrlpy_correction_{signal_integrity_threshold=}' for signal_integrity_threshold in signal_integrity_thresholds]
normalisations = ['lognorm']#,'sctransform']
layers = ['data']#,'scale_data']
references = ['matched_reference_combo','external_reference']
methods = ['rctd_class_aware']
levels = ['Level1','Level2.1']
colors = levels + ['sample']
plot_types = ['umap','facet_umap']

# contamination metrics params
use_precomputed = True
radius = 15
n_permutations = 30
n_splits = 5
n_repeats = 5
top_n = 20
scoring = 'precision'
cv_mode = 'spatial'
markers_modes = ['diffexpr']#,'common_markers'] #'/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/markers/cellmarker_cell_types_markers.json'
genes_dict = {
    'all':[],# default: use all genes if empty list
    # 'Xenium_NSCLC_5k_lung_chromium_common_genes':pd.read_csv(config['markers_dir']+'Xenium_NSCLC_5k_lung_chromium_common_genes.csv')['gene'].tolist(),
    # 'Xenium_hLung_v1_metadata':pd.read_csv(config['markers_dir']+'Xenium_hLung_v1_metadata.csv')['Gene'].tolist(),
    # 'CHUV_IO_340_panel':pd.read_csv(config['markers_dir']+'CHUV_IO_340_panel.csv')['Gene ID'].tolist(),
    # 'Xenium_hBreast_v1_metadata':pd.read_csv(config['markers_dir']+'Xenium_hBreast_v1_metadata.csv')['Gene'].tolist()
}
train_modes = ['multivariate']#,'univariate']

# scib metrics params
max_n_cells = 100_000

# allowed levels for each condition
CONDITIONS_LEVELS = {
    'melanoma':['sample','Level1'],
    'mesothelioma_pilot':['sample','Level2.1'],
    'NSCLC':['sample','Level2.1'],
    'breast':['sample','Level2.1'],
}

# allowed references for each condition
CONDITIONS_REFERENCES = {
    'melanoma':['external_reference'],
    'mesothelioma_pilot':['matched_reference_combo'],
    'NSCLC':['matched_reference_combo'],
    'breast':['matched_reference_combo'],
}