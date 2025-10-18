import dask

dask.config.set({"dataframe.query-planning": False})

import argparse
import pandas as pd
import scanpy as sc
import scipy
from pathlib import Path
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

# Set up argument parser
parser = argparse.ArgumentParser(description="Embed panel of Xenium donors.")
parser.add_argument("--panel", type=Path, help="Path to the panel file.")
parser.add_argument("--cell_type_annotation_dir", type=Path, help="Path to the cell type annotation dir.")
parser.add_argument("--out_file", type=str, help="Path to the output file.")
parser.add_argument("--normalisation", type=str, help="Normalisation method")
parser.add_argument("--layer", type=str, help="Name of saved layer of the seurat object, data or scale_data")
parser.add_argument("--reference", type=str, help="annotation reference")
parser.add_argument("--method", type=str, help="annotation method")
parser.add_argument("--level", type=str, help="annotation level")
parser.add_argument("--n_comps", type=int, help="Number of components.")
parser.add_argument("--min_counts", type=int, help="QC parameter from pipeline config")
parser.add_argument("--min_features", type=int, help="QC parameter from pipeline config")
parser.add_argument("--max_counts", type=float, help="QC parameter from pipeline config")
parser.add_argument("--max_features", type=float, help="QC parameter from pipeline config")
parser.add_argument("--min_cells", type=int, help="QC parameter from pipeline config")
parser.add_argument("--max_n_cells", type=int, help="Max number of cells to use.")
parser.add_argument("--num_samples", type=int, help="RESOLVI parameter from pipeline config")
parser.add_argument("--mixture_k", type=int, help="RESOLVI parameter from pipeline config")
parser.add_argument("--raw_corrected_counts", action="store_true", help="Use raw corrected counts")
parser.add_argument("--xenium_count_correction_dir", type=Path, help="xenium_count_correction_dir")
parser.add_argument("--results_dir", type=Path, help="results_dir")
parser.add_argument("--correction_method", type=str, help="correction_method")

args = parser.parse_args()

# Access the arguments
panel = args.panel
cell_type_annotation_dir = args.cell_type_annotation_dir
out_file = args.out_file
normalisation = args.normalisation
layer = args.layer
reference = args.reference
method = args.method
level = args.level
n_comps = args.n_comps
min_counts = args.min_counts
min_features = args.min_features
max_counts = args.max_counts
max_features = args.max_features
min_cells = args.min_cells
max_n_cells = args.max_n_cells
num_samples = args.num_samples
mixture_k = args.mixture_k
raw_corrected_counts = args.raw_corrected_counts
xenium_count_correction_dir = args.xenium_count_correction_dir
results_dir = args.results_dir
correction_method = args.correction_method

# variables
segmentation = panel.parents[1].stem
condition = panel.parents[0].stem

OBSM_KEY = "X_pca"
CT_KEY = (reference, method, level)
BATCH_KEY = "batch_key"
annotation_normalisation = "lognorm"  # fix this for now, even for sctransfrom
exclude_cell_type_containing = "malignant"

print(correction_method)
# read xenium samples
ads = {}
if raw_corrected_counts:
    for donor in (donors := panel.iterdir()):
        if "mixture_k" in donor.name or "lognorm" in donor.name or not donor.is_dir():
            continue
        for sample in (samples := donor.iterdir()):
            k = (segmentation, condition, panel.stem, donor.stem, sample.stem)
            name = "/".join(k)

            if correction_method == "split_fully_purified":
                name_corrected = f"{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/split_fully_purified/"
                sample_corrected_counts_path = xenium_count_correction_dir / f"{name_corrected}/corrected_counts.h5"

            else:
                if correction_method in ["resolvi", "resolvi_panel_use_batch=True", "resolvi_panel_use_batch=False"]:
                    name_corrected = f"{name}/{mixture_k=}/{num_samples=}/"
                elif correction_method in [
                    "resolvi_supervised",
                    "resolvi_panel_supervised_use_batch=True",
                    "resolvi_panel_supervised_use_batch=False",
                ]:
                    name_corrected = f"{name}/{normalisation}/reference_based/{reference}/{method}/{level}/{mixture_k=}/{num_samples=}"
                elif "ovrlpy" in correction_method:
                    name_corrected = f"{name}"

                sample_corrected_counts_path = results_dir / f"{correction_method}/{name_corrected}/corrected_counts.h5"

            ads[k] = sc.read_10x_h5(sample_corrected_counts_path)

            # read cell type annotation
            sample_annotation_dir = cell_type_annotation_dir / f"{name}/{annotation_normalisation}/reference_based"
            annot_file = sample_annotation_dir / f"{reference}/{method}/{level}/single_cell/labels.parquet"
            ads[k].obs[CT_KEY] = pd.read_parquet(annot_file).set_index("cell_id").iloc[:, 0]
    is_raw = True

else:
    for donor in (donors := panel.iterdir()):
        if "mixture_k" in donor.name or "lognorm" in donor.name or not donor.is_dir():
            continue
        for sample in (samples := donor.iterdir()):
            k = (
                segmentation,
                condition,
                panel.stem,
                donor.stem,
                sample.stem,
            )
            name = "/".join(k)

            sample_counts_path = sample / f"{normalisation}/normalised_counts/{layer}.parquet"
            sample_idx_path = sample / f"{normalisation}/normalised_counts/cells.parquet"

            ads[k] = sc.AnnData(pd.read_parquet(sample_counts_path))
            if layer != "scale_data":  # no need to sparsify scale_data which is dense
                ads[k].X = scipy.sparse.csr_matrix(ads[k].X)
            ads[k].obs_names = pd.read_parquet(sample_idx_path).iloc[:, 0]

            # read cell type annotation
            sample_annotation_dir = cell_type_annotation_dir / f"{name}/{annotation_normalisation}/reference_based"
            annot_file = sample_annotation_dir / f"{reference}/{method}/{level}/single_cell/labels.parquet"
            ads[k].obs[CT_KEY] = pd.read_parquet(annot_file).set_index("cell_id").iloc[:, 0]


# concatenate
xenium_levels = ["segmentation", "condition", "panel", "donor", "sample"]
for k in ads.keys():
    for i, lvl in enumerate(xenium_levels):
        ads[k].obs[lvl] = k[i]
ad_merge = sc.concat(ads)
ad_merge.obs[BATCH_KEY] = ad_merge.obs[xenium_levels].agg("_".join, axis=1)

# remove NaN  and exclude_cell_type_containing annotations
ad_merge = ad_merge[ad_merge.obs[CT_KEY].notna()]
ad_merge = ad_merge[~ad_merge.obs[CT_KEY].str.contains(exclude_cell_type_containing)].copy()

# subsample to reasonable size
if len(ad_merge) > max_n_cells:
    sc.pp.subsample(ad_merge, n_obs=max_n_cells)

# compute pca
sc.tl.pca(ad_merge, n_comps=n_comps)

# set up metrics
batchcor = BatchCorrection(
    silhouette_batch=True,
    ilisi_knn=True,
    kbet_per_label=True,
    graph_connectivity=True,
    pcr_comparison=True,
)

biocons = BioConservation(
    isolated_labels=True,
    nmi_ari_cluster_labels_leiden=True,
    nmi_ari_cluster_labels_kmeans=True,
    silhouette_label=True,
    clisi_knn=True,
)

# benchmark
bm = Benchmarker(
    ad_merge,
    batch_key=BATCH_KEY,
    label_key=CT_KEY,
    embedding_obsm_keys=[OBSM_KEY],
    pre_integrated_embedding_obsm_key=OBSM_KEY,
    bio_conservation_metrics=biocons,
    batch_correction_metrics=batchcor,
    n_jobs=-1,
)
bm.benchmark()

# save results
bm.get_results(min_max_scale=False).iloc[[0]].to_parquet(out_file)
