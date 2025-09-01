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
parser.add_argument("--max_n_cells", type=int, help="Max number of cells to use.")
parser.add_argument("--singlets", action="store_true", help="Use only singlets.")
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
max_n_cells = args.max_n_cells
singlets = args.singlets

# variables
segmentation = panel.parents[1].stem
condition = panel.parents[0].stem
OBSM_KEY = "X_pca"
CT_KEY = (reference, method, level)
BATCH_KEY = "batch_key"
annotation_normalisation = "lognorm"  # fix this for now, even for sctransfrom
exclude_cell_type_containing = "malignant"

# read xenium samples
ads = {}
for donor in (donors := panel.iterdir()):
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

        if singlets:
            # read spot class
            spot_class_file = (
                sample_annotation_dir / f"{reference}/{method}/{level}/single_cell/output/results_df.parquet"
            )

            ads[k].obs["spot_class"] = pd.read_parquet(spot_class_file, columns=["cell_id", "spot_class"]).set_index(
                "cell_id"
            )
            ads[k] = ads[k][ads[k].obs["spot_class"] == "singlet"]

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
