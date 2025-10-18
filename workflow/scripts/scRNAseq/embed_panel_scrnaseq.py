import dask

dask.config.set({"dataframe.query-planning": False})

from pathlib import Path
import argparse
import pandas as pd
import scanpy as sc
import sys

sys.path.append("workflow/scripts/")
import preprocessing

# Set up argument parser
parser = argparse.ArgumentParser(description="Embed panel of Xenium donors.")
parser.add_argument("--reference", type=Path, help="Path to the reference folder.")
parser.add_argument("--out_file", type=str, help="Path to the output file.")
parser.add_argument("--layer", type=str, help="Name of data layer to load.")
parser.add_argument("--n_comps", type=int, help="Number of components.")
parser.add_argument("--n_neighbors", type=int, help="Number of neighbors.")
parser.add_argument("--metric", type=str, help="Distance metric to use.")
parser.add_argument("--min_dist", type=float, help="Minimum distance parameter.")
parser.add_argument("--min_counts", type=int, help="QC parameter from pipeline config")
parser.add_argument("--min_features", type=int, help="QC parameter from pipeline config")
parser.add_argument("--max_counts", type=float, help="QC parameter from pipeline config")
parser.add_argument("--max_features", type=float, help="QC parameter from pipeline config")
parser.add_argument("--min_cells", type=int, help="QC parameter from pipeline config")
parser.add_argument("--genes", type=str, nargs="*", default=[], help="Restrict data to these genes for the UMAP.")
parser.add_argument("--samples", type=str, nargs="*", default=[], help="Restrict data to these samples for the UMAP.")

args = parser.parse_args()

# Access the arguments
reference = args.reference
out_file = args.out_file
layer = args.layer
n_comps = args.n_comps
n_neighbors = args.n_neighbors
metric = args.metric
min_dist = args.min_dist
min_counts = args.min_counts
min_features = args.min_features
max_counts = args.max_counts
max_features = args.max_features
min_cells = args.min_cells
genes = args.genes
samples = args.samples

print("Reading samples")
ad_merge = sc.read_10x_h5(reference / f"{layer}.h5")
ad_merge.obs = pd.read_parquet(reference / "metadata.parquet").set_index("cell_id")

# subset to samples
if len(samples):
    ad_merge = ad_merge[ad_merge.obs["donor"].isin(samples), :].copy()

# subset to genes
if len(genes):
    print("Subsetting")

    genes_found = [
        g
        for g in ad_merge.var_names
        if (g in genes) or (g.replace(".", "-") in genes)  # possible seurat renaming
    ]

    print(f"Found {len(genes_found)} out of {len(genes)} genes.")

    # read raw counts to reapply QC
    ad_merge_raw_counts = sc.read_10x_h5(reference / "RNA_counts.h5")
    ad_merge_raw_counts = ad_merge[:, genes_found].copy()

    # reapply QC to subset of genes
    preprocessing.preprocess(
        ad_merge_raw_counts,
        min_counts=min_counts,
        min_genes=min_features,
        max_counts=max_counts,
        max_genes=max_features,
        min_cells=min_cells,
        save_raw=False,
    )
    # subset
    ad_merge = ad_merge[ad_merge_raw_counts.obs_names, genes_found].copy()


print("Using", ad_merge.obs["donor"].nunique(), "samples and", ad_merge.n_vars, "genes")


if "counts" in layer:
    normalize = True
    log1p = True
else:
    normalize = False
    log1p = False

print("Computing PCA and UMAP")
# preprocess
preprocessing.preprocess(
    ad_merge,
    normalize=normalize,
    log1p=log1p,
    scale="none",
    n_comps=n_comps,
    metric=metric,
    min_dist=min_dist,
    n_neighbors=n_neighbors,
    pca=True,
    umap=True,
    save_raw=False,
    min_counts=None,
    min_genes=None,
    max_counts=None,
    max_genes=None,
    min_cells=None,
)

# save
df_umap = pd.DataFrame(ad_merge.obsm["X_umap"], index=ad_merge.obs_names, columns=["UMAP1", "UMAP2"])
df_umap.to_parquet(out_file)
