import dask

dask.config.set({"dataframe.query-planning": False})

import numpy as np
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import sys

sys.path.append("workflow/scripts/")
import _utils
import matplotlib as mpl

mpl.rcParams.update(
    {
        "pdf.fonttype": 42,  # embed TrueType fonts (keeps text as text)
        "ps.fonttype": 42,
        "svg.fonttype": "none",  # if exporting SVG
        "text.usetex": False,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "savefig.transparent": True,
    }
)

# Set up argument parser
parser = argparse.ArgumentParser(description="Embed panel of Xenium donors.")
parser.add_argument("--condition", type=str, help="condition name.")
parser.add_argument("--panel", type=str, help="panel name.")
parser.add_argument("--correction_methods", type=str, nargs="*", default=[], help="correction methods list.")
parser.add_argument("--std_seurat_analysis_dir", type=Path, help="Path to the std seurat analysis dir.")
parser.add_argument("--cell_type_annotation_dir", type=Path, help="Path to the cell type annotation dir.")
parser.add_argument("--scib_metrics_results_dir", type=Path, help="Path to the scib metrics dir.")
parser.add_argument("--out_file", type=str, help="Path to the output file.")
parser.add_argument("--normalisation", type=str, help="Normalisation method")
parser.add_argument("--layer", type=str, help="Name of saved layer of the seurat object, data or scale_data")
parser.add_argument("--reference", type=str, help="annotation reference")
parser.add_argument("--method", type=str, help="annotation method")
parser.add_argument("--level", type=str, help="annotation level")
parser.add_argument("--n_comps", type=int, help="Number of components.")
parser.add_argument("--max_n_cells", type=int, help="Max number of cells to use.")
parser.add_argument("--count_correction_palette", type=Path, help="Path to the count correction palette file.")
parser.add_argument("--dpi", type=int, help="Figure DPI.")
parser.add_argument("--score", type=str, help="conservation metric to plot")
args = parser.parse_args()

# Access the arguments
condition = args.condition
panel = args.panel
correction_methods = args.correction_methods
cell_type_annotation_dir = args.cell_type_annotation_dir
scib_metrics_results_dir = args.scib_metrics_results_dir
std_seurat_analysis_dir = args.std_seurat_analysis_dir
out_file = args.out_file
normalisation = args.normalisation
layer = args.layer
reference = args.reference
method = args.method
level = args.level
n_comps = args.n_comps
max_n_cells = args.max_n_cells
count_correction_palette = args.count_correction_palette
dpi = args.dpi
score = args.score

# variables
palette = pd.read_csv(count_correction_palette).set_index("correction_method").iloc[:, 0]
OBSM_KEY = "X_pca"
CT_KEY = (reference, method, level)
BATCH_KEY = "batch_key"
normalisation_annot = "lognorm"  # fix this for now, even for sctransfrom
exclude_cell_type_containing = "malignant"
xenium_levels = ["segmentation", "condition", "panel", "normalisation", "metric", "score"]
hue_segmentation = "segmentation"
hue_segmentation_order = [
    "MM 0µm",
    "MM",
    "MM 15µm",
    "0µm",
    "5µm",
    "15µm",
    "Baysor",
    "ProSeg",
    "ProSeg mode",
    "Segger",
]


hue_correction = "correction_method"
hue_correction_order = [
    "raw",
    "ResolVI",
    "ResolVI supervised",
    "ovrlpy 0.5",
    "ovrlpy 0.7",
    "SPLIT",
]


bio_score = "Bio conservation"
batch_score = "Batch correction"

biocons_metrics = [
    "cLISI",
    "Isolated labels",
    "KMeans NMI",
    "KMeans ARI",
    "Leiden NMI",
    "Leiden ARI",
    "Silhouette label",
    "Bio conservation",
]
batchcor_metrics = ["iLISI", "Graph connectivity", "KBET", "Silhouette batch", "Batch correction"]


# read results
df = {}
for correction_method in correction_methods:
    print(f"Reading {correction_method}")
    for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
        k = (segmentation.stem, condition, panel, normalisation)
        name = "/".join(k)
        scib_metrics_file = (
            scib_metrics_results_dir
            / f"{correction_method}/{name}/scib_metrics_{layer}_{reference}_{method}_{level}_{n_comps=}_{max_n_cells=}.parquet"
        )
        if scib_metrics_file.exists():
            df[correction_method, *k] = pd.read_parquet(scib_metrics_file).squeeze()
        else:
            print(f"File not found: {scib_metrics_file}")
df = pd.concat(df).reset_index()
df.columns = ["correction_method"] + xenium_levels
_utils.rename_methods(df)

# rename segmentations


u_segmentations = df["segmentation"].unique()
order = [h for h in hue_segmentation_order if h in u_segmentations]

# plotting params, palette
unique_labels = [c for c in hue_correction_order if c in np.unique(df[hue_correction].dropna())]
unique_labels = unique_labels + [c for c in np.unique(df[hue_correction].dropna()) if c not in unique_labels]
palette = {u: palette[u] for u in unique_labels}
legend_handles = [mpatches.Patch(color=color, label=label) for label, color in palette.items()]

# for score in biocons_metrics + batchcor_metrics:
df_score = df.query("metric == @score")

sns.set(style="ticks")
f = plt.figure(figsize=(7, 6))
ax = plt.subplot()
g = sns.barplot(
    df_score,
    x="segmentation",
    y="score",
    hue=hue_correction,
    hue_order=unique_labels,
    order=order,
    legend=False,
    palette=palette,
    ax=ax,
)

sns.despine()
ax.yaxis.grid(True)
ax.tick_params(axis="x", labelsize=14)
ax.tick_params(axis="y", labelsize=16)

# title = f"condition: {condition}, Panel: {panel}\n Reference: {reference}, Method: {method}, Level: {level}\n {score}"
# plt.suptitle(title)
# f.legend(
#     handles=legend_handles,
#     loc="center left",
#     bbox_to_anchor=(1, 0.5),
#     title=hue_correction,
#     frameon=False,
# )
# plt.tight_layout(rect=[0, 0, 1, 0.95])
df.to_csv(Path(out_file).with_suffix(".csv"))
plt.savefig(out_file, dpi=dpi, bbox_inches="tight")
# plt.show()
