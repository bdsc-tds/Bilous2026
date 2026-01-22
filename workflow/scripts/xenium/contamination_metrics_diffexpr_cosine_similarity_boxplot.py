# %%
import dask

dask.config.set({"dataframe.query-planning": False})

import argparse
import numpy as np
import pandas as pd
import sys
import matplotlib.patches as mpatches
import seaborn as sns
import scanpy as sc
import matplotlib.pyplot as plt
import gc
from pathlib import Path

sys.path.append("workflow/scripts/")
import _utils
import readwrite
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

cfg = readwrite.config()

sns.set_style("ticks")

# Set up argument parser
parser = argparse.ArgumentParser(description="")
parser.add_argument("--condition", type=str, help="condition name.")
parser.add_argument("--panel", type=str, help="panel name.")
parser.add_argument("--correction_methods", type=str, nargs="*", default=[], help="correction methods list.")
parser.add_argument("--results_dir", type=Path, help="Path to the results dir.")
parser.add_argument("--xenium_processed_data_dir", type=Path, help="Path to the xenium dir.")
parser.add_argument("--xenium_count_correction_dir", type=Path, help="Path to the count_correction dir.")
parser.add_argument("--scrnaseq_processed_data_dir", type=Path, help="Path to the scrnaseq_processed_data dir.")
parser.add_argument("--seurat_to_h5_dir", type=Path, help="Path to the seurat_to_h5 dir.")
parser.add_argument("--std_seurat_analysis_dir", type=Path, help="Path to the std seurat analysis dir.")
parser.add_argument("--cell_type_annotation_dir", type=Path, help="Path to the cell type annotation dir.")
parser.add_argument("--out_dir", type=Path, help="Path to the output dir.")
parser.add_argument("--normalisation", type=str, help="Normalisation method")
parser.add_argument("--layer", type=str, help="Name of saved layer of the seurat object, data or scale_data")
parser.add_argument("--reference", type=str, help="annotation reference")
parser.add_argument("--method", type=str, help="annotation method")
parser.add_argument("--level", type=str, help="annotation level")
parser.add_argument("--mixture_k", type=int, help="ResolVI parameter")
parser.add_argument("--num_samples", type=int, help="ResolVI parameter")
parser.add_argument("--use_precomputed", action="store_true", help="Use precomputed data.")
parser.add_argument("--count_correction_palette", type=Path, help="Path to the count correction palette file.")
parser.add_argument("--dpi", type=int, help="Figure DPI.")
parser.add_argument("--extension", type=str, help="conservation metric to plot")
args = parser.parse_args()

# Access the arguments
condition = args.condition
panel = args.panel
correction_methods = args.correction_methods
results_dir = args.results_dir
xenium_processed_data_dir = args.xenium_processed_data_dir
xenium_count_correction_dir = args.xenium_count_correction_dir
scrnaseq_processed_data_dir = args.scrnaseq_processed_data_dir
cell_type_annotation_dir = args.cell_type_annotation_dir
seurat_to_h5_dir = args.seurat_to_h5_dir
std_seurat_analysis_dir = args.std_seurat_analysis_dir
out_dir = args.out_dir
normalisation = args.normalisation
layer = args.layer
reference = args.reference
method = args.method
level = args.level
mixture_k = args.mixture_k
num_samples = args.num_samples
count_correction_palette = args.count_correction_palette
use_precomputed = args.use_precomputed
dpi = args.dpi
extension = args.extension


# Params
layer_scrnaseq = "RNA_counts"
xenium_levels = ["segmentation", "condition", "panel", "donor", "sample"]
order = ["breast", "chuvio", "lung", "5k"]
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


plot_metric = "cosine_similarity"

CONDITIONS_REFS = {
    "breast": "matched_combo_standard_breast_specific.rds",
    "melanoma": "external_melanoma.rds",
    "NSCLC": "matched_combo_standard_lung_specific.rds",
    "mesothelioma_pilot": "matched_combo_standard_lung_specific.rds",
}


# %% [markdown]
# # Load counts and corrected counts

xenium_paths = {}
xenium_annot_paths = {}

for correction_method in correction_methods:
    print(f"Processing {correction_method} for condition {condition} and panel {panel}")
    xenium_paths[correction_method] = {}
    xenium_annot_paths[correction_method] = {}

    for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
        if segmentation.stem in ["proseg_mode", "bats_expected", "bats_normalised"]:
            continue
        for condition_dir in (conditions := segmentation.iterdir()):
            if condition_dir.stem != condition:
                continue
            for panel_dir in (panels := condition_dir.iterdir()):
                if panel_dir.stem != panel:
                    continue
                for donor in (donors := panel_dir.iterdir()):
                    for sample in (samples := donor.iterdir()):
                        k = (segmentation.stem, condition_dir.stem, panel_dir.stem, donor.stem, sample.stem)
                        name = "/".join(k)

                        # raw samples
                        if "proseg" in segmentation.stem:
                            k_proseg = ("proseg", condition_dir.stem, panel_dir.stem, donor.stem, sample.stem)
                            name_proseg = "/".join(k_proseg)
                            sample_dir = xenium_processed_data_dir / f"{name_proseg}/raw_results"
                        else:
                            sample_dir = xenium_processed_data_dir / f"{name}/normalised_results/outs"

                        sample_annotation = (
                            cell_type_annotation_dir
                            / f"{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/labels.parquet"
                        )

                        if correction_method == "raw":
                            xenium_paths[correction_method][k] = sample_dir
                            xenium_annot_paths[correction_method][k] = sample_annotation

                        # corrected samples
                        else:
                            if correction_method == "split_fully_purified":
                                name_corrected = f"{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/split_fully_purified/"
                                sample_corrected_counts_path = (
                                    xenium_count_correction_dir / f"{name_corrected}/corrected_counts.h5"
                                )

                            else:
                                if correction_method in [
                                    "resolvi",
                                    "resolvi_panel_use_batch=True",
                                    "resolvi_panel_use_batch=False",
                                ]:
                                    name_corrected = f"{name}/{mixture_k=}/{num_samples=}/"
                                elif correction_method in [
                                    "resolvi_supervised",
                                    "resolvi_panel_supervised_use_batch=True",
                                    "resolvi_panel_supervised_use_batch=False",
                                ]:
                                    name_corrected = f"{name}/{normalisation}/reference_based/{reference}/{method}/{level}/{mixture_k=}/{num_samples=}"
                                elif "ovrlpy" in correction_method:
                                    name_corrected = f"{name}"

                                sample_corrected_counts_path = (
                                    results_dir / f"{correction_method}/{name_corrected}/corrected_counts.h5"
                                )
                            # sample_normalised_counts = (
                            #     std_seurat_analysis_dir / f"{name}/{normalisation}/normalised_counts/{layer}.parquet"
                            # )
                            # sample_idx = (
                            #     std_seurat_analysis_dir / f"{name}/{normalisation}/normalised_counts/cells.parquet"
                            # )

                            xenium_paths[correction_method][k] = sample_corrected_counts_path


ads = readwrite.read_count_correction_samples(xenium_paths, correction_methods[1:])
ads["raw"] = readwrite.read_xenium_samples(xenium_paths["raw"], anndata=True, transcripts=False, max_workers=6)

# fix obs names for proseg expected, load cell types, log normalize data
# filter out cells without labels (this will apply QC thresholds as well since annotation is done after QC)
pbs_xenium = {}
for correction_method in correction_methods:
    pbs_xenium[correction_method] = {}

    for k, ad in ads[correction_method].items():
        if ad is not None:
            if correction_method == "raw":
                if k[0] == "proseg_expected":
                    ad.obs_names = ad.obs_names.astype(str)
                    ad.obs_names = "proseg-" + ad.obs_names

                # filter cells and read labels for raw
                ad.obs[level] = pd.read_parquet(xenium_annot_paths["raw"][k]).set_index("cell_id").iloc[:, 0]

                ad = ad[ad.obs[level].notna()]
                if level == "Level2.1":
                    # for custom Level2.1, simplify subtypes
                    ad.obs.loc[ad.obs[level].str.contains("malignant"), level] = "malignant cell"
                    ad.obs.loc[ad.obs[level].str.contains("T cell"), level] = "T cell"

                # remove tissue from cell type name
                ad.obs[level] = ad.obs[level].str.replace(r" of .+", "", regex=True)

                ads["raw"][k] = ad

            # filter cells and add labels from raw
            if correction_method != "raw":
                ad.obs[level] = ads["raw"][k].obs[level]
                ad = ad[[c for c in ads["raw"][k].obs_names if c in ad.obs_names]]
                ads[correction_method][k] = ad

            pbs_xenium[correction_method][k] = _utils.pseudobulk(ads[correction_method][k], level)
            sc.pp.normalize_total(pbs_xenium[correction_method][k])
            sc.pp.log1p(pbs_xenium[correction_method][k])

# %% load scRNAseq data
ref = seurat_to_h5_dir / CONDITIONS_REFS[condition]
ref_name = ref.stem
ref_dir = seurat_to_h5_dir / ref_name

print("loading scRNAseq data", ref_name)

ad = sc.read_10x_h5(ref_dir / f"{layer_scrnaseq}.h5")
ad.obs[level] = pd.read_parquet(ref_dir / "metadata.parquet").set_index("cell_id")[level]
ad = ad[ad.obs[level].notna()]

if level == "Level2.1" and level in ad.obs.columns:
    # for custom Level2.1, simplify subtypes
    ad.obs.loc[ad.obs[level].str.contains("malignant"), level] = "malignant cell"
    ad.obs.loc[ad.obs[level].str.contains("T cell"), level] = "T cell"

# remove tissue from cell type name
ad.obs[level] = ad.obs[level].str.replace(r" of .+", "", regex=True)

# Prepare pseudo bulk data
pb_scrna = _utils.pseudobulk(ad, level)
sc.pp.normalize_total(pb_scrna)
sc.pp.log1p(pb_scrna)

del ad
gc.collect()

# %% [markdown]
# # Plot decontamination results diffexpr

# %%
df_count_correction_palette = pd.read_csv(count_correction_palette, index_col=0).iloc[:, 0]

# get cell identity score df
df_all = _utils.get_cosine_similarity_score(
    pbs_xenium,
    pb_scrna,
    level,
    correction_methods,
    columns=xenium_levels + ["correction_method", "cti", "cosine_similarity"],
)

# rename segmentations
_utils.rename_methods(df_all)


for cti in df_all["cti"].unique():
    df = df_all.query(f"panel == '{panel}' and cti == '{cti}'")
    cti_name = cti.replace(" ", "_")

    out_file = out_dir / f"{panel}_{cti_name}_{plot_metric}.{extension}"
    print(cti, out_file)

    # plotting params, palette
    title = f"Cell type identity score for {cti=}"
    unique_labels = [c for c in hue_correction_order if c in np.unique(df[hue_correction].dropna())]
    unique_labels = unique_labels + [c for c in np.unique(df[hue_correction].dropna()) if c not in unique_labels]
    palette = {u: df_count_correction_palette[u] for u in unique_labels}
    legend_handles = [mpatches.Patch(color=color, label=label) for label, color in palette.items()]

    ### hypergeometric pvalue boxplot
    f = plt.figure(figsize=(7, 6))
    ax = plt.subplot()
    g = sns.boxplot(
        df,
        x="segmentation",
        y=plot_metric,
        hue=hue_correction,
        hue_order=unique_labels,
        legend=False,
        palette=palette,
        ax=ax,
        order=[s for s in hue_segmentation_order if s in df["segmentation"].unique()],
        # flierprops={
        #     "marker": "o",
        #     "color": "black",
        #     "markersize": 1,
        #     "markerfacecolor": "w",
        # },
        boxprops={"alpha": 0.4},
        showfliers=False,
    )
    sns.stripplot(
        df,
        x="segmentation",
        y=plot_metric,
        hue=hue_correction,
        hue_order=unique_labels,
        legend=False,
        palette=palette,
        ax=ax,
        order=[s for s in hue_segmentation_order if s in df["segmentation"].unique()],
        dodge=True,
        jitter=True,
        s=3.5,
    )

    sns.despine()
    ax.yaxis.grid(True)
    ax.set_ylim(top=1)
    ax.set_yticks(ax.get_yticks().tolist() + [1] if 1 not in ax.get_yticks() else ax.get_yticks())
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=16)
    ax.set_ylabel(plot_metric, fontsize=14)

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
