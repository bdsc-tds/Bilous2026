# %%
import dask

dask.config.set({"dataframe.query-planning": False})

import argparse
import numpy as np
import pandas as pd
import sys
import matplotlib.patches as mpatches
import seaborn as sns
import matplotlib.pyplot as plt
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

sns.set_style("ticks")

# Set up argument parser
parser = argparse.ArgumentParser(description="")
parser.add_argument("--condition", type=str, help="condition name.")
parser.add_argument("--panel", type=str, help="panel name.")
parser.add_argument("--correction_methods", type=str, nargs="*", default=[], help="correction methods list.")
parser.add_argument("--results_dir", type=Path, help="Path to the results dir.")
parser.add_argument("--std_seurat_analysis_dir", type=Path, help="Path to the std seurat analysis dir.")
parser.add_argument("--cell_type_annotation_dir", type=Path, help="Path to the cell type annotation dir.")
parser.add_argument("--out_dir", type=Path, help="Path to the output dir.")
parser.add_argument("--normalisation", type=str, help="Normalisation method")
parser.add_argument("--layer", type=str, help="Name of saved layer of the seurat object, data or scale_data")
parser.add_argument("--reference", type=str, help="annotation reference")
parser.add_argument("--method", type=str, help="annotation method")
parser.add_argument("--level", type=str, help="annotation level")
parser.add_argument("--top_n", type=int, help="contamination diffexpr script parameter")
parser.add_argument("--mixture_k", type=int, help="ResolVI parameter")
parser.add_argument("--num_samples", type=int, help="ResolVI parameter")
parser.add_argument("--use_precomputed", action="store_true", help="Use precomputed data.")
parser.add_argument("--count_correction_palette", type=Path, help="Path to the count correction palette file.")
parser.add_argument("--radius", type=int, help="n° of neighbors to use to define the spatial graph")
parser.add_argument("--dpi", type=int, help="Figure DPI.")
parser.add_argument("--extension", type=str, help="conservation metric to plot")
args = parser.parse_args()

# Access the arguments
condition = args.condition
panel = args.panel
correction_methods = args.correction_methods
results_dir = args.results_dir
cell_type_annotation_dir = args.cell_type_annotation_dir
std_seurat_analysis_dir = args.std_seurat_analysis_dir
out_dir = args.out_dir
normalisation = args.normalisation
layer = args.layer
reference = args.reference
method = args.method
level = args.level
top_n = args.top_n
mixture_k = args.mixture_k
num_samples = args.num_samples
use_precomputed = args.use_precomputed
count_correction_palette = args.count_correction_palette
radius = args.radius
dpi = args.dpi
extension = args.extension


args = parser.parse_args()


# Params
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

rank_metrics = ["logfoldchanges", "-log10pvals_x_logfoldchanges", "-log10pvals_x_sign_logfoldchanges", "mean_zscore"]
plot_metrics = ["hypergeometric_pvalue", "NES", f"n_hits_{top_n=}", "mean_zscore_pvalue"]
labels_key = level

# %% [markdown]
# # Load results diffexpr

# %%
dfs = readwrite.read_contamination_metrics_results(
    results_dir,
    correction_methods,
    std_seurat_analysis_dir,
    reference,
    method,
    level,
    mixture_k,
    num_samples,
    normalisation,
    layer,
    radius=radius,
    top_n=top_n,
    markers_mode="diffexpr",
    evaluation="diffexpr",
    ref_condition=condition,
    ref_panel=panel,
)


# %% [markdown]
# # Plot decontamination results diffexpr

# %%

df_count_correction_palette = pd.read_csv(count_correction_palette, index_col=0).iloc[:, 0]

for rank_metric in rank_metrics:
    plot_metrics_ = plot_metrics[-1:] if rank_metric == "mean_zscore" else plot_metrics[:-1]
    for plot_metric in plot_metrics_:
        df = _utils.get_df_marker_rank_significance_plot(
            dfs["df_markers_rank_significance_diffexpr"],
            rank_metric=rank_metric,
            plot_metric=plot_metric,
            correction_methods=correction_methods,
            use_precomputed=use_precomputed,
        )

        # rename segmentations

        if plot_metric in ["hypergeometric_pvalue", "mean_zscore_pvalue"]:
            df["-log10pvalue"] = -np.log10(df[plot_metric].astype(float))
            plot_metric = "-log10pvalue"

        # df = df.query("condition == @condition and panel == @panel")

        for cti, ctj in df[["cti", "ctj"]].drop_duplicates().values:
            cti_name = cti.replace(" ", "_")
            ctj_name = ctj.replace(" ", "_")
            out_file = (
                out_dir / f"{panel}_{cti_name}_contaminated_by_{ctj_name}_{rank_metric}_{plot_metric}.{extension}"
            )

            df_plot = df.query("cti == @cti and ctj == @ctj")

            # plotting params, palette
            unique_labels = [c for c in hue_correction_order if c in np.unique(df_plot[hue_correction].dropna())]
            unique_labels = unique_labels + [
                c for c in np.unique(df_plot[hue_correction].dropna()) if c not in unique_labels
            ]
            palette = {u: df_count_correction_palette[u] for u in unique_labels}
            legend_handles = [mpatches.Patch(color=color, label=label) for label, color in palette.items()]

            ### hypergeometric pvalue boxplot
            f = plt.figure(figsize=(7, 6))
            ax = plt.subplot()
            g = sns.boxplot(
                df_plot,
                x="segmentation",
                y=plot_metric,
                hue=hue_correction,
                hue_order=unique_labels,
                legend=False,
                palette=palette,
                ax=ax,
                order=[s for s in hue_segmentation_order if s in df_plot["segmentation"].unique()],
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
                df_plot,
                x="segmentation",
                y=plot_metric,
                hue=hue_correction,
                hue_order=unique_labels,
                legend=False,
                palette=palette,
                ax=ax,
                order=[s for s in hue_segmentation_order if s in df_plot["segmentation"].unique()],
                dodge=True,
                jitter=True,
                s=3.5,
            )

            sns.despine()
            ax.yaxis.grid(True)
            ax.tick_params(axis="x", labelsize=14)
            ax.tick_params(axis="y", labelsize=16)
            if plot_metric == '"-log10pvalue"':
                ax.set.ylabel(r"$-\log_{10} \text{ p-value}$", fontsize=14)
            else:
                ax.set_ylabel(plot_metric, fontsize=14)

            # title = f"Condition: {condition}, Panel: {panel}, Reference: {reference}, Method: {method}, Level: {level} \n{cti} contaminated by {ctj}\n rank metric: {rank_metric}, plot metric: {plot_metric}"
            # plt.suptitle(title, y=1.05)
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
            plt.close()
