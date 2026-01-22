import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
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
parser = argparse.ArgumentParser(description="Plot panel of Xenium donors.")
parser.add_argument("--embed_file", type=str, help="Path to the embedding file.")
parser.add_argument("--reference", type=Path, help="annotation reference")
parser.add_argument("--color", type=str, help="annotation color")
parser.add_argument("--out_file", type=str, help="Path to the output file.")
parser.add_argument("--cell_type_palette", type=Path, help="Path to palette csv file")
parser.add_argument("--panel_palette", type=Path, help="Path to palette csv file")
parser.add_argument("--sample_palette", type=Path, help="Path to palette csv file")
parser.add_argument("--s", type=float, help="scatter point size")
parser.add_argument("--alpha", type=float, help="scatter alpha (transparency)")
parser.add_argument("--dpi", type=int, help="dpi of saved plot")
parser.add_argument(
    "--points_only",
    action="store_true",
    help="Remove axes, legend, title and only plot points",
)

args = parser.parse_args()

# Access the arguments
embed_file = args.embed_file
reference = args.reference
color = args.color
out_file = args.out_file
cell_type_palette = args.cell_type_palette
panel_palette = args.panel_palette
sample_palette = args.sample_palette
s = args.s
alpha = args.alpha
dpi = args.dpi
points_only = args.points_only

if color == "sample":
    palette = pd.read_csv(sample_palette, index_col=0).iloc[:, 0]
elif color == "panel":
    palette = pd.read_csv(panel_palette, index_col=0).iloc[:, 0]
else:
    if color == "Level2.1":
        palette_lvl2 = (
            pd.read_csv(cell_type_palette)[["Level2", "cols_Level2"]].drop_duplicates().set_index("Level2").squeeze()
        )
        palette = pd.read_csv(cell_type_palette)[[color, f"cols_{color}"]].drop_duplicates().set_index(color).squeeze()
        for k, v in palette_lvl2.items():
            if k not in palette.index:
                palette[k] = palette_lvl2[k]
    else:
        palette = pd.read_csv(cell_type_palette)[[color, f"cols_{color}"]].drop_duplicates().set_index(color).squeeze()

# load umap
df = pd.read_parquet(embed_file)
df["cell_id"] = df.index
annot = pd.read_parquet(reference / "metadata.parquet").set_index("cell_id")[color]
df = df.join(annot).dropna()

if color == "Level2.1":
    if reference.stem in ["external_lung", "matched_combo_standard_breast_specific", "matched_lung_standard"]:
        name_malignant = "malignant cell of lung"
    elif reference.stem in ["external_breast", "matched_combo_standard_lung_specific"]:
        name_malignant = "malignant cell of breast"
    else:
        name_malignant = "malignant cell"

    ct_to_replace = df[color][df[color].str.contains("malignant cell")].unique()
    replace_map = dict([[ct, name_malignant] for ct in ct_to_replace])
    df[color] = df[color].replace(replace_map)

# plotting color, palette
unique_labels = np.unique(df[color].dropna())
palette = {u: palette[u] for u in unique_labels}
legend_handles = [mpatches.Patch(color=color, label=label) for label, color in palette.items()]


# plot
figsize = (10, 10) if points_only else (12, 10)
f = plt.figure(figsize=figsize)
ax = plt.subplot()

sns.scatterplot(
    data=df,
    x="UMAP1",
    y="UMAP2",
    s=s,
    alpha=alpha,
    hue=color,
    ax=ax,
    palette=palette,
    legend=False,
)

if not points_only:
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    sns.despine()

    f.legend(
        handles=legend_handles,
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        title=color if isinstance(color, str) else ", ".join(color),
        frameon=False,
    )
    plt.tight_layout(rect=[0, 0, 0.85, 0.95])
else:
    ax.axis("off")

df.to_csv(Path(out_file).with_suffix(".csv"))
plt.savefig(out_file, dpi=dpi, bbox_inches="tight")
plt.savefig(Path(out_file).with_suffix(".png"), dpi=dpi, bbox_inches="tight")
plt.close()
