import dask

dask.config.set({"dataframe.query-planning": False})

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import argparse
import os
import sys

sys.path.append("workflow/scripts/")
import preprocessing
import readwrite

torch.set_float32_matmul_precision("medium")


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Run RESOLVI on a Xenium sample.")
    parser.add_argument("--panel", type=Path, help="Path to the panel file.")
    parser.add_argument("--xenium_processed_data_dir", type=Path, help="Path to the xenium processed data directories")
    parser.add_argument("--cell_type_annotation_dir", type=Path, help="Path to the cell_type_annotation_dir.")
    parser.add_argument("--annotation_mode", type=str, help="annotation mode")
    parser.add_argument("--normalisation", type=str, help="annotation normalisation method")
    parser.add_argument("--reference", type=str, help="annotation reference")
    parser.add_argument("--method", type=str, help="annotation method")
    parser.add_argument("--level", type=str, help="annotation level")
    parser.add_argument("--out_dir_resolvi_model", type=str, help="output directory with RESOLVI model weights")
    parser.add_argument("--min_counts", type=int, help="QC parameter from pipeline config")
    parser.add_argument("--min_features", type=int, help="QC parameter from pipeline config")
    parser.add_argument("--max_counts", type=float, help="QC parameter from pipeline config")
    parser.add_argument("--max_features", type=float, help="QC parameter from pipeline config")
    parser.add_argument("--min_cells", type=int, help="QC parameter from pipeline config")
    parser.add_argument(
        "--max_epochs",
        type=int,
        default=100,
        help="Maximum number of epochs to train the model.",
    )
    parser.add_argument("--mixture_k", type=int, help="mixture_k parameter for unsupervised RESOLVI")
    parser.add_argument("--use_batch", action="store_true", help="whether to use batch parameter for RESOLVI")
    parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )

    ret = parser.parse_args()

    return ret


if __name__ == "__main__":
    args = parse_args()

    print(args)

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    # process params
    if args.level is not None:
        print("Reading samples and cell type annotation")
        labels_key = "labels_key"
        semisupervised = True
    else:
        print("Reading samples")
        labels_key = None
        semisupervised = False

    if args.use_batch:
        batch_key = "sample"
        print(f"Using {batch_key=} for RESOLVI")
    else:
        batch_key = None

    condition = args.panel.parents[0].stem
    segmentation = args.panel.parents[1].stem

    # load data
    ads = {}
    for donor in (donors := args.panel.iterdir()):
        for sample in (samples_ := donor.iterdir()):
            print(donor.stem, sample.stem)

            if segmentation == "proseg_expected":
                k = ("proseg", condition, args.panel.stem, donor.stem, sample.stem)
                k_annot = (segmentation, condition, args.panel.stem, donor.stem, sample.stem)
                name = "/".join(k)
                name_annot = "/".join(k_annot)
                sample_dir = args.xenium_processed_data_dir / f"{name}/raw_results"
            else:
                k = (
                    segmentation.replace("proseg_mode", "proseg"),
                    condition,
                    args.panel.stem,
                    donor.stem,
                    sample.stem,
                )
                name = name_annot = "/".join(k)
                sample_dir = args.xenium_processed_data_dir / f"{name}/normalised_results/outs"

            ads[k] = readwrite.read_xenium_sample(sample_dir, anndata=True)
            ads[k].obs_names = ads[k].obs_names.astype(str)

            if segmentation == "proseg_expected":
                ads[k].obs_names = "proseg-" + ads[k].obs_names
                # need to round proseg expected counts for resolVI to run
                ads[k].X.data = ads[k].X.data.round()

            if args.level is not None:
                annot_file = (
                    args.cell_type_annotation_dir
                    / name_annot
                    / f"{args.normalisation}/{args.annotation_mode}/{args.reference}/{args.method}/{args.level}/single_cell/labels.parquet"
                )
                ads[k].obs[labels_key] = pd.read_parquet(annot_file).set_index("cell_id").iloc[:, 0]
                ads[k] = ads[k][ads[k].obs[labels_key].notna()].copy()

    print("Concatenating")
    # concatenate
    xenium_levels = ["segmentation", "condition", "panel", "donor", "sample"]
    for k in ads.keys():
        for i, lvl in enumerate(xenium_levels):
            ads[k].obs[lvl] = k[i]
    adata = sc.concat(ads)
    print("Done")

    adata.X.data = adata.X.data.astype(np.float32)

    # preprocess (QC filters only)
    # resolvi requires at least 5 counts in each cell
    preprocessing.preprocess(
        adata,
        normalize=False,
        log1p=False,
        scale="none",
        pca=False,
        umap=False,
        save_raw=False,
        min_counts=args.min_counts,
        min_genes=args.min_features,
        max_counts=args.max_counts,
        max_genes=args.max_features,
        min_cells=args.min_cells,
        backend="cpu",
    )

    scvi.external.RESOLVI.setup_anndata(
        adata, batch_key=batch_key, labels_key=labels_key, layer=None, prepare_data_kwargs={"spatial_rep": "spatial"}
    )
    resolvi = scvi.external.RESOLVI(adata, mixture_k=args.mixture_k, semisupervised=semisupervised)
    resolvi.train(max_epochs=args.max_epochs)
    resolvi.save(args.out_dir_resolvi_model, overwrite=True)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
