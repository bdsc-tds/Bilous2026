import dask

dask.config.set({"dataframe.query-planning": False})

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scvi
import scanpy as sc
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
    parser = argparse.ArgumentParser(description="Embed panel of Xenium samples.")
    parser.add_argument("--panel", type=Path, help="Path to the xenium sample file.")
    parser.add_argument("--xenium_processed_data_dir", type=Path, help="Path to the xenium processed data directories")
    parser.add_argument("--dir_resolvi_model", type=str, help="directory with saved RESOLVI model weights")
    parser.add_argument("--results_dir", type=Path, help="results directory")
    parser.add_argument("--cell_type_annotation_dir", type=Path, help="Path to the cell_type_annotation_dir.")
    parser.add_argument("--annotation_mode", type=str, help="annotation mode")
    parser.add_argument("--normalisation", type=str, help="annotation normalisation method")
    parser.add_argument("--reference", type=str, help="annotation reference")
    parser.add_argument("--method", type=str, help="annotation method")
    parser.add_argument("--level", type=str, help="annotation level")
    parser.add_argument("--min_counts", type=int, help="QC parameter from pipeline config")
    parser.add_argument("--min_features", type=int, help="QC parameter from pipeline config")
    parser.add_argument("--max_counts", type=float, help="QC parameter from pipeline config")
    parser.add_argument("--max_features", type=float, help="QC parameter from pipeline config")
    parser.add_argument("--min_cells", type=int, help="QC parameter from pipeline config")
    parser.add_argument(
        "--num_samples",
        type=int,
        help="Number of samples for RESOLVI generative model.",
    )
    parser.add_argument("--batch_size", type=int, help="batch size parameter")
    parser.add_argument("--mixture_k", type=str, help="mixture_k parameter for RESOLVI")
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
        out_dir = args.results_dir / f"resolvi_panel_supervised_use_batch={args.use_batch}"
    else:
        print("Reading samples")
        labels_key = None
        semisupervised = False
        out_dir = args.results_dir / f"resolvi_panel_use_batch={args.use_batch}"

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

    resolvi = scvi.external.RESOLVI.load(args.dir_resolvi_model, adata=adata)

    samples_corr = resolvi.sample_posterior(
        model=resolvi.module.model_corrected,
        return_sites=["obs"],
        return_observed=True,
        summary_fun={"post_sample_q50": np.median},
        num_samples=args.num_samples,
        batch_size=args.batch_size,
        summary_frequency=100,
    )
    samples_corr = pd.DataFrame(samples_corr).T

    samples = resolvi.sample_posterior(
        model=resolvi.module.model_residuals,
        return_sites=["mixture_proportions"],
        summary_fun={"post_sample_means": np.mean},
        num_samples=args.num_samples,
        batch_size=args.batch_size,
        summary_frequency=100,
    )
    samples_proportions = pd.DataFrame(samples).T

    ### save
    samples_corr = pd.DataFrame(
        samples_corr.loc["post_sample_q50", "obs"],
        index=adata.obs_names,
        columns=adata.var_names,
    )

    proportions_cols = ["true_proportion", "diffusion_proportion", "background_proportion"]
    samples_proportions = pd.DataFrame(
        samples_proportions.loc["post_sample_means", "mixture_proportions"],
        index=adata.obs_names,
        columns=proportions_cols,
    )

    # add expression and metadata to new anndata
    adata_out = ad.AnnData(samples_corr)
    adata_out.obs = adata.obs
    for c in proportions_cols:
        adata_out.obs[c] = samples_proportions[c].values

    # save
    for _, (seg, cond, pan, don, samp) in adata_out.obs[xenium_levels].drop_duplicates().iterrows():
        print('saving',seg, cond, pan, don, samp)

        if semisupervised:
            k = (
                seg.replace("proseg",segmentation),
                cond,
                pan,
                don,
                samp,
                args.normalisation,
                args.annotation_mode,
                args.reference,
                args.method,
                args.level,
                f"mixture_k={args.mixture_k}",
                f"num_samples={args.num_samples}",
            )
        else:
            k = (
                seg.replace("proseg",segmentation),
                cond,
                pan,
                don,
                samp,
                f"mixture_k={args.mixture_k}",
                f"num_samples={args.num_samples}",
            )
        name = "/".join(k)

        out_file_resolvi_corrected_counts = out_dir / f"{name}/corrected_counts.h5"
        out_file_resolvi_proportions = out_dir / f"{name}/proportions.parquet"

        print("saving to:",out_file_resolvi_corrected_counts,'\n',out_file_resolvi_proportions)
        # subset sample
        adata_out_sub = adata_out[(adata_out.obs[xenium_levels] == [seg, cond, pan, don, samp]).all(axis=1)].copy()

        # write
        out_file_resolvi_corrected_counts.parent.mkdir(exist_ok=True, parents=True)
        readwrite.write_10X_h5(adata_out_sub, out_file_resolvi_corrected_counts)
        adata_out_sub.obs[proportions_cols].to_parquet(out_file_resolvi_proportions)

        del adata_out_sub

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr
