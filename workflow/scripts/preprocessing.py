import scanpy as sc
import scvi
import torch
import os
import pathlib
import pandas as pd
import numpy as np
import importlib
import sklearn
import scipy
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import warnings
# import scarches as sca


if importlib.util.find_spec("rapids_singlecell") is not None:
    # Import gpu libraries, Initialize rmm and cupy
    import rapids_singlecell as rsc
    import cupy as cp
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator

    rmm.disable_logging()


def check_gpu_availability():
    import torch

    if torch.cuda.is_available():
        return True
    else:
        return False


def get_latent_keys(
    params_modes=["scvi", "ot_bw", "scanorama", "harmony"],
    params_hvg_modes=["union", "strict_intersection", "intersection"],
    params_n_latent=[10, 30, 50],
    params_scale=["all", "batch", "none"],
    integration_key="dataset_merge_id",
):
    LATENT_KEYS = []
    LATENT_KEYS_SIMPLE = []
    for mode in params_modes:
        for hvg_mode in params_hvg_modes:
            for n_latent in params_n_latent:
                for scale in params_scale:
                    if mode in ["scvi", "trvae"] and scale in ["all", "batch"]:
                        continue
                    else:
                        LATENT_KEY = f"X_mode={mode}_hvg_mode={hvg_mode}_n_latent={n_latent}_scale={scale}"
                        LATENT_KEY_SIMPLE = f"{mode}_{hvg_mode}_{n_latent}_{scale}"

                        LATENT_KEYS.append(LATENT_KEY)
                        LATENT_KEYS_SIMPLE.append(LATENT_KEY_SIMPLE)
    return LATENT_KEYS, LATENT_KEYS_SIMPLE


def split_batches(adata, batch, hvg=None, return_categories=False):
    """Split batches and preserve category information (taken from https://github.com/theislab/scib)

    :param adata:
    :param batch: name of column in ``adata.obs``. The data type of the column must be of ``Category``.
    :param hvg: list of highly variable genes
    :param return_categories: whether to return the categories object of ``batch``
    """
    split = []
    batch_categories = adata.obs[batch].cat.categories
    if hvg is not None:
        adata = adata[:, hvg]
    for i in batch_categories:
        split.append(adata[adata.obs[batch] == i].copy())
    if return_categories:
        return split, batch_categories
    return split


def merge_adata(*adata_list, **kwargs):
    """Merge adatas from list while removing duplicated ``obs`` and ``var`` columns

    :param adata_list: ``anndata`` objects to be concatenated
    :param kwargs: arguments to be passed to ``anndata.AnnData.concatenate``
    """

    if len(adata_list) == 1:
        return adata_list[0]

    # Make sure that adatas do not contain duplicate columns
    for _adata in adata_list:
        for attr in ("obs", "var"):
            df = getattr(_adata, attr)
            dup_mask = df.columns.duplicated()
            if dup_mask.any():
                print(f"Deleting duplicated keys `{list(df.columns[dup_mask].unique())}` from `adata.{attr}`.")
                setattr(_adata, attr, df.loc[:, ~dup_mask])

    return sc.AnnData.concatenate(*adata_list, **kwargs)


def scale_batch(adata, batch):
    """Batch-aware scaling of count matrix (taken from https://github.com/theislab/scib)

    Scaling counts to a mean of 0 and standard deviation of 1 using ``scanpy.pp.scale`` for each batch separately.

    :param adata: ``anndata`` object with normalised and log-transformed counts
    :param batch: ``adata.obs`` column
    """

    # Store layers for after merge (avoids vstack error in merge)
    adata_copy = adata.copy()
    tmp = dict()
    for lay in list(adata_copy.layers):
        tmp[lay] = adata_copy.layers[lay]
        del adata_copy.layers[lay]

    split = split_batches(adata_copy, batch)

    for i in split:
        sc.pp.scale(i)

    adata_scaled = merge_adata(*split, batch_key=batch, index_unique=None)
    # Reorder to original obs_name ordering
    adata_scaled = adata_scaled[adata.obs_names]

    # Add layers again
    for key in tmp:
        adata_scaled.layers[key] = tmp[key]

    del tmp
    del adata_copy

    return adata_scaled


def preprocess(
    adata,
    batch_key="dataset_merge_id",
    normalize=False,  # Normalize total counts
    log1p=False,  # Log1p transform
    pca=False,  # Perform PCA
    scale="none",  # Scale data
    umap=False,  # Perform UMAP
    n_comps=30,  # Number of PCA components
    n_neighbors=15,  # Number of neighbors for kNN
    metric="cosine",  # Metric for kNN
    min_dist=0.5,
    backend="gpu",  # "gpu" or "cpu"
    device=0,  # Device ID for GPU backend
    save_raw=False,  # Whether to save raw data
    verbose=True,  # Optional verbose output
    filter_empty=True,  # Whether to filter empty batches
    min_counts=None,
    min_genes=None,
    max_counts=None,
    max_genes=None,
    min_cells=None,
):
    """
    Preprocesses an AnnData object.

    Args:
        adata: AnnData object to preprocess.
        normalize (bool): Whether to normalize total counts.
        log1p (bool): Whether to log1p transform.
        pca (bool): Whether to perform PCA.
        scale (str): Whether to scale data or not. Can be "none", "batch" or "all".
        umap (bool): Whether to perform UMAP.
        n_comps (int): Number of PCA components.
        n_neighbors (int): Number of neighbors for kNN.
        metric (str): Metric for kNN.
        backend (str): "gpu" or "cpu".
        device (int): Device ID for GPU backend.
        save_raw (bool): Whether to save raw data in layers['counts']
        verbose (bool): Whether to print verbose output.

    Returns:
        None
    """
    if "preprocess" in adata.uns:
        print("Warning: preprocess key already found in adata.uns")

    if save_raw:
        if verbose:
            print("Saving raw counts in layers['counts']...")
        adata.layers["counts"] = adata.X

    n_cells_raw = adata.shape[0]
    n_genes_raw = adata.shape[1]

    if min_cells is not None:
        sc.pp.filter_genes(adata, min_cells=min_cells)
    if min_genes is not None:
        sc.pp.filter_cells(adata, min_genes=min_genes)
    if min_counts is not None:
        sc.pp.filter_cells(adata, min_counts=min_counts)
    if max_counts is not None:
        sc.pp.filter_cells(adata, max_counts=max_counts)
    if max_genes is not None:
        sc.pp.filter_cells(adata, max_genes=max_genes)

    if verbose:
        print("Removed", n_cells_raw - adata.shape[0], " cells...")
        print("Removed", n_genes_raw - adata.shape[1], " genes...")

    ### optional switch to GPU backend ###
    if backend == "gpu":
        if not check_gpu_availability():
            print("GPU not available. Switching to CPU backend...")
            xsc = sc
            backend = "cpu"
        else:
            # allow memory oversubscription, transfer data to GPU
            rmm.reinitialize(managed_memory=True, devices=device)
            cp.cuda.set_allocator(rmm_cupy_allocator)
            if verbose:
                print("Transferring data to GPU...")
            rsc.get.anndata_to_GPU(adata)
            xsc = rsc
    else:
        xsc = sc

    ### preprocessing ###
    if normalize:
        if verbose:
            print("Normalizing total counts...")
        xsc.pp.normalize_total(adata, target_sum=1e4)
    if log1p:
        if verbose:
            print("Applying log1p transformation...")
        xsc.pp.log1p(adata)
    if scale == "all":
        if verbose:
            print("Scaling data...")
        adata = scale_batch(adata, batch_key)
    elif scale == "batch":
        if verbose:
            print("Batch-aware scaling of data...")
        xsc.pp.scale(adata)
    if pca:
        if verbose:
            print("Performing PCA...")
        xsc.tl.pca(adata, n_comps=n_comps)
    if umap:
        if verbose:
            print("Performing UMAP...")
        xsc.pp.neighbors(adata, n_pcs=n_comps, n_neighbors=n_neighbors, metric=metric)
        xsc.tl.umap(adata, min_dist=min_dist)

    # Transfer data back to CPU if using GPU backend
    if backend == "gpu":
        if verbose:
            print("Transferring data back to CPU...")
        rsc.get.anndata_to_CPU(adata)
    adata.uns["preprocess"] = dict(
        normalize=normalize,
        log1p=log1p,
        pca=pca,
        scale=scale,
        umap=umap,
        n_comps=n_comps,
        n_neighbors=n_neighbors,
        metric=metric,
        backend=backend,
        device=device,
    )


def prepare_adatas_hvg_split(ads, path=None, label="dataset_merge_id", overwrite=False):
    """
    Prepare data by finding HVGs for each AnnData object and
    creating a joint AnnData with union, intersection, and strict_intersection of HVGs.

    Args:
        ads (dict): Dictionary of AnnData objects.
        path (str): Path to save the processed data.
        label (str): Label column of dataset ids in the joint AnnData object.
    Returns:
        None
    """

    # Get common genes across all AnnData objects
    common_genes = list(set.intersection(*map(set, [a.var_names for a in ads.values()])))
    assert len(common_genes) > 0

    # Concatenate AnnData objects and save as ad_all
    ad_all = sc.concat(ads, label=label, join="outer")[:, common_genes]
    ad_all.obs["obs_names"] = ad_all.obs_names
    ad_all.obs_names_make_unique()
    ad_all.write(f"{path}_all.h5ad")

    # Loop over HVG modes
    for hvg_mode in ["union", "strict_intersection", "intersection"]:
        print(hvg_mode)

        # Check if file exists and skip if so
        if (pathlib.Path(f"{path}_{hvg_mode}.h5ad").is_file()) and (not overwrite):
            print("\tfound saved file")
            continue
        else:
            hvgs = []  # List to store HVGs for each AnnData object

            # Find HVGs for each AnnData object
            for k in ads.keys():
                ad = ads[k][:, common_genes]
                sc.pp.highly_variable_genes(ad, flavor="seurat_v3", n_top_genes=2000)
                hvgs.append(ad.var["highly_variable"][ad.var["highly_variable"]].index)

            # Create joint HVG and save
            if hvg_mode == "union":
                union_hvgs = list(set.union(*map(set, hvgs)))
                ad_joint = ad_all[:, union_hvgs].copy()
            elif hvg_mode == "strict_intersection":
                intersection_hvgs = list(set.intersection(*map(set, hvgs)))
                ad_joint = ad_all[:, intersection_hvgs].copy()
            elif hvg_mode == "intersection":
                sc.pp.highly_variable_genes(
                    ad_all,
                    flavor="seurat_v3",
                    batch_key="dataset_merge_id",
                    n_top_genes=2000,
                )
                ad_joint = ad_all[
                    :,
                    ad_all.var["highly_variable"][ad_all.var["highly_variable"]].index,
                ].copy()

            if path is None:
                return ad_joint
            else:
                # Save joint HVG
                ad_joint.write(f"{path}_{hvg_mode}.h5ad")
                del ad_joint


def prepare_adatas_hvg(ads, path=None, label="dataset_merge_id"):
    """
    Prepare data by finding HVGs for each AnnData object and
    creating a joint AnnData.

    Args:
        ads (dict): Dictionary of AnnData objects.
        path (str): Path to save the processed data.
        label (str): Label column of dataset ids in the joint AnnData object.
    Returns:
        None
    """

    # Find common genes across all AnnData objects
    common_genes = list(set.intersection(*map(set, [a.var_names for a in ads.values()])))
    assert len(common_genes) > 0

    # Concatenate AnnData objects and save as ad_all
    ad_all = sc.concat(ads, label=label, join="outer")[:, common_genes]
    ad_all.obs["obs_names"] = ad_all.obs_names
    ad_all.obs_names_make_unique()
    ad_all.layers["counts"] = ad_all.X

    # Loop over HVG modes
    for hvg_mode in ["union", "strict_intersection", "intersection"]:
        print(hvg_mode)

        # Find HVGs for each AnnData object
        hvgs_batch = []
        for k in ads.keys():
            ad = ads[k][:, common_genes]
            sc.pp.highly_variable_genes(ad, flavor="seurat_v3", n_top_genes=2000)
            hvgs_batch.append(ad.var["highly_variable"][ad.var["highly_variable"]].index)

        # Create joint HVG and save
        if hvg_mode == "union":
            hvgs = list(set.union(*map(set, hvgs_batch)))
        elif hvg_mode == "strict_intersection":
            hvgs = list(set.intersection(*map(set, hvgs_batch)))
        elif hvg_mode == "intersection":
            hvgs = (
                sc.pp.highly_variable_genes(
                    ad_all,
                    flavor="seurat_v3",
                    batch_key="dataset_merge_id",
                    n_top_genes=2000,
                    inplace=False,
                )
                .query("highly_variable==True")
                .index
            )

        # Add joint HVG to ad_all and save
        ad_all.var[f"highly_variable_{hvg_mode}"] = False
        ad_all.var.loc[hvgs, f"highly_variable_{hvg_mode}"] = True

    if path is None:
        return ad_all
    else:
        # Save joint HVG
        ad_all.write(path)
        del ad_all


def embed_adata_cellxgene_scvi(
    adata,
    batch_key="donor_name",
    cxg_scvi_retrain=False,
    n_latent=None,
    model_filename="/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/jbac/projects/data/cellxgene/models/2023-12-15-scvi-homo-sapiens/scvi.model",
):
    """
    Embeds adata using scvi and cellxgene retrained model.

    Args:
        adata : AnnData
            Annotated data matrix.
        cxg_scvi_retrain : bool, optional
            Whether to use the cellxgene retrained model, by default False.
        model_filename : str, optional
            File path to the scvi model
    Returns:
        np.ndarray
            Latent representation of the input data.
    """
    adata_ = sc.AnnData(adata.layers["counts"], obs=adata.obs, var=adata.var)
    adata_.var_names = adata.var_names
    adata_.obs_names = adata.obs_names
    adata_.var["ensembl_id"] = adata_.var.index
    adata_.obs["n_counts"] = adata_.X.sum(axis=1)
    adata_.obs["joinid"] = list(range(adata_.n_obs))
    adata_.obs["batch"] = adata_.obs[batch_key]

    scvi.model.SCVI.prepare_query_anndata(adata_, model_filename)
    if cxg_scvi_retrain:
        model = scvi.model.SCVI.load_query_data(adata_, model_filename, freeze_dropout=True)
        model.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))
    else:
        model = scvi.model.SCVI.load_query_data(adata_, model_filename)
        # This allows for a simple forward pass
        model.is_trained = True
    return model.get_latent_representation()


def ot_mapping(xs, xt, mode="ot_bw"):
    """
    Compute optimal transport mapping between source and target distributions.

    Args:
        xs: numpy array, source distribution
        xt: numpy array, target distribution
        mode: string, optional, mode of optimal transport, default is 'ot_bw'
    Returns:
        xst: numpy array, optimal transport mapping from source to target distribution
    """
    import ot

    if mode == "ot_bw":
        Ae, be = ot.gaussian.empirical_bures_wasserstein_mapping(xs, xt)
        xst = xs.dot(Ae) + be
    elif mode == "ot_gw":
        Ae, be = ot.gaussian.empirical_gaussian_gromov_wasserstein_mapping(xs, xt)
        xst = xs.dot(Ae) + be
    elif mode == "ot_emd":
        xst = ot.da.EMDTransport().fit_transform(Xs=xs, Xt=xt)
    return xst


def integrate(
    adata,
    integration_key,
    mode,
    ref=None,
    n_latent=10,
    compute_embeddings=["mde", "umap"],
    model_save_path=None,
):
    """Integrates the joint dataset according to integration_key using mode as integration method.

    Args:
        adata (anndata.AnnData):
            adata object with joint data
        integration_key (str):
            adata.obs column to use for integration
        mode (str): 'ot_bw','ot_gw','ot_emd','scvi', 'cxg_scvi', 'cxg_scvi_retrain', 'scanorama', 'harmony',
            Integration method to use
            obsm['X_pca'] is used as input to OT methods
            layers['counts'] is used as raw counts input to scvi models
        ref (str, optional):
            For OT based integration only, the dataset to use as the reference dataset. Defaults to None.
        compute_embeddings (bool, optional):
            Whether to compute UMAP and MDE on the integrated data. Defaults to True.

    Returns:
        new keys in adata.obs:
            LATENT_KEY = f"X_{mode}"
            MDE_KEY = f"{LATENT_KEY}_MDE"
            UMAP_KEY = f"{LATENT_KEY}_UMAP"
    """
    import scarches as sca

    if "counts" in adata.layers:
        adata.layers["X"] = adata.X
        adata.X = adata.layers["counts"]
    else:
        print('WARNING: adata.layers["counts"] does not exist. Using adata.X with assumption it contains raw counts')

    LATENT_KEY = f"X_{mode}"
    MDE_KEY = f"{LATENT_KEY}_mde"
    UMAP_KEY = f"{LATENT_KEY}_umap"

    if mode in ["ot_bw", "ot_gw", "ot_emd"]:
        if ref is None:
            ref = adata.obs[integration_key].value_counts().idxmax()
            print(
                "ref parameter is None. Using dataset with most cells as reference:",
                ref,
            )
        ref_idx = adata.obs[integration_key] == ref

        adata.obsm[LATENT_KEY] = adata.obsm["X_pca"].copy()
        for donor in adata.obs[integration_key].unique():
            if donor != ref:
                xs_idx = adata.obs[integration_key] == donor
                adata.obsm[LATENT_KEY][xs_idx] = ot_mapping(
                    adata.obsm[LATENT_KEY][xs_idx],
                    adata.obsm["X_pca"][ref_idx],
                    mode=mode,
                )

    elif mode == "scanorama":
        sc.external.pp.scanorama_integrate(adata, key=integration_key, basis="X_pca")
    elif mode == "harmony":
        sc.external.pp.harmony_integrate(adata, key=integration_key, basis="X_pca", adjusted_basis=LATENT_KEY)
    elif mode == "scvi":
        sca.models.SCVI.setup_anndata(adata, batch_key=integration_key)
        model = sca.models.SCVI(adata, n_hidden=128, n_layers=2, n_latent=n_latent, gene_likelihood="nb")
        model.train()
        if model_save_path is not None:
            model.save(model_save_path, overwrite=True)
        adata.obsm[LATENT_KEY] = model.get_latent_representation()

    elif mode == "trvae":
        adata.X = adata.X.astype(np.float32)
        model = sca.models.TRVAE(
            adata=adata,
            condition_key=integration_key,
            conditions=adata.obs[integration_key].unique().tolist(),
            latent_dim=n_latent,
            hidden_layer_sizes=[128, 128],
            recon_loss="nb",
        )
        model.train()
        if model_save_path is not None:
            model.save(model_save_path, overwrite=True)
        adata.obsm[LATENT_KEY] = model.get_latent()

    elif mode == "cxg_scvi":
        adata.obsm[LATENT_KEY] = embed_adata_cellxgene_scvi(adata, cxg_scvi_retrain=False)

    elif mode == "cxg_scvi_retrain":
        adata.obsm[LATENT_KEY] = embed_adata_cellxgene_scvi(adata, cxg_scvi_retrain=True)

    if "mde" in compute_embeddings:
        adata.obsm[MDE_KEY] = scvi.model.mde(adata.obsm[LATENT_KEY])
    if "umap" in compute_embeddings:
        rsc.get.anndata_to_GPU(adata)
        rsc.pp.neighbors(adata, metric="cosine", use_rep=LATENT_KEY)
        adata.obsm[UMAP_KEY] = rsc.tl.umap(adata, copy=True, min_dist=0.2).obsm["X_umap"]
        rsc.get.anndata_to_CPU(adata)

    if "counts" in adata.layers:
        adata.X = adata.layers["X"]
        del adata.layers["X"]


def preprocess_and_integrate(
    input_h5ad,
    output_dir=None,
    hvg_modes=["union", "strict_intersection", "intersection"],
    modes=["scvi", "cxg_scvi", "trvae", "ot_bw", "ot_emd", "scanorama", "harmony"],
    integration_key="dataset_merge_id",
    scale=False,
    n_latent=10,
    overwrite=False,
):
    # Read input
    name = pathlib.Path(input_h5ad).stem

    ad_all = sc.read(input_h5ad)

    # Create output dir
    if output_dir is not None:
        os.makedirs(os.path.join(output_dir, "models"), exist_ok=True)

    # Preprocess
    if "preprocess" not in ad_all.uns:
        preprocess(
            ad_all,
            normalize=True,
            log1p=True,
            scale=scale,
            pca=False,
            n_comps=n_latent,
            umap=False,
            save_raw=True,
            verbose=True,
        )

    obsm = {}
    # Loop over HVG modes
    for hvg_mode in hvg_modes:
        ad_joint = ad_all[:, ad_all.var[f"highly_variable_{hvg_mode}"]].copy()

        # Loop over integration modes
        for mode in modes:
            LATENT_KEY = f"X_mode={mode}_hvg_mode={hvg_mode}_n_latent={n_latent}_scale={scale}"
            print(LATENT_KEY)

            if not overwrite or (output_dir is not None and os.path.exists(f"{output_dir}/{name}_{LATENT_KEY}.csv")):
                print("\t\tAlready computed", mode, "skipping...")
                continue
            else:
                # Compute PCA on HVG
                if scale in ["all", "batch"]:
                    filter_empty = False
                else:
                    filter_empty = True
                preprocess(
                    ad_joint,
                    filter_empty=filter_empty,
                    normalize=False,
                    log1p=False,
                    scale=False,
                    pca=True,
                    n_comps=n_latent,
                    umap=False,
                    save_raw=False,
                    verbose=True,
                )

                integrate(
                    ad_joint,
                    integration_key=integration_key,
                    mode=mode,
                    ref=ad_joint.obs[integration_key][0],
                    n_latent=n_latent,
                    compute_embeddings=[],
                    model_save_path=f"{output_dir}/models/{name}_{LATENT_KEY}",
                )
                obsm[LATENT_KEY] = ad_joint.obsm[f"X_{mode}"]
                if output_dir is not None:
                    pd.DataFrame(obsm[LATENT_KEY], index=ad_joint.obs_names).to_csv(
                        f"{output_dir}/{name}_{LATENT_KEY}.csv",
                    )
                    print("\t\tSaved output at", f"{output_dir}/{name}_{LATENT_KEY}.csv")

        del ad_joint
        return obsm


def transfer_labels_scarches(
    ad_ref,
    ad_query=None,
    latent_key=None,
    label_key=None,
    n_neighbors=50,
    knn_model=None,
):
    """
    Transfer labels from a reference dataset to a query dataset using scanvi weighted k-nearest neighbors (KNN) utility function.

    Parameters:
        ad_ref (AnnData): The reference dataset.
        ad_query (AnnData): The query dataset.
        latent_key (str): The location of the joint embedding.
        n_neighbors (int, optional): The number of nearest neighbors to consider. Defaults to 30.

    Returns:
        tuple: A tuple containing the transferred labels and their uncertainties.
    """
    import scarches as sca

    if knn_model is None:
        knn_model = sca.utils.weighted_knn_trainer(
            train_adata=ad_ref,
            train_adata_emb=latent_key,  # location of our joint embedding
            n_neighbors=n_neighbors,
        )
    if ad_query is not None:
        pred_labels, pred_uncertainty = sca.utils.weighted_knn_transfer(
            query_adata=ad_query,
            query_adata_emb=latent_key,  # location of our embedding, query_adata.X in this case
            label_keys=label_key,  # (start of) obs column name(s) for which to transfer labels
            knn_model=knn_model,
            ref_adata_obs=ad_ref.obs,
        )
        return pred_labels, pred_uncertainty, knn_model
    else:
        return knn_model


def transfer_labels(
    ad_ref,
    ad_query=None,
    latent_key=None,
    label_keys=None,
    n_neighbors=50,
    knn_model=None,
    use_gpu=True,
    n_jobs=-1,
    weights="distance",
):
    """
    Transfer labels from a reference dataset to a query dataset using a k-nearest neighbors (KNN) classifier.

    Parameters:
        ad_ref (AnnData): The reference dataset.
        ad_query (AnnData, optional): The query dataset.
        latent_key (str): The location of the joint embedding.
        label_keys (str or list): The observation column name(s) for labels to transfer.
        n_neighbors (int, optional): The number of nearest neighbors to consider. Defaults to 50.
        knn_model (KNeighborsClassifier, optional): A pre-trained KNN model.
        use_gpu (bool, optional): Whether to use GPU acceleration with cuML. Defaults to False.
        weights (str, optional): Weight function used in prediction. Defaults to "distance".

    Returns:
        tuple: A tuple containing the transferred labels, their uncertainties, and the knn_model.
    """

    if isinstance(label_keys, str):
        label_keys = [label_keys]

    # Check if we need to initialize a new KNN model
    if knn_model is None:
        # Select the appropriate backend for KNN
        if use_gpu:
            if torch.cuda.is_available():
                try:
                    from cuml.neighbors import (
                        KNeighborsClassifier as GPUKNeighborsClassifier,
                    )
                except ImportError:
                    GPUKNeighborsClassifier = None

                if GPUKNeighborsClassifier is None:
                    print("Failed to import cuml.neighbors.KNeighborsClassifier... Using CPU instead of GPU")
            else:
                GPUKNeighborsClassifier = None

        if use_gpu and GPUKNeighborsClassifier is not None:
            print("Only uniform weighing strategy is supported for GPU.")
            knn_model = GPUKNeighborsClassifier(n_neighbors=n_neighbors, weights="uniform")

            le = sklearn.preprocessing.LabelEncoder().fit(ad_ref.obs[label_keys].values.flat)
            encoded_labels = ad_ref.obs[label_keys].apply(le.transform)

            # Fit the model on the reference data embedding and labels
            knn_model.fit(ad_ref.obsm[latent_key], encoded_labels)

        else:
            knn_model = sklearn.neighbors.KNeighborsClassifier(n_neighbors=n_neighbors, weights=weights, n_jobs=n_jobs)

            # Fit the model on the reference data embedding and labels
            knn_model.fit(ad_ref.obsm[latent_key], ad_ref.obs[label_keys])

    # If query data is provided, make predictions
    if ad_query is not None:
        pred_labels = knn_model.predict(ad_query.obsm[latent_key])
        pred_probability = knn_model.predict_proba(ad_query.obsm[latent_key])

        if "cuml" in str(type(knn_model)):
            le = sklearn.preprocessing.LabelEncoder().fit(ad_ref.obs[label_keys].values.flat)
            pred_labels = np.apply_along_axis(le.inverse_transform, 0, pred_labels)

        # Calculate prediction uncertainty
        if isinstance(label_keys, list) and len(label_keys) > 1:
            pred_uncertainty = [1 - np.max(probas, axis=1) for probas in pred_probability]
        else:
            pred_uncertainty = 1 - np.max(pred_probability, axis=1)

        return pred_labels, pred_uncertainty, knn_model
    else:
        # Return the trained model if no query data is provided
        return knn_model


def plot_confusion_matrix(
    ytrue,
    ypred,
    row_order=None,
    col_order=None,
    figsize=(8, 8),
    ax=None,
    color="Purples",
    tp_color="Greens",
    annot=False,
    highlight_tp=False,
    vmin=0,
    vmax=100,
    return_dfs=False,
    cbar=True,
):
    import scarches as sca

    df_confusion_counts = sca.classifiers.scHPL.evaluate.confusion_matrix(ytrue, ypred)
    df_confusion = df_confusion_counts.div(df_confusion_counts.sum(axis=1), axis=0)

    if row_order == "entropy":
        row_order = scipy.stats.entropy(df_confusion, axis=1).argsort()
    elif row_order == "sorted":
        row_order = np.argsort(df_confusion.index)
    elif row_order is None:
        row_order = range(len(df_confusion))

    if col_order == "entropy":
        col_order = scipy.stats.entropy(df_confusion, axis=0).argsort()
    elif col_order == "sorted":
        col_order = np.argsort(df_confusion.columns)
    elif col_order is None:
        col_order = range(len(df_confusion))

    df_confusion = df_confusion.iloc[row_order].iloc[:, col_order]
    df_confusion *= 100
    if ax is None:
        f = plt.figure(figsize=figsize)
        ax = plt.subplot()

    if highlight_tp:
        cmap_non_diag = sns.color_palette(color, as_cmap=True)
        cmap_diag = sns.color_palette(tp_color, as_cmap=True)
        mask_diag = np.zeros_like(df_confusion, dtype=bool)
        for i, i_name in enumerate(df_confusion.index):
            for j, j_name in enumerate(df_confusion.columns):
                if i_name == j_name:
                    mask_diag[i, j] = True

        sns.heatmap(
            df_confusion,
            mask=mask_diag,
            cmap=cmap_non_diag,
            cbar=cbar,
            linewidths=0.8,
            linecolor="darkgrey",
            xticklabels=True,
            yticklabels=True,
            ax=ax,
            vmin=vmin,
            vmax=vmax,
        )

        if cbar:
            if tp_color == color:
                cbar = False
            else:
                cbar = True
        sns.heatmap(
            df_confusion.where(mask_diag),
            mask=~mask_diag,
            cmap=cmap_diag,
            cbar=cbar,
            linewidths=0.8,
            linecolor="darkgrey",
            xticklabels=True,
            yticklabels=True,
            ax=ax,
            vmin=vmin,
            vmax=vmax,
            annot=annot,
            fmt=".0f",
            annot_kws={"size": 35 / np.sqrt(len(df_confusion))},
        )

    else:
        sns.heatmap(
            df_confusion,  # mask=(df_confusion==0),
            linewidths=0.8,
            linecolor="darkgrey",
            cbar=cbar,
            cmap=color,
            xticklabels=True,
            yticklabels=True,
            ax=ax,
            vmin=vmin,
            vmax=vmax,
            annot=annot,
            fmt=".0f",
            annot_kws={"size": 35 / np.sqrt(len(df_confusion))},
        )
    ax.grid(False)

    if return_dfs:
        return df_confusion_counts, df_confusion


def _compute_nmi_ari_cluster_labels(adata, labels: np.ndarray, resolution: float = 1.0) -> tuple[float, float]:
    if f"leiden_{resolution}" not in adata.obs:
        sc.tl.leiden(
            adata,
            resolution,
            key_added=f"leiden_{resolution}",
            flavor="igraph",
            n_iterations=2,
        )
    labels_pred = adata.obs[f"leiden_{resolution}"]
    nmi = sklearn.metrics.cluster.normalized_mutual_info_score(labels, labels_pred)
    ari = sklearn.metrics.cluster.adjusted_rand_score(labels, labels_pred)
    return nmi, ari


def leiden_best_nmi_ari(
    adata,
    labels: np.ndarray,
    optimize_resolution: bool = True,
    n=10,
    resolutions=None,
    n_jobs: int = 1,
) -> dict[str, float]:
    """Compute nmi and ari between leiden clusters and labels.
    Modified from https://github.com/YosefLab/scib-metrics

    This deviates from the original implementation in scib by using leiden instead of
    louvain clustering. Installing joblib allows for parallelization of the leiden
    resoution optimization.

    Parameters
    ----------
    X
        A :class:`~scib_metrics.utils.nearest_neighbors.NeighborsResults` object.
    labels
        Array of shape (n_cells,) representing label values
    optimize_resolution
        Whether to optimize the resolution parameter of leiden clustering by searching over
        10 values
    resolution
        Resolution parameter of leiden clustering. Only used if optimize_resolution is False.
    n_jobs
        Number of jobs for parallelizing resolution optimization via joblib. If -1, all CPUs
        are used.

    Returns
    -------
    nmi
        Normalized mutual information score
    ari
        Adjusted rand index score
    """
    if resolutions is None:
        resolutions = np.array([2 * x / n for x in range(1, n + 1)])
    try:
        from joblib import Parallel, delayed

        out = Parallel(n_jobs=n_jobs)(delayed(_compute_nmi_ari_cluster_labels)(adata, labels, r) for r in resolutions)
    except ImportError:
        warnings.warn("Using for loop over clustering resolutions. `pip install joblib` for parallelization.")
        out = [_compute_nmi_ari_cluster_labels(adata, labels, r) for r in resolutions]

    nmi_ari = np.array(out)
    nmi_ind = np.argmax(nmi_ari[:, 0])
    nmi, ari = nmi_ari[nmi_ind, :]
    resolution = resolutions[nmi_ind]
    return nmi_ari, {
        "nmi": nmi,
        "ari": ari,
        "resolution": resolution,
    }


def plot_transfer_labels(adata, UMAP_KEY, BATCH_KEY, CT_KEYS):
    # Define shared xlim, ylim based on the entire dataset
    xlim = (adata.obsm[UMAP_KEY][:, 0].min(), adata.obsm[UMAP_KEY][:, 0].max())
    ylim = (adata.obsm[UMAP_KEY][:, 1].min(), adata.obsm[UMAP_KEY][:, 1].max())

    u_batches = adata.obs[BATCH_KEY].unique()
    for CT_KEY in CT_KEYS:
        # Set a consistent color palette for all cell type predictions
        u_ct = adata.obs[CT_KEY].unique()
        ncol = 1

        if len(u_ct) < 20:
            color_palette = sc.pl.palettes.default_20
        elif len(u_ct) < 28:
            color_palette = sc.pl.palettes.default_28
        else:
            color_palette = sc.pl.palettes.default_102
            ncol = 2

        color_palette = dict(zip(u_ct, color_palette))

        # Plot each batch
        fig, axes = plt.subplots(1, len(u_batches), figsize=(25, 7), sharex=True, sharey=True)
        for i, batch in enumerate(u_batches):
            ad_batch = adata[adata.obs[BATCH_KEY] == batch]
            ax = axes[i]

            if i < 2:
                legend = False
            else:
                legend = False

            sns.scatterplot(
                x=ad_batch.obsm[UMAP_KEY][:, 0],
                y=ad_batch.obsm[UMAP_KEY][:, 1],
                hue=ad_batch.obs[CT_KEY],
                palette=color_palette,
                ax=ax,
                s=5,
                alpha=0.7,
                legend=legend,
            )

            # Set axis limits
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            # Set titles
            ax.set_title(f"{batch} - {CT_KEY}")

            # Add legend
            handles = [
                Line2D([0], [0], marker="o", color="w", markerfacecolor=color, label=label)
                for label, color in color_palette.items()
            ]
            plt.legend(
                handles=handles,
                title=CT_KEY,
                loc="upper left",
                bbox_to_anchor=(1.05, 1),
                prop={"size": 13},
                markerscale=2,
                frameon=False,
                ncol=ncol,
            )

        # Adjust layout for clarity
        plt.tight_layout()
        plt.show()


def subsample(adata, obs_key, n_obs, random_state=0, copy=True):
    """
    subsample each class to same cell numbers (N). Classes are given by obs_key pointing to categorical in adata.obs.
    """
    rs = np.random.RandomState(random_state)
    counts = adata.obs[obs_key].value_counts()

    indices = []
    for group, count in counts.items():
        idx = np.where(adata.obs[obs_key] == group)[0]

        if count <= n_obs:
            indices.extend(idx)
        else:
            indices.extend(rs.choice(idx, size=n_obs, replace=False))

    if copy:
        return adata[indices].copy()
    else:
        return adata[indices]
