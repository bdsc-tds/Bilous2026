import yaml
import pandas as pd
import os
import json
import h5py
import numpy as np
import scipy
import geopandas as gpd
import dask

dask.config.set({"dataframe.query-planning": False})
import dask.dataframe as dd
import spatialdata
import spatialdata_io
import warnings
import anndata as ad
import scanpy as sc
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from collections import defaultdict
from tqdm import tqdm
from types import MappingProxyType
from spatialdata.models import (
    # Image2DModel,
    # Labels2DModel,
    # Labels3DModel,
    ShapesModel,
    PointsModel,
)
from functools import lru_cache


@lru_cache(maxsize=None)
def print_once(msg):
    print(msg)


try:
    import msgspec
except ImportError:
    pass

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "../../config/config.yml")


def config(path=config_path):
    """
    Read the configuration file and return a dictionary of config values.

    Parameters
    ----------
    path : str
        The path to the configuration file. Defaults to the value of
        `config_path` if not provided.

    Returns
    -------
    cfg : dict
        A dictionary of configuration values. All values are strings and
        have been converted to absolute paths by prepending the value of
        `cfg["base_dir"]`.
    """
    with open(path, "r") as ymlfile:
        cfg = yaml.safe_load(ymlfile)
        for k in cfg.keys():
            if "umap_" in k:
                continue
            cfg[k] = os.path.join(cfg["base_dir"], cfg[k])
    return cfg


def soma_to_anndata(soma_uri, measurement_name, X_layer_name, return_experiment=False, **kwargs):
    """
    Export a SOMA experiment to an anndata object.

    Parameters
    ----------
    soma_uri (str): The URI (path) of the SOMA experiment.
    measurement_name (str): The measurement name (e.g., 'RNA') to extract.
    X_layer_name (str): The layer name in the X matrix to extract (e.g., 'counts').
    return_experiment (bool): Whether to return the SOMA experiment.

    Returns
    -------
    anndata.AnnData: The anndata object.
    """
    import tiledbsoma as soma
    import tiledbsoma.io

    # Open the SOMA experiment
    experiment = soma.open(soma_uri)
    # Export the SOMA experiment to anndata format
    ad = soma.io.to_anndata(
        experiment=experiment,
        measurement_name=measurement_name,
        X_layer_name=X_layer_name,
        **kwargs,
    )

    if return_experiment:
        return ad, experiment
    else:
        return ad


def xenium_specs(path):
    path = Path(path)
    with open(path / "experiment.xenium") as f:
        specs = json.load(f)
    return specs


######### Xenium readers
def xenium_samples_files(dir_segmentation_condition, segmentation=None, samples=None):
    """
    Get a dictionary of files for each sample in a Xenium segmentation run.

    Parameters
    ----------
    dir_segmentation_condition (str): The directory path of the segmentation run.
    segmentation (str): The segmentation name, e.g., 'default', '10x_5um', 'baysor'.
    samples (list): The sample names to include. If None, include all samples.

    Returns
    -------
    dict: A dictionary of files for each sample.
    """
    files = {}
    for sample_path in Path(dir_segmentation_condition).iterdir():
        for sample_path in sample_path.iterdir():
            sample_name = sample_path.stem

            if samples is not None and sample_name not in samples:
                continue
            elif "corrupted" in sample_name:
                continue
            else:
                if segmentation != "default":
                    files[sample_name] = sample_path / "normalised_results/outs"
                else:
                    files[sample_name] = sample_path

    return files


def read_xenium_specs(xenium_specs_file):
    xenium_specs_file = Path(xenium_specs_file)
    with open(xenium_specs_file) as f:
        specs = json.load(f)
    return specs


def xenium_proseg(
    path,
    cells_boundaries=True,
    cells_boundaries_layers=True,
    nucleus_boundaries=False,
    cells_as_circles=False,
    cells_labels=False,
    nucleus_labels=False,
    transcripts=True,
    morphology_mip=False,
    morphology_focus=False,
    aligned_images=False,
    cells_table=True,
    n_jobs=1,
    imread_kwargs=MappingProxyType({}),
    image_models_kwargs=MappingProxyType({}),
    labels_models_kwargs=MappingProxyType({}),
    cells_metadata=True,
    xeniumranger_dir=None,
    xenium_specs=True,
    pandas_engine="pyarrow",
    verbose=False,
):
    """
    Reads a Xenium segmentation run and returns a `SpatialData` object.

    Parameters
    ----------
    path : str
        The directory path of the segmentation run.
    cells_boundaries : bool or str, optional
        The file path of the cell boundaries GeoJSON file.
    cells_boundaries_layers : bool or str, optional
        The file path of the cell boundaries layers GeoJSON file.
    nucleus_boundaries : bool or str, optional
        Whether to read nucleus boundaries. Not implemented.
    cells_as_circles : bool or str, optional
        Whether to read cells as circles.
    cells_labels : bool or str, optional
        Whether to read cells labels. Not implemented.
    nucleus_labels : bool or str, optional
        Whether to read nucleus labels. Not implemented.
    transcripts : bool or str, optional
        The file path of the transcripts CSV file.
    morphology_mip : bool, optional
        Whether to read morphology MIP images.
    morphology_focus : bool, optional
        Whether to read morphology focus images.
    aligned_images : bool, optional
        Whether to read aligned images.
    cells_table : bool or str, optional
        The file path of the cells table CSV file.
    n_jobs : int, optional
        The number of jobs to use. Not implemented.
    imread_kwargs : dict, optional
        Keyword arguments to pass to `imread`.
    image_models_kwargs : dict, optional
        Keyword arguments to pass to `ImageModel`.
    labels_models_kwargs : dict, optional
        Keyword arguments to pass to `LabelsModel`.
    cells_metadata : bool or str, optional
        The file path of the cells metadata CSV file.
    xeniumranger_dir : str, optional
        The directory path of the XeniumRanger run.
    xenium_specs : bool or str, optional
        The file path of the Xenium specs file.
    pandas_engine : str, optional
        The pandas engine to use when reading CSV files.

    Returns
    -------
    sdata : SpatialData
        The `SpatialData` object.
    """
    path = Path(path)

    # unsupported options compared to spatialdata_io.xenium
    if nucleus_boundaries:
        raise ValueError("reading nucleus_boundaries not implemented for proseg")
    if cells_labels:
        raise ValueError("reading cells_labels not implemented for proseg")
    if nucleus_labels:
        raise ValueError("reading nucleus_labels not implemented for proseg")
    if n_jobs > 1:
        raise ValueError("n_jobs>1 not supported")

    # default expected file paths
    def parse_arg(arg, default):
        if arg:
            return default
        elif isinstance(arg, str):
            return Path(arg)
        else:
            return arg

    cells_metadata = parse_arg(cells_metadata, path / "cell-metadata.csv.gz")
    cells_boundaries = parse_arg(cells_boundaries, path / "cell-polygons.geojson.gz")
    cells_boundaries_layers = parse_arg(cells_boundaries_layers, path / "cell-polygons-layers.geojson.gz")
    transcripts = parse_arg(transcripts, path / "transcript-metadata.csv.gz")
    cells_table = parse_arg(cells_table, path / "expected-counts.csv.gz")
    xeniumranger_dir = parse_arg(xeniumranger_dir, None)
    if xeniumranger_dir is not None or isinstance(xenium_specs, str):
        xenium_specs = parse_arg(xenium_specs, xeniumranger_dir / "experiment.xenium")

    ### images
    if morphology_mip or morphology_focus or aligned_images:
        if verbose:
            print("Reading images...")
        sdata_images = spatialdata_io.xenium(
            xeniumranger_dir,
            cells_table=False,
            cells_as_circles=False,
            cells_boundaries=False,
            nucleus_boundaries=False,
            cells_labels=False,
            nucleus_labels=False,
            transcripts=False,
            morphology_mip=morphology_mip,
            morphology_focus=morphology_focus,
            aligned_images=aligned_images,
            imread_kwargs=imread_kwargs,
            image_models_kwargs=image_models_kwargs,
            labels_models_kwargs=labels_models_kwargs,
        )

        images = sdata_images.images
    else:
        images = {}

    ### tables
    region = "cell_polygons"
    region_key = "region"
    instance_key = "cell_id"

    # flag columns not corresponding to genes
    if isinstance(cells_table, Path):
        if verbose:
            print("Reading cells table...")
        df_table = pd.read_csv(cells_table, engine=pandas_engine)

        control_columns = df_table.columns.str.contains("|".join(["BLANK_", "Codeword", "NegControl"]))

        table = ad.AnnData(
            df_table.iloc[:, ~control_columns],
            uns={
                "spatialdata_attrs": {
                    "region": region,
                    "region_key": region_key,
                    "instance_key": instance_key,
                }
            },
        )

        if isinstance(cells_metadata, Path):
            if verbose:
                print("Reading cells metadata...")

            df_cells_metadata = pd.read_csv(cells_metadata, engine=pandas_engine).rename(columns={"cell": "cell_id"})
            table.obs = pd.concat((df_cells_metadata, df_table.iloc[:, control_columns]), axis=1)
            table.obsm["spatial"] = table.obs[["centroid_x", "centroid_y"]].values
            table.obs = table.obs
        table.obs[region_key] = region

        # sparsify .X
        table.X = scipy.sparse.csr_matrix(table.X)
        tables = {"table": table}
    else:
        tables = {}

    ### labels
    # not implemented
    labels = {}

    ### points
    if isinstance(transcripts, Path):
        if verbose:
            print("Reading transcripts...")

        df_transcripts = dd.read_csv(transcripts, blocksize=None).rename(
            columns={"gene": "feature_name", "assignment": "cell_id"}
        )
        points = {"transcripts": PointsModel.parse(df_transcripts)}
    else:
        points = {}

    ### shapes
    shapes = {}

    # read specs
    if isinstance(xenium_specs, Path):
        if verbose:
            print("Reading specs...")
        specs = read_xenium_specs(xenium_specs)

        # get xenium pixel size
        scale = spatialdata.transformations.Scale(
            [1.0 / specs["pixel_size"], 1.0 / specs["pixel_size"]], axes=("x", "y")
        )
        transformations = {"global": scale}
    else:
        transformations = None
        if isinstance(cells_boundaries, Path) or isinstance(cells_boundaries_layers, Path):
            warnings.warn(
                """
                Couldn't load xenium specs file with pixel size. 
                Not applying scale transformations to shapes.
                Please specify xeniumranger_dir or xenium_specs
                """
            )

    # read cells boundaries
    if isinstance(cells_boundaries, Path):
        if verbose:
            print("Reading cells boundaries...")

        df_cells_boundaries = gpd.read_file("gzip://" + cells_boundaries.as_posix()).rename(columns={"cell": "cell_id"})
        shapes["cells_boundaries"] = ShapesModel.parse(df_cells_boundaries, transformations=transformations)

    # read cells boundaries layers
    if isinstance(cells_boundaries_layers, Path):
        if verbose:
            print("Reading cells boundaries layers...")

        df_cells_boundaries_layers = gpd.read_file("gzip://" + cells_boundaries_layers.as_posix()).rename(
            columns={"cell": "cell_id"}
        )
        shapes["cells_boundaries_layers"] = ShapesModel.parse(
            df_cells_boundaries_layers, transformations=transformations
        )

    # convert cells boundaries to circles
    if cells_as_circles:
        if verbose:
            print("Converting cells boundaries to circle...")

        shapes["cells_boundaries_circles"] = spatialdata.to_circles(shapes["cells_boundaries"])

    ### sdata
    sdata = spatialdata.SpatialData(images=images, labels=labels, points=points, shapes=shapes, tables=tables)

    return sdata


def read_xenium_sample(
    path,
    cells_as_circles=False,
    cells_boundaries=False,
    cells_boundaries_layers=False,
    nucleus_boundaries=False,
    cells_labels=False,
    nucleus_labels=False,
    transcripts=False,
    morphology_mip=False,
    morphology_focus=False,
    aligned_images=False,
    cells_table=True,
    anndata=False,
    xeniumranger_dir=None,
    sample_name=None,
):
    """
    Reads a xenium sample from a directory path.

    Parameters
    ----------
    path (str): The directory path of the segmentation run.
    cells_as_circles (bool): Whether to include cell polygons as circles or not.
    cells_boundaries (bool): Whether to include cell boundaries or not.
    nucleus_boundaries (bool): Whether to include nucleus boundaries or not.
    cells_labels (bool): Whether to include cell labels or not.
    nucleus_labels (bool): Whether to include nucleus labels or not.
    transcripts (bool): Whether to include transcript locations or not.
    morphology_mip (bool): Whether to include morphology MIP or not.
    morphology_focus (bool): Whether to include morphology focus or not.
    aligned_images (bool): Whether to include aligned images or not.
    cells_table (bool): Whether to include cells table or not.
    anndata (bool): Whether to return only the anndata object or the full spatialdata object.
    xeniumranger_dir (str): Path to xeniumranger output directory (for proseg raw only)
    sample_name (str): The sample name.

    Returns
    -------
    If anndata, returns a tuple of the sample name and anndata object.
    Otherwise, returns a tuple of the sample name and spatialdata object.

    If sample_name is None, sample_name is not returned
    """
    path = Path(path)
    kwargs = dict(
        cells_as_circles=cells_as_circles,
        cells_boundaries=cells_boundaries,
        nucleus_boundaries=nucleus_boundaries,
        cells_labels=cells_labels,
        nucleus_labels=nucleus_labels,
        transcripts=transcripts,
        morphology_mip=morphology_mip,
        morphology_focus=morphology_focus,
        aligned_images=aligned_images,
        cells_table=cells_table,
    )

    # automatically check whether path is a folder with proseg raw outputs or in xeniumranger format
    if (path / "expected-counts.csv.gz").exists():
        reader = xenium_proseg
        kwargs["cells_boundaries_layers"] = cells_boundaries_layers
        kwargs["xeniumranger_dir"] = xeniumranger_dir
    else:
        reader = spatialdata_io.xenium

    sdata = reader(path, **kwargs)

    adata = sdata["table"]
    adata.obs_names = adata.obs["cell_id"].values

    metrics_path = path / "metrics_summary.csv"
    if metrics_path.exists():
        adata.uns["metrics_summary"] = pd.read_csv(metrics_path)
    else:
        print("metrics_summary.csv not found at:", metrics_path)

    if sample_name is None:
        if anndata:
            return adata
        else:
            return sdata
    else:
        if anndata:
            return sample_name, adata
        else:
            return sample_name, sdata


def read_xenium_samples(
    data_dirs,
    cells_as_circles=False,
    cells_boundaries=False,
    cells_boundaries_layers=False,
    nucleus_boundaries=False,
    cells_labels=False,
    nucleus_labels=False,
    transcripts=False,
    morphology_mip=False,
    morphology_focus=False,
    aligned_images=False,
    cells_table=True,
    anndata=False,
    sample_name_as_key=True,
    xeniumranger_dir=None,
    max_workers=None,
    pool_mode="thread",
):
    """
    Reads in a dictionary of sample directories and returns a dictionary of
    AnnData objects or spatialdata objects depending on the anndata flag.

    Parameters
    ----------
    data_dirs : dict or list
        A dictionary of sample directories or a list of paths to sample directories.
    cells_as_circles : bool, optional
        Whether to include cell boundary data as circles, by default False
    cells_boundaries : bool, optional
        Whether to include cell boundary data, by default False
    cells_boundaries : bool, optional
        Whether to include cell boundary layers data (for proseg raw only), by default False
    nucleus_boundaries : bool, optional
        Whether to include nucleus boundary data, by default False
    cells_labels : bool, optional
        Whether to include cell labels, by default False
    nucleus_labels : bool, optional
        Whether to include nucleus labels, by default False
    transcripts : bool, optional
        Whether to include transcript data, by default False
    morphology_mip : bool, optional
        Whether to include morphology data at the maximum intensity projection, by default False
    morphology_focus : bool, optional
        Whether to include morphology data at the focus, by default False
    aligned_images : bool, optional
        Whether to include aligned images, by default False
    cells_table (bool):
        Whether to include cells table or not, by default True
    anndata : bool, optional
        Whether to only return an AnnData object, by default False
    sample_name_as_key: bool, optional
        Whether to use the sample name as the key in the return dictionary,
        otherwise returns full path as key
    xeniumranger_dir: str, optional
        Path to xeniumranger output dir (for proseg raw only)
    max_workers : int, optional
        Maximum number of workers to use for parallel processing, by default None
    pool_mode : str, optional
        Pool mode for parallel processing, "thread" or "process", by default "thread"

    Returns
    -------
    dict
        A dictionary of sample names mapped to AnnData objects or spatialdata objects.
    """
    if isinstance(data_dirs, list):
        sample_names = [Path(path).stem if sample_name_as_key else path for path in data_dirs]
        data_dirs = {sample_name: path for sample_name, path in zip(sample_names, data_dirs)}

    # Parallel processing
    if pool_mode == "process":
        pool = ProcessPoolExecutor
    elif pool_mode == "thread":
        pool = ThreadPoolExecutor

    sdatas = {}
    with pool(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                read_xenium_sample,
                path,
                cells_as_circles,
                cells_boundaries,
                cells_boundaries_layers,
                nucleus_boundaries,
                cells_labels,
                nucleus_labels,
                transcripts,
                morphology_mip,
                morphology_focus,
                aligned_images,
                cells_table,
                anndata,
                xeniumranger_dir,
                sample_name,
            )
            for sample_name, path in data_dirs.items()
        ]

        for future in as_completed(futures):
            try:
                sample_name, result = future.result()
                sdatas[sample_name] = result
            except Exception as e:
                print(f"Error processing {e}")

    return sdatas


######### 10x writers


def write_10X_h5(adata, file):
    """Writes adata to a 10X-formatted h5 file.
    taken from https://github.com/scverse/anndata/issues/595

    Note that this function is not fully tested and may not work for all cases.
    It will not write the following keys to the h5 file compared to 10X:
    '_all_tag_keys', 'pattern', 'read', 'sequence'

    Args:
        adata (AnnData object): AnnData object to be written.
        file (str): File name to be written to. If no extension is given, '.h5' is appended.

    Raises:
        FileExistsError: If file already exists.

    Returns:
        None
    """

    if isinstance(file, Path):
        file = file.as_posix()

    if ".h5" not in file:
        file = f"{file}.h5"
    if Path(file).exists():
        raise FileExistsError(f"There already is a file `{file}`.")

    def int_max(x):
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)

    def str_max(x):
        return max([len(i) for i in x])

    if not scipy.sparse.issparse(adata.X):
        adata.X = scipy.sparse.csr_matrix(adata.X)
    if "genome" not in adata.var:
        adata.var["genome"] = "undefined"
    if "feature_types" not in adata.var:
        adata.var["feature_types"] = "Gene Expression"
    if "gene_ids" not in adata.var:
        adata.var["gene_ids"] = adata.var_names

    w = h5py.File(file, "w")
    grp = w.create_group("matrix")
    grp.create_dataset(
        "barcodes",
        data=np.array(adata.obs_names, dtype=f"|S{str_max(adata.obs_names)}"),
    )
    grp.create_dataset("data", data=np.array(adata.X.data, dtype=f"<i{int_max(adata.X.data)}"))
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset(
        "feature_type",
        data=np.array(adata.var.feature_types, dtype=f"|S{str_max(adata.var.feature_types)}"),
    )
    ftrs.create_dataset(
        "genome",
        data=np.array(adata.var.genome, dtype=f"|S{str_max(adata.var.genome)}"),
    )
    ftrs.create_dataset(
        "id",
        data=np.array(adata.var.gene_ids, dtype=f"|S{str_max(adata.var.gene_ids)}"),
    )
    ftrs.create_dataset("name", data=np.array(adata.var.index, dtype=f"|S{str_max(adata.var.index)}"))
    grp.create_dataset("indices", data=np.array(adata.X.indices, dtype=f"<i{int_max(adata.X.indices)}"))
    grp.create_dataset("indptr", data=np.array(adata.X.indptr, dtype=f"<i{int_max(adata.X.indptr)}"))
    grp.create_dataset(
        "shape",
        data=np.array(list(adata.X.shape)[::-1], dtype=f"<i{int_max(adata.X.shape)}"),
    )


######### RCTD readers


def read_json(file_path):
    with open(file_path, "r") as file:
        return json.load(file)


def read_json_msgspec(file_path):
    with open(file_path, "rb") as file:
        return msgspec.json.decode(file.read())


def _rds2py_dict_to_df(r_obj_df, mode="results_df"):
    if mode == "results_df":
        r_obj_df_columns = r_obj_df["attributes"]["names"]["data"]
        r_obj_df_index = r_obj_df["attributes"]["row.names"]["data"]
        pandas_df = pd.DataFrame(
            [r_obj_df["data"][i]["data"] for i in range(len(r_obj_df["data"]))],
            index=r_obj_df_columns,
            columns=r_obj_df_index,
        ).T
    elif mode == "weights":
        r_obj_df_columns = r_obj_df["attributes"]["dimnames"]["data"][1]["data"]
        r_obj_df_index = r_obj_df["attributes"]["dimnames"]["data"][0]["data"]
        r_obj_df["data"] = r_obj_df["data"].reshape(r_obj_df["attributes"]["dim"]["data"], order="F")
        pandas_df = pd.DataFrame(r_obj_df["data"], index=r_obj_df_index, columns=r_obj_df_columns)

    return pandas_df


def read_rctd_sample(sample_name, rctd_results_path):
    """
    Reads RCTD results from a single sample and returns a dictionary containing:

    - results_df: a pandas DataFrame with columns to be added to the anndata object's obs
    - weights: a pandas Series with the weights for each cell for the given reference
    - weights_doublet: (not implemented) a pandas Series with the weights for each cell for doublets for the given reference
    - singlet_scores: (not implemented) a pandas Series with the singlet scores for each cell for the given reference

    Parameters
    ----------
    sample_name: str
        The name of the sample
    rctd_results_path : str
        The path to the sample RCTD results
    rsuffix : str, optional
        The suffix to append to the reference name when storing to the anndata objects

    Returns
    -------
    A tuple containing the sample name and the results dictionary
    """
    from rds2py import read_rds

    r_obj = read_rds(rctd_results_path)

    results = r_obj["attributes"]["results"]
    results_keys = results["attributes"]["names"]["data"]
    results_keys_idx = {k: results_keys.index(k) for k in results_keys}

    pandas_results = {}
    for k in ["results_df", "weights"]:
        pandas_results[k] = _rds2py_dict_to_df(results["data"][results_keys_idx[k]], mode=k)

    return sample_name, pandas_results


def read_rctd_samples(ads, rctd_results_paths, prefix=""):
    """
    Read RCTD results into anndata objects in parallel using ProcessPoolExecutor.

    Parameters
    ----------
    ads : dict of anndata.AnnData
        The anndata objects to be updated.
    rctd_results_paths : str
        The directory containing the RCTD results.
    prefix : str, optional
        The prefix to append to the reference name when storing to the anndata objects.

    Returns
    -------
    None
    """

    # Use ProcessPoolExecutor for CPU-bound tasks
    with ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(
                read_rctd_sample,
                sample_name,
                rctd_results_paths[sample_name],
            ): sample_name
            for sample_name in ads.keys()
        }

        # Update anndata objects in the parent process
        for future in as_completed(futures):
            try:
                sample_name, results = future.result()
                if results:
                    ad = ads[sample_name]
                    ad.obs = ad.obs.join(results["results_df"].add_prefix(prefix))
                    ad.uns[f"{prefix}_weights"] = results["weights"]

            except Exception as e:
                print(f"Error processing sample {futures[future]}: {e}")


###### coexpression files readers
def read_coexpression_file(k, method, target_count, results_dir):
    """
    Worker function to read the coexpression and positivity rate parquet for a single sample.

    Parameters
    ----------
    k : tuple
        The sample name tuple (segmentation, condition, panel, sample, sample).
    method : str
        The coexpression method.
    target_count : int
        The target count of the coexpression method.
    results_dir : Path
        The directory containing the coexpression results.

    Returns
    -------
    method : str
        The coexpression method.
    target_count : int
        The target count of the coexpression method.
    cc : pd.DataFrame
        The coexpression matrix.
    pos_rate : pd.Series
        The positivity rate.
    """
    out_file_coexpr = results_dir / f"{'/'.join(k)}/coexpression_{method}_{target_count}.parquet"
    out_file_pos_rate = results_dir / f"{'/'.join(k)}/positivity_rate_{method}_{target_count}.parquet"

    cc = pd.read_parquet(out_file_coexpr)
    pos_rate = pd.read_parquet(out_file_pos_rate)[0]
    return method, target_count, cc, pos_rate


def read_coexpression_files(cc_paths, results_dir):
    """
    Reads coexpression parquet files for multiple methods and target counts in parallel using ThreadPoolExecutor.

    Parameters
    ----------
    cc_paths : list of tuples
        A list of tuples containing the key `k`, i.e., a sample name tuple (segmentation, condition, panel, sample, sample)
        and the method and target count to read.
    results_dir : str
        The directory containing the coexpression results.

    Returns
    -------
    CC : dict
        A dictionary with the coexpression matrices for each method and target count.
    pos_rate : dict
        A dictionary with the positivity rates for each method and target count.
    """
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(read_coexpression_file, k, method, target_count, results_dir)
            for k, method, target_count in cc_paths
        ]

        CC = {}
        pos_rate = {}
        for future in as_completed(futures):
            method, target_count, cc, pr = future.result()
            k = cc_paths[futures.index(future)][0]  # Retrieve the `k` corresponding to this future

            if k not in CC:
                CC[k] = {}
            if k not in pos_rate:
                pos_rate[k] = {}

            CC[k][method, target_count] = cc
            pos_rate[k][method, target_count] = pr
    return CC, pos_rate


def get_gene_panel_info(path):
    with open(path, "r") as f:
        gene_panel = json.load(f)["payload"]["targets"]

    gene_panel_info = pd.DataFrame(columns=["codewords"])
    for i, g in enumerate(gene_panel):
        gene_panel_info.at[i, "gene_coverage"] = g["info"]["gene_coverage"]
        gene_panel_info.at[i, "id"] = g["type"]["data"].get("id")
        gene_panel_info.at[i, "name"] = g["type"]["data"]["name"]
        gene_panel_info.at[i, "codewords"] = g["codewords"]
        gene_panel_info.at[i, "source_category"] = g["source"]["category"]
        gene_panel_info.at[i, "source_design_id"] = g["source"]["identity"]["design_id"]
        gene_panel_info.at[i, "source_name"] = g["source"]["identity"]["name"]
        gene_panel_info.at[i, "source_version"] = g["source"]["identity"].get("version")
    return gene_panel_info


def _df_split(df):
    """
    Utility function for read_anndata_folder to parse a dataframe into separate .obsm and .varm field
    Splits a DataFrame into multiple DataFrames based on a common prefix in the column names.

    Given a DataFrame with columns like 'gene1', 'gene2', 'control1', 'control2', etc., this function
    will split the DataFrame into two DataFrames, one with the 'gene' columns and one with the 'control'
    columns.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame to split.

    Returns
    -------
    dfs : dict
        A dictionary where the keys are the prefixes and the values are the DataFrames.
    """

    def _split_text_from_numbers(s):
        head = s.rstrip("0123456789")
        tail = s[len(head) :]
        return head, tail

    columns = df.columns
    columns_split = [_split_text_from_numbers(col) for col in columns]

    # Group columns by prefix
    column_groups = defaultdict(list)
    for prefix, suffix in columns_split:
        column_groups[prefix].append((prefix, suffix))

    # Create DataFrames
    dfs = defaultdict(pd.DataFrame)
    for prefix, columns in column_groups.items():
        # Sort columns by suffix
        columns.sort(key=lambda x: int(x[1]))
        columns = ["".join(c) for c in columns]
        # Create a DataFrame with sorted columns
        dfs[prefix][columns] = df[columns]

    return dfs


def write_anndata_folder(adata, adata_dir):
    """
    Writes an AnnData object to a folder as a series of csv files and .h5 files for X and layers.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object to write.
    adata_dir : str
        The path to the folder where the AnnData object will be written.

    Notes
    -----
    This function writes `adata.obs`, `adata.var`, `adata.uns` and `adata.X` to separate csv files.
    `adata.layers` are written as separate .h5 files in the same folder.
    """
    adata = adata.copy()
    del adata.uns
    adata.write_csvs(adata_dir)
    write_10X_h5(adata, f"{adata_dir}/X.h5")

    for layer in adata.layers:
        adata.X = adata.layers["pvals"]
        write_10X_h5(adata, f"{adata_dir}/{layer}.h5")


def read_anndata_folder(adata_dir):
    """
    Reads an anndata folder saved with write_anndata_folder and returns an AnnData object.

    Parameters
    ----------
    adata_dir : str
        The path to the anndata folder.

    Returns
    -------
    adata : anndata.AnnData
        The AnnData object.
    """

    # read X
    adata = sc.read_10x_h5(f"{adata_dir}/X.h5")

    # read layers
    for f in Path(adata_dir).glob("*.h5"):
        if f.name == "X.h5":
            continue
        adata.layers[f.stem] = sc.read_10x_h5(f).X

    # read obs, obsm, var, varm
    for attr in ["obs", "obsm", "var", "varm"]:
        if attr in ["obs", "var"]:
            try:
                df = pd.read_csv(f"{adata_dir}/{attr}.csv", index_col=0)
                setattr(adata, attr, df)
            except pd.errors.EmptyDataError:
                pass
        elif attr == "obsm":
            try:
                df = pd.read_csv(f"{adata_dir}/{attr}.csv")
                dfs = _df_split(df)
                for k, v in dfs.items():
                    adata.obsm[k] = v.values
            except pd.errors.EmptyDataError:
                pass
        elif attr == "varm":
            try:
                df = pd.read_csv(f"{adata_dir}/{attr}.csv")
                dfs = _df_split(df)
                for k, v in dfs.items():
                    adata.varm[k] = v
            except pd.errors.EmptyDataError:
                pass
    return adata


def _read_count_correction_sample(sample_name, corrected_counts_path):
    """Reads a 10x h5 file using scanpy."""
    try:
        adata = sc.read_10x_h5(corrected_counts_path)
        return sample_name, adata
    except Exception as e:
        print(f"Error reading {corrected_counts_path}: {e}")
        return sample_name, None  # Return None in case of an error


def read_count_correction_samples(xenium_paths, correction_methods):
    """
    Reads corrected count samples in parallel using ThreadPoolExecutor.

    Args:
        xenium_paths (dict): A dictionary where keys are correction methods and values are dictionaries
                            mapping sample names to corrected counts file paths.  Assumes `xenium_paths[correction_method]`
                            is a dictionary with keys as sample_name and values as path to the .h5 file.
        correction_methods (list): A list of correction methods.
    Returns:
        dict: A dictionary where keys are correction methods, and values are dictionaries mapping sample names
              to AnnData objects (or None if reading failed).
    """

    xenium_corrected_counts = {}

    for correction_method in correction_methods:  # Skip the first correction method
        xenium_corrected_counts[correction_method] = {}

        with ThreadPoolExecutor() as executor:
            futures = {
                executor.submit(_read_count_correction_sample, sample_name, xenium_corr_path): (
                    correction_method,
                    sample_name,
                )
                for sample_name, xenium_corr_path in xenium_paths[correction_method].items()
            }

            # Progress bar with total number of samples
            for future in tqdm(as_completed(futures), total=len(futures), desc=f"Processing {correction_method}"):
                try:
                    sample_name, adata = future.result()
                    if adata is not None:
                        xenium_corrected_counts[correction_method][sample_name] = adata
                    else:
                        xenium_corrected_counts[correction_method][sample_name] = None
                except Exception as e:
                    correction_method, sample_name = futures[future]
                    xenium_corrected_counts[correction_method][sample_name] = None  # Store None in case of error

    return xenium_corrected_counts


def _read_contamination_metrics_results_sample(
    results_dir: Path,
    reference: str,
    method: str,
    level: str,
    mixture_k: int,
    num_samples: int,
    correction_method: str,
    segmentation: Path,
    condition: Path,
    panel: Path,
    donor: Path,
    sample: Path,
    normalisation: str,
    layer: str,
    radius=10,
    n_splits=5,
    n_permutations=30,
    n_repeats=5,
    top_n=20,
    markers_mode="diffexpr",
    cv_mode="spatial",
    scoring="precision",
    evaluation: str = "diffexpr",
    genes_name="all",
    train_mode: str = "multivariate",
):
    """
    Reads contamination metrics results for a single sample

    Args:
        results_dir (Path): Path to the results directory.
        reference (str): Reference string.
        method (str): Method string.
        level (str): Level string.
        mixture_k (int): Mixture K value.
        num_samples (int): Number of samples.
        correction_method (str): The correction method used.
        segmentation (Path): Path to the segmentation directory.
        condition (Path): Path to the condition directory.
        panel (Path): Path to the panel directory.
        donor (Path): Path to the donor directory.
        sample (Path): Path to the sample directory.
        normalisation (str): The normalisation method used.
        layer (str): The layer used.
        radius (int): The radius used.
        n_permutations (int): The number of permutations used.
        n_repeats (int): The number of repeats used.
        top_n (int): The top N used.
        markers_mode (str): The markers mode used.
        cv_mode (str): The cross-validation mode used.
        scoring (str): The scoring metric used.
        evaluation (str): evaluation test to load, 'diffexpr' or 'logreg'.
        genes_name (str): Load results computed on a marker genes subset with folder name 'genes_name' or all genes ('all').
        train_mode (str): The train mode used.

    Returns:
        tuple or None:
            - If all required files exist and are successfully read:
                A tuple containing:
                    - correction_method (str): The correction method used.
                    - k (tuple): A tuple of keys (segmentation, condition, panel, donor, sample).
                    - loaded_data (tuple): A tuple containing the loaded DataFrames and summary stats based on the boolean flags.
            - If any required file does not exist or an error occurs during reading:
                None.
    """

    def _read(out_dict, key, path):
        if path.exists():
            if path.suffix == ".parquet":
                res = pd.read_parquet(path, engine="pyarrow")
            elif path.suffix == ".json":
                with open(path, "r") as f:
                    res = json.load(f)
            if len(res):
                out_dict[key] = res
            else:
                print(f"File is empty: {path}")
        else:
            print(f"File does not exist: {path}")

    if condition.stem == "melanoma":
        if level == "Level2.1":
            level = "Level1"
        if reference == "matched_reference_combo":
            reference = "external_reference"

    k = (segmentation.stem, condition.stem, panel.stem, donor.stem, sample.stem)
    name = "/".join(k)
    name_params_diffexpr = f"{markers_mode}_{radius=}_{top_n=}"
    name_params_logreg = (
        f"{markers_mode}_{radius=}_{n_permutations=}_{n_splits=}_{top_n=}_{scoring}_{cv_mode}_{train_mode}"
    )

    # folder name
    folder_diffexpr = f"contamination_metrics_{name_params_diffexpr}"
    folder_logreg = f"contamination_metrics_{name_params_logreg}_logreg"

    if genes_name is not None:
        folder_diffexpr += f"/{genes_name}"
        folder_logreg += f"/{genes_name}"

    if evaluation == "logreg":
        folder = folder_logreg
    elif evaluation == "diffexpr":
        folder = folder_diffexpr

    # full path prefix
    if correction_method == "raw":
        name = f"{correction_method}/{name}"
    elif correction_method in ["resolvi", "resolvi_panel_use_batch=True", "resolvi_panel_use_batch=False"]:
        name = f"{correction_method}/{name}/{mixture_k=}/{num_samples=}/"
    elif correction_method in [
        "resolvi_supervised",
        "resolvi_panel_supervised_use_batch=True",
        "resolvi_panel_supervised_use_batch=False",
    ]:
        name = f"{correction_method}/{name}/{normalisation}/reference_based/{reference}/{method}/{level}/{mixture_k=}/{num_samples=}"
    elif "ovrlpy" in correction_method:
        name = f"{correction_method}/{name}"
    elif correction_method == "split_fully_purified":
        name = f"{correction_method}/{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/split_fully_purified/"

    prefix = (results_dir / f"{folder}/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_").as_posix()

    # always get it from diffexpr folder even for logreg (which does not output it)
    prefix_ctj_marker_genes = (
        results_dir / f"{folder_diffexpr}/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_"
    ).as_posix()

    # Load data
    loaded_data = {}

    # markers files
    out_file_df_ctj_marker_genes = Path(prefix_ctj_marker_genes + "marker_genes.parquet")

    if out_file_df_ctj_marker_genes.exists():
        loaded_data["df_ctj_marker_genes"] = pd.read_parquet(out_file_df_ctj_marker_genes, engine="pyarrow")
    else:
        print(f"File does not exist: {out_file_df_ctj_marker_genes}")

    if evaluation == "diffexpr":
        # diffexpr files
        out_file_df_diffexpr = Path(prefix + "diffexpr.parquet")
        out_file_df_markers_rank_significance_diffexpr = Path(prefix + "markers_rank_significance_diffexpr.parquet")
        out_file_summary_stats = Path(prefix + "summary_stats.json")

        _read(loaded_data, "df_diffexpr", out_file_df_diffexpr)
        _read(loaded_data, "df_markers_rank_significance_diffexpr", out_file_df_markers_rank_significance_diffexpr)
        _read(loaded_data, "summary_stats", out_file_summary_stats)

    elif evaluation == "logreg":
        # logreg files
        out_file_df_permutations_logreg = Path(prefix + "permutations_logreg.parquet")
        out_file_df_importances_logreg = Path(prefix + "importances_logreg.parquet")
        out_file_df_markers_rank_significance_logreg = Path(prefix + "markers_rank_significance_logreg.parquet")

        _read(loaded_data, "df_permutations_logreg", out_file_df_permutations_logreg)
        _read(loaded_data, "df_importances_logreg", out_file_df_importances_logreg)
        _read(loaded_data, "df_markers_rank_significance_logreg", out_file_df_markers_rank_significance_logreg)

    else:
        raise ValueError(f"Unknown evaluation: {evaluation}")

    return correction_method, k, loaded_data


def read_contamination_metrics_results(
    results_dir: Path,
    correction_methods: list[str],
    xenium_std_seurat_analysis_dir: Path,
    reference: str,
    method: str,
    level: str,
    mixture_k: int,
    num_samples: int,
    normalisation: str,
    layer: str,
    radius=10,
    n_splits=5,
    n_permutations=30,
    n_repeats=5,
    top_n=20,
    markers_mode="diffexpr",
    cv_mode="spatial",
    scoring="precision",
    ref_condition: str = None,
    ref_panel: str = None,
    evaluation: str = "diffexpr",
    genes_name: str = "all",
    train_mode: str = "multivariate",
):
    """
    Reads contamination metrics results for multiple samples in parallel.

    Args:
        results_dir (Path): Path to the results directory.
        correction_methods (list): List of correction methods to process.
        xenium_std_seurat_analysis_dir (Path): Path to the analysis directory.
        reference (str): Reference string.
        method (str): Method string.
        level (str): Level string.
        mixture_k (int): Mixture K value.
        num_samples (int): Number of samples.
        normalisation (str): The normalisation method used.
        layer (str): The layer used.
        radius (int): The radius used.
        n_splits (int): The number of splits used.
        n_permutations (int): The number of permutations used.
        n_repeats (int): The number of repeats used.
        top_n (int): The top N used.
        markers_mode (str): The markers mode used.
        cv_mode (str): The cross-validation mode used.
        scoring (str): The scoring metric used.
        evaluation (str): evaluation test to load, 'diffexpr' or 'logreg'.
        genes_name (None or str): Load results computed on a marker genes subset with folder name genes_name rather than all genes.
        train_mode (str): train mode to load for evaluation='logreg', 'multivariate' or 'univariate'.

    Returns:
        dict: A dictionary where keys are correction methods, and values are dictionaries mapping sample keys
              to DataFrames.
    """

    dfs = {}
    with ThreadPoolExecutor() as executor:
        futures = []
        for correction_method in correction_methods:
            for segmentation in xenium_std_seurat_analysis_dir.iterdir():
                if segmentation.stem in ["proseg_mode", "bats_normalised", "bats_expected"]:
                    continue
                for condition in segmentation.iterdir():
                    if ref_condition is not None and condition.name != ref_condition:
                        continue
                    if condition.stem == "melanoma":
                        if level == "Level2.1":
                            print_once("Using Level1 for melanoma\t")
                        if reference == "matched_reference_combo":
                            print_once("Using external_reference for melanoma\t")

                    for panel in condition.iterdir():
                        if ref_panel is not None and panel.name != ref_panel:
                            continue
                        for donor in panel.iterdir():
                            for sample in donor.iterdir():
                                futures.append(
                                    executor.submit(
                                        _read_contamination_metrics_results_sample,
                                        results_dir,
                                        reference,
                                        method,
                                        level,
                                        mixture_k,
                                        num_samples,
                                        correction_method,
                                        segmentation,
                                        condition,
                                        panel,
                                        donor,
                                        sample,
                                        normalisation,
                                        layer,
                                        radius,
                                        n_splits,
                                        n_permutations,
                                        n_repeats,
                                        top_n,
                                        markers_mode,
                                        cv_mode,
                                        scoring,
                                        evaluation,
                                        genes_name,
                                        train_mode,
                                    )
                                )

        with tqdm(total=len(futures), desc="Processing futures") as pbar:
            for future in as_completed(futures):
                result = future.result()
                if result:
                    correction_method, k, dfs_ = result

                    if correction_method not in dfs:
                        dfs[correction_method] = {}

                    dfs[correction_method][k] = dfs_
                pbar.update(1)  # Update progress bar

    # reorganize to have innermost keys (df types) as top level keys
    new_dfs = {}

    for correction_method, nested_dict in dfs.items():
        for k, sub_dict in nested_dict.items():
            for df_key, df_value in sub_dict.items():
                if df_key not in new_dfs:
                    new_dfs[df_key] = {}
                if correction_method not in new_dfs[df_key]:
                    new_dfs[df_key][correction_method] = {}

                new_dfs[df_key][correction_method][k] = df_value

    return new_dfs
