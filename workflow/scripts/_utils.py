import sklearn
import numpy as np
import pandas as pd
import scipy
import gseapy
import anndata
from sklearn.model_selection import train_test_split, permutation_test_score, StratifiedKFold, GroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, BisectingKMeans
from sklearn.utils import check_random_state
from sklearn.utils.validation import check_array
from typing import Dict, List
from joblib import Parallel, delayed
import warnings

xenium_levels = ["segmentation", "condition", "panel", "donor", "sample"]


def get_knn_labels(
    adata, label_key, obsm=None, knnidx=None, n_neighbors=None, radius=None, n_jobs=-1, return_sparse_neighbors=False
):
    """
    Compute the k-nearest neighbors (KNN) labels for each cell in the anndata object.

    Parameters
    ----------
    adata : AnnData
        The anndata object containing the data to be labeled.
    label_key : str
        The key in adata.obs containing the true labels.
    obsm : str
        The key in adata.obsm containing the KNN data.
    knnidx : array-like, optional
        The indices of the nearest neighbors. If not provided, it is computed from the knn_key.
    n_neighbors : int, optional
        The number of nearest neighbors to consider. If not provided, it is computed from the knn_key.
    radius : float, optional
        The radius of the ball to consider neighbors in. If not provided, it is computed from the knn_key.
    n_jobs : int, optional
        The number of jobs to run in parallel. If not provided, it is set to -1, which means all available cores.
    return_sparse_neighbors : bool, optional
        Return neighbors matrix as sparse in addition to knndis, knnidx list

    Returns
    -------
    knnlabels: pd.DataFrame
        A pandas DataFrame containing the KNN labels for each cell in the anndata object.
    knndis, knnidx: np.array
        Distance and indices of nearest neighbors
    knn_graph:
        if return_neighbors_sparse, additional return of sparse kNN matrix
    """
    if knnidx is None:
        nn = sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors, radius=radius, n_jobs=n_jobs).fit(
            adata.obsm[obsm]
        )
        if n_neighbors is not None:
            knndis, knnidx = nn.kneighbors(adata.obsm[obsm])
        elif radius is not None:
            knndis, knnidx = nn.radius_neighbors(adata.obsm[obsm])

    df_dummies = pd.get_dummies(adata.obs[label_key])

    # Convert df_dummies to a numpy array for efficient indexing
    if not isinstance(df_dummies, np.ndarray):
        dummy_array = np.array(df_dummies)

    if knnidx.ndim == 2:
        # knnidx contains same length list of neighbor indices for each donor
        knnlabels = dummy_array[knnidx].sum(1)
    else:
        # Initialize an empty list to store the summed labels
        knnlabels = []

        # Loop over each row in knnidx
        for neighbors in knnidx:
            # Get the one-hot encoded labels for the current neighbors
            neighbor_labels = dummy_array[neighbors]

            # Sum the labels across the neighbors (axis=0 sums column-wise)
            summed_labels = neighbor_labels.sum(axis=0)

            # Append the summed labels to the list
            knnlabels.append(summed_labels)

        # Convert the list back to a numpy array (optional)
        knnlabels = np.array(knnlabels)

    knnlabels = pd.DataFrame(knnlabels, index=adata.obs.index, columns=df_dummies.columns)

    if return_sparse_neighbors:
        if n_neighbors is not None:
            knn_graph = nn.kneighbors_graph()
        elif radius is not None:
            knn_graph = nn.radius_neighbors_graph()

        return knnlabels, knndis, knnidx, knn_graph
    else:
        return knnlabels, knndis, knnidx


def logit_pvalue(model, x):
    """Calculate z-scores for scikit-learn LogisticRegression.
    parameters:
        model: fitted sklearn.linear_model.LogisticRegression with intercept and large C
        x:     matrix on which the model was fit
    This function uses asymtptics for maximum likelihood estimates.
    """
    p = model.predict_proba(x)
    n = len(p)
    m = len(model.coef_[0]) + 1
    coefs = np.concatenate([model.intercept_, model.coef_[0]])
    x_full = np.matrix(np.insert(np.array(x), 0, 1, axis=1))
    ans = np.zeros((m, m))
    for i in range(n):
        ans = ans + np.dot(np.transpose(x_full[i, :]), x_full[i, :]) * p[i, 1] * p[i, 0]
    vcov = np.linalg.inv(np.matrix(ans))
    se = np.sqrt(np.diag(vcov))
    t = coefs / se
    p = (1 - scipy.stats.norm.cdf(abs(t))) * 2
    return p


def get_marker_rank_significance(rnk, gene_set, top_n=None):
    """
    Calculate the significance of marker ranks using pre-ranked GSEA and hypergeometric test.
    rnk must contain all features, its length is assumed to be the total number N for hypergoemetric testing

    Parameters:
    - rnk (pd.Series): index with feature names and values representing feature scores with decreasing ranks.
    - adata (AnnData): Anndata object containing the gene expression data.
    - gene_set (list): List of genes to test for enrichment.
    - top_n (optional, int): Number of top-ranked genes to consider for the hypergeometric test.

    Returns:
    - markers_rank_significance (pd.DataFrame): DataFrame containing the significance results.
    """

    if np.sum(rnk > 0.0) == 0:
        print("Warning: all zero or all negative feature score vector. Returning empty dataframe.")
        return pd.DataFrame(
            np.nan,
            index=[0],
            columns=[
                "Name",
                "Term",
                "ES",
                "NES",
                "NOM p-val",
                "FDR q-val",
                "FWER p-val",
                "Tag %",
                "Gene %",
                "Lead_genes",
                "hypergeometric_pvalue",
            ],
        )
    # Calculate marker rank significance from pre-ranked GSEA
    markers_rank_significance = gseapy.prerank(
        rnk=rnk,
        gene_sets=[{"gene_set": gene_set}],
        min_size=0,
    ).res2d

    # Calculate marker rank significance from hypergeometric test
    if top_n is not None:
        N = len(rnk)  # Total genes in ranked list
        K = len(gene_set)  # Genes in the pathway/set of interest
        x = np.isin(rnk.index[:top_n], gene_set).sum()  # Overlapping genes in top n

        # Add hypergeometric p-value to the results
        markers_rank_significance["hypergeometric_pvalue"] = scipy.stats.hypergeom.sf(x - 1, N, K, top_n)
        markers_rank_significance[f"n_hits_{top_n=}"] = x

    return markers_rank_significance


def _logreg_univariate(X, y, feature_names, compute_pvalue=True, **init_params):
    def _fit_univariate(x_i, y, feature_name, compute_pvalue=True, **init_params):
        model = LogisticRegression(penalty=None, **init_params)
        model.fit(x_i, y)
        if compute_pvalue:
            try:
                pvalue = logit_pvalue(model, x_i)[1]
            except:
                pvalue = np.nan
        output = {
            "feature_name": feature_name,
            "importances": model.coef_[0][0],
            "intercepts": model.intercept_[0],
            "pvalues": pvalue,
        }
        return output

    df_importances = pd.DataFrame(
        Parallel(n_jobs=-1)(
            delayed(_fit_univariate)(
                X[:, [i]],
                y,
                feature_names[i],
                compute_pvalue=compute_pvalue,
                **init_params,
            )
            for i in range(X.shape[1])
        )
    ).set_index("feature_name")

    return df_importances


def logreg(
    X,
    y,
    feature_names=None,
    scoring="precision",
    test_size=0.2,
    n_splits=5,
    n_permutations=30,
    n_repeats=5,
    max_iter=100,
    random_state=0,
    importance_mode="coef",
    class_weight="balanced",
    cv_mode="spatial",
    spatial_coords=None,
    accept_partial_cv=False,
    scale=True,
    train_mode="multivariate",
):
    """
    Perform logistic regression with permutation test and compute feature importances.

    Parameters:
    - X (array-like): Input data for model training/testing.
    - y (vector-like): Input vector of labels for model training/testing.
    - feature_names (vector-like): Names of X features (optional).
    - scoring (str): Scoring metric for the permutation test (e.g., 'f1', 'accuracy').
    - test_size (float): Proportion of data to use for testing.
    - max_iter (int): Maximum number of iterations for the logistic regression model.
    - n_splits (int): Number of splits for cross-validation.
    - n_permutations (int): Number of permutations for the permutation test.
    - n_repeats (int): Number of repeats for the permutation importance calculation.
    - random_state (int): Random seed for reproducibility.
    - importance_mode (str): Mode for feature importance calculation ('permutation' or 'coef').
    - class_weight (str): Class weight for the logistic regression model ('balanced' or None or dict).
    - cv_mode (str): Cross-validation mode ('stratified' or 'spatial').
    - spatial_coords (array-like): Spatial coordinates for spatial cross-validation.
    - accept_partial_cv (bool): Whether to run the function despite at least one cross validation split only having one class (for spatial cv).
    - scale (bool): Whether to scale the input data.
    - train_mode (str): Training mode ('univariate' or 'multivariate').

    Returns:
    - df_permutations (pd.DataFrame): Summary of permutation test results. For train_mode='univariate', df_permutations will be empty.
    - df_importances (pd.DataFrame): Feature importances from permutation importance.
    """

    if feature_names is None:
        feature_names = np.arange(X.shape[1])

    if scipy.sparse.issparse(X):
        X = X.toarray()
    if scale:
        X = StandardScaler().fit_transform(X)

    # Split data into cross validation sets
    if cv_mode == "stratified":
        cv = StratifiedKFold(n_splits=n_splits, shuffle=False)
    elif cv_mode == "spatial":
        if spatial_coords is None:
            raise ValueError("spatial_coords must be provided when cv_mode is 'spatial'")
        cv = list(SpatialClusterGroupKFold(algorithm="bisectingkmeans", n_splits=n_splits).split(spatial_coords, y))

        # check that all splits have more than one class
        single_class_splits = [(len(np.unique(y[train])) == 1 or len(np.unique(y[test])) == 1) for (train, test) in cv]
        if any(single_class_splits):
            if accept_partial_cv:
                cv = [cv[i] for i in range(len(cv)) if not single_class_splits[i]]
                if len(cv) == 0:
                    warnings.warn("All splits have only one class. Aborting.")
                    return None, None
                else:
                    print("Some splits have only one class. Running on", len(cv), "splits.")
            else:
                warnings.warn("Some splits have only one class. Aborting.")
            return None, None
    else:
        raise ValueError("cv_mode must be 'stratified' or 'spatial'")

    if train_mode == "multivariate":
        # Initialize logistic regression model
        model = LogisticRegression(max_iter=max_iter, class_weight=class_weight)

        # Empirical p-value calculation using permutation test
        try:
            score, perm_scores, p_value = permutation_test_score(
                model, X, y, scoring=scoring, n_permutations=n_permutations, cv=cv, n_jobs=-1, verbose=1
            )
        except ValueError:
            score = np.nan
            perm_scores = np.array([np.nan])
            p_value = np.nan

        # Summarize permutation test results
        df_permutations = pd.DataFrame(
            [[score, perm_scores.mean(), perm_scores.std(), p_value]],
            columns=[f"{scoring}_score", f"perm_mean{scoring}_score", f"perm_std{scoring}_score", "p_value"],
        )
        df_permutations["effect_size"] = (
            df_permutations[f"{scoring}_score"] - df_permutations[f"perm_mean{scoring}_score"]
        ) / df_permutations[f"perm_std{scoring}_score"]

        # Fit the model and compute feature importances from permutations
        if importance_mode == "permutation":
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=test_size, stratify=y, random_state=random_state
            )
            model.fit(X_train, y_train)
            importances = permutation_importance(
                model,
                pd.DataFrame.sparse.from_spmatrix(X_test),
                y_test,
                scoring=scoring,
                n_repeats=n_repeats,
                n_jobs=-1,
            )
            importances.pop("importances")
            df_importances = pd.DataFrame(importances, index=feature_names).sort_values(
                "importances_mean", ascending=False
            )

        elif importance_mode == "coef":
            model.fit(X, y)
            # Feature importances from model coefs
            # cv_results = cross_validate(model,X,y,return_estimator=True, scoring=scoring, n_jobs=-1)
            # importances = np.std(X, axis=0) * np.vstack([m.coef_[0] for m in cv_results["estimator"]])
            # importances = StandardScaler(with_mean=False).fit(X).scale_ * model.coef_[0]
            importances = model.coef_[0]
            df_importances = pd.DataFrame(importances, index=feature_names, columns=["importances"]).sort_values(
                "importances", ascending=False
            )

            # coef pvalues from formula
            # df_importances["pvalues"] = logit_pvalue(model, X.toarray())[1:]
        else:
            raise ValueError("Importance mode must be 'permutation' or 'coef'")

    elif train_mode == "univariate":
        df_importances = _logreg_univariate(
            X, y, feature_names, compute_pvalue=True, max_iter=max_iter, class_weight=class_weight
        ).sort_values("importances", ascending=False)
        df_permutations = pd.DataFrame()

    return df_permutations, df_importances


def get_mean_cell_identity_score_markers(
    ads: Dict[str, Dict[str, anndata.AnnData]],
    df_markers: pd.DataFrame,
    correction_methods: List[str],
    labels_key: str,
    columns=None,
) -> pd.DataFrame:
    """
    Calculates the mean number of genes (`n_genes`) for all unique cell types (`cti`) found
    in the `labels_key` column across different AnnData objects and correction methods.

    Args:
        ads (Dict[str, Dict[str, anndata.AnnData]]): A nested dictionary structure containing
            AnnData objects. The outer dictionary's keys are correction methods, and the
            inner dictionary's keys are sample identifiers. Each inner dictionary value
            is an AnnData object. Assumes `ads[correction_method][k]` is an AnnData object.
        df_markers (pd.DataFrame): A DataFrame containing cell type markers.
        correction_methods (List[str]): A list of correction methods to iterate through.
        labels_key (str): The key in `ad.obs` that contains cell type labels.
        columns (list): optional column names for the returned DataFrame
    Returns:
        pd.DataFrame: A DataFrame where each row represents a unique combination of sample
            identifier, correction method, and cell type. The columns include sample identifiers,
            correction method, cell type (`cti`), and the calculated mean number of genes
            ('n_genes').
    """

    data = []
    for correction_method in correction_methods:
        for k, ad in ads[correction_method].items():
            if ad is not None:
                unique_ctis = ad.obs[labels_key].unique()

                for cti in unique_ctis:
                    cti_markers = df_markers.query("cell_type == @cti")["gene_name"].tolist()
                    cti_markers = [g for g in cti_markers if g in ad.var_names]
                    if len(cti_markers):
                        n_genes_cti_markers = (ad[ad.obs[labels_key] == cti, cti_markers].X > 0).sum(1).A1
                        mean_n_genes = np.mean(n_genes_cti_markers)
                        data.append((*k, correction_method, cti, mean_n_genes))  # Append cell type to the data

    # Create the DataFrame from the collected data
    df = pd.DataFrame(data, columns=columns)

    if columns is not None:
        df.columns = columns
    return df


def get_mean_cell_identity_score_scrna(
    ads: Dict[str, Dict[str, anndata.AnnData]],
    ads_scrnaseq: Dict[str, anndata.AnnData],
    df_markers: pd.DataFrame,
    correction_methods: List[str],
    labels_key: str,
    columns=None,
) -> pd.DataFrame:
    """
    Calculates the mean number of genes (`n_genes`) for all unique cell types (`cti`) found
    in the `labels_key` column across different AnnData objects and correction methods.

    Args:
        ads (Dict[str, Dict[str, anndata.AnnData]]): A nested dictionary structure containing
            AnnData objects. The outer dictionary's keys are correction methods, and the
            inner dictionary's keys are sample identifiers. Each inner dictionary value
            is an AnnData object. Assumes `ads[correction_method][k]` is an AnnData object.
        ads_scrnaseq Dict[str, anndata.AnnData]: A dictionary structure containing
            AnnData objects. The dictionary's keys are tissue identifiers corresponding to keys in `ads`.
        correction_methods (List[str]): A list of correction methods to iterate through.
        labels_key (str): The key in `ad.obs` that contains cell type labels.
        columns (list): optional column names for the returned DataFrame
    Returns:
        pd.DataFrame: A DataFrame where each row represents a unique combination of sample
            identifier, correction method, and cell type. The columns include sample identifiers,
            correction method, cell type (`cti`), and the calculated mean number of genes
            ('n_genes').
    """

    data = []
    for correction_method in correction_methods:
        for k, ad in ads[correction_method].items():
            if ad is not None:
                unique_ctis = ad.obs[labels_key].unique()

                for cti in unique_ctis:
                    cti_markers = df_markers.query("cell_type == @cti")["gene_name"].tolist()
                    cti_markers = [g for g in cti_markers if g in ad.var_names]
                    if len(cti_markers):
                        n_genes_cti_markers = (ad[ad.obs[labels_key] == cti, cti_markers].X > 0).sum(1).A1
                        mean_n_genes = np.mean(n_genes_cti_markers)
                        data.append((*k, correction_method, cti, mean_n_genes))  # Append cell type to the data

    # Create the DataFrame from the collected data
    df = pd.DataFrame(data, columns=columns)

    if columns is not None:
        df.columns = columns
    return df


def get_cosine_similarity_score(pbs_xenium, pb_scrna, labels_key, correction_methods, columns=None):
    """
    Calculates the cosine similarity score between each cell type's expression profile and its
    corresponding pseudobulk expression profile.

    Args:
        pbs_xenium (Dict[str, Dict[str, anndata.AnnData]]): A nested dictionary structure containing
            AnnData objects. The outer dictionary's keys are correction methods, and the
            inner dictionary's keys are sample identifiers. Each inner dictionary value
            is an AnnData object. Assumes `ads[correction_method][k]` is an AnnData object.
        pb_scrna (anndata.AnnData): Anndata containing pseudobulk expression
            profiles.
        labels_key (str): The key in `pb_xenium.obs` that contains cell type labels.
        correction_methods (List[str]): A list of correction methods to iterate through.
        columns (list): optional column names for the returned DataFrame

    Returns:
        pd.DataFrame: A DataFrame where each row represents a unique combination of sample
            identifier, correction method, and cell type. The columns include sample identifiers,
            correction method, cell type (`cti`), and the calculated cosine similarity score
            ('cosine_sim').

    Notes:
        If a cell type is not found in the pseudobulk expression profile, a warning is printed.
    """
    data = []
    for correction_method in correction_methods:
        for k, pb_xenium in pbs_xenium[correction_method].items():
            print(correction_method, k)
            if pb_xenium is not None:
                common_genes = np.intersect1d(pb_xenium.var_names, pb_scrna.var_names)
                unique_ctis = pb_xenium.obs_names.unique()

                for cti in unique_ctis:
                    if cti in pb_scrna.obs_names:
                        x = pb_xenium[pb_xenium.obs_names == cti, common_genes].X.toarray().squeeze()
                        y = pb_scrna[pb_scrna.obs_names == cti, common_genes].X.toarray().squeeze()

                        cosine_sim = 1 - scipy.spatial.distance.cosine(x, y)
                        data.append((*k, correction_method, cti, cosine_sim))  # Append cell type to the data
                    else:
                        print("Warning: cell type", cti, "not found in pseudobulk")

    # Create the DataFrame from the collected data
    df = pd.DataFrame(data, columns=columns)

    if columns is not None:
        df.columns = columns

    return df


def rename_methods(
    df,
    segmentation_column="segmentation",
    correction_column="correction_method",
    rename_segmentation={
        "10x_mm_0um": "MM 0µm",
        "10x_mm_5um": "MM",
        "10x_mm_15um": "MM 15µm",
        "10x_0um": "0µm",
        "10x_5um": "5µm",
        "10x_15um": "15µm",
        "baysor": "Baysor",
        "proseg_expected": "ProSeg",
        "proseg_mode": "ProSeg mode",
        "segger": "Segger",
    },
    rename_correction={
        "resolvi": "ResolVI",
        "resolvi_panel_use_batch=True": "ResolVI",
        "resolvi_panel_use_batch=False": "ResolVI",
        "resolvi_supervised": "ResolVI supervised",
        "resolvi_panel_supervised_use_batch=True": "ResolVI supervised",
        "resolvi_panel_supervised_use_batch=False": "ResolVI supervised",
        "ovrlpy_correction_signal_integrity_threshold=0.5": "ovrlpy 0.5",
        "ovrlpy_correction_signal_integrity_threshold=0.7": "ovrlpy 0.7",
        "ovrlpy_0.5": "ovrlpy 0.5",
        "ovrlpy_0.7": "ovrlpy 0.7",
        "split_fully_purified": "SPLIT",
        "split": "SPLIT",
    },
):
    if segmentation_column is not None:
        df[segmentation_column] = df[segmentation_column].replace(rename_segmentation)
    if correction_column is not None:
        df[correction_column] = df[correction_column].replace(rename_correction)


def extract_info_cell_type_pair(series, cti, ctj, flag):
    """Add info columns for number of cti cells being neighbor or not to ctj cells"""
    s = series.map(
        lambda df_: df_.query(f"label_key == '{cti}' and variable == 'has_{ctj}_neighbor' and value == {flag}")[
            "count"
        ].squeeze()
    )
    s = s.map(lambda x: 0 if isinstance(x, pd.Series) else x)  # replace empty series with 0

    return s


def get_df_summary_stats_plot(dfs, plot_metric="n_cells"):
    """
    Generate a summary DataFrame from a dictionary of summary statistics.

    Parameters:
    - dfs (dict): A dictionary containing summary statistics with correction methods as keys, as read by readwrite.read_diffexpr_results_samples.
    - plot_metric (str, optional): The metric to plot. Defaults to "n_cells". Can also be "df_has_neighbor_counts"
      or any other metric present in the summary statistics.

    Returns:
    - pd.DataFrame: A DataFrame containing the summary statistics, with columns based on the `xenium_levels`,
      correction method, and the specified plot metric. If the plot metric is a dictionary, additional columns are
      included for the dictionary keys and values.
    """

    data = []  # List to store data for the DataFrame

    for correction_method, k_values in dfs["summary_stats"].items():  # Iterate through correction methods
        for k, values in k_values.items():  # Iterate through keys k
            v = values[plot_metric]
            if plot_metric == "df_has_neighbor_counts":
                v = pd.DataFrame(v)
                row = list(k) + [correction_method, v]
                data.append(row)

            elif isinstance(v, dict):
                for k1, values1 in v.items():
                    row = list(k) + [correction_method, k1, values1]
                    data.append(row)

            else:
                row = list(k) + [correction_method, v]  # Create a row of data
                data.append(row)

    if isinstance(v, dict):
        columns = xenium_levels + ["correction_method", plot_metric + "_key", plot_metric + "_value"]
    else:
        columns = xenium_levels + ["correction_method", plot_metric]

    df = pd.DataFrame(data, columns=columns)
    rename_methods(df)
    return df


def get_df_permutations_logreg_plot(df_permutations_logreg, correction_methods):
    df = {}
    for correction_method in correction_methods:
        for k, v in df_permutations_logreg[correction_method].items():
            for cti, ctj, _ in v.index:
                df[(correction_method, *k, cti, ctj)] = v.loc[(cti, ctj, 0)]
    df = pd.DataFrame(df).T.reset_index()
    df.columns = ["correction_method"] + xenium_levels + ["cti", "ctj"] + df.columns[-5:].tolist()
    rename_methods(df)
    return df


def get_df_marker_rank_significance_plot(
    dfs_marker_rank_significance,
    rank_metric,
    plot_metric,
    correction_methods,
    use_precomputed,
    n=20,
):
    """
    Generate a DataFrame from a dictionary of DataFrames, to be used for plotting.

    Parameters:
    - dfs_marker_rank_significance (dict): A dictionary containing DataFrames with correction methods as keys, as read by readwrite.read_diffexpr_results_samples.
    - rank_metric (str): The ranking metric to use for the plot.
    - plot_metric (str): The metric to plot.
    - correction_methods (list): A list of correction methods to include in the plot.
    - use_precomputed (bool): Whether to use precomputed results (True) or not (False).
    - n (int): The number of ctj markers used for score computation.

    Returns:
    - pd.DataFrame: A DataFrame containing the data for the plot, with columns for the correction method, xenium levels, cell types, and plot metric.
    """
    df = {}
    for correction_method in correction_methods:
        if correction_method not in dfs_marker_rank_significance.keys():
            continue
        for k, v in dfs_marker_rank_significance[correction_method].items():
            if use_precomputed:
                rank_metric_ = (
                    f"{rank_metric}_{n=}" if correction_method == "raw" else f"{rank_metric}_precomputed_{n=}"
                )
            else:
                rank_metric_ = f"{rank_metric}_{n=}"
            df[(correction_method, *k)] = v.loc[rank_metric_, v.columns.get_level_values(2) == plot_metric]
    df = pd.concat(df).reset_index()
    df.columns = ["correction_method"] + xenium_levels + ["cti", "ctj", "plot_metric", plot_metric]
    rename_methods(df)
    return df


def get_df_ctj_marker_genes(
    dfs_ctj_marker_genes,
    correction_methods,
):
    df = {}
    for correction_method in correction_methods:
        for k, v in dfs_ctj_marker_genes[correction_method].items():
            for ctj in v.columns:
                df[(correction_method, *k, ctj)] = v[ctj]

    df = pd.concat(df).reset_index().drop("level_7", axis=1)
    df.columns = ["correction_method"] + xenium_levels + ["cell_type", "gene"]
    rename_methods(df)
    return df


def get_df_diffexpr_cti_ctj(dfs_diffexpr, cti, ctj, ref_panel, correction_methods, rank_metric):
    df = {}
    for correction_method in correction_methods:
        for k, v in dfs_diffexpr[correction_method].items():
            if k[2] != ref_panel:
                continue

            if (cti, ctj) in v.index:
                df[(correction_method, *k, cti, ctj)] = v.loc[(cti, ctj)].set_index("names").loc[:, rank_metric]

    df = pd.DataFrame(df).T.reset_index()
    df.columns = ["correction_method"] + xenium_levels + ["cti", "ctj"] + df.columns[8:].tolist()
    rename_methods(df)
    return df


def get_df_importances_cti_ctj(
    df_importances_logreg, cti, ctj, ref_panel, correction_methods, rank_metric="importances"
):
    df = {}
    for correction_method in correction_methods:
        for k, v in df_importances_logreg[correction_method].items():
            if k[2] != ref_panel:
                continue
            if (cti, ctj) in v.index:
                df[(correction_method, *k, cti, ctj)] = v.loc[(cti, ctj), rank_metric]

    df = pd.DataFrame(df).T.reset_index()
    df.columns = ["correction_method"] + xenium_levels + ["cti", "ctj"] + df.columns[8:].tolist()
    rename_methods(df)
    return df


def pseudobulk(ad, key, mode="sum"):
    """
    Generate pseudo-bulk RNA-seq data from single-cell RNA-seq datasets.

    Parameters
    ----------
    ad : Anndata
    key : str
        The key in `ad.obs` used to identify cell type labels.
    mode : str, optional
        Whether to sum or mean the expression values across cells of the same cell type.
        Options are "sum" (default) or "mean".

    Returns
    -------
    Anndata
        Contains pseudo-bulk expression profiles for each cell type,
        with cell types as columns and genes as rows.
    """

    def _agg_obs(x):
        val = x.mode()
        if len(val):
            return val[0]
        else:
            return None

    agg = {}
    for c in ad.obs[key].unique():
        if pd.isna(c):
            idx = ad.obs[key].isna()
        else:
            idx = ad.obs[key] == c
        if mode == "sum":
            agg[c] = np.asarray(ad[idx].X.sum(0)).squeeze()
        elif mode == "mean":
            agg[c] = np.asarray(ad[idx].X.mean(0)).squeeze()
    ad_states = anndata.AnnData(pd.DataFrame(agg).T)
    ad_states.var_names = ad.var_names
    ad_states.obs = ad.obs.join(ad.obs.groupby(key).agg(_agg_obs))
    return ad_states


def get_expression_percent_per_celltype(
    adata: anndata.AnnData,
    label_key: str,
    layer: str = None,
) -> pd.DataFrame:
    """
    Calculates the percentage of cells expressing each gene within each cell type.

    Expression is defined as having a count > 0 in the specified layer
    (or adata.X if layer is None).

    Args:
        adata: An AnnData object containing gene expression data and cell labels.
        label_key: The key in adata.obs where cell type labels are stored
                    (e.g., 'cell_type', 'leiden').
        layer: Optional name of the layer in adata to use for expression counts.
               If None, uses adata.X.
        gene_symbols: Optional key in adata.var where gene symbols are stored.
                      If None, uses adata.var_names (often Ensembl IDs).

    Returns:
        A pandas DataFrame where rows are genes (using gene_symbols if provided,
        otherwise var_names), columns are cell types, and values are the
        percentage [0-100] of cells expressing the gene within that cell type.

    """

    # Get gene identifiers
    gene_names = adata.var_names.tolist()
    gene_index_name = "gene"  # Default index name

    # Get unique cell type labels
    cell_types = adata.obs[label_key].unique()
    # Filter out any potential NaN/missing labels if they exist
    cell_types = [ct for ct in cell_types if pd.notna(ct)]

    # Choose the data matrix
    if layer:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers.")
        expression_matrix = adata.layers[layer]
    else:
        expression_matrix = adata.X

    # Initialize dictionary to store results (percentage vectors for each cell type)
    results = {}

    # Iterate through each cell type
    for cell_type in cell_types:
        cell_mask = (adata.obs[label_key] == cell_type).values

        # Get the number of cells in this type
        n_cells_in_type = np.sum(cell_mask)  # More direct than subsetting shape

        if n_cells_in_type == 0:
            # Handle cases where a category might exist but have 0 cells
            # Assign 0% expression for all genes for this type.
            percent_expressing = np.zeros(adata.shape[1])
        else:
            # Subset the expression matrix rows (cells) for the current cell type
            subset_matrix = expression_matrix[cell_mask, :]

            # Count cells with expression > 0 for each gene (column)
            if scipy.sparse.issparse(subset_matrix):
                expression_counts = subset_matrix.getnnz(axis=0)
            else:
                expression_counts = np.sum(np.asarray(subset_matrix) > 0, axis=0)
                expression_counts = np.asarray(expression_counts).flatten()

            # Calculate percentage
            percent_expressing = (expression_counts / n_cells_in_type) * 100

        # Store results for this cell type
        results[cell_type] = percent_expressing

    # Create the final DataFrame
    percent_df = pd.DataFrame(results, index=gene_names)
    percent_df.index.name = gene_index_name
    percent_df.columns.name = label_key  # Set the name of the columns axis

    print("Calculation complete.")
    return percent_df


class SpatialClusterGroupKFold:
    """
    Spatial K-Folds cross-validator using coordinate clustering and GroupKFold.

    Creates folds by first clustering spatial coordinates (e.g., lon/lat)
    using KMeans or BisectingKMeans. The resulting cluster labels are then
    used as group IDs for sklearn's GroupKFold, ensuring that all points
    within the same spatial cluster stay in the same fold (train or test).

    Parameters
    ----------
    n_splits : int, default=5
        Number of folds. Must be at least 2. This will also be the number
        of clusters generated.

    algorithm : {'kmeans', 'bisectingkmeans'}, default='bisectingkmeans'
        The clustering algorithm to use on the coordinates.

    random_state : int, RandomState instance or None, default=None
        Controls the randomness of the clustering. Pass an int
        for reproducible output. Note: GroupKFold itself is deterministic.

    **kwargs : dict
        Additional keyword arguments passed directly to the underlying
        clustering algorithm (`KMeans` or `BisectingKMeans`).
    """

    def __init__(self, n_splits=5, *, algorithm="bisectingkmeans", random_state=None, **kwargs):
        if n_splits < 2:
            raise ValueError("n_splits must be at least 2.")
        if algorithm not in ["kmeans", "bisectingkmeans"]:
            raise ValueError("algorithm must be 'kmeans' or 'bisectingkmeans'")

        self.n_splits = n_splits
        self.algorithm = algorithm
        self.random_state = random_state
        self.kwargs = kwargs

    def _get_cluster_labels(self, X_coords):
        """Performs clustering and returns cluster labels."""
        X_coords = check_array(X_coords, accept_sparse=False, ensure_2d=True, dtype="numeric")
        rng = check_random_state(self.random_state)
        # Pass random_state only if it's not None for reproducibility
        clustering_random_state = rng if self.random_state is not None else None

        if self.algorithm == "kmeans":
            # Use n_init='auto' to avoid future warning with default
            n_init = self.kwargs.pop("n_init", "auto")
            clustering_model = KMeans(
                n_clusters=self.n_splits, random_state=clustering_random_state, n_init=n_init, **self.kwargs
            )
        elif self.algorithm == "bisectingkmeans":
            clustering_model = BisectingKMeans(
                n_clusters=self.n_splits, random_state=clustering_random_state, **self.kwargs
            )
        else:
            raise ValueError("algorithm must be 'kmeans' or 'bisectingkmeans'")
        try:
            # Fit and predict cluster labels
            cluster_labels = clustering_model.fit_predict(X_coords)
            n_clusters_found = len(np.unique(cluster_labels))
            if n_clusters_found < self.n_splits:
                warnings.warn(
                    f"Clustering found only {n_clusters_found} clusters, "
                    f"less than n_splits={self.n_splits}. GroupKFold might behave "
                    "unexpectedly or yield fewer folds if some group IDs are missing "
                    "or if a cluster is empty.",
                    UserWarning,
                )
                # Note: GroupKFold might raise an error if a fold ends up empty.
            return cluster_labels
        except Exception as e:
            raise RuntimeError(f"Clustering failed: {e}") from e

    def split(self, X, y=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features_coords)
            Spatial coordinates used for clustering to determine groups.
            **Important**: This method expects coordinates here.

        y : object, optional
            Always ignored in this implementation, exists for compatibility.

        Yields
        ------
        train : ndarray
            The training set indices for that split.
        test : ndarray
            The testing set indices for that split.
        """
        X_coords = X  # Assume X passed is the coordinates
        self.cluster_labels_ = self._get_cluster_labels(X_coords)
        gkf = GroupKFold(n_splits=self.n_splits)
        return gkf.split(X_coords, y, groups=self.cluster_labels_)
