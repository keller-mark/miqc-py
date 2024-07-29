from stepmix.stepmix import StepMix
import numpy as np


def calculate_miqc(
    adata,
    model_type="linear",
    detected="total_genes_by_counts",
    subsets_mito_percent="pct_counts_mito",
    key_added="prob_compromised",
    inplace=True,
    random_state=None,
    n_init=10,
    max_iter=10000,
):
    """
    Function to fit a two-distribution mixture model on an AnnData object.

    :param adata: The AnnData object.
    :type adata: anndata.AnnData
    :param str model_type: The type of model to fit. Currently only 'linear' is supported.
    :param str detected: The key in adata.obs that contains the number of unique genes detected per cell.
    :param str subsets_mito_percent: The key in adata.obs that contains the percentage of mitochondrial genes per cell.
    :param str key_added: The key to add to adata.obs that will contain the posterior probability of each cell being compromised.
    :param bool inplace: Whether to modify adata in place. If False, a copy of adata.obs with a new column is returned.
    :param random_state: A random seed to use. By default, None.
    :type random_state: int or None
    :param int n_init: The number of initializations to run. By default, 10.
    :param int max_iter: The maximum number of iterations to run. By default, 10000.

    :returns: None if inplace else the new obs dataframe
    :rtype: None or pd.DataFrame

    """
    assert model_type == "linear", "Only model_type 'linear' is currently supported."

    Y = adata.obs[subsets_mito_percent].values
    X = np.expand_dims(adata.obs[detected].values, axis=1)

    model = StepMix(
        n_components=2,
        measurement="binary",
        structural="continuous",
        assignment="soft",
        n_init=n_init,
        max_iter=max_iter,
        random_state=random_state,
        progress_bar=0,
    )
    model.fit(X, Y)

    parameters = model.get_parameters()
    intercept0 = parameters["structural"]["means"][0, 0]
    intercept1 = parameters["structural"]["means"][1, 0]
    posterior_cols = model.predict_proba_class(X, Y)

    if intercept0 > intercept1:
        compromised_dist = 0
    else:
        compromised_dist = 1

    prob_compromised = posterior_cols[:, compromised_dist]
    
    if inplace:
        adata.obs[key_added] = prob_compromised
        return None
    else:
        obs_df = adata.obs.copy()
        obs_df[key_added] = prob_compromised
        return obs_df


def filter_miqc(adata, posterior_cutoff=0.75, key_added="prob_compromised", inplace=True):
    """
    Filter cells based on the posterior probability of being compromised.

    :param adata: The AnnData object.
    :type adata: anndata.AnnData
    :param float posterior_cutoff: The posterior probability cutoff. Cells with a probability lower than this will be kept. By default, 0.75.
    :param str key_added: The key in adata.obs that contains the posterior probability of each cell being compromised.
    :param bool inplace: Whether to modify adata in place. If False, a boolean mask (with True representing cells to keep) is returned.

    :returns: None if inplace else the boolean mask
    :rtype: None or np.ndarray
    """
    assert key_added in adata.obs.columns, f"Key '{key_added}' not found in adata.obs."
    
    prob_compromised = adata.obs[key_added].values
    keep = (prob_compromised <= posterior_cutoff)

    if inplace:
        adata._inplace_subset_obs(keep)
        return None
    else:
        return keep