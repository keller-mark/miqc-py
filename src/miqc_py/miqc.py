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
    assert key_added in adata.obs.columns, f"Key '{key_added}' not found in adata.obs."
    
    prob_compromised = adata.obs[key_added].values
    keep = (prob_compromised <= posterior_cutoff)

    if inplace:
        adata._inplace_subset_obs(keep)
        return None
    else:
        return keep