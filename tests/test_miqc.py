import pytest

import numpy as np
import pandas as pd
from os.path import join
from anndata import AnnData

import miqc_py

def test_results_match_r():
    # Results should match R implementation
    posterior_cutoff = 0.75

    r_df = pd.read_csv(join("tests", "data", "zeisel_brain.r_results.csv"), index_col=0)
    r_prob_compromised = r_df["posterior_1"]
    num_keep_r = np.sum(r_prob_compromised <= posterior_cutoff)
    ids_keep_r = r_df.index[r_prob_compromised <= posterior_cutoff]

    adata = AnnData(X=None, obs=r_df)
    miqc_py.calculate_miqc(adata, detected="detected", subsets_mito_percent="subsets_mito_percent", inplace=True, random_state=123)

    num_keep_python = np.sum(miqc_py.filter_miqc(adata, posterior_cutoff=posterior_cutoff, inplace=False))
    ids_keep_python = adata.obs.index[adata.obs["prob_compromised"] <= posterior_cutoff]

    # We will tolerate 2% difference in the number of cells kept
    tolerance = 0.02 * adata.shape[0]
    assert abs(num_keep_python - num_keep_r) <= tolerance

    jaccard_similarity = (
        len(set(ids_keep_r).intersection(set(ids_keep_python)))
        / len(set(ids_keep_r).union(set(ids_keep_python)))
    )
    assert jaccard_similarity >= 0.97
    