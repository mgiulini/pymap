"""Library to perform entropy-related tasks."""
import numpy as np
from scipy.stats import entropy
import pandas as pd


def calculate_entropies(cg_clust):
    """Calculate resolution and relevance from CG clustering."""
    if "records" not in cg_clust.columns:
        raise ValueError("cg_clust does not have a 'records' column")
    if (cg_clust["records"] < 0).any() is True:
        raise ValueError("'records' col in cg_clust contains negative values")
    # resolution
    hs = entropy(cg_clust["records"])
    # relevance
    ks = np.unique(cg_clust["records"], return_counts=True)
    pk = np.multiply(ks[0], ks[1])
    hk = entropy(pk)
    return hs, hk


def calculate_pbar(at_clust, cg_clust, nobs, mapping):
    """Calculate smeared pdf pbar."""
    omega_1 = at_clust.groupby(at_clust.columns[mapping].tolist()).size().reset_index().rename(columns={0: 'omega_1'})
    # smeared probability distribution
    p_bar_r = cg_clust["records"]/nobs/omega_1["omega_1"]
    # keep track of all the cg configurations
    new_p_bar = pd.concat([omega_1, p_bar_r], axis=1)
    new_p_bar.columns = pd.concat(
        [
            pd.Series(at_clust.columns[mapping]),
            pd.Series("omega_1"),
            pd.Series("p_bar")
        ]
    )
    return new_p_bar


def calculate_smap_inf(nat, ncg, hs_at, hs_cg, V):
    """Calculate the infinite-sampling mapping entropy."""
    if ncg > nat:
        raise ValueError(f"n {nat} < N {ncg}")
    if hs_at < hs_cg:
        raise ValueError(f"hs_at {hs_at} < hs_cg {hs_cg}")
    delta_s_conf = (nat - ncg)*np.log(V) - hs_at + hs_cg
    return delta_s_conf


def calculate_smap(at_clust, mapping, pr, new_p_bar):
    """Calculate the mapping entropy in its standard definition."""
    # state-wise mapping entropy
    mapping_entropy = []
    for n in range(at_clust.shape[0]):
        s = np.where((at_clust.iloc[n, mapping] == new_p_bar.iloc[:, :-2]).all(1) == True)[0][0]
        mapping_entropy.append(pr[n]*np.log(pr[n]/new_p_bar.iloc[s, -1]))
    tot_smap = sum(mapping_entropy)
    # infinite sampling mapping entropy
    return tot_smap
