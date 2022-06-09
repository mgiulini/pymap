"""Library to perform entropy-related tasks."""
import numpy as np
import pandas as pd
from scipy.stats import entropy


def calculate_entropies(cg_clust):
    """Calculate resolution and relevance from CG clustering."""
    if "records" not in cg_clust.columns:
        raise ValueError("cg_clust does not have a 'records' column")
    if (cg_clust["records"] < 0).any():
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
    gby = at_clust.groupby(at_clust.columns[mapping].tolist())
    omega_1 = gby.size().reset_index().rename(columns={0: 'omega_1'})
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

def calculate_pbar_indices(at_clust, cg_clust, nobs, mapping):
    #adding a column "omega_1
    at_clust.insert(at_clust.columns.size, "omega_1", 0, True)
    at_clust = at_clust.reset_index()
    mapping = mapping + 1
    
    #grouping the atomistic configs, writing the number of atomistic configs inside every cluster under the "omega_1" column and their indices under the "index" column
    omega_1 = at_clust.groupby(by = at_clust.columns[mapping].tolist()).agg({'omega_1': 'size', 'index': lambda x: tuple(x)}).reset_index()
    
    #calculating p_bar
    p_bar_r = cg_clust["records"]/nobs/omega_1["omega_1"]
    
    new_p_bar = pd.concat([omega_1,p_bar_r],axis=1) # to keep track of all the cg configurations
    new_p_bar.columns = pd.concat([pd.Series(at_clust.columns[mapping]), pd.Series("omega_1"), pd.Series("index"),pd.Series("p_bar")])
    
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
        eq_at_pbar = (at_clust.iloc[n, mapping] == new_p_bar.iloc[:, :-2])
        npw = np.where(eq_at_pbar.all(1))[0][0]
        mapping_entropy.append(pr[n]*np.log(pr[n]/new_p_bar.iloc[npw, -1]))
    tot_smap = sum(mapping_entropy)
    # infinite sampling mapping entropy
    return tot_smap

def calculate_smap_fast(mapping, pr, p_bar):
    """
    faster calculation of smap
    """
    mapping_entropy = []
    
    #for loop on the CG configs contained in "p_bar"
    indices = p_bar['index']
    for ii in range(p_bar.shape[0]):
        #for loop on the indices contained in "p_bar['index']"
        for n in indices[ii]:
            mapping_entropy.append(pr[n] * np.log(pr[n] / p_bar.iloc[ii,-1]))
    tot_smap = sum(mapping_entropy)
    
    return tot_smap
