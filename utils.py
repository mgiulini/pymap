import numpy as np 
import pandas as pd

def check_volume(dataframe,ncols):
    """
    calculate the volume associated to each degree of freedom
    """
    V = np.unique([len(np.unique(dataframe.iloc[:,n])) for n in range(ncols-1)], return_counts=True)
    if len(V[1]) != 1:
        print("non-constant volume detected, V = ", V, "choosing the most likely one")
        agmax = np.argmax(V[1])
        V = V[0][agmax]
    else:
        print("constant volume detected, V = ", V)
        V = V[0][0]
    print("V = ", V)
    return V

def validate_clust(clust):
    """
    do some validation of the at/cg dataframes
    """
    if clust.iloc[:,:-1].duplicated().any():
        raise Exception("Duplicate row detected in dataframe")

def validate_dataframe(dataframe):
    """
    do some checks on the original dataframe
        - you expect some duplicates, to give rise to non-trivial probabilities
    """
    if dataframe.duplicated().any() == False:
        raise Warning("No duplicate row detected")

def get_clust(dataframe,mapping):
    """
    compute the clustering
    """
    clust = dataframe.groupby(dataframe.columns[mapping].tolist()).size().reset_index().rename(columns={0:'records'})
    clust.columns = pd.concat([pd.Series(dataframe.columns[mapping]), pd.Series("records")])
    return clust

def calculate_entropies(cg_clust):
    """
    calculate resolution and relevance from the coarse-grained clustering
    """
    hs = entropy(cg_clust["records"])
    ks = np.unique(cg_clust["records"], return_counts=True)
    pk = np.multiply(ks[0], ks[1])
    hk = entropy(pk)
    return hs,hk

def calculate_smap(at_clust, cg_clust, V):
    """
    calculate the mapping entropy in its two definitions
    """
    omega_1 = at_clust.groupby(at_clust.columns[mapping].tolist()).size().reset_index().rename(columns={0:'records'})
    # smeared probability distribution
    p_bar_r = cg_clust["records"]/df.shape[0]/omega_1["records"]
    new_p_bar = pd.concat([omega_1,p_bar_r],axis=1) # to keep track of all the cg configurations
    # state-wise mapping entropy
    mapping_entropy = []
    for n in range(at_clust.shape[0]):
        s = np.where((at_clust.iloc[n,mapping] == new_p_bar.iloc[:,:-2]).all(1) == True)[0][0]
        mapping_entropy.append(pr[n]*np.log(pr[n]/new_p_bar.iloc[s,-1]))
    tot_smap = sum(mapping_entropy)
    # infinite sampling mapping entropy
    delta_s_conf = (len(at_clust.columns) - 1 - len(mapping))*np.log(V) - entropy(at_clust["records"]) + hs
    return smap, delta_s_conf
