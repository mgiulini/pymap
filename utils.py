from typing import Optional
from pathlib import Path
import numpy as np 
import pandas as pd
from numpy.random import default_rng
from scipy.stats import entropy

def check_output_file(cleaned_pars):
    """If the output_filename already exists, block the execution"""
    output_path = Path(cleaned_pars["output_filename"])
    if output_path.is_file():
        raise Exception(f"Output file {cleaned_pars['output_filename']} already existing. Aborting")
    return

def system_parameters_setup(parfile):
    """Sets up the parameters."""
    pars = read_parfile(parfile)
    check_mandatory_parameters(pars)
    cleaned_pars = check_optional_parameters(pars)
    print(f"Cleaned Parameters {cleaned_pars}")
    check_output_file(cleaned_pars)
    return cleaned_pars

def read_parfile(parfile_string):
    # receives the name of the parameters file in input
    # parses it and gives back a dictionary
    ###################################################
    parfile = Path(parfile_string)
    print("reading parameters file ", parfile)
    if parfile.is_file == False:
        raise Exception("Parameter file not existing. Aborting.")
    parameters = {}
    with open(parfile, "r") as pars:
        for ln in pars:
            if ln.startswith("#") == False:
                split_list = ln.split()
                if len(split_list) != 2:
                    raise Exception(f"badly formatted parameter line\n{ln}")
                else:
                    par_name = split_list[0]
                    par_value = split_list[1]
                    parameters[par_name] = par_value
                    print("parameter ", par_name, " = ", par_value)
    return parameters

def check_mandatory_parameters(parameters):
    """Check the existence of mandatory parameters."""
    mandatory_keys = [
        "input_filename",
        "output_filename"
    ]

    observed_pars = parameters.keys()

    for par in mandatory_keys:
        if par not in observed_pars:
            raise Exception(f"missing parameter {par}")
    
    return

def check_optional_parameters(parameters):
    """Check the existence of optional parameters. Add their default value if absent"""
    optional_keys = {
        "max_binom" : ["integer", 100000],
    }
    
    observed_pars = parameters.keys()

    for optk in optional_keys.keys():
        if optk not in observed_pars:
            parameters[optk] = optional_keys[optk][1]
        else:
            if optional_keys[optk][0] == "integer":
                parameters[optk] = int(parameters[optk])
            if optional_keys[optk][0] == "float":
                parameters[optk] = float(parameters[optk])
        
    return parameters


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
    validate_clust(clust)
    return clust

def calculate_entropies(cg_clust):
    """
    calculate resolution and relevance from the coarse-grained clustering
    """
    if "records" not in cg_clust.columns:
        raise ValueError("cg_clust does not have a 'records' column")
    if (cg_clust["records"] < 0).any() == True:
        raise ValueError("'records' column in cg_clust contains negative values")
    # resolution        
    hs = entropy(cg_clust["records"])
    # relevance
    ks = np.unique(cg_clust["records"], return_counts=True)
    pk = np.multiply(ks[0], ks[1])
    hk = entropy(pk)
    return hs,hk

def calculate_pbar(at_clust, cg_clust, nobs, mapping):
    omega_1 = at_clust.groupby(at_clust.columns[mapping].tolist()).size().reset_index().rename(columns={0:'omega_1'})
    # smeared probability distribution
    p_bar_r = cg_clust["records"]/nobs/omega_1["omega_1"]
    #print("p_bar_r\n",p_bar_r)
    new_p_bar = pd.concat([omega_1,p_bar_r],axis=1) # to keep track of all the cg configurations
    new_p_bar.columns = pd.concat([pd.Series(at_clust.columns[mapping]), pd.Series("omega_1"), pd.Series("p_bar")])
    return new_p_bar

def calculate_smap_inf(nat, ncg, hs_at, hs_cg, V):
    """
    calculate the infinite-sampling mapping entropy
    """
    if ncg > nat:
        raise ValueError(f"n {nat} < N {ncg}")
    if hs_at < hs_cg:
        raise ValueError(f"hs_at {hs_at} < hs_cg {hs_cg}")
    delta_s_conf = (nat - ncg)*np.log(V) - hs_at + hs_cg
    return delta_s_conf

def calculate_smap(at_clust, mapping, pr, new_p_bar):
    """
    calculate the mapping entropy in its standard definition
    """
    # state-wise mapping entropy
    mapping_entropy = []
    for n in range(at_clust.shape[0]):
        s = np.where((at_clust.iloc[n,mapping] == new_p_bar.iloc[:,:-2]).all(1) == True)[0][0]
        mapping_entropy.append(pr[n]*np.log(pr[n]/new_p_bar.iloc[s,-1]))
    tot_smap = sum(mapping_entropy)
    # infinite sampling mapping entropy
    return tot_smap

# def calculate_smap_faster(at_clust):
#     """
#     faster calculation of smap
#     """

class cg_mappings():
    """cg mappings object."""

    def __init__(
            self,
            ncg,
            n_mappings
            ):
        self.ncg = ncg
        self.n_mappings = n_mappings
        self.hs_vect = np.zeros(self.n_mappings)
        self.hk_vect = np.zeros(self.n_mappings)
        self.smap_vect = np.zeros(self.n_mappings)
        self.smap_inf_vect = np.zeros(self.n_mappings)
        self.mappings_order = []
    
    def add_mapping()
    