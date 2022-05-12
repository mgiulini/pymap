import os 
from pathlib import Path
import pandas as pd 
import numpy as np 
from scipy.stats import entropy
from scipy.special import binom
import math
import sys
import time
import argparse
# local modules
from utils import (
    check_optional_parameters,
    check_volume,
    get_clust,
    system_parameters_setup,
    validate_clust,
    calculate_pbar,
    calculate_smap,
    calculate_smap_inf,
    calculate_entropies,
    system_parameters_setup,
    check_mandatory_parameters,
    check_optional_parameters
)

def parse_arguments():
    """
    parse and check the command-line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--parameters", help="dat input parameter file")
    parser.add_argument("-v","--verbose", help="increase output verbosity",action="store_true")
    args = parser.parse_args()
    print(args.__dict__)
    # checks
    if not args.parameters:
        raise Exception("Parameter file is mandatory")
    if args.verbose:
        print("verbosity turned on")
    return args

def main():
    """
    main function
    """
    start_time = time.time()

    args = parse_arguments()

    # read_data
    cleaned_pars = system_parameters_setup(args.parameters)   

    with open(cleaned_pars["input_filename"], "r") as f:
        ncols = len(f.readline().split(','))
    print("number of columns in the dataset = ", ncols)
    
    df = pd.read_csv(cleaned_pars["input_filename"], sep = ",",usecols=range(1,ncols))
    print("df shape", df.shape)
    print("df.columns", df.columns)
    n_at = df.shape[1]
    
    # creating fully detailed mapping
    print("total number of mappings = ", math.pow(2,n_at) - 1)
    at_mapping = np.array(range(n_at))
    print("df.columns[at_mapping]",df.columns[at_mapping])

    # atomistic clustering: the original dataframe is divided in microstates
    at_clust = get_clust(df,at_mapping)
    #at_clust = df.groupby(df.columns.tolist()).size().reset_index().rename(columns={0:'records'})
    print("at_clust", at_clust)
    print("atomistic columns", at_clust.columns)
    print("at_clust shape", at_clust.shape)
    pr = at_clust["records"]/df.shape[0] # detailed, atomistic probability distribution
    # check volume
    V = check_volume(df,ncols)

    # atomistic quantities
    hs_at, hk_at = calculate_entropies(at_clust)
    print("at_clust.columns[at_mapping]",at_clust.columns[at_mapping])
    print("atomistic resolution ", hs_at) # computing fully atomistic resolution
    print("atomistic relevance ", hk_at) # computing fully atomistic resolution
    
    cg_mappings = dict()
    cg_mappings_vect = []
    # going through the levels of coarse-graining
    for ncg in range(1,n_at+1):
        print("ncg = ", ncg, ", elapsed time (seconds) = %8.6lf" % (time.time() - start_time))
        cg_count = int(binom(n_at,ncg))
        print("cg_count", cg_count)
        k = 0
        max_range = min(cg_count,cleaned_pars["max_binom"])
        print(f"max_range = {max_range}")
        fixed_n_mappings = []
        
        while k < max_range:
            mapping = np.random.choice(at_mapping, ncg, replace=False)
            mapping.sort()
            key = tuple(mapping)
            if key not in cg_mappings.keys():
                k += 1
                if args.verbose:
                    print("adding key", key, " k = ", k)
                cg_clust = get_clust(df,mapping)
                hs, hk = calculate_entropies(cg_clust)
                smap_inf = calculate_smap_inf(n_at, ncg, hs_at, hs, V)
                p_bar = calculate_pbar(at_clust, cg_clust, df.shape[0], mapping)
                smap = calculate_smap(at_clust, mapping, pr, p_bar)
                cg_mappings[key] = len(mapping),mapping,list(at_clust.columns[mapping]),hs,hk,smap,smap_inf
                fixed_n_mappings.append(key)
        print(f"fixed_n_mappings {fixed_n_mappings}")
        fixed_n_mappings.sort()
        print(f"sorted fixed_n_mappings {fixed_n_mappings}")


    
    output_df = pd.DataFrame(cg_mappings.values())
    output_df.columns = ["N","mapping","trans_mapping","hs","hk","smap","smap_inf"]
    output_df.to_csv(cleaned_pars["output_filename"], sep=",", float_format = "%8.6lf")
    
    print("Total execution time (seconds) %8.6lf" % (time.time() - start_time))

# running main
main()