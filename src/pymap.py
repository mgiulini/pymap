"""pymap main file."""

import math
import time

import numpy as np
import pandas as pd

# local modules
from libs.libclust import check_volume, get_clust
from libs.libentropy import calculate_entropies
from libs.libio import parse_arguments, system_parameters_setup
from measure import MEASURE
from optimize import OPTIMIZE


def main():
    """Define main function."""
    start_time = time.time()

    args = parse_arguments()

    # read_data
    cleaned_pars = system_parameters_setup(args.parameters, args.task)

    with open(cleaned_pars["input_filename"], "r") as f:
        ncols = len(f.readline().split(','))
    print("number of columns in the dataset = ", ncols)

    df = pd.read_csv(
        cleaned_pars["input_filename"],
        sep=",",
        usecols=range(1, ncols)
    )
    print("df shape", df.shape)
    print("df.columns", df.columns)
    n_at = df.shape[1]

    # creating fully detailed mapping
    print("total number of mappings = ", math.pow(2, n_at) - 1)
    at_mapping = np.array(range(n_at))
    print("df.columns[at_mapping]", df.columns[at_mapping])

    # atomistic clustering: the original dataframe is divided in microstates
    at_clust = get_clust(df, at_mapping)
    print("at_clust", at_clust)
    print("atomistic columns", at_clust.columns)
    print("at_clust shape", at_clust.shape)
    pr = at_clust["counts"]/df.shape[0]  # at. probability distribution
    # check volume
    V = check_volume(df, ncols)

    # atomistic quantities
    hs_at, hk_at = calculate_entropies(at_clust)
    print("at_clust.columns[at_mapping]", at_clust.columns[at_mapping])
    print(f"atomistic resolution {hs_at:.6f}")  # computing fully at. resolution
    print(f"atomistic relevance {hk_at:.6f}")  # computing fully at. resolution

    if args.task == "measure":
        measure_obj = MEASURE(
                    n_at = n_at,
                    at_clust = at_clust,
                    hs_at = hs_at,
                    V = V,
                    pr = pr,
                    at_mapping = at_mapping,
                    max_binom = cleaned_pars["max_binom"],
                    args = args
        )
        measure_obj.run(df)
        measure_obj.output_mappings(cleaned_pars["output_filename"])

    elif args.task == "optimize":
        optimize_obj = OPTIMIZE(
                    n_at,
                    cleaned_pars["ncg"],
                    at_clust,
                    hs_at,
                    V,
                    pr,
                    cleaned_pars["nsteps"]
                    )
        optimize_obj.run(df)
        optimize_obj.output_mappings(cleaned_pars["output_filename"])

    print("Total execution time (seconds) %8.6lf" % (time.time() - start_time))


# running main
main()
