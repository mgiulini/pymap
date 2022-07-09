"""Library to perform clustering-related tasks."""
import numpy as np
import pandas as pd


def check_volume(dataframe, ncols):
    """Calculate the volume associated to each degree of freedom."""
    nique_list = [len(np.unique(dataframe.iloc[:, n])) for n in range(ncols-1)]
    V = np.unique(nique_list, return_counts=True)
    if len(V[1]) != 1:
        print(f"non-constant volume detected, V = {V}")
        print("choosing the most likely one...")
        agmax = np.argmax(V[1])
        V = V[0][agmax]
    else:
        print("constant volume detected, V = ", V)
        V = V[0][0]
    print("V = ", V)
    return V


def validate_clust(clust):
    """Do some validation of the at/cg dataframes."""
    if clust.iloc[:, :-1].duplicated().any():
        raise Exception("Duplicate row detected in dataframe")


def validate_dataframe(dataframe):
    """Do some checks on the original dataframe."""
    # you expect some duplicates, to give rise to non-trivial probabilities
    if dataframe.duplicated().any() is False:
        raise Warning("No duplicate row detected")


def get_clust(dataframe, mapping):
    """Compute the clustering."""
    gby = dataframe.groupby(dataframe.columns[mapping].tolist())
    clust = gby.size().reset_index().rename(columns={0: 'counts'})
    clust.columns = pd.concat(
        [
            pd.Series(dataframe.columns[mapping]),
            pd.Series("counts")
        ]
    )
    validate_clust(clust)
    return clust
