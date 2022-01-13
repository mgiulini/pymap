import os 
import pandas as pd 
import numpy as np 
from scipy.stats import entropy
from scipy.special import binom
import math
import sys
import time

start_time = time.time()

def compute_entropies(at_clust, df,mapping):
    """
    starting from original and atomistic data frames, resolution, relevance, and mapping entropy are computed
    """
    # coarse-grained clustering
    cg_clust = df.groupby(df.columns[mapping].tolist()).size().reset_index().rename(columns={0:'records'})
    cg_clust.columns = pd.concat([pd.Series(df.columns[mapping]), pd.Series("records")])
    # resolution with respect to the full sampling
    hs = entropy(cg_clust["records"])
    print("hs [", str(mapping),  "] = ", hs)
    ks = np.unique(cg_clust["records"], return_counts=True)
    pk = np.multiply(ks[0], ks[1])
    hk = entropy(pk)
    print("hk [",str(mapping) ,"] = ", hk)
    # marginalised coarse-grained distribution
    cg_marginalised = at_clust.groupby(at_clust.columns[mapping].tolist()).size().reset_index().rename(columns={0:'records'})
    # smeared probability distribution
    p_bar_r = cg_clust["records"]/sum(cg_clust["records"])/cg_marginalised["records"]
    new_p_bar = pd.concat([cg_marginalised,p_bar_r],axis=1)
    #print("new_p_bar\n", new_p_bar)
    # detailed, atomistic proba
    pr = at_clust["records"]/sum(cg_clust["records"])
    print("sliced p_bar", new_p_bar.iloc[:,:-2])
    mapping_entropy = []
    for n in range(at_clust.shape[0]):
        s = np.where((at_clust.iloc[n,mapping] == new_p_bar.iloc[:,:-2]).all(1) == True)[0][0]
        #print(s)
        mapping_entropy.append(pr[n]*np.log(pr[n]/new_p_bar.iloc[s,-1]))
    tot_smap = sum(mapping_entropy)
    print("mapping entropy", tot_smap)

# read_data
datafile = "./data/" + sys.argv[1] + ".csv"
with open(datafile) as f:
    ncols = len(f.readline().split(','))

print("number of columns in the dataset = ", ncols)
df = pd.read_csv(datafile, sep = ",",usecols=range(1,ncols))
print("df shape", df.shape)
print("df.columns", df.columns)
n_at = df.shape[1]
# atomistic clustering: the original dataframe is divided in microstates
at_clust = df.groupby(df.columns.tolist()).size().reset_index().rename(columns={0:'records'})
print("at_clust", at_clust)
print("atomistic columns", at_clust.columns)
print("at_clust shape", at_clust.shape)

print("total number of mappings = ", math.pow(2,n_at) - 1)
at_mapping = np.array(range(n_at))
print("df.columns[at_mapping]",df.columns[at_mapping])
print("at_clust.columns[at_mapping]",at_clust.columns[at_mapping])
# compute fully atomistic resolution
print("atomistic resolution ", entropy(at_clust["records"]))

cg_mappings = dict()
#cg_mappings_at = dict()
for ncg in range(1,3):
#for ncg in range(1,n_at+1):
    print("ncg = ", ncg)
    cg_count = int(binom(n_at,ncg))
    print("cg_count", cg_count)
    k = 0
    #while k < cg_count :
    for s in range(1):
        mapping = np.random.choice(at_mapping, ncg, replace=False)
        print("unsorted mapping", mapping)
        mapping.sort()
        print("sorted mapping", mapping)
        key = tuple(mapping)
        if key not in cg_mappings.keys():
            k += 1
            #print("adding key", key, " k = ", k)
            cg_mappings[key] = compute_entropies(at_clust, df, mapping)
            ## calculating entropies wrt the atomistic clust

print(cg_mappings)

print("Total execution time (seconds) ", time.time() - start_time)