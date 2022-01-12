import os 
import pandas as pd 
import numpy as np 
from scipy.stats import entropy
from scipy.special import binom
import math
import sys
import time

start_time = time.time()

def compute_res_rel(dict_cg_states):
    #print("values ", list(dict_cg_states.values()))
    data = list(dict_cg_states.values())
    freqs = [len(el) for el in data]
    #print(freqs)
    hs = entropy(freqs)
    print("hs = ", hs)
    ks = np.unique(freqs, return_counts=True)
    pk = np.multiply(ks[0], ks[1])
    print(pk)
    hk = entropy(pk)
    print("hk" , hk)
    return hs, hk

def compute_at_res_rel(dict_cg_states):
    #print("values ", list(dict_cg_states.values()))
    data = list(dict_cg_states.values())
    freqs = [len(el) for el in data]
    #print(freqs)
    hs = entropy(freqs)
    print("hs = ", hs)
    ks = np.unique(freqs, return_counts=True)
    pk = np.multiply(ks[0], ks[1])
    print(pk)
    hk = entropy(pk)
    print("hk" , hk)
    return hs, hk

def compute_entropies(at_clust, cg_clust,mapping):
    """
    starting from AT and CG data frames, mapping entropy is computed
    """
    cg_states = {}
    cg_at_freqs = {} # for each microstate in the CG state, I save the number of sampled confs
    for n in range(at_clust.shape[0]):
        #print("==\n", (at_clust.iloc[n,mapping] == cg_clust.iloc[:,:-1]))
        mapped_state = np.where((at_clust.iloc[n,mapping] == cg_clust.iloc[:,:-1]).all(1) == True)[0][0]
        #print("mapped_state = ", mapped_state)
        if mapped_state in cg_states.keys():
            cg_states[mapped_state].append(n)
            cg_at_freqs[mapped_state] += at_clust.iloc[n,-1]  # -1 column should be "records"
        else:
            cg_states[mapped_state] = [n]
            cg_at_freqs[mapped_state] = at_clust.iloc[n,-1]
    #print("cg states", cg_states)
    #print("cg_at_freqs ",cg_at_freqs)
    #print("cg_at_freqs values ",cg_at_freqs.values())
    hs , hk = compute_res_rel(cg_states)
    cg_at_freqs_values = list(cg_at_freqs.values())
    print("cg_at_freqs_values\n",cg_at_freqs_values)
    at_hs = entropy(cg_at_freqs_values)
    print("at hs = ", at_hs)
    ks = np.unique(cg_at_freqs_values, return_counts=True)
    pk = np.multiply(ks[0], ks[1])
    at_hk = entropy(pk)
    print("at hk = ", at_hk)
    #hs_at, hk_at = compute_at_res_rel(cg_states, at_clust)
    mapping_entropy = []
    for n in range(at_clust.shape[0]):
        mapped_state = np.where((at_clust.iloc[n,mapping] == cg_clust.iloc[:,:-1]).all(1) == True)[0][0]
        p_bar_r = cg_clust.iloc[mapped_state,-1]/(df.shape[0]*len(cg_states[mapped_state]))
        pr = at_clust.iloc[n,-1]/df.shape[0]
        mapping_entropy.append(pr*np.log(pr/p_bar_r))
    tot_smap = sum(mapping_entropy)
    print("mapping entropy", tot_smap)
    return hs, hk, tot_smap, at_hs, at_hk

def compute_entropies_fast(at_clust, cg_clust,mapping):
    """
    starting from AT and CG data frames, mapping entropy is computed
    """
    at_hs = entropy(cg_clust["records"])
    print("at hs = ", at_hs)
    ks = np.unique(cg_clust["records"], return_counts=True)
    pk = np.multiply(ks[0], ks[1])
    at_hk = entropy(pk)
    print("at hk = ", at_hk)
    cg_marginalised = at_clust.groupby(at_clust.columns[mapping].tolist()).size().reset_index().rename(columns={0:'records'})
    #todo: sub sum(cg_clust["records"]) with value
    p_bar_r = cg_clust["records"]/sum(cg_clust["records"])/cg_marginalised["records"]
    new_p_bar = pd.concat([cg_marginalised,p_bar_r],axis=1)
    print("new_p_bar\n", new_p_bar)
    pr = at_clust["records"]/sum(cg_clust["records"])
    print("sliced p_bar", new_p_bar.iloc[:,:-2])
    print("at_clust record", at_clust.iloc[1000,mapping])
    print("where", np.where((at_clust.iloc[1000,mapping] == new_p_bar.iloc[:,:-2]).all(1) == True))
    mapping_entropy = []
    for n in range(at_clust.shape[0]):
        s = np.where((at_clust.iloc[n,mapping] == new_p_bar.iloc[:,:-2]).all(1) == True)[0][0]
        #print(s)
        mapping_entropy.append(pr[n]*np.log(pr[n]/new_p_bar.iloc[s,-1]))
    tot_smap = sum(mapping_entropy)
    print("mapping entropy", tot_smap)
    #    print("n", n, "where", np.where((new_p_bar.iloc[:,:-1] == at_clust.iloc[n,mapping]).all(1) == True))
    #mapping_entropy = [pr[n]*np.log(pr/new_p_bar.where[at_clust.iloc[n,mapping]]) ]
    #:
    #    print("at_clust[mapping]", at_clust.iloc[:,mapping])


    


datafile = "./data/" + sys.argv[1] + ".csv"
with open(datafile) as f:
    ncols = len(f.readline().split(','))

print("ncols = ", ncols)
df = pd.read_csv(datafile, sep = ",",usecols=range(1,ncols))
print("df\n", df)
print("df shape", df.shape)
print("df.columns", df.columns)
n_at = df.shape[1]
## here usecols defines the dimensionality of the HR space
#print("groupby\n", df.groupby(df.columns.tolist()).size().reset_index().rename(columns={0:'records'}))
at_clust = df.groupby(df.columns.tolist()).size().reset_index().rename(columns={0:'records'})
#at_clust = df.groupby(df.columns[1:].tolist()).size().reset_index().\
#    rename(columns={0:'records'})
print("at_clust", at_clust)
print("atomistic columns", at_clust.columns)
print("at_clust shape", at_clust.shape)

print("total number of mappings = ", math.pow(2,n_at) - 1)
at_mapping = np.array(range(n_at))
print("df.columns[at_mapping]",df.columns[at_mapping])
print("at_clust.columns[at_mapping]",at_clust.columns[at_mapping])
# compute fully atomistic resolution
print("atomistic resolution ", entropy(at_clust["records"]))

# I here comment everything
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
            cg_clust = df.groupby(df.columns[mapping].tolist()).size().reset_index().rename(columns={0:'records'})
            cg_clust.columns = pd.concat([pd.Series(df.columns[mapping]), pd.Series("records")])
            print("cg_clust", cg_clust)
            #print("cg_clust.columns", cg_clust.columns)
            cg_mappings[key] = compute_entropies_fast(at_clust, cg_clust, mapping)
            ## calculating entropies wrt the atomistic clust

print(cg_mappings)

print("Total execution time (seconds) ", time.time() - start_time)