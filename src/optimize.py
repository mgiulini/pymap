import numpy as np

from libs.libclust import check_volume, get_clust
from libs.libentropy import (calculate_entropies, calculate_pbar,
                             calculate_smap, calculate_smap_inf)


def convert_mapping(mapping):
    print("conv_mapping: ",np.nonzero(mapping)[0])


def make_a_move(mapping, np_retained, np_not_retained, n_at, ncg):
    """
    routine that creates a new mapping by changing a site
    """
    at1_idx = np.random.randint(0,ncg)
    at1 = np_retained[at1_idx]
    at2_idx = np.random.randint(0,n_at-ncg)
    at2 = np_not_retained[at2_idx]
    mapping[at1] = 0
    mapping[at2] = 1
    print(f"move made: at1 {at1} at1_idx {at1_idx}, at2 {at2} at2_idx {at2_idx}")
    return mapping,at1_idx,at2_idx,at1,at2

def optimize(df, parameters, n_at, hs_at, V, at_clust, pr): # TODO: find a better wat to pass objects
    """Simulated Annealing optimisation."""
    cg_mappings = dict()
    cg_mappings_order = []
    mapping = np.zeros(n_at, dtype=int)
    # define initial mapping
    rd_sel = np.random.choice(np.arange(0,n_at), size=parameters["ncg"], replace=False)
    for el in rd_sel:
        mapping[el] = 1
    #mapping = np.random.choice(at_mapping, parameters["ncg"], replace=False)
    print(f"first mapping {mapping}")
    #convert_mapping(mapping)
    atom_ret = np.nonzero(mapping)[0]
    #print(f"atom_ret {atom_ret}")
    #print("shape atom_ret", atom_ret.shape)
    atom_nnret = np.nonzero(mapping == 0)[0]
    #atom_nnret = np.array([el for el in range(n_at) if el not in mapping])
    print(f"atom_nnret {atom_nnret}")
    mapping_prime = np.array(mapping, copy=True)
    # calculating first mapping obs
    cg_clust = get_clust(df, atom_ret)
    hs, hk = calculate_entropies(cg_clust)
    smap_inf = calculate_smap_inf(n_at,
                                  parameters["ncg"],
                                  hs_at,
                                  hs,
                                  V)
    p_bar = calculate_pbar(at_clust,
                           cg_clust,
                           df.shape[0],
                           atom_ret)
    smap = calculate_smap(atom_ret, pr, p_bar)
    cg_mappings[tuple(atom_ret)] = (parameters["ncg"],
                        atom_ret,
                        list(at_clust.columns[mapping]),
                        hs,
                        hk,
                        smap,
                        smap_inf)
    
    print("beginning optimization")
    t_zero = 1.0
    decay_time = 10.0
    cg_mappings_order = [tuple(atom_ret)]

    smaps = [smap]
    for n in range(1, parameters["nsteps"]):
        mapping_prime, at1_idx, at2_idx, at1, at2 = make_a_move(mapping_prime,atom_ret,atom_nnret,n_at,parameters["ncg"])
        at1, at2 =  atom_ret[at1_idx], atom_nnret[at2_idx]
        atom_ret[at1_idx] = atom_ret[-1]
        atom_nnret[at2_idx] = atom_nnret[-1]
        atom_ret[-1] = at2
        atom_nnret[-1] = at1
        # calculating stuff
        cg_clust = get_clust(df, atom_ret)
        hs_prime, hk_prime = calculate_entropies(cg_clust)
        smap_inf_prime = calculate_smap_inf(n_at,
                                      parameters["ncg"],
                                      hs_at,
                                      hs_prime,
                                      V)
        p_bar_prime = calculate_pbar(at_clust,
                                   cg_clust,
                                   df.shape[0],
                                   atom_ret)
        smap_prime = calculate_smap(atom_ret, pr, p_bar_prime)
        # metropolis rule
        print(f"smap {smap} vs smap_prime {smap_prime}")
        if (smap_prime < smap):
            print("move accepted")
            accepted = True
        else:
            temp = t_zero * np.exp(-n/ decay_time)
            p = np.exp((smap - smap_prime)/temp)
            r = np.random.uniform()
            print(f"temp p r {temp} {p} {r}")
            if(r < p): #move accepted, update mapping and stuff
                print("p > r : move accepted")
                accepted = True
            else:
                print("move rejected")
                accepted = False
        # update observable
        if accepted:
            smap, hs, hk, smap_inf = smap_prime, hs_prime, hk_prime, smap_inf_prime
            mapping[at1]=0
            mapping[at2]=1
            # saving stuff
            atom_ret.sort()
            key = tuple(atom_ret)
            smaps.append(smap)
            cg_mappings_order.append(key)
        else:
            mapping_prime[at1]=1
            mapping_prime[at2]=0
            atom_ret[-1] = at1
            atom_nnret[-1] = at2

        # print/saving stuff
        print(f"atom_ret {atom_ret}")
        print(f"atom_nnret {atom_nnret}")
        cg_mappings[tuple(atom_ret)] = (parameters["ncg"],
                            atom_ret,
                            list(at_clust.columns[atom_ret]),
                            hs,
                            hk,
                            smap,
                            smap_inf)
    print(cg_mappings_order)
    print(smaps)
    return cg_mappings, cg_mappings_order