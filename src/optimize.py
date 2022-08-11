import numpy as np
import pandas as pd

from libs.libclust import get_clust
from libs.libentropy import (calculate_entropies, calculate_pbar,
                             calculate_smap, calculate_smap_inf)


class OPTIMIZE:
    """OPTIMIZE class."""

    def __init__(
            self,
            n_at,
            ncg,
            at_clust,
            hs_at,
            V,
            pr,
            nsteps
            ):
        """
        Initialize the class.

        Parameters
        ----------
        identificator : str
        """
        self.cg_mappings = []
        self.n_at = n_at
        self.ncg = ncg
        self.at_clust = at_clust
        self.hs_at = hs_at
        self.V = V
        self.pr = pr
        self.t_zero = 1.0
        self.decay_time = 10.0
        self.mapping = np.zeros(n_at, dtype=int)
        self.nsteps = nsteps
        self.step = 0 # before starting
        # define initial mapping
        rd_sel = np.random.choice(np.arange(0,self.n_at), size=self.ncg, replace=False)
        for el in rd_sel:
            self.mapping[el] = 1
        print(f"first mapping {self.mapping}")
        self.atom_ret = np.nonzero(self.mapping)[0]
        self.atom_nnret = np.nonzero(self.mapping == 0)[0]
        print(f"atom_nnret {self.atom_nnret}")


    def convert_mapping(self):
        print("conv_mapping: ",np.nonzero(self.mapping)[0])

    def make_a_move(self):
        """
        routine that creates a new mapping by changing a site
        """
        at1_idx = np.random.randint(0,self.ncg)
        at1 = self.atom_ret[at1_idx]
        at2_idx = np.random.randint(0,self.n_at-self.ncg)
        at2 = self.atom_nnret[at2_idx]
        self.mapping[at1] = 0
        self.mapping[at2] = 1
        print(f"move made: at1 {at1} at1_idx {at1_idx}, at2 {at2} at2_idx {at2_idx}")
        self.atom_ret[at1_idx] = self.atom_ret[-1]
        self.atom_nnret[at2_idx] = self.atom_nnret[-1]
        self.atom_ret[-1] = at2
        self.atom_nnret[-1] = at1
        return at1,at2

    def metropolis(self, smap, smap_prime):
        """routine that applies the Metropolis rule."""
        print(f"smap {smap:.6f} vs smap_prime {smap_prime:.6f}")
        if (smap_prime < smap):
            print("move accepted")
            accepted = True
        else:
            temp = self.t_zero * np.exp(-self.step/ self.decay_time)
            p = np.exp((smap - smap_prime)/temp)
            r = np.random.uniform()
            print(f"temp p r {temp:.3f} {p:.3f} {r:.3f}")
            if(r < p): #move accepted, update mapping and stuff
                print("p > r : move accepted")
                accepted = True
            else:
                print("move rejected")
                accepted = False
        return accepted

    def calculate_observables(self, df):
        """calculates the observables."""
        cg_clust = get_clust(df, self.atom_ret)
        hs, hk = calculate_entropies(cg_clust)
        smap_inf = calculate_smap_inf(self.n_at,
                                      self.ncg,
                                      self.hs_at,
                                      hs,
                                      self.V)
        p_bar = calculate_pbar(self.at_clust,
                               cg_clust,
                               df.shape[0],
                               self.atom_ret)
        smap = calculate_smap(self.atom_ret, self.pr, p_bar)
        return hs,hk,smap,smap_inf

    def append_quantities(self,hs,hk,smap,smap_inf):
        """Append the relevant quantities to the self.cg_mappings."""
        self.cg_mappings.append(
                        [self.ncg,
                        self.atom_ret.copy(),
                        list(self.at_clust.columns[self.atom_ret]),
                        hs,
                        hk,
                        smap,
                        smap_inf]
                        )

    def output_mappings(self, output_filename):
        """output self.cg_mappings to the desired output filename."""
        output_df = pd.DataFrame(self.cg_mappings)
        output_df.columns = ["N","mapping","trans_mapping","hs","hk","smap","smap_inf"]
        output_df.to_csv(output_filename, sep = "\t", float_format="%.6f", index=None)


    def run(self, df):
        """Run the Simulated Annealing optimisation."""
        # calculating first mapping obs
        hs,hk,smap,smap_inf = self.calculate_observables(df)
        self.append_quantities(hs,hk,smap,smap_inf)
        print("beginning optimization")
        for n in range(self.nsteps):
            at1, at2 = self.make_a_move()
            hs_prime,hk_prime,smap_prime,smap_inf_prime = self.calculate_observables(df)
            accepted = self.metropolis(smap, smap_prime)
            if accepted:
                smap, hs, hk, smap_inf = smap_prime, hs_prime, hk_prime, smap_inf_prime
                self.atom_ret.sort()
                self.append_quantities(hs,hk,smap,smap_inf)
            else:
                # revert move
                self.mapping[at1]=1
                self.mapping[at2]=0
                self.atom_ret[-1] = at1
                self.atom_nnret[-1] = at2
            print(f"self.atom_ret {self.atom_ret}")
            self.step += 1