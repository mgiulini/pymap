"""task measure file."""
import os
import time

import numpy as np
from scipy.special import binom

from libs.libclust import get_clust
from libs.libentropy import (calculate_entropies, calculate_pbar,
                             calculate_smap, calculate_smap_inf)


class MEASURE:
    """MEASURE class."""

    def __init__(
            self,
            n_at,
            at_clust,
            hs_at,
            V,
            pr,
            max_binom,
            at_mapping,
            args,
            ):
        """Initialize module."""
        self.cg_mappings = dict()
        self.cg_mappings_order = []
        self.n_at = n_at
        self.at_clust = at_clust
        self.hs_at = hs_at
        self.V = V
        self.pr = pr
        self.max_binom = max_binom
        self.at_mapping = at_mapping
        self.args = args
        self.start_time = time.time()
        print(self.args)

    def output_mappings(self, output_filename):
        """
        Functions that outputs the mappings to a file.

        Parameters
        ----------
        mapping_dict : dict
            dictionary of mappings

        mapping_order : list
            list of ordered mappings

        output_filename : str
        """
        header = "N\tmapping\ttrans_mapping\ths\thk\tsmap\tsmap_inf"
        header += os.linesep
        with open(output_filename, "w") as wfile:
            # write header
            wfile.write(header)
            for ord_map in self.cg_mappings_order:
                output_str = []
                for elem in self.cg_mappings[ord_map]:
                    if isinstance(elem, float):
                        output_str.append(f"{elem:.6f}")
                    else:
                        output_str.append(f"{elem}")
                wfile.write("\t".join(output_str) + os.linesep)

    def run(self, df):
        """Run measure module."""
        for ncg in range(1, self.n_at+1):
            cg_count = int(binom(self.n_at, ncg))
            print("cg_count", cg_count)
            k = 0
            max_range = min(cg_count, self.max_binom)
            print(f"max_range = {max_range}")
            fixed_n_mappings = []
            while k < max_range:
                mapping = np.random.choice(self.at_mapping, ncg, replace=False)
                mapping.sort()
                key = tuple(mapping)
                if key not in self.cg_mappings.keys():
                    k += 1
                    if self.args.verbose:
                        print("adding key", key, " k = ", k)
                    cg_clust = get_clust(df, mapping)
                    hs, hk = calculate_entropies(cg_clust)
                    smap_inf = calculate_smap_inf(self.n_at,
                                                  ncg,
                                                  self.hs_at,
                                                  hs,
                                                  self.V)
                    p_bar = calculate_pbar(self.at_clust,
                                           cg_clust,
                                           df.shape[0],
                                           mapping)
                    smap = calculate_smap(mapping, self.pr, p_bar)
                    trans_mapping = list(self.at_clust.columns[mapping])
                    self.cg_mappings[key] = (len(mapping),
                                             mapping,
                                             trans_mapping,
                                             hs,
                                             hk,
                                             smap,
                                             smap_inf)
                    fixed_n_mappings.append(key)
            fixed_n_mappings.sort()
            # extending the original list
            self.cg_mappings_order.extend(fixed_n_mappings)
            elap_time = time.time() - self.start_time
            print(f"ncg = {ncg} took  {elap_time:.6f} seconds")
