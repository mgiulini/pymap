This folder contains the data frames produced by running pymap.

More specifically, *results_m1.csv* and *results_m2.csv* refer to the two reduced models of the Nasdaq stock market.

Each Data Frame contains the following columns:
- *N* : number of retained degrees of freedom
- *mapping* : the index of the retained degrees of freedom
- *trans_mapping* : the names of the retained degrees of freedom, that is, the column names in the original Dataframe 
- *hs* : the coarse-grained resolution
- *hk* : the coarse-grained relevance
- *smap* : the mapping entropy
- *smap_inf* : the infinite-sampling mapping entropy, calculated as
 
![equ](https://latex.codecogs.com/gif.latex?S_{map}^{\infty}&space;=&space;\sum_{\phi}p(\phi)&space;\ln\left(p(\phi)&space;\right)-\sum_{\Psi}P(\Psi)\ln\left(P(\Psi)\right)+(n-N)\ln3)

![equ](https://latex.codecogs.com/gif.latex?S_{map}^{\infty}&space;=H_s^{\phi}-H_s^{\Psi}+(n-N)\ln3)

where the first two terms are the fine and coarse-grained resolutions.

## results_6d93

This folder includes the following files:
- *6d93_opt_smaps.txt*, with the optimal values of mapping entropy obtained through 48 independent optimisations;
- *6d93_opt_mappings.txt*, with the corresponding mappings;
- *6d93_opt_probs.pdb*, the PDB of tamapin in which the atom-wise probabilities of being conserved by the optimisations are insterted in the *beta factor* column;
- *6d93_smaps.txt*, values of mapping entropy calculated for all the mappings present [here](https://github.com/CIML-VARIAMOLS/GRAWL/blob/master/dataset/6d93_mappings_def.txt). The values differ from those showed [there](6d93_smaps_def_scaled.txt), as here the calculation is carried out directly using the Kullback-Leibler definition of the mapping entropy.
