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
