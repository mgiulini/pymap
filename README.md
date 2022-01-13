# pymap

A program for describing how different selections of degrees of freedom (mappings) affect the amount of information retained about a full data set composed by discrete degrees of freedom.

Three quantities are calculated for each coarse-grained mappings, namely the mapping entropy:

![equ](https://latex.codecogs.com/gif.latex?S_{map}&space;=&space;\sum_{\phi}p(\phi)&space;\ln\left(\frac{p(\phi)}{\overline{p(\phi)}}&space;\right))

the resolution:

![equ](https://latex.codecogs.com/gif.latex?H_{s}&space;=&space;-\sum_{\phi}p(\phi)&space;\ln\left(p(\phi)\right))

and the relevance:

![equ](https://latex.codecogs.com/gif.latex?H_{k}&space;=&space;-\sum_{K}p(k)\ln\left(p(k)\right).)

where K is the set of unique frequencies observed in the sample.

If you use pymap please cite the following publications:



# Setup

A minimal conda environment to run the calculations can be constructed from the .yml file pymap.yml.

If you don't have conda installed, please see [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

```
conda env create --file pymap.yml
```

Once the environment is correctly created, it must be activated via:

```
conda activate pymap
```

# Usage

The program must be provided with the name of a data set contained in the *data* folder. The second, optional parameter, *max_binom*, that can be given to pymap is the maximum number of mappings that must be generated for each degree of coarse-graining.

## non-interacting spin system

```
python3 pymap spins
```

## financial market

To obtain the full results of [this paper]() one can run

```
python3 pymap m1
```

and 

```
python3 pymap m2
```

In the latter case, the mapping space starts to be quite big, and it is possible to explore just a portion of it in few minutes using *max_binom*

```
python3 pymap m2 5
```
