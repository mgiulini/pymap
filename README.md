# pymap

A program for describing how different selections of *N* out of *n* degrees of freedom (mappings) affect the amount of information retained about a full data set.

Three quantities are calculated for each low-resolution representation, namely the mapping entropy:

![equ](https://latex.codecogs.com/gif.latex?S_{map}&space;=&space;\sum_{\phi}p(\phi)&space;\ln\left(\frac{p(\phi)}{\overline{p(\phi)}}&space;\right))

the resolution:

![equ](https://latex.codecogs.com/gif.latex?H_{s}&space;=&space;-\sum_{\phi}p(\phi)&space;\ln\left(p(\phi)\right))

and the relevance:

![equ](https://latex.codecogs.com/gif.latex?H_{k}&space;=&space;-\sum_{K}p(k)\ln\left(p(k)\right).)

where K is the set of unique frequencies observed in the sample.

If you use pymap please cite [this paper](https://arxiv.org/abs/2203.00100).


# Setup

A minimal conda environment (see [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)) to run the calculations can be generated from the .yml file pymap.yml using the following command:

```
conda env create --file pymap.yml
```

Once the environment is correctly created, it must be activated via:

```
conda activate pymap
```

# Usage

The program must be provided with the name of a data set contained in the *data* folder. The second, optional parameter, *max_binom*, that can be given to pymap is the maximum number of mappings that must be generated for each degree of coarse-graining. The default choice is to generate all the coarse-grained mappings for each *N*, a task that becomes prohibitive when *n > 15*.

## non-interacting spin system

The first data set described in [this article](https://arxiv.org/abs/2203.00100) contains 20 non-interacting spins. The variables of interest can be calculated with the following command

```
python3 pymap spins
```

In this context, the mapping space is quite big, and *max_binom* allows one to explore just a portion of it in few minutes: 

```
python3 pymap spins 5
```

## financial market

To obtain the full results for the simple model of the Nasdaq stock market reported [here](https://arxiv.org/abs/2203.00100) one can use the following command:

```
python3 pymap m1
```

and 

```
python3 pymap m2
```
