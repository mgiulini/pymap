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

# Testing

Pytest is employed to test the correct installation of pymap. In order to do so, run the following command from the main directory:

```
python -m pytest tests
```

Or directly run *pytest* inside the *tests* folder:

```
cd tests
pytest
```

# Contributing

If you like to add a contribution to this open-source project, follow these steps:

**1.** create an issue with a brief explanation of the contribution;

**2.** add a reasonable label to the issue, or create a new one;

**3.** create a new branch *entirely dedicated to this contribution* either on this repo on your fork;

**4.** develop the code

**5.** use *tox* to test the code. In particular, you should run the following commands:

```
tox -e py310
tox -e lint  
``` 

The first command tests the code with a standard python 3.10 environment, while the second checks the code-style.
    
**6.** open a new Pull-Request on this page, correctly linking the issue. Ask for a review from anyone of the contributors 

Enjoy!

# Usage

The program must be provided with one single command line argument (-p), namely a (relative) path to a parameter file, containing the parameters to be employed. A list of the accepted parameters is provided here:

| Parameter | Description | Type | Mandatory |
| ----------- | ----------- | ---- | ------- |
| *input_filename* | relative path to the input data | str | yes |
| *output_filename* | relative path to the desired output file | str | yes |
| *max_binom* *| max number of mappings that must be generated for each degree of coarse-graining | int | no |


*The default choice is to generate all the coarse-grained mappings for each *N*, a task that becomes prohibitive when *n > 15*. 

Verbosity can be turned on with the *-v* (*--verbose*) flag.

In general, running

```
python pymap -h
```

shows the available command line arguments.

## non-interacting spin system

The first data set described in [this article](https://arxiv.org/abs/2203.00100) contains 20 non-interacting spins. The variables of interest can be calculated with the following command

```
python3 pymap.py -p parameters/parameters_spins.dat
```

In this context, the mapping space is quite big, and *max_binom* allows one to explore just a portion of it in few minutes: 

```
python3 pymap.py -p parameters/parameters_spins_test.dat
```

## financial market

To obtain the full results for the simple model of the Nasdaq stock market reported [here](https://arxiv.org/abs/2203.00100) one can use the following command:

```
python3 pymap.py -p parameters/parameters_m1.dat
```

and 

```
python3 pymap.py -p parameters/parameters_m2.dat
```
