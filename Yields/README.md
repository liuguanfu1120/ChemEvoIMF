# A brief introduction

It is the collected and processed yields from the literature.
Please feel free to use it.

# How to process the original yields

The yields from the literature are processed by the following procedure:
- Convert the yields in mass to the yields in number. Namely, if the yield is given in mass, it is divided by the mass of the isotope to get the yield in number.
- Sum the yields of different isotopes of the same element.


# The considered literature

| Literature             | Short Name |
|------------------------|------------|
| [Nomoto et al. (2013)](https://www.annualreviews.org/doi/10.1146/annurev-astro-082812-140956)   | N13        |
| [Iwamoto et al. (1999)](https://iopscience.iop.org/article/10.1086/313278) | I99        |
| [Sukhbold et al. (2016)](https://iopscience.iop.org/article/10.3847/0004-637X/821/1/38) | S16        |
|[Karakas (2010)](https://academic.oup.com/mnras/article-lookup/doi/10.1111/j.1365-2966.2009.16198.x) | K10        |

The processed yields are named something like `S16W18-0.0200-SNcc-total.csv`:
- `S16` is the short name of the literature, i.e., [Sukhbold et al. (2016)](https://iopscience.iop.org/article/10.3847/0004-637X/821/1/38).
- `W18` is one of the models implemented in S16
- `0.0200` is the initial metallicity, which is a four-digit float
- `SNcc` is the type of yield, which is the core-collapse supernovae
- `total` means that it is the total yields, namely, the ejecta or the wind, while the `net` is the net yields, i.e., the ejecta minus the initial abundances.
  
As for Type Ia supernova, the initial metallicity in the names of yields files is "X.XXXX".

All the Python scripts used to get the processed yields are provided here, such as the Python script `Nomoto2013-Yields-Conversion`.ipynb` is used to process the yields from [Nomoto et al. (2013)](https://www.annualreviews.org/doi/10.1146/annurev-astro-082812-140956).

# Directory structure

The original yields are stored in the `Original` directory.
There are several sub-directories in the `Yields` directory except for the `Original` directory, i.e., `AGB`, `SNIa`, `SNcc`, and `HN`, the definition of which reads

| Directory | Meaning                             |
|-----------|-------------------------------------|
| AGB       | Asymptotic Giant Branch (AGB) stars |
| SNcc      | Core-collapse Supernovae (SNcc)     |
| SNIa      | Type Ia Supernovae (SNIa)           |
| HN        | Hyper Novae                         |


# Explanation of the processed yields files

All the processed yields are stored in the format of `csv` files. The digits of each processed yield are the same as in the original yields.
An example of the processed yields is shown below:

| element | 0.9      | 1        | 1.2      | 1.5      |
|---------|----------|----------|----------|----------|
| H       | 6.99E-02 | 1.04E-01 | 2.29E-01 | 3.55E-01 |
| He      | 7.84E-03 | 1.29E-02 | 2.88E-02 | 4.49E-02 |
| Li      | 1.07E-11 | 2.83E-11 | 5.27E-11 | 6.59E-11 |
| Be      | 0.00E+00 | 0.00E+00 | 0.00E+00 | 0.00E+00 |

where 0.9, 1, 1.2, and 1.5 are the initial mass in the unit of solar mass.
The unit of the number 6.99E-02, and 1.04E-01 is:

$$
\frac{1 M_{\odot}}{1 \mathrm{Da}},
$$

where $\mathrm{Da}$ is the Dalton or unified atomic mass unit (symbols: Da or u), see [Dalton (unit)](https://en.wikipedia.org/wiki/Dalton_(unit)).å

$$
1\mathrm{Da} = 1.660\times 10^{-27} ~ \mathrm{kg}.
$$

The exact value of $\frac{1 M_{\odot}}{1 \mathrm{Da}}$ is not important if only the ratio of the yield of different elements is considered.

An example of the processed yields of the Type Ia Supernovae (SNIa) is shown below:
| element | W7       | W70      | WDD1     | WDD2     |
|---------|----------|----------|----------|----------|
| C       | 4.83E-02 | 5.08E-02 | 5.42E-03 | 8.99E-03 |
| N       | 1.16E-06 | 4.46E-07 | 2.85E-04 | 2.70E-04 |
| O       | 1.43E-01 | 1.33E-01 | 8.82E-02 | 6.58E-02 |

Since the initial mass of the progenitor of the SNIa is not important, there is no initial mass in the header of the table but the name of the SNIa model, e.g., W7, W70, WDD1, and WDD2.
The number like 4.83E-02, and 5.08E-02 have the same unit, which is 

$$
\frac{1 M_{\odot}}{1 \mathrm{Da}}.
$$
