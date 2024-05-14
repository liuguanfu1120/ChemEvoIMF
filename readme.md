# Structure of the repository

## `./inputs/`
It contains the essential input files for the calculation of chemical evolution with variant or invariant IMF.
The essential input files are:
- Stellar yields
- Star formation history

### `./inputs/K10_N13_S16/`

The pre-processed stellar yields from [Karakas et al. (2010)](https://academic.oup.com/mnras/article-lookup/doi/10.1111/j.1365-2966.2009.16198.x), [Nomoto et al. (2013)](https://www.annualreviews.org/doi/10.1146/annurev-astro-082812-140956), and [Sukhbold et al. (2016)](https://iopscience.iop.org/article/10.3847/0004-637X/821/1/38) are stored in this directory.

A single sub-directory `./inputs/K10_N13_S16/` contains a single set of stellar yields, which should contain `yields1.h5` and `yields2.h5`.
There may be some other files but all the information is stored in these two files.

### `./inputs/NuPyCEE/`

The pre-processed stellar yields from [NuPyCEE](https://github.com/NuGrid/NuPyCEE)

A single sub-directory `./inputs/NuPyCEE/` contains a single set of stellar yields, which should contain `yields1.h5` and `yields2.h5`.
There may be some other files but all the information is stored in these two files.

### `./inputs/SFH.h5`, `./inputs/SFH-1delta.txt`, and `./inputs/SFH-1square.txt`

The star formation history (SFH) files.

## `./outputs/`

It contains the output files of the calculation of chemical evolution with variant or invariant IMF.
The output files are named `YYYY-MM-DD-HH-MM-SS.`h5`.
The name convention can be changed according to the user's preference.

## `./primordial_gas.py`

Define the primordial gas composition.

## `./IMF.py`

Define the IMF.

## `./SupernoaveIa.py`

Define the delay-time distribution of SNe Ia.

## `./utils.py`

The utility functions, like the mass lifetime relation.

## `./constants.py`

Define some important constants, like the solar abundance table, atomic weights, solar metallicity, the minimum and maximum initial mass of stars, etc.


# Some conventions

## Conventions of the stellar yields file

There are two kinds of stellar yield files, which are named `yields1.h5` and `yields2.h5`.
They are just different formats but contain the same information.

Take `yields1.h5` as an example, the stellar yields of different initial metallicities are stored in different groups, like `Z=0.0001`, `Z=0.0003`, etc.
The name of groups cannot contain special characters like `=` and `.` when using `pandas.to_hdf()` to save the data frame to the hdf5 file.
Therefore, the `=` and `.` in the group name are replaced by `_` in `yields2.h`.

The are several datasets in each group, like `Interpolated`, `Original`, and `MassLifetime`.
ONLY the dataset `Interpolated` is used in the calculation of chemical evolution.


<figure style="text-align: center;">
  <img src="./readme-fig/Yields.png" alt="Fig1" style="display: block; margin: 0 auto; width: 120%; height: auto;">
  <figcaption>

   The yields files `./inputs/NuPyCEE/agb_and_massive_stars_C15_LC18_R_mix/yields1.h5` and `./inputs/NuPyCEE/agb_and_massive_stars_C15_LC18_R_mix/yields2.h5`
  </figcaption>
</figure>



