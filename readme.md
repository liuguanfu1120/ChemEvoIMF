# A quick start

1. Download the code
2. Run `ChemEvo_example.ipynb`
3. The output file will be saved in `./outputs/`
4. Run `read_results_example.ipynb` to read and analyse the output file

I recommends you to gain some basic concepts about how HDF5 file works before you start.
If you need some specific post-processing of the output file, you may write your own code to read and analyse the output file.

# Set the input

Here, more details about the input parameters are provided.

## Set the star formation history (`SFH`)

The star formation history (`SFH`) is stored in a dict. The keys of the dict are `Age`, and `SFR`, both of which are case sensitive.

- `SFH['Age']` should be a numpy array in the shape of `(N,)`, where `N` is the number of ages. The ages are in the unit of `yr.
- `SFH['SFR']` should be a numpy array in the shape of `(N,)`, where `N` is the number of ages. The star formation rate is in the unit of $M_{\odot}$/yr.

`ChemEvo.py` will check if the shape of `SFH['Age']` and `SFH['SFR']` are the same, if there is any negative value in `SFH['SFR']`, and if `SFH['Age']` is in ascending order. However, please ensure that you input an **correct** `SFH`.

`SFH['SFR'][i]` is the star formation rate in the age interval from `SFH['SFR'][i]` to `SFH['SFR'][i+1]`. Therefore, the star formation rate at the last age won't take any effect.

## Set the star formation efficiency (`SFE`)

The star formation efficiency (`SFE`) is a float that should be within $(0,0.5]$. The typical value of SFE is 0.1, and an enhanced value is 0.3. A SFE larger than 0.5 is unphysical, and hence the `ChemEvo.py` will raise an error and stop.
The total mass of the initial gas is determined by
$$
\text{The total mass of the initial gas} = \frac{\text{The cumulative stellar mass}}{\text{SFE}}.
$$

## Set the initial gas composition

Note that the code does not distinguish the "initial" gas and "primordial" gas.
The code just needs to initilize the gas composition at the beginning of the calculation.
The initial gas composition could be the primordial gas or any other gas composition you want.

### set the input gas composition (`input_primordial_gas`)

The default value of `input_primordial_gas` is `None`, which means the initial gas composition is set to be same as that described in the table II of Cyburt et al. (2016).

If you want to set the initial gas composition by yourself, you can set `input_primordial_gas` to be a numpy array of a shape of (32,).
`input_primordial_gas` does not need to be normalized, as the code will normalize it automatically.

### set the mass fraction of metals except Li (`Z_0`) and solar abundance table `solar_set`

Table II of Cyburt et al. (2016) provides the mass ratio of He/H and Li/H.
If you provide the mass fraction of metals except Li (`Z_0`, a float), the code will calculate the mass fraction of H, He, Li, and other metals.
The mass of other metals will be redistributed according to the solar abundance table you provide.

`solar_set` is a string, whose possible values can be found in `constants.py`.

The `input_primordial_gas` and `Z_0` are mutually exclusive. If you set `input_primordial_gas`, the `Z_0` will be ignored.

The recommended set is `input_primordial_gas = None`, `Z_0 = 0`, and `solar_set = "Default"`.
For more details, please refer to `primordial_gas.py`.

## Set the yield tables (`yield_files`)

It is dict like

```py
yield_files = {
                "AGB+SNcc": "./inputs/NuPyCEE/agb_and_massive_stars_C15_N13_0_5_HNe/yields2.h5",
                "SNIa": "./inputs/NuPyCEE/sn1a_i99_W7/yields2.h5",
               }
```

The keys should exactly be "AGB+SNcc" and "SNIa". The yield tables of AGB and SNcc are stored in one HDF5 file, while the yield table(s) of SNIa is in the other HDF5 file.

If you want to generate the HDF file of your own yield tables, please see `./Yields-Example/Generate-Yields.ipynb`.

The provided jupyter notebook file will tell you how to generate the yield files step-by-step. There are two kinds of yield files, `yields1.h5` and `yields2.h5`, which are equivalent but of different format.

## Set the SNIa

### Set the calibration factor delay time distribution (DTD) of SNIa (`p_preset`)

Maoz & Mannucci (2012) gives an observational constraint on the SNIa rate, which goes as follows
$$
\Psi(t) =
4\times 10^{-13} ~ \mathrm{SN}~\mathrm{yr}^{-1}~M_{\odot}^{-1}\left( \frac{t}{1~\mathrm{Gyr}} \right)^{-1},
$$
where $\Psi(t)$ is the specific SNIa rate at time $t$ after the star formation. 
    $4\times 10^{-13}$ is used to ensure the integral time between 10 Myr and 10 Gyr is $2.2\times 10^{-3} M_{\odot}^{-1}$.
However, such a SNIa rate is based on dite-Salpeter IMF, which should be calibrated to the IMF you use.
The code will calculate the calibration factor automatically, but you can also set it manually by setting `p_preset`.

Unless you have a good reason and clear understanding, I recommend you to **use the default value** of `p_preset`, that is `None`.

### Trun on or off the SNIa (`SNIaOn`)

`SNIaOn` is a boolean variable, which is `True` by default. If you want to turn off the SNIa, you can set it to be `False`.

## Set the IMF

## Set the evolution function of IMF (`imf_evolve`)

A simpe example goes as follows,

```py

def imf_evolve(Z_gas, age):
    """
    The IMF could evolve with the metallicity of the gas or with the age.
    The first argument is the metallicity of the gas, while the second argument is the age of the gas.
    There must be two arguments in the function even if one of them is not used (which is often the case).
    Here, I provide a very simple case where the IMF is top-heavy if the gas metallicity is below 0.02,
    but it is Salpeter-like otherwise.

    Parameters
    ----------
    Z_gas : float
        The metallicity of the gas.
    age : float
        The age in yr.
    
    Returns
    -------
    function
        The IMF function.
    """
    if Z_gas < 0.02:
        return lambda m: m**(-1.3) if (m>=constants.Mstar_min and m<=constants.Mstar_max) else 0
    else:
        return lambda m: m**(-2.35) if (m>=constants.Mstar_min and m<=constants.Mstar_max) else 0
```

The `imf_evolve` should return a function of the IMF, which takes the initial mass of the star as the argument and returns the IMF at that mass.
If you do not want the IMF to evolve, you can set `imf_evolve` to be `None`.

### Set the IMF by `imf_dict` (pre-determined IMF)

If you set the IMF by `imf_evolve`, the exact form of the IMF is determined by evolution, namely, it is pre-determined by the user.
`imf_dict` is a dict, which could take the following form

1. Invariant IMF

```py
imf_dict = {"Salpeter": None}
imf_dict = {"Kroupa": None}
imf_dict = {"Charier": None}
imf_dict = {"DietSalpeter": None}
# np.array([2.35]) is just an example, you can change it to any other value.
imf_dict = {"PowerLaw": np.array([2.35])}
# np.array([[0.08, 100], [1, 10], [10, 1], [150, 0.1]]) is just an example, you can change it to any other value.
imf_dict = {"Custom": np.array([[0.08, 100], [1, 10], [10, 1], [150, 0.1]])}
```

"Invariant" means the IMF of different star formation epochs is the same.

- For the Salpeter, Kroupa, Chabrier, and diet-Salpeter IMF, the value of `imf_dict` is set to be `None`, as the value is not needed.
- For the power-law IMF, the value of `imf_dict` is a numpy array of the shape of `(1,)`, which is the power index of the IMF.
- For the custom IMF, the value of `imf_dict` is a numpy array of the shape of `(N, 2)`, where `N` is the number of mass points. The first column is the initial mass of the star, while the second column is the IMF at that mass. The code will interpolate the input array of the shape of `(n, 2)` to recover the IMF.
  
2. Variant IMF

"Variant" means the IMF of different star formation epochs could be different.
The difference between the variant IMF set by `imf_dict` and the IMF set by `imf_evolve` is that the former is pre-determined, while the latter is determined by the evolution.

- Input the power law indices of the IMF of different star formation epochs like

```py
# np.array([2.35, 2.3, 2.4, 2.5]) is just an example, you can change it to any other value.
imf_dict = {"PowerLaw": np.array([2.35, 2.3, 2.4, 2.5])}
```

The number of the power law indices should be no less than the number of the star formation epochs (positive star formation rate).
For example, if there are 10 ages with positive star formation rate, the length of the numpy array should be no less than 10 (extra indices will be ignored but will raise a warning.).

- Input a list of the custom IMF of different star formation epochs like

```py
# Just an example, you can change it to any other value.
imf_dict = {"Custom": [np.array([[0.08, 100], [1, 20], [10, 2], [150, 0.21]]), 
                      np.array([[0.08, 120], [1, 30], [10, 2], [150, 0.31]]), 
                      np.array([[0.08, 50], [1, 40], [10, 2], [150, 0.41]]), 
                      np.array([[0.08, 60], [1, 25], [10, 2], [150, 0.51]])]
            }
```

`imf_evolve` is prioritized over `imf_dict`. If you set `imf_evolve`, the `imf_dict` will be ignored and a warning will be raised.
I recommend you to set `imf_dict` to be `None` if you set `imf_evolve`.

## Set the elements you consider (`ElemNotice`)

`ElemNotice` is a list whose elements are element names, like

```py
ElemNotice = ["H", "He", "O", "Other"]
```

**ChemEvoIMF will only calculate the evolution of the elements shown in `ElemNotice`**
The evolution of the elements not shown in `ElemNotice` will be "added up" and treated as "Other".
For the above `ElemNotice`, the code will calculate the evolution of H, He, O, and all the other elements in "Other".

The computation time is proportional to the number of elements you consider, and hence you may not want to consider all the elements.
`ChemEvoIMF` only consider the first 30 elements, and hence you cannot the element whose atomic number is larger than 30.
`ElemNotice` should contain "H", "He", and "Other".
Actually, if any of "H", "He", and "Other" does not appear in `ElemNotice`, the code will add it automatically.

## Set the interpolation method (`interp_kind`)

`interp_kind` is a string, whose possible values are "linear-linear", "linear-log", "log-linear", "log-log", and "TNG-like" ("TNG-like" is equivalent to "linear-log").
"linear-log" means the interpolation is done by linearly interpolating the yield with respecto to the **log** of metallicity **Z**.

The gas metallicity will be affected by the number of elements you consider when using "log-log" or "log-linear" interpolation. This is unphysical and these two interpolation methods are not recommended and will be removed in the future.

## Set the output (`output_dir`, `out_file`, `comments`)

- `output_dir` is the directory where the output file is saved.
  `output_dir` should be a string, like `output_dir = "./outputs"`.
  It should not end with "/".
- `out_file` is the name of the output file.
  `out_file` should be a string, like `out_file = "N13.h5"`.
- `comments` is a string, which is used to describe the calculation.
  I recommend you to write some important information about the calculation in `comments`, which will be saved in the output file.

# Post-processing

After running the code, you will get an output file, which is saved in the directory you set by `output_dir`.
`read_results_example.ipynb` provides an example of how to read and analyse the output file.

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


## Units

- Mass: solar mass $M_{\odot}$.
- Time: yr.

## The interpolation/extrapolation of the stellar yields

Both the mass range of the IMF and the stellar yield tables is restricted to be $[0.08,150]~M_{\odot}$

   - If the minimum mass $M_{\rm min}$ of the IMF is larger than $0.08~M_{\odot}$, then the IMF in the range $[0.08,M_{\rm min}]~M_{\odot}$ is set to be zero.
   - As is often the case, the $M_{\rm min}$ in the stellar yield tables is larger than $0.08~M_{\odot}$, so the yield in the range $[0.08, M_{\rm min}]~M_{\odot}$ is set to be zero, and the $M_{\rm max}$ in the stellar yield tables is smaller than $150~M_{\odot}$, so the yield in the range $[M_{\rm max},150]~M_{\odot}$ is set to be the same as that at $M_{\rm max}$. As for the $M_{\rm rem}$, it is set to be $M_{\rm ini}$ in the range $[0.08,M_{\rm min}]~M_{\odot}$ and we keep $M_{\rm ini}-M_{\rm rem}$ a constant in the range $[M_{\rm max},150]~M_{\odot}$.
   - The number of sampled mass points in the stellar yield tables is often just a few, so we will interpolate the yield tables to get the yield at any mass in the range $[0.08,150]~M_{\odot}$.

## Solar abundance table

- The solar abundance tables are taken from [SPEX](https://spex-xray.github.io/spex-help/reference/commands/abundance.html) but converted to the ratio of the mass of each element to the mass of H.
- The shape of a single solar abundance table `abund_table` is `(32,)`
 `abund_table[0]` is zero and should be ignored, which is used to make the index of the element the same as the atomic number.
  `abund_table[1:31]` is the mass ratio of each element to the mass of H. We only consider the first 30 elements.
  `abund_table[31]` is also zero, which stands for the mass ratio of the elements heavier than Zn to the mass of H.

