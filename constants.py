# -*- coding: utf-8 -*-
"""
Created on April 24, 2024
Last modified on April 24, 2024
@Author: Guan-Fu Liu

Description: Define some constants for the chemical evolution model.
It contains
1. the element names,
2. the atomic weights of the elements,
3. the solar abundances of the elements.

"""
import numpy as np
# If you import numpy as np, there will be something like constants.np
# which is not a good practice.
# The underscore prefix (_) is widely recognized in the Python community as an indicator of private or internal elements.



###### The basic atomic data and the solar abundances ######

# The atomic masses of the elements are taken from https://github.com/jzuhone/soxs/blob/main/soxs/constants.py
elem_names = [
    "",  # Just to make the index start from 1.
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Other"
]

atomic_weights = np.array(
    [
        0.0,   # Just to make the index start from 1.
        1.00794,
        4.00262,
        6.941,
        9.012182,
        10.811,
        12.0107,
        14.0067,
        15.9994,
        18.9984,
        20.1797,
        22.9898,
        24.3050,
        26.9815,
        28.0855,
        30.9738,
        32.0650,
        35.4530,
        39.9480,
        39.0983,
        40.0780,
        44.9559,
        47.8670,
        50.9415,
        51.9961,
        54.9380,
        55.8450,
        58.9332,
        58.6934,
        63.5460,
        65.3800,
        0.0  # The atomic weight of the other elements is set to ZERO.
    ]
)

# Taken from https://spex-xray.github.io/spex-help/reference/commands/abundance.html
# which is converted to the mass fraction.
# It refers to how much mass of each element per unit mass of hydrogen in the sun.
#  The abbrevation and corresponding references are as follows:
# "Default": Lodders et al. (2009), https://ui.adsabs.harvard.edu/abs/2009LanB...4B..712L/abstract
# "AG": Anders & Grevesse (1989), https://ui.adsabs.harvard.edu/abs/1989GeCoA..53..197A/abstract
# "Allen": Allen (1973), https://ui.adsabs.harvard.edu/abs/1973asqu.book.....A/abstract
# "RA": Ross & Aller (1976), https://ui.adsabs.harvard.edu/abs/1976Sci...191.1223R/abstract
# "Grevesse": Grevesse et al. (1992), https://ui.adsabs.harvard.edu/abs/1992ESASP.348..305G/abstract
# "GS": Grevesse & Sauval (1998), https://ui.adsabs.harvard.edu/abs/1998SSRv...85..161G/abstract
# "Lodders": Lodders proto-Solar Lodders, https://ui.adsabs.harvard.edu/abs/2003ApJ...591.1220L/abstract
# "solar": Lodders Solar photospheric Lodders, https://ui.adsabs.harvard.edu/abs/2003ApJ...591.1220L/abstract
abund_tables = {
    "Default": np.array(
        [
            0.0,  # Just to make the index start from 1.
            1.000e+00,
            3.854e-01,
            1.476e-08,
            2.111e-10,
            7.770e-09,
            3.305e-03,
            1.135e-03,
            9.609e-03,
            5.838e-07,
            2.538e-03,
            5.071e-05,
            9.578e-04,
            8.722e-05,
            1.074e-03,
            9.830e-06,
            5.159e-04,
            7.002e-06,
            1.416e-04,
            5.620e-06,
            9.257e-05,
            5.920e-08,
            4.525e-06,
            5.567e-07,
            2.603e-05,
            1.938e-05,
            1.809e-03,
            5.296e-06,
            1.099e-04,
            1.314e-06,
            3.251e-06,
            0.0 # The mass ratio of the other elements to H is set to ZERO.
        ]
    ),
    "AG": np.array(
        [
            0.0,  # Just to make the index start from 1.
            1.000e+00,
            3.881e-01,
            9.954e-11,
            1.263e-10,
            4.270e-09,
            4.326e-03,
            1.559e-03,
            1.351e-02,
            6.844e-07,
            2.463e-03,
            4.876e-05,
            9.168e-04,
            7.900e-05,
            9.887e-04,
            8.661e-06,
            5.159e-04,
            1.112e-05,
            1.439e-04,
            5.114e-06,
            9.109e-05,
            5.615e-10,
            4.641e-06,
            5.054e-07,
            2.413e-05,
            1.338e-05,
            2.591e-03,
            4.863e-06,
            1.036e-04,
            1.022e-06,
            2.582e-06,
            0.0  # The mass ratio of the other elements to H is set to ZERO.
        ]
    ),
    "Allen": np.array(
        [
            0.0,  # Just to make the index start from 1.
            1.000e+00,
            3.380e-01,
            3.451e-11,
            1.126e-10,
            1.073e-08,
            3.946e-03,
            1.267e-03,
            1.049e-02,
            7.504e-07,
            1.665e-03,
            4.056e-05,
            6.343e-04,
            6.571e-05,
            9.227e-04,
            1.018e-05,
            5.042e-04,
            1.400e-05,
            3.148e-04,
            3.457e-06,
            7.934e-05,
            7.402e-10,
            6.406e-06,
            1.270e-06,
            3.652e-05,
            1.369e-05,
            2.206e-03,
            7.361e-06,
            1.162e-04,
            1.994e-06,
            1.028e-06,
            0.0  # The mass ratio of the other elements to H is set to ZERO.
        ]
    ),
    "RA": np.array(
        [
            0.0,  # Just to make the index start from 1.
            1.000e+00,
            2.506e-01,
            6.886e-11,
            1.263e-10,
            1.382e-09,
            4.967e-03,
            1.210e-03,
            1.098e-02,
            6.844e-07,
            5.909e-04,
            4.346e-05,
            9.381e-04,
            8.864e-05,
            1.245e-03,
            9.718e-06,
            5.042e-04,
            1.112e-05,
            4.056e-05,
            5.607e-06,
            8.902e-05,
            4.890e-10,
            5.328e-06,
            5.292e-07,
            2.646e-05,
            1.434e-05,
            1.752e-03,
            4.644e-06,
            1.110e-04,
            7.239e-07,
            1.828e-06,
            0.0  # The mass ratio of the other elements to H is set to ZERO.
        ]
    ),
    "Grevesse": np.array(
        [
            0.0,  # Just to make the index start from 1.
            1.000e+00,
            3.706e-01,
            9.954e-11,
            1.263e-10,
            4.270e-09,
            4.228e-03,
            1.297e-03,
            1.177e-02,
            6.844e-07,
            2.407e-03,
            4.876e-05,
            9.168e-04,
            7.900e-05,
            9.887e-04,
            8.661e-06,
            5.159e-04,
            1.112e-05,
            1.312e-04,
            5.114e-06,
            9.109e-05,
            7.069e-08,
            4.973e-06,
            5.054e-07,
            2.413e-05,
            1.338e-05,
            1.793e-03,
            4.863e-06,
            1.036e-04,
            1.022e-06,
            2.582e-06,
            0.0  # The mass ratio of the other elements to H is set to ZERO.
        ]
    ),
    "GS": np.array(
        [
            0.0,  # Just to make the index start from 1.
            1.000e+00,
            3.881e-01,
            1.406e-08,
            2.352e-10,
            6.613e-09,
            3.946e-03,
            1.156e-03,
            1.073e-02,
            5.692e-07,
            2.407e-03,
            4.765e-05,
            9.168e-04,
            8.272e-05,
            1.012e-03,
            1.116e-05,
            5.042e-04,
            6.702e-06,
            9.955e-05,
            5.233e-06,
            8.902e-05,
            5.615e-08,
            4.136e-06,
            5.292e-07,
            2.527e-05,
            1.847e-05,
            1.752e-03,
            4.753e-06,
            1.036e-04,
            1.229e-06,
            3.034e-06,
            0.0   # The mass ratio of the other elements to H is set to ZERO.
        ]
    ),
    "Lodders": np.array(
        [
            0.0,  # Just to make the index start from 1.
            1.000e+00,
            3.792e-01,
            1.542e-08,
            2.700e-10,
            7.593e-09,
            3.437e-03,
            1.104e-03,
            9.134e-03,
            6.387e-07,
            1.784e-03,
            5.347e-05,
            1.005e-03,
            9.282e-05,
            1.135e-03,
            1.066e-05,
            5.789e-04,
            7.520e-06,
            1.652e-04,
            5.871e-06,
            1.022e-04,
            6.300e-08,
            4.749e-06,
            5.938e-07,
            2.707e-05,
            2.072e-05,
            1.921e-03,
            5.584e-06,
            1.135e-04,
            1.379e-06,
            3.251e-06,
            0.0  # The mass ratio of the other elements to H is set to ZERO.
        ]
    ),
    "solar": np.array(
        [
            0.0,  # Just to make the index start from 1.
            1.000e+00,
            3.154e-01,
            1.312e-08,
            2.298e-10,
            6.463e-09,
            2.925e-03,
            9.395e-04,
            7.774e-03,
            5.436e-07,
            1.484e-03,
            4.551e-05,
            8.556e-04,
            7.720e-05,
            9.662e-04,
            8.863e-06,
            4.927e-04,
            6.401e-06,
            1.406e-04,
            4.997e-06,
            8.699e-05,
            5.240e-08,
            3.950e-06,
            5.054e-07,
            2.304e-05,
            1.724e-05,
            1.635e-03,
            4.753e-06,
            9.664e-05,
            1.147e-06,
            2.767e-06,
            0.0  # The mass ratio of the other elements to H is set to ZERO.
        ]
    ),
}

Z_sun = { }
for key, arr in abund_tables.items():
    Z_sun[key] = arr[3:].sum() / arr[1:].sum()  # The solar metallicity.

Mstar_min = 0.08  # The minimum mass of a star in solar mass.
Mstar_max = 150.0  # The maximum mass of a star in solar mass.

__all__ = [
    'elem_names',
    'atomic_weights',
    'abund_tables',
    'Z_sun',
    'Mstar_min',
    'Mstar_max'
]
# The solar abundance tables have been checked on April 24, 2024.