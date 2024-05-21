# -*- coding: utf-8 -*-
"""
Created on April 24, 2024
Last modified on April 25, 2024
@Author: Guan-Fu Liu

Description: Define some utilities for the chemical evolution model, which includes
- 1. The primordial gas class.
- 2. The mass lifetime relation.
- 3. The IMF class
- 4. The SNIa class
"""

import constants
import numpy as np
from scipy import interpolate
from scipy.integrate import quad


###### The primordial gas mass class ######
class Primordial_gas:
    """
    Define the primordial gas class.
    """
    def __init__(self, mass, Z_0=None):
        """
        Initialize the primordial gas class.

        Parameters
        ----------
        mass: float
            The total mass of the primordial gas, in the unit of solar mass.
            It is usually determined by the following equation:
            the total mass of primordial gas = the total initial mass of formed star * star formation efficiency.
        Z_0: float or None
            The primordial metallicity but without lithium. 
            If it is set to be None, the primordial gas is set to be the default primordial gas, i.e., primor_gas.
            You can also set it to be a float, which should be a much small number like 1e-6.
            However, you can also set it to be larger than 1e-4, it will raise a warning.
        
        Returns
        -------
        None

        self.primor_gas will contain the following keys:
        'H': float
            The mass of H in the primordial gas in the unit of solar mass.
        'He': float
            The mass of He in the primordial gas in the unit of solar mass.
        'Li': float
            The mass of Li in the primordial gas in the unit of solar mass.
        'Z_0': float
            The mass of metals without lithium in the primordial gas.
        'Z': float
            The total mass of metals (with lithium) in the primordial gas.
        """
        # From the table II in Cyburt et al. (2016), whose DOI is 10.1103/RevModPhys.88.015004,
        # Y_p=0.2470, D/H=2.579e-5, 3He/H=0.9996e-5, 7Li/H=4.648e-10, 6Li/H=1.288e-14.
        # Y_p is the primordial helium mass fraction,
        # D/H is the ratio of primordial deuterium mass fraction to hydrogen mass fraction,
        # 3He/H is the ratio of primordial helium-3 mass fraction to primordial hydrogen mass fraction, 
        # 7Li/H is the ratio of primordial lithium-7 mass fraction to primordial hydrogen mass fraction,
        # and 6Li/H is the ratio of primordial lithium-6 mass fraction to primordial hydrogen mass fraction.
        # The primordial metallicity Z_0 is set to be 10^-6 (references to be found).
        # Obsviouly, X+Y_p+X(D/H+3He/H+7Li/H+6Li/H)+Z_0=1, where X is the primordial hydrogen mass fraction.
        primor_gas = { }
        primor_gas['Y_p'] = 0.2470
        primor_gas['D/H'] = 2.579e-5
        primor_gas['3He/H'] = 0.9996e-5
        primor_gas['7Li/H'] = 4.648e-10
        primor_gas['6Li/H'] = 1.288e-14
        primor_gas['Z_0'] = 1e-6  # The primordial metallicity but without lithium, you can also set it to be 0.
        primor_gas['X'] = (1 - primor_gas['Y_p'] - primor_gas['Z_0'])/\
                        (1 + primor_gas['D/H'] + primor_gas['3He/H'] + primor_gas['7Li/H'] + primor_gas['6Li/H'])
        # The difference between the primordial hydrogen mass fraction primor_gas['X'] and the rough value 
        # (1 - primor_gas['Y_p'] - primor_gas['Z_0']) is of order 1e-5, which is negligible but considered here.
        primor_gas['H'] = primor_gas['X'] + primor_gas['D/H']*primor_gas['X']  
        # Add the mass fraction of deuterium to hydrogen.
        primor_gas['He'] = primor_gas['Y_p'] + primor_gas['3He/H']*primor_gas['X']  
        # Add the mass fraction of helium-3 to helium.
        primor_gas['Li'] = primor_gas['7Li/H']*primor_gas['X'] + primor_gas['6Li/H']*primor_gas['X']  
        # Add the mass fraction of lithium-7 and lithium-6 to lithium.
        primor_gas['Z'] = primor_gas['Z_0'] + primor_gas['Li']  # The primordial metallicity, you can also set it to be 0.
        self.primor_gas = primor_gas.copy()
        if Z_0 is None:
            # Z_0 is None, use the default primordial gas.
            pass
        else:
            # Z_0 is not None and pre-set, re-calculate the primordial gas.
            self.primor_gas['Z_0'] = Z_0
            self.primor_gas['X'] = (1 - self.primor_gas['Y_p'] - self.primor_gas['Z_0'])/\
                    (1 + self.primor_gas['D/H'] + self.primor_gas['3He/H'] + self.primor_gas['7Li/H'] + self.primor_gas['6Li/H'])
            self.primor_gas['H'] = self.primor_gas['X'] + self.primor_gas['D/H']*self.primor_gas['X']
            self.primor_gas['He'] = self.primor_gas['Y_p'] + self.primor_gas['3He/H']*self.primor_gas['X']
            self.primor_gas['Li'] = self.primor_gas['7Li/H']*self.primor_gas['X'] + self.primor_gas['6Li/H']*self.primor_gas['X']
            self.primor_gas['Z'] = self.primor_gas['Z_0'] + self.primor_gas['Li']
            self.primor_gas.pop('X')
            self.primor_gas.pop('D/H')
            self.primor_gas.pop('Y_p')
            self.primor_gas.pop('3He/H')
            self.primor_gas.pop('7Li/H')
            self.primor_gas.pop('6Li/H')
        if self.primor_gas['Z_0'] > 1e-4:
            print("Warning: The primordial metallicity is larger than 1e-4, which is not recommended.")
            print("It will move with such a large primordial metallicity.")
        else:
            pass
        for key in self.primor_gas.keys():
            self.primor_gas[key] = self.primor_gas[key] * mass
        
        primor_gas.pop('X')  
        # Remove the primordial hydrogen mass fraction X from the dictionary.
        primor_gas.pop('D/H')  
        # Remove the ratio of primordial deuterium mass fraction to hydrogen mass fraction from the dictionary.
        primor_gas.pop('Y_p')  
        # Remove the primordial helium mass fraction from the dictionary.
        primor_gas.pop('3He/H')  
        # Remove the ratio of primordial helium-3 mass fraction to primordial hydrogen mass fraction from the dictionary.
        primor_gas.pop('7Li/H')  
        # Remove the ratio of primordial lithium-7 mass fraction to primordial hydrogen mass fraction from the dictionary.
        primor_gas.pop('6Li/H')  
        # Remove the ratio of primordial lithium-6 mass fraction to primordial hydrogen mass fraction from the dictionary.


    def add_metals(self, abund_table=constants.abund_tables['Default']):
        """
        Calculate the composition of the metals in the primordial gas.

        Parameters
        ----------
        abund_table: numpy array
        It should be an array with the shape (31, ), which contains the solar abundances of the first 30 elements.
        The first element in the array should be ZERO, which is just to make the index start from 1.
        For more details, please refer to the module constants.py.
        """
        for i, elem in enumerate(constants.elem_names[4:]):
            self.primor_gas[elem] = abund_table[i+4]/abund_table[4:].sum() * self.primor_gas['Z_0']
            # Suppose the composition of the metals except lithium in the primordial gas can be scaled by the solar abundance.
            # Namely, the relative abundance between the metals except is the same as the solar abundance,
            # but the total mass fraction of the metals is Z_0.
            # Such as O/Fe is the same as the solar abundance, but the total mass fraction of the metals is Z_0.
        self.primor_gas.pop('Z_0')  # Remove the primordial metallicity Z_0 without lithium from the dictionary.
        self.primor_gas['Gas'] = np.zeros(len(constants.elem_names))  # The mass of the primordial gas in the unit of solar mass.
        for i in range(1, 31):
            self.primor_gas['Gas'][i] = self.primor_gas[constants.elem_names[i]]
            # The index starts from 1, so the first element is ZERO.
        return self.primor_gas


### This part is tested on April 25, 2024.
###### The primordial gas mass class ######


###### The mass lifetime relation ######
# The following data is taken from Table 14 in Portinari et al. (1998), https://ui.adsabs.harvard.edu/abs/1998A&A...334..505P.
lifetime = { }
# The initial mass, first column in Table 14 in Portinari et al. (1998).
lifetime['Mini'] = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 
                 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 12.0, 15.0, 20.0, 30.0, 40.0, 60.0, 100.0, 120.0])
# The lifetime of the stars with metallicity Z=0.0004, second column in Table 14 in Portinari et al. (1998).
lifetime['Z=0.0004'] = np.array([4.28e10, 2.37e10, 1.41e10, 8.97e9, 6.03e9, 4.23e9, 3.08e9, 2.34e9, 
                                 1.92e9, 1.66e9, 1.39e9, 1.18e9, 1.11e9, 9.66e8, 8.33e8, 4.64e8, 3.03e8,
                                 1.61e8, 1.01e8, 7.15e7, 5.33e7, 3.42e7, 2.13e7, 1.54e7, 1.06e7, 6.90e6, 
                                 5.45e6, 4.20e6, 3.32e6, 3.11e6])
# The lifetime of the stars with metallicity Z=0.004, third column in Table 14 in Portinari et al. (1998).
lifetime['Z=0.004'] = np.array([5.35e10, 2.95e10, 1.73e10, 1.09e10, 7.13e9, 4.93e9, 3.52e9, 2.64e9,
                                2.39e9, 1.95e9, 1.63e9, 1.28e9, 1.25e9, 1.23e9, 1.08e9, 5.98e8, 3.67e8, 
                                1.82e8, 1.11e8, 7.62e7, 5.61e7, 3.51e7, 2.14e7, 1.52e7, 1.05e7, 6.85e6,
                                5.44e6, 4.19e6, 3.38e6, 3.23e6])
# The lifetime of the stars with metallicity Z=0.008, fourth column in Table 14 in Portinari et al. (1998).
lifetime['Z=0.008'] = np.array([6.47e10, 3.54e10, 2.09e10, 1.30e10, 8.46e9, 5.72e9, 4.12e9, 2.92e9,
                                2.36e9, 2.18e9, 1.82e9, 1.58e9, 1.41e9, 1.25e9, 1.23e9, 6.86e8, 4.12e8,
                                1.93e8, 1.15e8, 7.71e7, 5.59e7, 3.44e7, 2.10e7, 1.49e7, 1.01e7, 6.65e6, 
                                5.30e6, 4.15e6, 3.44e6, 3.32e6])
# The lifetime of the stars with metallicity Z=0.02, fifth column in Table 14 in Portinari et al. (1998).
lifetime['Z=0.02'] = np.array([7.92e10, 4.45e10, 2.61e10, 1.59e10, 1.03e10, 6.89e9, 4.73e9, 3.59e9,
                               2.87e9, 2.64e9, 2.18e9, 1.84e9, 1.59e9, 1.38e9, 1.21e9, 7.64e8, 4.56e8, 
                               2.03e8, 1.15e8, 7.45e7, 5.31e7, 3.17e7, 1.89e7, 1.33e7, 9.15e6, 6.13e6,
                               5.12e6, 4.12e6, 3.39e6, 3.23e6])
# The lifetime of the stars with metallicity Z=0.05, sixth column in Table 14 in Portinari et al. (1998).
lifetime['Z=0.05'] = np.array([7.18e10, 4.00e10, 2.33e10, 1.42e10, 8.88e9, 5.95e9, 4.39e9, 3.37e9,
                               3.10e9, 2.51e9, 2.06e9, 1.76e9, 1.51e9, 1.34e9, 1.24e9, 6.58e8, 3.81e8, 
                               1.64e8, 8.91e7, 5.67e7, 3.97e7, 2.33e7, 1.39e7, 9.95e6, 6.99e6, 5.15e6, 
                               4.34e6, 3.62e6, 3.11e6, 3.10e6])
# In the original table, the last lifetime is still 3.11e6, which is the same as the lifetime of the star with 100 solar mass.
# But it is not allowed by interpolate.interp1d, we change it to 3.10e6.


# According to Ritter et al. (2018), https://iopscience.iop.org/article/10.3847/1538-4365/aad691, 
# the lifetime of the star with a certain  mass and a certain metallicity is given by
# linear interpolation of log t and log M.


def mass_to_lifetime(mass, Z):
    """
    Get the lifetime of the star with a certain mass and a certain metallicity.

    Parameters
    ----------
    mass: float
        The mass of the star in the unit of solar mass.
    Z: float
        The initial metallicity of the star.
    
    Returns
    -------
    t: float
        The lifetime of the star in the unit of year.
    """
    # The keys of lifetime dict (except the key "Mini")
    Z_list = [key for key in lifetime.keys() if 'Z=' in key]
    Z_float = np.array([float(a.replace("Z=", "")) for a in Z_list])
    # Get the closest metallicity of the metallicity of the selecet stellar population
    Z_index = np.abs(Z-Z_float).argmin()
    Z_key = Z_list[Z_index]
    # Linear interpolation or extrapolation of log10 t versus log10 Mini
    interp = interpolate.interp1d(np.log10(lifetime['Mini']), np.log10(lifetime[Z_key]), 
                                  kind=1, fill_value='extrapolate')
    t = interp(np.log10(mass))
    t = 10**t
    return t


mass_to_lifetime = np.vectorize(mass_to_lifetime, otypes=[np.float64])


def lifetime_to_mass(t, Z):
    """
    Get the mass of the star with a certain lifetime and a certain metallicity.
    
    Parameters
    ----------
    t: float
        The lifetime of the star in the unit of year.
    Z: float
        The initial metallicity of the star.
    """
    # The keys of lifetime dict (except the key "Mini")
    Z_list = [key for key in lifetime.keys() if 'Z=' in key]
    Z_float = np.array([float(a.replace("Z=", "")) for a in Z_list])
    # Get the closest metallicity of the metallicity of the selecet stellar population
    Z_index = np.abs(Z-Z_float).argmin()
    Z_key = Z_list[Z_index]
    # Linear interpolation or extrapolation of log10 M versus log10 t
    interp = interpolate.interp1d(np.log10(lifetime[Z_key]), np.log10(lifetime['Mini']), 
                                  kind=1, fill_value='extrapolate')
    mass = interp(np.log10(t))
    mass = 10**mass
    return mass


lifetime_to_mass = np.vectorize(lifetime_to_mass, otypes=[np.float64])

### This part is tested on April 27, 2024.
###### The mass lifetime relation ######


###### The IMF class ######

class IMF:
    """
    Define the IMF class.
    """
    def __init__(self, IMF_type="Salpeter", IMF_arr=None, power_index=None):
        """
        Initialize the IMF class.

        Parameters
        ----------
        IMF_type: str
            The type of the IMF. It should be in ['Salpeter', 'Chabrier', 'Kroupa', 'Custom'].
            'Salpeter': 
                The Salpeter IMF from Salpeter (1955), https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S
            'Chabrier': 
                The Chabrier IMF from Chabrier (2003), https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C
            'Kroupa': 
                The Kroupa IMF from Kroupa (2001), https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K
            'PowerLaw':
                The power-law IMF, which is set by the user. You should provide the power-law index of the IMF.
            'DietSalpeter':
                The diet Salpeter IMF from Bell & de Jong 2001, https://ui.adsabs.harvard.edu/abs/2001ApJ...550..212B/abstract
                It is almost the same as the Salpeter IMF, but the slope of mass below 0.35 solar mass is set to be 1.
            'Custom': 
                The custom IMF, which is set by the user. You should provide the IMF by inputing a numpy array with 
                the shape (n, 2), where n is the number of the mass bins and the first column is the mass of the stars in the
                unit of solar mass and the second column is the IMF of the stars. The implemented IMF will be the interpolation
                of the input IMF array. Therefore, n should be as large as possible to get a smooth IMF.
                The method of interpolation is just simple linear interpolation. We do not recommend using too complicated 
                interpolation methods, but recommend using a large n to get a smooth IMF.
            
            The IMF should be normalized to 1, i.e., the integral of the IMF over the mass range should be 1.
            The mass range is set to be [0.08, 150] solar mass, which is the common range for the IMF.
            The IMF outside the mass range is set to be ZERO.
        IMF_arr: numpy array or None
            The IMF array for the custom IMF. It should be an array with the shape (n, 2).
            If it is not None, the IMF_type should be set to be 'Custom'.
        power_index: float or None
            The power-law index of the power-law IMF. It should be a float.
            The IMF is set to be m^(-power_index), where m is the mass of the star.
            If it is not None, the IMF_type should be set to be 'PowerLaw'.
        """
        self.IMF_type = IMF_type
        self.IMF_arr = IMF_arr
        self.power_index = power_index
        # The normalization factor of IMFs
        self.A = { }
        self.A['Salpeter'] = quad(self._Salpeter, constants.Mstar_min, constants.Mstar_max)[0]
        self.A['Chabrier'] = quad(self._Chabrier, constants.Mstar_min, constants.Mstar_max)[0]
        self.A['Kroupa'] = quad(self._Kroupa, constants.Mstar_min, constants.Mstar_max)[0]
        self.A['DietSalpeter'] = quad(self._DietSalpeter, constants.Mstar_min, constants.Mstar_max)[0]
        # Check the input.
        if self.IMF_type not in ['Salpeter', 'Chabrier', 'Kroupa', 'Custom', 'PowerLaw', 'DietSalpeter']:
            print("Error: The IMF type is not supported.")
            print("Please set the IMF type to be in ['Salpeter', 'Chabrier', 'Kroupa', 'Custom', 'PowerLaw', 'DietSalpeter'].")
        if self.IMF_type == 'Custom' and self.IMF_arr is None:
            print("Error: The IMF type is set to be 'Custom', but the IMF array is not provided.")
            print("Please provide the IMF array.")
        if self.IMF_type == 'PowerLaw' and self.power_index is None:
            print("Error: The IMF type is set to be 'PowerLaw', but the power-law index is not provided.")
            print("Please provide the power-law index.")
        if self.IMF_type != 'Custom' and self.IMF_arr is not None:
            print("Warning: The IMF type is not 'Custom', but the IMF array is provided.")
            print("The IMF array will be ignored.")
        if self.IMF_type != 'PowerLaw' and self.power_index is not None:
            print("Warning: The IMF type is not 'PowerLaw', but the power-law index is provided.")
            print("The power-law index will be ignored.")
        # self.imf is the one of the IMF functions that corresponds to the IMF type.
        # However, you can still use self.Salpeter, self.Chabrier, self.Kroupa.
        if self.IMF_type == 'Salpeter':
            self.imf = self.Salpeter
        elif self.IMF_type == 'Chabrier':
            self.imf = self.Chabrier
        elif self.IMF_type == 'Kroupa':
            self.imf = self.Kroupa
        elif self.IMF_type == 'Custom':
            self.A['Custom'] = quad(self._Custom, constants.Mstar_min, constants.Mstar_max)[0]
            self.imf = self.Custom
        elif self.IMF_type == 'PowerLaw':
            self.A['PowerLaw'] = quad(self._PowerLaw, constants.Mstar_min, constants.Mstar_max)[0]
            self.imf = self.PowerLaw
        elif self.IMF_type == 'DietSalpeter':
            self.imf = self.DietSalpeter
        else:
            pass


    # The Salpeter IMF from Salpeter (1955), https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S.

    
    def _Salpeter(self, m):
        """
        The Salpeter IMF from Salpeter (1955), https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S.

        Note that it is not normalized to 1 in the mass range [0.08, 150] solar mass.

        Parameters
        ----------
        m: float
            The mass of the star
        """
        if m<0.1 or m>100:
            # The mass range of Salpeter IMF is [0.1, 100] solar mass.
            return 0.0
            # Do not return 0, or np.vectorize will assume the return value is an integer
        else:
            return m**(-2.35)
    

    def Salpeter(self, m):
        """
        The Salpeter IMF from Salpeter (1955), https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S.

        Note that it is NORMALIZED to 1 in the mass range [0.08, 150] solar mass.

        Parameters
        ----------
        m: float
            The mass of the star
        
        Returns
        -------
        xi: float
            The Salpeter IMF of the star
        """
        # imf = np.vectorize(self._Salpeter, otypes=[np.float64])
        # Vectorization will speed down the calculation of quad, so we do not use it here.
        xi = self._Salpeter(m)/self.A['Salpeter']
        return xi


    # The Chabrier IMF from Chabrier (2003), https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C.

    
    def _Chabrier(self, m):
        """
        The Chabrier IMF from Chabrier (2003), https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C.

        Note that it is not normalized to 1 in the mass range [0.08, 150] solar mass.

        Parameters
        ----------
        m: float
            The mass of the star
        """
        A = np.exp((np.log10(0.08)**2/2/0.69**2))
        # To make the IMF continous at m=1.
        if m<constants.Mstar_min or m>constants.Mstar_max:
            return 0.0
        else:
            if m<1:
                return A*np.exp(-(np.log10(m)-np.log10(0.08))**2/2/0.69**2)/(m*np.log(10))
            else:
                return m**(-2.3)/np.log(10)
    

    def Chabrier(self, m):
        """
        The Chabrier IMF from Chabrier (2003), https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C.

        Note that it is NORMALIZED to 1 in the mass range [0.08, 150] solar mass.

        Parameters
        ----------
        m: float
            The mass of the star
        
        Returns
        -------
        xi: float
            The Chabrier IMF of the star
        """
        xi = self._Chabrier(m)/self.A['Chabrier']
        return xi
    

    # The Kroupa IMF from Kroupa (2001), https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K.


    def _Kroupa(self, m):
        """
        The Kroupa IMF from Kroupa (2001), https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K.

        Note that it is not normalized to 1 in the mass range [0.08, 150] solar mass.

        Parameters
        ----------
        m: float
            The mass of the star
        """
        if m<constants.Mstar_min or m>constants.Mstar_max:
            return 0.0
        else:
            if m<0.08:
                # 1/0.04 is used to insure the IMF is continous at m=0.08.
                return 1/0.04*m**(-0.3)
            elif m<0.5:
                # 1/0.5 is used to insure the IMF is continous at m=0.5.
                return 1/0.5*m**(-1.3)
            else:
                return m**(-2.3)
    

    def Kroupa(self, m):
        """
        The Kroupa IMF from Kroupa (2001), https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K.

        Note that it is NORMALIZED to 1 in the mass range [0.08, 150] solar mass.

        Parameters
        ----------
        m: float
            The mass of the star
        
        Returns
        -------
        xi: float
            The Kroupa IMF of the star
        """
        xi = self._Kroupa(m)/self.A['Kroupa']
        return xi
    

    # The custom IMF, which is set by the user.
    
    def _Custom(self, m):
        """
        The custom IMF, which is set by the user.

        Note that it is not normalized to 1 in the mass range [0.08, 150] solar mass.

        Parameters
        ----------
        m: float
            The mass of the star
        """
        n = self.IMF_arr.shape[0]
        if n<500:
            print("Warning: The number of the mass bins is less than 500, which is too few.")
            print("The program will move on but a large n is recommended to get a smooth IMF.")
        if m<constants.Mstar_min or m>constants.Mstar_max:
            return 0.0
        else:
            interp = interpolate.interp1d(self.IMF_arr[:, 0], self.IMF_arr[:, 1], kind=1, fill_value='extrapolate')
            return interp(m)
    

    def Custom(self, m):
        """
        The custom IMF, which is set by the user.

        Note that it is NORMALIZED to 1 in the mass range [0.08, 150] solar mass.

        Parameters
        ----------
        m: float
            The mass of the star
        
        Returns
        -------
        xi: float
            The custom IMF of the star
        """
        xi = self._Custom(m)/self.A['Custom']
        return xi
    

    # The power-law IMF.


    def _PowerLaw(self, m):
        """
        The power-law IMF.

        Note that it is not normalized to 1 in the mass range [0.08, 150] solar mass.

        Parameters
        ----------
        m: float
            The mass of the star
        """
        if m<constants.Mstar_min or m>constants.Mstar_max:
            return 0.0
        else:
            return m**(-self.power_index)
    

    def PowerLaw(self, m):
        """
        The power-law IMF.

        Note that it is NORMALIZED to 1 in the mass range [0.08, 150] solar mass.

        Parameters
        ----------
        m: float
            The mass of the star
        
        Returns
        -------
        xi: float
            The power-law IMF of the star
        """
        xi = self._PowerLaw(m)/self.A['PowerLaw']
        return xi
    

    # The diet Salpeter IMF.


    def _DietSalpeter(self, m):
        """
        The diet Salpeter IMF.
        The slope of mass below 0.35 solar mass is set to be 1, 
        while the slope of mass above 0.35 solar mass is set to be the same as the Salpeter IMF (2.35).

        Note that it is not normalized to 1 in the mass range [0.08, 150] solar mass.

        Parameters
        ----------
        m: float
            The mass of the star
        """
        if m<constants.Mstar_min or m>constants.Mstar_max:
            return 0.0
        elif m<0.35:
            return 0.35**(-1.35)*m**(-1)
        else:
            return m**(-2.35)
    

    def DietSalpeter(self, m):
        """
        The diet Salpeter IMF.

        Note that it is NORMALIZED to 1 in the mass range [0.08, 150] solar mass.

        Parameters
        ----------
        m: float
            The mass of the star
        
        Returns
        -------
        xi: float
            The diet Salpeter IMF of the star
        """
        imf = np.vectorize(self._DietSalpeter, otypes=[np.float64])
        xi = imf(m)/self.A['DietSalpeter']
        return xi
    

###### The IMF class ######


###### Extrapolate and interpolate the yields ######
# The miminum and maximum mass of the star are 0.08 and 150 solar mass, respectively.

# However, the miminum mass available in the yield table is uaually larger than 0.08 solar mass, like 1.0 solar mass.
# The remnant mass of the stars with mass less than the miminum mass available in the yield table is set to be 
# its initial mass. Namely, they are assumed to experience no mass loss.
# For consistency, their yields are set to be ZERO.

# For the stars with mass larger than the maximum mass available in the yield table,
# we keep the Mini - Mrem the same as that one of the larest mass available in the yield table.
# Namely, the additional mass is all locked in the remnant.
# For consistency, their yields are set to be the SAME as that of the larest mass available in the yield table.

def extra_interpolate_yields(x1, y1, yields_type, x):
    """
    Extrapolate and interpolate the yields.

    Parameters
    ----------
    x1: numpy array
        The initial mass of the stars in the unit of solar mass.
    y1: numpy array
        The yields of the stars.
    yields_type: str
        The type of the yields. It should be 'Mrem' or an element name.
    x: float
        The initial mass of the stars in the unit of solar mass at which the yields are evaluated.
    Returns
    -------
    y: float
        Extra and interpolated yields of the stars.
    """
    interp = interpolate.interp1d(x1, y1, kind=1, fill_value='extrapolate')
    if yields_type == 'Mrem':
        # The remnant mass of the stars with mass less than the miminum mass available in the yield table is set to be 
        # its initial mass. Namely, they are assumed to experience no mass loss.
        # For consistency, their yields are set to be ZERO.
        if (x>x1.min()) and (x<x1.max()):
            y = interp(x)
        elif x <= x1.min():
            y = x
        else:
            y = x - x1.max() + y1[x1.argmax()]
    else:
        # For the stars with mass larger than the maximum mass available in the yield table,
        # we keep the Mini - Mrem the same as that one of the larest mass available in the yield table.
        # Namely, the additional mass is all locked in the remnant.
        # For consistency, their yields are set to be the SAME as that of the larest mass available in the yield table.
        if (x>x1.min()) and (x<x1.max()):
            y = interp(x)
        elif x <= x1.min():
            y = 0
        else:
            y = y1[x1.argmax()]
    return y


### This part is tested on April 30, 2024.
###### Extrapolate and interpolate the yields ######

###### Get the rate of SNIa per stellar mass ######
# The rate of SNIa per stellar mass can be prescribed as 
# the Equation 13 in Maoz & Mannucci (2012), https://ui.adsabs.harvard.edu/abs/2012PASA...29..447M/abstract
# which is a constraint from observation.
# Note that it assumes a diet Salpeter IMF, which should be calibrated to the adotped IMF.
# The calibration factor is the Equation 4 in Yan et al. (2019), https://www.aanda.org/10.1051/0004-6361/201936029

class SNIa:
    """
    The class of the rate of SNIa per stellar mass.
    """


    def __init__(self, imf):
        """
        Initialize the class of the rate
        
        Parameters
        ----------
        imf: function
            The IMF function.
        """
        self.imf = imf
        imf_diet = IMF(IMF_type='DietSalpeter').imf  # The IMF function of the diet Salpeter IMF.
        p_denominator = (quad(imf_diet, 1.5, 8)[0])**2/\
                        (quad(imf_diet, 0.08, 150)[0]*quad(lambda m: m*imf_diet(m), 0.08, 150)[0])
        p_numerator = (quad(imf, 1.5, 8)[0])**2/\
                        (quad(imf, 0.08, 150)[0]*quad(lambda m: m*imf(m), 0.08, 150)[0])
        self.p = p_numerator/p_denominator
        Denominator=4.819420e-04
        Numerator=1.575560e-04
        # self.p1 is taken from Yan's code
        self.p1 = Numerator/Denominator
    def _rate(self, t):
        """
        The original constraint on the rate of SNIa per stellar mass from Maoz & Mannucci (2012),
        https://ui.adsabs.harvard.edu/abs/2012PASA...29..447M/abstract
        Maoz & Mannucci (2012) adoptedthe diet Salpeter IMF.
        Within this function, it has not been calibrated to the adotped IMF.
        The calibration will be done in rate().

        Parameters
        ----------
        t: float
            The age of the stellar population in the unit of year.
        
        Returns
        -------
        rate: float
            The rate of SNIa per stellar mass.
            The unit is 1/yr/solar mass.
        """
        if t<4e7:
            # t<40 Myr
            y = 0.0
        else:
            y = 4e-13*(t/1e9)**(-1)
        # 4e-13 is used to normalize the intergral from 40 Myr to 10 Gyr to be 2.2e-3 Msun^{-1}.
        # Note that there no need to implement a cutoff when t>10 Gyr, 
        # for the reason that the uncertainty of the observational result is 50% (see Table 1 of Maoz & Mannucci 2012).
        # The uncertainty due the cutoff is relatively small.
        # For simplicity, we do not implement the cutoff when t>10 Gyr.
        # There is also no cutoff in the code of Yan et al. (2019), https://www.aanda.org/10.1051/0004-6361/201936029.
        return y


    def rate(self, t):
        """
        The calibrated constraint on the rate of SNIa per stellar mass from Maoz & Mannucci (2012),
        https://ui.adsabs.harvard.edu/abs/2012PASA...29..447M/abstract

        Parameters
        ----------
        t: float
            The age of the stellar population in the unit of year.
        imf: function
            The IMF function.

        Returns
        -------
        rate: float
            The rate of SNIa per stellar mass.
            The unit is 1/yr/solar mass.
        """
        # y = self._rate(t)*self.p
        # self.p1 is taken from Yan's code
        y = self._rate(t)*self.p1
        return y
    

    def number(self, Mstar, tmin, tmax):
        """
        The number of SNIa in the time interval [tmin, tmax] for the stellar population with the total mass Mstar.

        Parameters
        ----------
        Mstar: float
            The total mass of the stellar population in the unit of solar mass.
        tmin: float
            The minimum age of the stellar population in the unit of year.
        tmax: float
            The maximum age of the stellar population in the unit of year.

        Returns
        -------
        N: float
            The number of SNIa in the time interval [tmin, tmax] for the stellar population with the total mass Mstar.
        """
        N = Mstar*quad(self.rate, tmin, tmax)[0]
        return N


###### Get the rate of SNIa per stellar mass ######