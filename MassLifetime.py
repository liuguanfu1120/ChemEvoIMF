# -*- coding: utf-8 -*-
"""
Created on May 20, 2024
Last modified on May 20, 2024
@Author: Guan-Fu Liu

Define the MassLifetime class.
"""
import numpy as np
from scipy import interpolate


class MassLifetime:
    """
    To define the MassLifetime class.
    """
    def __init__(self, lifetime=None):
        # Set the default lifetime.
        # The following data is taken from Table 14 in Portinari et al. (1998), 
        # https://ui.adsabs.harvard.edu/abs/1998A&A...334..505P.
        default_lifetime = { }
        default_lifetime['Z=0.0004'] = { }
        # The initial mass, first column in Table 14 in Portinari et al. (1998).
        default_lifetime['Z=0.0004']['Mini'] = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                                        1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
                                                        1.8, 1.9, 2.0, 2.5, 3.0, 4.0,
                                                        5.0, 6.0, 7.0, 9.0, 12.0, 15.0,
                                                        20.0, 30.0, 40.0, 60.0, 100.0, 120.0])
        # The lifetime of the stars with metallicity Z=0.0004, second column in Table 14 in Portinari et al. (1998).
        default_lifetime['Z=0.0004']['lifetime'] = np.array([4.28e10, 2.37e10, 1.41e10, 8.97e9, 6.03e9, 4.23e9,
                                                            3.08e9, 2.34e9, 1.92e9, 1.66e9, 1.39e9, 1.18e9,
                                                            1.11e9, 9.66e8, 8.33e8, 4.64e8, 3.03e8, 1.61e8,
                                                            1.01e8, 7.15e7, 5.33e7, 3.42e7, 2.13e7, 1.54e7,
                                                            1.06e7, 6.90e6, 5.45e6, 4.20e6, 3.32e6, 3.11e6])
        default_lifetime['Z=0.004'] = { }
        # The initial mass, first column in Table 14 in Portinari et al. (1998).
        default_lifetime['Z=0.004']['Mini'] = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 
                                                        1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 
                                                        1.8, 1.9, 2.0, 2.5, 3.0, 4.0, 
                                                        5.0, 6.0, 7.0, 9.0, 12.0, 15.0, 
                                                        20.0, 30.0, 40.0, 60.0, 100.0, 120.0])
        # The lifetime of the stars with metallicity Z=0.004, third column in Table 14 in Portinari et al. (1998).
        default_lifetime['Z=0.004']['lifetime'] = np.array([5.35e10, 2.95e10, 1.73e10, 1.09e10, 7.13e9, 4.93e9,
                                                            3.52e9, 2.64e9, 2.39e9, 1.95e9, 1.63e9, 1.28e9,
                                                            1.25e9, 1.23e9, 1.08e9, 5.98e8, 3.67e8, 1.82e8,
                                                            1.11e8, 7.62e7, 5.61e7, 3.51e7, 2.14e7, 1.52e7,
                                                            1.05e7, 6.85e6, 5.44e6, 4.19e6, 3.38e6, 3.23e6])

        default_lifetime['Z=0.008'] = { }
        # The initial mass, first column in Table 14 in Portinari et al. (1998).
        default_lifetime['Z=0.008']['Mini'] = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                                        1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
                                                        1.8, 1.9, 2.0, 2.5, 3.0, 4.0,
                                                        5.0, 6.0, 7.0, 9.0, 12.0, 15.0,
                                                        20.0, 30.0, 40.0, 60.0, 100.0, 120.0])
        # The lifetime of the stars with metallicity Z=0.008, fourth column in Table 14 in Portinari et al. (1998).
        default_lifetime['Z=0.008']['lifetime'] = np.array([6.47e10, 3.54e10, 2.09e10, 1.30e10, 8.46e9, 5.72e9,
                                                            4.12e9, 2.92e9, 2.36e9, 2.18e9, 1.82e9, 1.58e9,
                                                            1.41e9, 1.25e9, 1.23e9, 6.86e8, 4.12e8, 1.93e8,
                                                            1.15e8, 7.71e7, 5.59e7, 3.44e7, 2.10e7, 1.49e7,
                                                            1.01e7, 6.65e6, 5.30e6, 4.15e6, 3.44e6, 3.32e6])
        default_lifetime['Z=0.02'] = { }
        # The initial mass, first column in Table 14 in Portinari et al. (1998).
        default_lifetime['Z=0.02']['Mini'] = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                                        1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
                                                        1.8, 1.9, 2.0, 2.5, 3.0, 4.0,
                                                        5.0, 6.0, 7.0, 9.0, 12.0, 15.0,
                                                        20.0, 30.0, 40.0, 60.0, 100.0, 120.0])
        # The lifetime of the stars with metallicity Z=0.02, fifth column in Table 14 in Portinari et al. (1998).
        default_lifetime['Z=0.02']['lifetime'] = np.array([7.92e10, 4.45e10, 2.61e10, 1.59e10, 1.03e10, 6.89e9,
                                            4.73e9, 3.59e9, 2.87e9, 2.64e9, 2.18e9, 1.84e9, 
                                            1.59e9, 1.38e9, 1.21e9, 7.64e8, 4.56e8, 2.03e8,
                                            1.15e8, 7.45e7, 5.31e7, 3.17e7, 1.89e7, 1.33e7,
                                            9.15e6, 6.13e6, 5.12e6, 4.12e6, 3.39e6, 3.23e6])
        default_lifetime['Z=0.05'] = { }
        # The initial mass, first column in Table 14 in Portinari et al. (1998).
        default_lifetime['Z=0.05']['Mini'] = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                                        1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
                                                        1.8, 1.9, 2.0, 2.5, 3.0, 4.0,
                                                        5.0, 6.0, 7.0, 9.0, 12.0, 15.0,
                                                        20.0, 30.0, 40.0, 60.0, 100.0, 120.0])
        # The lifetime of the stars with metallicity Z=0.05, sixth column in Table 14 in Portinari et al. (1998).
        default_lifetime['Z=0.05']['lifetime'] = np.array([7.18e10, 4.00e10, 2.33e10, 1.42e10, 8.88e9, 5.95e9,
                                            4.39e9, 3.37e9, 3.10e9, 2.51e9, 2.06e9, 1.76e9,
                                            1.51e9, 1.34e9, 1.24e9, 6.58e8, 3.81e8, 1.64e8,
                                            8.91e7, 5.67e7, 3.97e7, 2.33e7, 1.39e7, 9.95e6,
                                            6.99e6, 5.15e6, 4.34e6, 3.62e6, 3.11e6, 3.10e6])
        # If the lifetime is given, we use the given lifetime.
        # The following code is used to check the format of the lifetime.
        flag = True
        if not isinstance(lifetime, dict):
            flag = False
            print("The lifetime should be a dictionary.")
        else:
            for key in lifetime.keys():
                if key[:2] != 'Z=':
                    flag = False
                    print("The key of the lifetime dictionary should start with 'Z='.")
                    break
                try:
                    float(key[2:])
                except:
                    flag = False
                    print("The metallicity should be a float number.")
                    break
                if not isinstance(lifetime[key], dict):
                    flag = False
                    print("The value of the lifetime dictionary should be a dictionary.")
                    break
                if 'Mini' not in lifetime[key].keys():
                    flag = False
                    print("The key 'Mini' should be in the dictionary life[%s]."%key)
                    break
                if 'lifetime' not in lifetime[key].keys():
                    flag = False
                    print("The key 'lifetime' should be in the dictionary life[%s]."%key)
                    break
                argsort = lifetime[key]['Mini'].argsort()
                lifetime[key]['lifetime'] = lifetime[key]['lifetime'][argsort]
                lifetime[key]['Mini'] = lifetime[key]['Mini'][argsort]
                if (np.diff(lifetime[key]['lifetime'])).all()>0:
                    # The lifetime of less massive stars should be longer.
                    continue
                else:
                    flag = False
                    print("The lifetime of less massive stars should be longer.")
                    print("Please check the mass lifetime relation of %s" % key)

        if lifetime is None:
            self.lifetime = default_lifetime
        else:
            if flag:
            # flag is True if the format of the lifetime is correct.
                self.lifetime = lifetime
            else:
            # flag is False if the format of the lifetime is not correct. Use the default lifetime instead.
                print("The format of the input lifetime is not correct. Use the default lifetime instead.")
                self.lifetime = default_lifetime
        self.mass_to_lifetime = np.vectorize(self._mass_to_lifetime, otypes=[np.float64])
        self.lifetime_to_mass = np.vectorize(self._lifetime_to_mass, otypes=[np.float64])


    def _mass_to_lifetime(self, mass, Z, kind='linear-log'):
        """
        Get the lifetime of the star with a certain mass and a certain metallicity.

        Parameters
        ----------
        mass: float
            The mass of the star in the unit of solar mass.
        Z: float
            The initial metallicity of the star.
        kind: str, optional
            The kind of interpolation. It should be in 
            ['linear-linear', 'linear-log', 'log-linear', 'log-log', 'nearest', 'TNG-like'].
            The default is 'linear-log'.
            For more details see InterpolateYields.py        
        Returns
        -------
        t: float
            The lifetime of the star in the unit of year.
        """
        # The keys of lifetime dict (except the key "Mini")
        Z_list = [key for key in self.lifetime.keys() if 'Z=' in key]
        Z_float = np.array([float(a.replace("Z=", "")) for a in Z_list])
        Z_index = { }

        if Z<=Z_float.min():
            Z_index['low'] = Z_float.argmin()
            Z_index['high'] = Z_index['low']
        elif Z>=Z_float.max():
            Z_index['low'] = Z_float.argmax()
            Z_index['high'] = Z_index['low']
        else:
            Z_index['low'] = np.where(Z_float==Z_float[Z_float<=Z].max())[0][0]
            Z_index['high'] = np.where(Z_float==Z_float[Z_float>=Z].min())[0][0]

        Z_key = { }
        Z_key['low'] = Z_list[Z_index['low']]
        Z_key['high'] = Z_list[Z_index['high']]

        if kind == 'nearest':
            Z_index['low'] = np.abs(Z_float-Z).argmin()
            Z_index['high'] = Z_index['low']
            interp = interpolate.interp1d(np.log10(self.lifetime[Z_key['low']]['Mini']),
                                          np.log10(self.lifetime[Z_key['low']]['lifetime']),
                                          kind=1, fill_value='extrapolate')
            t = interp(np.log10(mass))
            t = 10**t
            return t
        else:
            pass

        if Z_key['low'] == Z_key['high']:
            interp = interpolate.interp1d(np.log10(self.lifetime[Z_key['low']]['Mini']),
                                          np.log10(self.lifetime[Z_key['low']]['lifetime']),
                                          kind=1, fill_value='extrapolate')
            t = interp(np.log10(mass))
            t = 10**t
            return t
        else:
            interp_low = interpolate.interp1d(np.log10(self.lifetime[Z_key['low']]['Mini']),
                                              np.log10(self.lifetime[Z_key['low']]['lifetime']),
                                              kind=1, fill_value='extrapolate')
            interp_high = interpolate.interp1d(np.log10(self.lifetime[Z_key['high']]['Mini']),
                                               np.log10(self.lifetime[Z_key['high']]['lifetime']),
                                               kind=1, fill_value='extrapolate')
            t_low = interp_low(np.log10(mass))
            t_high = interp_high(np.log10(mass))
            t_low = 10**t_low
            t_high = 10**t_high
        
        Z_low = Z_float[Z_index['low']]
        Z_high = Z_float[Z_index['high']]
        # In case of Z_low is zero, we set Z_low = 1e-10.
        Z_low = np.maximum(Z_low, 1e-10)
        # In case of Z_high is zero, we set Z_high = 1e-10.
        Z_high = np.maximum(Z_high, 1e-10)
        
        if kind == 'linear-linear':
            t = t_low + (t_high-t_low)/(Z_high-Z_low)*(Z-Z_low)
        elif kind == 'linear-log' or kind == 'TNG-like':
            t = t_low + (t_high-t_low)/(np.log10(Z_high)-np.log10(Z_low))\
                      * (np.log10(Z)-np.log10(Z_low))
        elif kind == 'log-linear':
            t = t_low**((Z_high-Z)/(Z_high-Z_low))*t_high**((Z-Z_low)/(Z_high-Z_low))
        elif kind == 'log-log':
            t = t_low**((np.log10(Z_high)-np.log10(Z))/(np.log10(Z_high)-np.log10(Z_low)))*\
                t_high**((np.log10(Z)-np.log10(Z_low))/(np.log10(Z_high)-np.log10(Z_low)))
        else:
            print("The kind of interpolation is not supported.")
            print("Please choose from ['linear-linear', 'linear-log', 'log-linear', 'log-log', 'nearest', 'TNG-like']")
            return None
        return t

    def _lifetime_to_mass(self, t, Z, kind='linear-log'):
        """
        Get the mass of the star with a certain lifetime and a certain metallicity.
        
        Parameters
        ----------
        t: float
            The lifetime of the star in the unit of year.
        Z: float
            The initial metallicity of the star.
        kind: str, optional
            The kind of interpolation. It should be in 
            ['linear-linear', 'linear-log', 'log-linear', 'log-log', 'nearest', 'TNG-like'].
            The default is 'nearest'.
            For more details see InterpolateYields.py
        """
        # The keys of lifetime dict (except the key "Mini")
        Z_list = [key for key in self.lifetime.keys() if 'Z=' in key]
        Z_float = np.array([float(a.replace("Z=", "")) for a in Z_list])
        Z_index = { }

        if Z<=Z_float.min():
            Z_index['low'] = Z_float.argmin()
            Z_index['high'] = Z_index['low']
        elif Z>=Z_float.max():
            Z_index['low'] = Z_float.argmax()
            Z_index['high'] = Z_index['low']
        else:
            Z_index['low'] = np.where(Z_float==Z_float[Z_float<=Z].max())[0][0]
            Z_index['high'] = np.where(Z_float==Z_float[Z_float>=Z].min())[0][0]

        Z_key = { }
        Z_key['low'] = Z_list[Z_index['low']]
        Z_key['high'] = Z_list[Z_index['high']]

        if kind == 'nearest':
            Z_index['low'] = np.abs(Z_float-Z).argmin()
            Z_index['high'] = Z_index['low']
            interp = interpolate.interp1d(np.log10(self.lifetime[Z_key['low']]['lifetime']),
                                          np.log10(self.lifetime[Z_key['low']]['Mini']),
                                          kind=1, fill_value='extrapolate')
            mass = interp(np.log10(t))
            mass = 10**mass
            return mass
        else:
            pass

        if Z_key['low'] == Z_key['high']:
            interp = interpolate.interp1d(np.log10(self.lifetime[Z_key['low']]['lifetime']),
                                          np.log10(self.lifetime[Z_key['low']]['Mini']),
                                          kind=1, fill_value='extrapolate')
            mass = interp(np.log10(t))
            mass = 10**mass
            return mass
        else:
            interp_low = interpolate.interp1d(np.log10(self.lifetime[Z_key['low']]['lifetime']),
                                              np.log10(self.lifetime[Z_key['low']]['Mini']),
                                              kind=1, fill_value='extrapolate')
            interp_high = interpolate.interp1d(np.log10(self.lifetime[Z_key['high']]['lifetime']),
                                               np.log10(self.lifetime[Z_key['high']]['Mini']),
                                               kind=1, fill_value='extrapolate')
            mass_low = interp_low(np.log10(t))
            mass_high = interp_high(np.log10(t))
            mass_low = 10**mass_low
            mass_high = 10**mass_high

        Z_low = Z_float[Z_index['low']]
        Z_high = Z_float[Z_index['high']]
        # In case of Z_low is zero, we set Z_low = 1e-10.
        Z_low = np.maximum(Z_low, 1e-10)
        # In case of Z_high is zero, we set Z_high = 1e-10.
        Z_high = np.maximum(Z_high, 1e-10)
        
        if kind == 'linear-linear':
            mass = mass_low + (mass_high-mass_low)/(Z_high-Z_low)*(Z-Z_low)
        elif kind == 'linear-log' or kind == 'TNG-like':
            mass = mass_low + (mass_high-mass_low)/(np.log10(Z_high)-np.log10(Z_low))\
                   *(np.log10(Z)-np.log10(Z_low))
        elif kind == 'log-linear':
            mass = mass_low**((Z_high-Z)/(Z_high-Z_low))*mass_high**((Z-Z_low)/(Z_high-Z_low))
        elif kind == 'log-log':
            mass = mass_low**((np.log10(Z_high)-np.log10(Z))/(np.log10(Z_high)-np.log10(Z_low)))*\
                   mass_high**((np.log10(Z)-np.log10(Z_low))/(np.log10(Z_high)-np.log10(Z_low)))
        else:
            print("The kind of interpolation is not supported.")
            print("Please choose from ['linear-linear', 'linear-log', 'log-linear', 'log-log', 'nearest', 'TNG-like']")
            return None
        
        return mass


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