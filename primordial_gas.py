# -*- coding: utf-8 -*-
"""
Created on May 11, 2024
Last modified on May 14, 2024
@Author: Guan-Fu Liu

Define the Primordial_gas class
"""

import constants
import numpy as np


###### The primordial gas mass class ######
class Primordial_gas:
    """
    Define the primordial gas class.
    """
    def __init__(self, mass, Z_0=None, input_array=None):
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
        input_array: numpy array or None
            If it is None, use the default primordial gas from the table II in Cyburt et al. (2016).
            If it is a numpy array, it shuld have a shape of (len(constants.elem_names), ).
            input_array[0]: should be zero, which is just to make the index start from 1 and should NOT be used.
            input_array[1:len(constants.elem_names)-1]: the mass fraction of H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl,
            Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, and Zn. Namely, input_array[i] is the ith element mass fraction.
            input_array[len(constants.elem_names)-1]: the mass fraction of the rest elements ((len(constants.elem_names)-1)th and after) 
            in the primordial gas, should be a small number.
            input_array will be normalized by input_array[1:] = input_array[1:]/input_array[1:].sum().

            Z_0==None and input_array==None: use the default primordial gas.
            Z_0!=None and input_array==None: use the default primordial gas but with the pre-set primordial metallicity
              (without lithium).
            Z_0==None and input_array!=None: use the input_array as the primordial gas.
            Z_0!=None and input_array!=None: use the input_array as the primordial gas but raise a warning, which is exactly
                the same as the case of Z_0==None and input_array!=None but not recommended.
        Returns
        -------
        self.primor_gas: dictionary

        self.primor_gas will contain the following keys:
        'Gas': numpy array
            The mass of different elements in the primordial gas in the unit of solar mass.
            The shape is (len(constants.elem_names), ), where the first element is zero, which is just to make the index start from 1.
        'H': float
            The mass of H in the primordial gas in the unit of solar mass.
        'He': float
            The mass of He in the primordial gas in the unit of solar mass.
        'Li': float
            The mass of Li in the primordial gas in the unit of solar mass.
        ......
        'Zn': float
            The mass of Zn in the primordial gas in the unit of solar mass.
        'Other': float
            The mass of the rest elements ((len(constants.elem_names)-1)th and after) in the primordial gas in the unit of solar mass.
        
        Examples
        --------
        >>> import primordial_gas as pg
        >>> import constants
        >>> solar_set = "Default"  # The solar abundance set
        >>> mass = 2000  # The mass of the gas in the unit of Msun
        >>> pr_gas = pg.Primordial_gas(mass=mass, Z_0=0).add_metals(abund_table=constants.abund_tables[solar_set])
        >>> print(pr_gas)
        {'H': 1505.9849458627261, 'He': 494.0150534372907, 'Li': 6.999831473577275e-07, 'Z': 6.999831473577275e-07, 
        'Be': 0.0, 'B': 0.0, 'C': 0.0, 'N': 0.0, 'O': 0.0, 'F': 0.0, 'Ne': 0.0, 'Na': 0.0, 'Mg': 0.0, 
        'Al': 0.0, 'Si': 0.0, 'P': 0.0, 'S': 0.0, 'Cl': 0.0, 'Ar': 0.0, 'K': 0.0, 'Ca': 0.0, 'Sc': 0.0, 
        'Ti': 0.0, 'V': 0.0, 'Cr': 0.0, 'Mn': 0.0, 'Fe': 0.0, 'Co': 0.0, 'Ni': 0.0, 'Cu': 0.0, 'Zn': 0.0, 
        'Other': 0.0, 'Gas': array([0.00000000e+00, 1.50598495e+03, 4.94015053e+02, 6.99983147e-07,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])}
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

        self.Z_0 = Z_0
        self.input_array = input_array

        if input_array is None:
            # input_array is None, use the default primordial gas.
            self.primor_gas = primor_gas.copy()
            if self.primor_gas['Z_0'] > 1e-4:
                print("Warning: The primordial metallicity is larger than 1e-4, which is not recommended.")
                print("It will move with such a large primordial metallicity.")
            else:
                pass
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
            for key in self.primor_gas.keys():
                self.primor_gas[key] = self.primor_gas[key] * mass

            # Remove the items that are not needed.
            self.primor_gas.pop('X')
            self.primor_gas.pop('D/H')
            self.primor_gas.pop('Y_p')
            self.primor_gas.pop('3He/H')
            self.primor_gas.pop('7Li/H')
            self.primor_gas.pop('6Li/H')
            self.primor_gas.pop('Z')

        else:
            if len(input_array) != len(constants.elem_names):
                raise ValueError("The length of the input_array should be %d, but it is not."%len(constants.elem_names))
            # input_array is not None and pre-set, use the input_array as the primordial gas.
            input_array[0] = 0  # The first element should be zero, which is just to make the index start from 1.
            input_array = input_array/input_array[1:].sum()
            self.primor_gas = { }
            self.primor_gas['Gas'] = np.zeros(len(constants.elem_names))
            for i, elem in enumerate(constants.elem_names[1:]):
                self.primor_gas[elem] = input_array[i+1] * mass
                self.primor_gas['Gas'][i+1] = self.primor_gas[elem]
            
            if Z_0 is not None:
                print("Warning: The primordial metallicity is pre-set, which is not recommended.")
                print("The primordial metallicity will be ignored and only the input array will be used.")
    

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
        It should be an array with the shape (len(constants.elem_names), ), which contains the solar abundances of the first 30 elements.
        The first element in the array should be ZERO, which is just to make the index start from 1.
        For more details, please refer to the module constants.py.
        """
        if self.input_array is not None:
            print("Warning: The input array is pre-set, which is not recommended.")
            print("Nothing will be done in this function.")
            print("Just use the input array as the primordial gas.")
            return self.primor_gas
        if self.Z_0 is None:
            raise ValueError("The primordial metallicity is should be pre-set, but it is not.")
        for i, elem in enumerate(constants.elem_names[4:]):
            self.primor_gas[elem] = abund_table[i+4]/abund_table[4:].sum() * self.primor_gas['Z_0']
            # Suppose the composition of the metals except lithium in the primordial gas can be scaled by the solar abundance.
            # Namely, the relative abundance between the metals except is the same as the solar abundance,
            # but the total mass fraction of the metals is Z_0.
            # Such as O/Fe is the same as the solar abundance, but the total mass fraction of the metals is Z_0.
        self.primor_gas.pop('Z_0')  # Remove the primordial metallicity Z_0 without lithium from the dictionary.
        self.primor_gas['Gas'] = np.zeros(len(constants.elem_names))  # The mass of the primordial gas in the unit of solar mass.
        for i in range(1, len(constants.elem_names)-1):
            self.primor_gas['Gas'][i] = self.primor_gas[constants.elem_names[i]]
            # The index starts from 1, so the first element is ZERO.
        return self.primor_gas


### This part is tested on May 14, 2024.
###### The primordial gas mass class ######