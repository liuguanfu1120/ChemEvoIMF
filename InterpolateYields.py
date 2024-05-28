# -*- coding: utf-8 -*-
"""
Created on May 27, 2024
Last modified on May 27, 2024
@Author: Guan-Fu Liu

To define the class of interpolating the yields.
"""
from scipy import interpolate
import numpy as np
class InterpolateYields:
    """
    Define the class of interpolating the yields.
    For AGB and SNcc, the yield of the element i is denoted as Y_{i}(Z, m).
    Y_{i}(Z, m) is the yield of the element i for a star with mass m and metallicty Z.
    The yield table only gives a few values of Z and m, and hence some interpolation is needed.
    For a given metallicity Z and a given element, we interpolate the yield of the element with respect to the mass of the star.
    For a given element and mass, we interpolate the yield with respect to the metallicity, which could be serval methods 
    (for details, see below).
    For SNIa, the yield of the element i is denoted as Y_{i}(Z).
    Therefore, the interpolation is only conducted with respect to the metallicity.
    The interpolation with respect to the mass is linear interpolation.
    """
    def __init__(self, dfs, kind='linear-linear'):
        """
        Initialize the class that enables the interpolation of the yields.

        Parameters
        ----------
        dfs : dict
            The dictionary containing the yield tables.
            Its keys are dict_keys(['AGB+SNcc', 'SNIa']), and the values are something like:
            dfs['AGB+SNcc']: dict_keys(['Z=0.0001', 'Z=0.0003', 'Z=0.001', 'Z=0.02', 'Z=0.05']), and
            dfs['SNIa']: dict_keys(['Z=0.0002', 'Z=0.002', 'Z=0.01', 'Z=0.02'])).
            The value of the metallicities may vary according to the yield table.
            dfs['AGB+SNcc']['Z=0.0001'] is a DataFrame
        kind : str
            The kind of interpolation. It should be in ['linear-linear', 'linear-log', 'log-linear', 'log-log', 'nearest', 'TNG-like'].
            The default is 'linear-linear'.
            For linear-linear interpolation:
                Y = Y_low + (Y_high - Y_low)/(Z_high - Z_low)*(ZGas - Z_low)
                Y and Z are linearly interpolated.
            For linear-log interpolation:
                Y = Y_low + (Y_high - Y_low)/(log(Z_high) - log(Z_low))*(log(ZGas) - log(Z_low))
                Y and log Z are linearly interpolated.
                In case of Z_low is zero, we set Z_low = 1e-10.
                This is the same as what is done in TNG simulation.
            For log-linear interpolation:
                log Y = log Y_low + (log Y_high - log Y_low)/(Z_high - Z_low)*(ZGas - Z_low)
                Y = 10**(log10(Y_low) + (log10(Y_high) - log10(Y_low))/(Z_high - Z_low)*(ZGas - Z_low))
                log Y and Z are linearly interpolated.
                In case of Y_low or Y_high is zero, we re-write the above equation as:
                Y = Y_low ** ((Z_high-ZGas)/(Z_high-Z_low)) * Y_high ** ((ZGas-Z_low)/(Z_high-Z_low))
            For log-log interpolation:
                log Y = log Y_low + (log Y_high - log Y_low)/(log Z_high - log Z_low)*(log ZGas - log Z_low)
                Y = 10**(log10(Y_low) + (log10(Y_high) - log10(Y_low))/(log(Z_high) - log(Z_low))*(log(ZGas) - log(Z_low))
                In case of Y_low or Y_high is zero, we re-write the above equation as:
                Y = Y_low ** ((log(Z_high)-log(ZGas))/(log(Z_high)-log(Z_low))) * 
                    Y_high ** ((log(ZGas)-log(Z_low))/(log(Z_high)-log(Z_low)))
                In case of Z_low is zero, we set Z_low = 1e-10.
            For nearest interpolation:
                Y = Y_nearest
            For TNG-like interpolation:
                The interpolation is the same as the linear-log interpolation.

            If ZGas is less than the lowest metallicity in the yield table, the lowest metallicity is used.
            If ZGas is larger than the highest metallicity in the yield table, the highest metallicity is used.

            If there is only one metallicity in the yield table (as is often the case for SNIa), the interpolation is not needed
            We will use the nearest interpolation.
        
        Returns
        -------
        self.interps : dict
            The dictionary containing the interpolation functions.
            self.interps['AGB+SNcc'] is the interpolation function for the yield of AGB+SNcc.
            self.interps['SNIa'] is the interpolation function for the yield of SNIa.
        
        Examples
        --------
        >>> interp_Y = InterpolateYields(dfs, kind=kind)
        >>> interp = interp_Y.interps['AGB+SNcc'](ZGas[0], elem)
        >>> dEjectElement[j, ElemIndex] = Nstar[0]*(quad(lambda m:
                                        IMF.imf(m)*interp(m),
                                        mass_bounds[j+1], mass_bounds[j],
                                        epsrel=1e-5, limit=80, full_output=1)[0])
    
        ###### interp = interp_Y.interps['AGB+SNcc'](ZGas[0], elem) is a MUST! ######
        # You SHOULD NOT use the following code:
        # dEjectElement[j, ElemIndex] = Nstar[0]*(quad(lambda m:
        #                                 IMF.imf(m)*interp_Y.interps['AGB+SNcc'](ZGas[0], elem)(m),
        #                                 mass_bounds[j+1], mass_bounds[j],
        #                                 epsrel=1e-5, limit=80, full_output=1)[0])
        # It will slow down the calculation.
        # quad will call the interp_Y.interps['AGB+SNcc'](ZGas[0], elem) again and again, which is not necessary
        # and slow down the calculation.
        ###### interp = interp_Y.interps['AGB+SNcc'](ZGas[0], elem) is a MUST! ######
        
        """
        self.dfs = dfs
        groups = {key: list(dfs[key].keys()) for key in dfs.keys()}
        if "Z_" in groups['AGB+SNcc'][0]:
            # for yields2.h5, the group names cannot contain "." and "_".
            # Therefore their group names are like "Z_0_0001", "Z_0_0003", etc.
            groups_Z = {key: [group.replace("Z_","Z=").replace("_",".") for group in groups[key]] for key in groups.keys()}
        else:
            groups_Z = groups
        # groups_Z is a dictionary like:
        # {'AGB+SNcc': ['Z=0.0001', 'Z=0.0003', 'Z=0.001', 'Z=0.002',],
        #  'SNIa': ['Z=0.0002', 'Z=0.002', 'Z=0.01', 'Z=0.02']}
        Zyield = {key: np.array([float(group.split("=")[1]) for group in groups_Z[key]]) for key in groups.keys()}
        # Zyield is a dictionary like:
        # {'AGB+SNcc': array([0.0001, 0.0003, 0.001 , 0.002 ]),
        #  'SNIa': array([0.0002, 0.002 , 0.01  , 0.02  ])}
        # Sort the metallicity values in ascending order
        self.groups = {key: [groups[key][a] for a in Zyield[key].argsort()] for key in groups.keys()}
        self.groups_Z = {key: [groups_Z[key][a] for a in Zyield[key].argsort()] for key in groups_Z.keys()}
        # The available metallicity values in yield table
        self.Zyield = {key: Zyield[key][Zyield[key].argsort()] for key in Zyield.keys()}
        self.kind = kind
        self.Zindex = { }
        self.Zindex['SNIa'] = { }
        self.Zindex['AGB+SNcc'] = { }
        self.interps = { }

        if len(self.Zyield['SNIa']) == 1:
            self.Zindex['SNIa']['low'] = 0
            self.Zindex['SNIa']['high'] = 0
            self.interps['SNIa'] = lambda ZGas, elem: dfs['SNIa'][self.groups['SNIa'][0]].loc[elem].to_numpy().astype(np.float64)
        else:
            if kind == 'nearest':
                self.interps['SNIa'] = lambda ZGas, elem: self._interp_Z_SNIa_nearest(ZGas, elem)
            elif kind == 'linear-linear':
                self.interps['SNIa'] = lambda ZGas, elem: self._interp_Z_SNIa_linear_linear(ZGas, elem)
            elif kind == 'linear-log':
                self.interps['SNIa'] = lambda ZGas, elem: self._interp_Z_SNIa_linear_log(ZGas, elem)
            elif kind == 'log-linear':
                self.interps['SNIa'] = lambda ZGas, elem: self._interp_Z_SNIa_log_linear(ZGas, elem)
            elif kind == 'log-log':
                self.interps['SNIa'] = lambda ZGas, elem: self._interp_Z_SNIa_log_log(ZGas, elem)
            else:
                print("The kind should be in ['linear-linear', 'linear-log', 'log-linear', 'log-log', 'nearest]")
                print("Return None!")
                self.interps['SNIa'] = None
        if len(self.Zyield['AGB+SNcc']) == 1:
            self.Zindex['AGB+SNcc']['low'] = 0
            self.Zindex['AGB+SNcc']['high'] = 0
            self.interps['AGB+SNcc'] = lambda ZGas, elem: self.interpolate_m(0, elem)
            print("There is only one metallicity in the yield table of AGB+SNcc!")
            print("Please check your yield tables!")
        else:
            if kind == 'nearest':
                self.interps['AGB+SNcc'] = lambda ZGas, elem: self._interp_Z_AGB_SNcc_nearest(ZGas, elem)
            elif kind == 'linear-linear':
                self.interps['AGB+SNcc'] = lambda ZGas, elem: self._interp_Z_AGB_SNcc_linear_linear(ZGas, elem)
            elif kind == 'linear-log':
                self.interps['AGB+SNcc'] = lambda ZGas, elem: self._interp_Z_AGB_SNcc_linear_log(ZGas, elem)
            elif kind == 'log-linear':
                self.interps['AGB+SNcc'] = lambda ZGas, elem: self._interp_Z_AGB_SNcc_log_linear(ZGas, elem)
            elif kind == 'log-log':
                self.interps['AGB+SNcc'] = lambda ZGas, elem: self._interp_Z_AGB_SNcc_log_log(ZGas, elem)
            else:
                print("The kind should be in ['linear-linear', 'linear-log', 'log-linear', 'log-log', 'nearest]")
                print("Return None!")
                self.interps['AGB+SNcc'] = None


    def interpolate_m(self, Zi, elem):
        """
        Interpolate the yield of the element with respect to the mass of the star.
        It means that one yield table is selected and the interpolation is done with respect to the mass.

        Parameters
        ----------
        Zi : int
            The index of the metallicity value in the yield table.
        elem : str
            The element name.
        """
        dfs = self.dfs
        groups = self.groups
        m = dfs['AGB+SNcc'][groups['AGB+SNcc'][Zi]].columns.to_numpy().astype(float)
        y = dfs['AGB+SNcc'][groups['AGB+SNcc'][Zi]].loc[elem].to_numpy()
        return interpolate.interp1d(m, y, kind='linear', fill_value='extrapolate')


####### How interpolate the yields of SNIa with respect to Z     ######


    def _interp_Z_SNIa_linear_linear(self, ZGas, elem):
        """
        Interpolate the yield of the element from SNIa with respect to the metallicity.
        The interpolation is done with the linear-linear method.

        Parameters
        ----------
        elem : str
            The element name.
        ZGas : float
            The metallicity of the gas.
        """
        dfs = self.dfs
        Zindex = {}
        Zindex['SNIa'] = { }
        Zyield = self.Zyield
        groups = self.groups
        if ZGas <= Zyield['SNIa'].min():
            Zindex['SNIa']['low'] = Zyield['SNIa'].argmin()
            Zindex['SNIa']['high'] = Zindex['SNIa']['low']
        elif ZGas >= Zyield['SNIa'].max():
            Zindex['SNIa']['high'] = Zyield['SNIa'].argmax()
            Zindex['SNIa']['low'] = Zindex['SNIa']['high']
        else:
            Zindex['SNIa']['high'] = np.where(Zyield['SNIa']==Zyield['SNIa'][Zyield['SNIa']>ZGas].min(),)[0][0]
            Zindex['SNIa']['low']  = np.where(Zyield['SNIa']==Zyield['SNIa'][Zyield['SNIa']<ZGas].max(),)[0][0]
        Z_low = Zyield['SNIa'][Zindex['SNIa']['low']]
        Z_high = Zyield['SNIa'][Zindex['SNIa']['high']]
        # In case of Z_low is zero, we set Z_low = 1e-10.
        Z_low = np.maximum(Z_low, 1e-10)
        # In case of Z_high is zero, we set Z_high = 1e-10.
        Z_high = np.maximum(Z_high, 1e-10)
        self.Zindex['SNIa'] = Zindex['SNIa']
        if Z_low == Z_high:
            return dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['low']]].loc[elem, groups['SNIa'][Zindex['SNIa']['low']]]
        else:
            Y_low = dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['low']]].loc[elem, groups['SNIa'][Zindex['SNIa']['low']]]
            Y_high = dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['high']]].loc[elem, groups['SNIa'][Zindex['SNIa']['high']]]
            Y = Y_low + (Y_high-Y_low)/(Z_high-Z_low)*(ZGas-Z_low)
            return Y
        # SNIa, linear-linear interpolation

    
    def _interp_Z_SNIa_linear_log(self, ZGas, elem):
        """
        Interpolate the yield of the element from SNIa with respect to the metallicity.
        The interpolation is done with the linear-log method.

        Parameters
        ----------
        elem : str
            The element name.
        ZGas : float
            The metallicity of the gas.
        """
        dfs = self.dfs
        Zindex = {}
        Zindex['SNIa'] = { }
        Zyield = self.Zyield
        groups = self.groups
        if ZGas <= Zyield['SNIa'].min():
            Zindex['SNIa']['low'] = Zyield['SNIa'].argmin()
            Zindex['SNIa']['high'] = Zindex['SNIa']['low']
        elif ZGas >= Zyield['SNIa'].max():
            Zindex['SNIa']['high'] = Zyield['SNIa'].argmax()
            Zindex['SNIa']['low'] = Zindex['SNIa']['high']
        else:
            Zindex['SNIa']['high'] = np.where(Zyield['SNIa']==Zyield['SNIa'][Zyield['SNIa']>ZGas].min(),)[0][0]
            Zindex['SNIa']['low']  = np.where(Zyield['SNIa']==Zyield['SNIa'][Zyield['SNIa']<ZGas].max(),)[0][0]
        Z_low = Zyield['SNIa'][Zindex['SNIa']['low']]
        Z_high = Zyield['SNIa'][Zindex['SNIa']['high']]
        # In case of Z_low is zero, we set Z_low = 1e-10.
        Z_low = np.maximum(Z_low, 1e-10)
        # In case of Z_high is zero, we set Z_high = 1e-10.
        Z_high = np.maximum(Z_high, 1e-10)
        self.Zindex['SNIa'] = Zindex['SNIa']
        if Z_low == Z_high:
            return dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['low']]].loc[elem, groups['SNIa'][Zindex['SNIa']['low']]]
        else:
            Y_low = dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['low']]].loc[elem, groups['SNIa'][Zindex['SNIa']['low']]]
            Y_high = dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['high']]].loc[elem, groups['SNIa'][Zindex['SNIa']['high']]]
            Y = Y_low + (Y_high-Y_low)/(np.log10(Z_high)-np.log10(Z_low))*(np.log10(ZGas)-np.log10(Z_low))
            return Y
        # SNIa, linear-log interpolation
    

    def _interp_Z_SNIa_log_linear(self, ZGas, elem):
        """
        Interpolate the yield of the element from SNIa with respect to the metallicity.
        The interpolation is done with the log-linear method.

        Parameters
        ----------
        elem : str
            The element name.
        ZGas : float
            The metallicity of the gas.
        """
        dfs = self.dfs
        Zindex = {}
        Zindex['SNIa'] = { }
        Zyield = self.Zyield
        groups = self.groups
        if ZGas <= Zyield['SNIa'].min():
            Zindex['SNIa']['low'] = Zyield['SNIa'].argmin()
            Zindex['SNIa']['high'] = Zindex['SNIa']['low']
        elif ZGas >= Zyield['SNIa'].max():
            Zindex['SNIa']['high'] = Zyield['SNIa'].argmax()
            Zindex['SNIa']['low'] = Zindex['SNIa']['high']
        else:
            Zindex['SNIa']['high'] = np.where(Zyield['SNIa']==Zyield['SNIa'][Zyield['SNIa']>ZGas].min(),)[0][0]
            Zindex['SNIa']['low']  = np.where(Zyield['SNIa']==Zyield['SNIa'][Zyield['SNIa']<ZGas].max(),)[0][0]
        Z_low = Zyield['SNIa'][Zindex['SNIa']['low']]
        Z_high = Zyield['SNIa'][Zindex['SNIa']['high']]
        # In case of Z_low is zero, we set Z_low = 1e-10.
        Z_low = np.maximum(Z_low, 1e-10)
        # In case of Z_high is zero, we set Z_high = 1e-10.
        Z_high = np.maximum(Z_high, 1e-10)
        self.Zindex['SNIa'] = Zindex['SNIa']
        if Z_low == Z_high:
            return dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['low']]].loc[elem, groups['SNIa'][Zindex['SNIa']['low']]]
        else:
            Y_low = dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['low']]].loc[elem, groups['SNIa'][Zindex['SNIa']['low']]]
            Y_high = dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['high']]].loc[elem, groups['SNIa'][Zindex['SNIa']['high']]]
            Y = Y_low ** ((Z_high-ZGas)/(Z_high-Z_low)) * Y_high ** ((ZGas-Z_low)/(Z_high-Z_low))
            return Y
        # SNIa, log-linear interpolation


    def _interp_Z_SNIa_log_log(self, ZGas, elem):
        """
        Interpolate the yield of the element from SNIa with respect to the metallicity.
        The interpolation is done with the log-log method.

        Parameters
        ----------
        elem : str
            The element name.
        ZGas : float
            The metallicity of the gas.
        """
        dfs = self.dfs
        Zindex = {}
        Zindex['SNIa'] = { }
        Zyield = self.Zyield
        groups = self.groups
        if ZGas <= Zyield['SNIa'].min():
            Zindex['SNIa']['low'] = Zyield['SNIa'].argmin()
            Zindex['SNIa']['high'] = Zindex['SNIa']['low']
        elif ZGas >= Zyield['SNIa'].max():
            Zindex['SNIa']['high'] = Zyield['SNIa'].argmax()
            Zindex['SNIa']['low'] = Zindex['SNIa']['high']
        else:
            Zindex['SNIa']['high'] = np.where(Zyield['SNIa']==Zyield['SNIa'][Zyield['SNIa']>ZGas].min(),)[0][0]
            Zindex['SNIa']['low']  = np.where(Zyield['SNIa']==Zyield['SNIa'][Zyield['SNIa']<ZGas].max(),)[0][0]
        Z_low = Zyield['SNIa'][Zindex['SNIa']['low']]
        Z_high = Zyield['SNIa'][Zindex['SNIa']['high']]
        # In case of Z_low is zero, we set Z_low = 1e-10.
        Z_low = np.maximum(Z_low, 1e-10)
        # In case of Z_high is zero, we set Z_high = 1e-10.
        Z_high = np.maximum(Z_high, 1e-10)
        self.Zindex['SNIa'] = Zindex['SNIa']
        if Z_low == Z_high:
            return dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['low']]].loc[elem, groups['SNIa'][Zindex['SNIa']['low']]]
        else:
            Y_low = dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['low']]].loc[elem, groups['SNIa'][Zindex['SNIa']['low']]]
            Y_high = dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['high']]].loc[elem, groups['SNIa'][Zindex['SNIa']['high']]]
            Y = Y_low ** ((np.log10(Z_high)-np.log10(ZGas))/(np.log10(Z_high)-np.log10(Z_low))) * \
                Y_high ** ((np.log10(ZGas)-np.log10(Z_low))/(np.log10(Z_high)-np.log10(Z_low)))
            return Y
        # SNIa, log-log interpolation


    def _interp_Z_SNIa_nearest(self, ZGas, elem):
        """
        Interpolate the yield of the element from SNIa with respect to the metallicity.
        The interpolation is done with the nearest method.

        Parameters
        ----------
        elem : str
            The element name.
        ZGas : float
            The metallicity of the gas.
        """
        dfs = self.dfs
        Zindex = {}
        Zindex['SNIa'] = { }
        Zyield = self.Zyield
        groups = self.groups
        Zindex['SNIa']['low'] = np.abs(Zyield['SNIa']-ZGas).argmin()
        Zindex['SNIa']['high'] = Zindex['SNIa']['low']
        Y = dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['low']]].loc[elem, groups['SNIa'][Zindex['SNIa']['low']]]
        self.Zindex['SNIa'] = Zindex['SNIa']
        return Y
        # SNIa, nearest interpolation


    def interpolate_Z_SNIa(self, ZGas, elem):
        """
        Interpolate the yield of the element from SNIa with respect to the metallicity.

        All the interpolation methods are implemented in this function.
        However, I still keep the interpolation methods in separate functions to avoid the repeated "if kind == ..." check,
        which may slow down the calculation.

        Parameters
        ----------
        elem : str
            The element name.
        ZGas : float
            The metallicity of the gas.
        """
        kind = self.kind
        dfs = self.dfs
        Zyield = self.Zyield
        Zindex = {}
        Zindex['SNIa'] = { }
        Zyield = self.Zyield
        groups = self.groups
        if len(Zyield['SNIa']) == 1:
            Zindex['SNIa']['low'] = 0
            Zindex['SNIa']['high'] = 0
            return dfs['SNIa'][self.groups['SNIa'][0]].loc[elem].to_numpy()

        if kind == 'TNG-like':
            # TNG uses the linear-log interpolation
            kind = 'linear-log'
        if kind == 'nearest':
            Zindex['SNIa']['low'] = np.abs(Zyield['SNIa']-ZGas).argmin()
            Zindex['SNIa']['high'] = Zindex['SNIa']['low']
            return dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['low']]].loc[elem, groups['SNIa'][Zindex['SNIa']['low']]]
        elif kind in ['linear-linear', 'linear-log', 'log-linear', 'log-log']:
            if ZGas <= Zyield['SNIa'].min():
                Zindex['SNIa']['low'] = Zyield['SNIa'].argmin()
                Zindex['SNIa']['high'] = Zindex['SNIa']['low']
            elif ZGas >= Zyield['SNIa'].max():
                Zindex['SNIa']['high'] = Zyield['SNIa'].argmax()
                Zindex['SNIa']['low'] = Zindex['SNIa']['high']
            else:
                Zindex['SNIa']['high'] = np.where(Zyield['SNIa']==Zyield['SNIa'][Zyield['SNIa']>ZGas].min(),)[0][0]
                Zindex['SNIa']['low']  = np.where(Zyield['SNIa']==Zyield['SNIa'][Zyield['SNIa']<ZGas].max(),)[0][0]
        else:
            print("The kind should be in ['linear-linear', 'linear-log', 'log-linear', 'log-log', 'nearest]")
            print("Return None!")
            return None
        Z_low = Zyield['SNIa'][Zindex['SNIa']['low']]
        Z_high = Zyield['SNIa'][Zindex['SNIa']['high']]
        # In case of Z_low is zero, we set Z_low = 1e-10.
        Z_low = np.maximum(Z_low, 1e-10)
        # In case of Z_high is zero, we set Z_high = 1e-10.
        Z_high = np.maximum(Z_high, 1e-10)
        self.Zindex['SNIa'] = Zindex['SNIa']
        if Z_low == Z_high:
            return dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['low']]].loc[elem, groups['SNIa'][Zindex['SNIa']['low']]]
        else:
            Y_low = dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['low']]].loc[elem, groups['SNIa'][Zindex['SNIa']['low']]]
            Y_high = dfs['SNIa'][groups['SNIa'][Zindex['SNIa']['high']]].loc[elem, groups['SNIa'][Zindex['SNIa']['high']]]
            if kind == 'linear-linear':
                Y = Y_low + (Y_high-Y_low)/(Z_high-Z_low)*(ZGas-Z_low)
                return Y
            elif kind == 'linear-log':
                Y = Y_low + (Y_high-Y_low)/(np.log10(Z_high)-np.log10(Z_low))*(np.log10(ZGas)-np.log10(Z_low))
                return Y
            elif kind == 'log-linear':
                Y = Y_low ** ((Z_high-ZGas)/(Z_high-Z_low)) * Y_high ** ((ZGas-Z_low)/(Z_high-Z_low))
                return Y
            else:
                Y = Y_low ** ((np.log10(Z_high)-np.log10(ZGas))/(np.log10(Z_high)-np.log10(Z_low))) * \
                    Y_high ** ((np.log10(ZGas)-np.log10(Z_low))/(np.log10(Z_high)-np.log10(Z_low)))
                return Y
            # SNIa, all the interpolation methods


####### How interpolate the yields of SNIa with respect to Z     ######

####### How interpolate the yields of AGB+SNcc with respect to Z     ######


    def _interp_Z_AGB_SNcc_linear_linear(self, ZGas, elem):
        """
        Interpolate the yield of the element from AGB stars and SNcc with respect to the metallicity.
        The interpolation is done with the linear-linear method.

        Parameters
        ----------
        elem : str
            The element name.
        ZGas : float
            The metallicity of the gas.
        """
        Zindex = {}
        Zindex['AGB+SNcc'] = { }
        Zyield = self.Zyield
        if ZGas <= Zyield['AGB+SNcc'].min():
            Zindex['AGB+SNcc']['low'] = Zyield['AGB+SNcc'].argmin()
            Zindex['AGB+SNcc']['high'] = Zindex['AGB+SNcc']['low']
        elif ZGas >= Zyield['AGB+SNcc'].max():
            Zindex['AGB+SNcc']['high'] = Zyield['AGB+SNcc'].argmax()
            Zindex['AGB+SNcc']['low'] = Zindex['AGB+SNcc']['high']
        else:
            Zindex['AGB+SNcc']['high'] = np.where(Zyield['AGB+SNcc']==Zyield['AGB+SNcc'][Zyield['AGB+SNcc']>ZGas].min(),)[0][0]
            Zindex['AGB+SNcc']['low']  = np.where(Zyield['AGB+SNcc']==Zyield['AGB+SNcc'][Zyield['AGB+SNcc']<ZGas].max(),)[0][0]
        Z_low = Zyield['AGB+SNcc'][Zindex['AGB+SNcc']['low']]
        Z_high = Zyield['AGB+SNcc'][Zindex['AGB+SNcc']['high']]
        # In case of Z_low is zero, we set Z_low = 1e-10.
        Z_low = np.maximum(Z_low, 1e-10)
        # In case of Z_high is zero, we set Z_high = 1e-10.
        Z_high = np.maximum(Z_high, 1e-10)
        if Z_low == Z_high:
            interp = self.interpolate_m(Zindex['AGB+SNcc']['low'], elem)
        else:
            interp_low = self.interpolate_m(Zindex['AGB+SNcc']['low'], elem)
            interp_high = self.interpolate_m(Zindex['AGB+SNcc']['high'], elem)
            interp = lambda x: interp_low(x) + (interp_high(x)-interp_low(x))/(Z_high-Z_low)*(ZGas-Z_low)
        self.Zindex['AGB+SNcc'] = Zindex['AGB+SNcc']
        return interp
        # AGB+SNcc, linear-linear interpolation
    
    def _interp_Z_AGB_SNcc_linear_log(self, ZGas, elem):
        """
        Interpolate the yield of the element from AGB stars and SNcc with respect to the metallicity.
        The interpolation is done with the linear-log method.

        Parameters
        ----------
        elem : str
            The element name.
        ZGas : float
            The metallicity of the gas.
        """
        Zindex = {}
        Zindex['AGB+SNcc'] = { }
        Zyield = self.Zyield
        if ZGas <= Zyield['AGB+SNcc'].min():
            Zindex['AGB+SNcc']['low'] = Zyield['AGB+SNcc'].argmin()
            Zindex['AGB+SNcc']['high'] = Zindex['AGB+SNcc']['low']
        elif ZGas >= Zyield['AGB+SNcc'].max():
            Zindex['AGB+SNcc']['high'] = Zyield['AGB+SNcc'].argmax()
            Zindex['AGB+SNcc']['low'] = Zindex['AGB+SNcc']['high']
        else:
            Zindex['AGB+SNcc']['high'] = np.where(Zyield['AGB+SNcc']==Zyield['AGB+SNcc'][Zyield['AGB+SNcc']>ZGas].min(),)[0][0]
            Zindex['AGB+SNcc']['low']  = np.where(Zyield['AGB+SNcc']==Zyield['AGB+SNcc'][Zyield['AGB+SNcc']<ZGas].max(),)[0][0]
        Z_low = Zyield['AGB+SNcc'][Zindex['AGB+SNcc']['low']]
        Z_high = Zyield['AGB+SNcc'][Zindex['AGB+SNcc']['high']]
        # In case of Z_low is zero, we set Z_low = 1e-10.
        Z_low = np.maximum(Z_low, 1e-10)
        # In case of Z_high is zero, we set Z_high = 1e-10.
        Z_high = np.maximum(Z_high, 1e-10)
        if Z_low == Z_high:
            interp = self.interpolate_m(Zindex['AGB+SNcc']['low'], elem)
        else:
            interp_low = self.interpolate_m(Zindex['AGB+SNcc']['low'], elem)
            interp_high = self.interpolate_m(Zindex['AGB+SNcc']['high'], elem)
            interp = lambda x: interp_low(x) + (interp_high(x)-interp_low(x))/\
                    (np.log10(Z_high)-np.log10(Z_low))*(np.log10(ZGas)-np.log10(Z_low))
        self.Zindex['AGB+SNcc'] = Zindex['AGB+SNcc']
        return interp
        # AGB+SNcc, linear-log interpolation


    def _interp_Z_AGB_SNcc_log_linear(self, ZGas, elem):
        """
        Interpolate the yield of the element from AGB stars and SNcc with respect to the metallicity.
        The interpolation is done with the log-linear method.

        Parameters
        ----------
        elem : str
            The element name.
        ZGas : float
            The metallicity of the gas.
        """
        Zindex = {}
        Zindex['AGB+SNcc'] = { }
        Zyield = self.Zyield
        if ZGas <= Zyield['AGB+SNcc'].min():
            Zindex['AGB+SNcc']['low'] = Zyield['AGB+SNcc'].argmin()
            Zindex['AGB+SNcc']['high'] = Zindex['AGB+SNcc']['low']
        elif ZGas >= Zyield['AGB+SNcc'].max():
            Zindex['AGB+SNcc']['high'] = Zyield['AGB+SNcc'].argmax()
            Zindex['AGB+SNcc']['low'] = Zindex['AGB+SNcc']['high']
        else:
            Zindex['AGB+SNcc']['high'] = np.where(Zyield['AGB+SNcc']==Zyield['AGB+SNcc'][Zyield['AGB+SNcc']>ZGas].min(),)[0][0]
            Zindex['AGB+SNcc']['low']  = np.where(Zyield['AGB+SNcc']==Zyield['AGB+SNcc'][Zyield['AGB+SNcc']<ZGas].max(),)[0][0]
        Z_low = Zyield['AGB+SNcc'][Zindex['AGB+SNcc']['low']]
        Z_high = Zyield['AGB+SNcc'][Zindex['AGB+SNcc']['high']]
        # In case of Z_low is zero, we set Z_low = 1e-10.
        Z_low = np.maximum(Z_low, 1e-10)
        # In case of Z_high is zero, we set Z_high = 1e-10.
        Z_high = np.maximum(Z_high, 1e-10)
        if Z_low == Z_high:
            interp = self.interpolate_m(Zindex['AGB+SNcc']['low'], elem)
        else:
            interp_low = self.interpolate_m(Zindex['AGB+SNcc']['low'], elem)
            interp_high = self.interpolate_m(Zindex['AGB+SNcc']['high'], elem)
            interp = lambda x: interp_low(x) ** ((Z_high-ZGas)/(Z_high-Z_low)) * \
                        interp_high(x) ** ((ZGas-Z_low)/(Z_high-Z_low))
        self.Zindex['AGB+SNcc'] = Zindex['AGB+SNcc']
        return interp
        # AGB+SNcc, log-linear interpolation


    def _interp_Z_AGB_SNcc_log_log(self, ZGas, elem):
        """
        Interpolate the yield of the element from AGB stars and SNcc with respect to the metallicity.
        The interpolation is done with the log-log method.

        Parameters
        ----------
        elem : str
            The element name.
        ZGas : float
            The metallicity of the gas.
        """
        Zindex = {}
        Zindex['AGB+SNcc'] = { }
        Zyield = self.Zyield
        if ZGas <= Zyield['AGB+SNcc'].min():
            Zindex['AGB+SNcc']['low'] = Zyield['AGB+SNcc'].argmin()
            Zindex['AGB+SNcc']['high'] = Zindex['AGB+SNcc']['low']
        elif ZGas >= Zyield['AGB+SNcc'].max():
            Zindex['AGB+SNcc']['high'] = Zyield['AGB+SNcc'].argmax()
            Zindex['AGB+SNcc']['low'] = Zindex['AGB+SNcc']['high']
        else:
            Zindex['AGB+SNcc']['high'] = np.where(Zyield['AGB+SNcc']==Zyield['AGB+SNcc'][Zyield['AGB+SNcc']>ZGas].min(),)[0][0]
            Zindex['AGB+SNcc']['low']  = np.where(Zyield['AGB+SNcc']==Zyield['AGB+SNcc'][Zyield['AGB+SNcc']<ZGas].max(),)[0][0]
        Z_low = Zyield['AGB+SNcc'][Zindex['AGB+SNcc']['low']]
        Z_high = Zyield['AGB+SNcc'][Zindex['AGB+SNcc']['high']]
        # In case of Z_low is zero, we set Z_low = 1e-10.
        Z_low = np.maximum(Z_low, 1e-10)
        # In case of Z_high is zero, we set Z_high = 1e-10.
        Z_high = np.maximum(Z_high, 1e-10)
        if Z_low == Z_high:
            interp = self.interpolate_m(Zindex['AGB+SNcc']['low'], elem)
        else:
            interp_low = self.interpolate_m(Zindex['AGB+SNcc']['low'], elem)
            interp_high = self.interpolate_m(Zindex['AGB+SNcc']['high'], elem)
            interp = lambda x: interp_low(x) ** ((np.log10(Z_high)-np.log10(ZGas))/(np.log10(Z_high)-np.log10(Z_low))) * \
                        interp_high(x) ** ((np.log10(ZGas)-np.log10(Z_low))/(np.log10(Z_high)-np.log10(Z_low)))
        self.Zindex['AGB+SNcc'] = Zindex['AGB+SNcc']
        return interp
        # AGB+SNcc, log-log interpolation

    
    def _interp_Z_AGB_SNcc_nearest(self, ZGas, elem):
        """
        Interpolate the yield of the element from AGB stars and SNcc with respect to the metallicity.
        The interpolation is done with the nearest method.

        Parameters
        ----------
        elem : str
            The element name.
        ZGas : float
            The metallicity of the gas.
        """
        Zindex = {}
        Zindex['AGB+SNcc'] = {}
        Zyield = self.Zyield
        Zindex['AGB+SNcc']['low'] = np.abs(Zyield['AGB+SNcc']-ZGas).argmin()
        Zindex['AGB+SNcc']['high'] = Zindex['AGB+SNcc']['low']
        interp = self.interpolate_m(Zindex['AGB+SNcc']['low'], elem)
        self.Zindex['AGB+SNcc'] = Zindex['AGB+SNcc']
        return interp
        # AGB+SNcc, nearest interpolation


    def interpolate_Z_AGB_SNcc(self, ZGas, elem):
        """
        Interpolate the yield of the element from AGB stars and SNcc with respect to the metallicity.
        All the interpolation methods are implemented in this function.
        However, I still keep the interpolation methods in separate functions to avoid the repeated "if kind == ..." check,
        which may slow down the calculation.

        Parameters
        ----------
        elem : str
            The element name.
        ZGas : float
            The metallicity of the gas.
        """
        kind = self.kind
        Zindex = {}
        Zindex['AGB+SNcc'] = {}
        Zyield = self.Zyield
        if len(self.Zyield['AGB+SNcc']) == 1:
            interp = lambda ZGas, elem: self.interpolate_m(0, elem)
            Zindex['SNIa']['low'] = 0
            Zindex['SNIa']['high'] = 0
            return interp

        if kind == 'TNG-like':
            # TNG uses the linear-log interpolation
            kind = 'linear-log'
        if kind == 'nearest':
            Zindex['AGB+SNcc']['low'] = np.abs(Zyield['AGB+SNcc']-ZGas).argmin()
            Zindex['AGB+SNcc']['high'] = Zindex['AGB+SNcc']['low']
            interp = self.interpolate_m(Zindex['AGB+SNcc']['low'], elem)
            self.Zindex['AGB+SNcc'] = Zindex['AGB+SNcc']
            return interp
        elif kind in ['linear-linear', 'linear-log', 'log-linear', 'log-log']:
            if ZGas <= Zyield['AGB+SNcc'].min():
                Zindex['AGB+SNcc']['low'] = Zyield['AGB+SNcc'].argmin()
                Zindex['AGB+SNcc']['high'] = Zindex['AGB+SNcc']['low']
            elif ZGas >= Zyield['AGB+SNcc'].max():
                Zindex['AGB+SNcc']['high'] = Zyield['AGB+SNcc'].argmax()
                Zindex['AGB+SNcc']['low'] = Zindex['AGB+SNcc']['high']
            else:
                Zindex['AGB+SNcc']['high'] = np.where(Zyield['AGB+SNcc']==Zyield['AGB+SNcc'][Zyield['AGB+SNcc']>ZGas].min(),)[0][0]
                Zindex['AGB+SNcc']['low']  = np.where(Zyield['AGB+SNcc']==Zyield['AGB+SNcc'][Zyield['AGB+SNcc']<ZGas].max(),)[0][0]
        else:
            print("The kind should be in ['linear-linear', 'linear-log', 'log-linear', 'log-log', 'nearest]")
            print("Return None!")
            return None
        Z_low = Zyield['AGB+SNcc'][Zindex['AGB+SNcc']['low']]
        Z_high = Zyield['AGB+SNcc'][Zindex['AGB+SNcc']['high']]
        # In case of Z_low is zero, we set Z_low = 1e-10.
        Z_low = np.maximum(Z_low, 1e-10)
        # In case of Z_high is zero, we set Z_high = 1e-10.
        Z_high = np.maximum(Z_high, 1e-10)
        if Z_low == Z_high:
            interp = self.interpolate_m(Zindex['AGB+SNcc']['low'], elem)
        else:
            interp_low = self.interpolate_m(Zindex['AGB+SNcc']['low'], elem)
            interp_high = self.interpolate_m(Zindex['AGB+SNcc']['high'], elem)
            if kind == 'linear-linear':
                interp = lambda x: interp_low(x) + (interp_high(x)-interp_low(x))/(Z_high-Z_low)*(ZGas-Z_low)
            elif kind == 'linear-log':
                interp = lambda x: interp_low(x) + (interp_high(x)-interp_low(x))/\
                        (np.log10(Z_high)-np.log10(Z_low))*(np.log10(ZGas)-np.log10(Z_low))
            elif kind == 'log-linear':
                interp = lambda x: interp_low(x) ** ((Z_high-ZGas)/(Z_high-Z_low)) * \
                        interp_high(x) ** ((ZGas-Z_low)/(Z_high-Z_low))
            else:
                interp = lambda x: interp_low(x) ** ((np.log10(Z_high)-np.log10(ZGas))/(np.log10(Z_high)-np.log10(Z_low))) * \
                        interp_high(x) ** ((np.log10(ZGas)-np.log10(Z_low))/(np.log10(Z_high)-np.log10(Z_low)))
        self.Zindex['AGB+SNcc'] = Zindex['AGB+SNcc']
        return interp


####### How interpolate the yields of AGB+SNcc with respect to Z     ######