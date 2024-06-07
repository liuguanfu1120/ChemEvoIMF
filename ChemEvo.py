# -*- coding: utf-8 -*-
"""
Created on May 27, 2024
Last modified on May 27, 2024
@Author: Guan-Fu Liu

To calculate how the abundance evolves with time with variable IMF
To reproduce the Figure 10 in Nomoto et al. (2013).

Update: 
1. Separate utils.py into IMF.py, primordial_gas.py, SupernovaeIa.py, and MassLifetime.py
2. Interpolate the stellar yields of different initial metallicities.
3. Add stellar mass
4. Add SNcc events
5. Add AGB events
6. Add the abundances of the ejecta from SNcc
7. Add the abundances of the ejecta from AGB
8. Add the abundances of the ejecta from SNIa
"""
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, './')
from scipy.integrate import quad
import h5py
from tqdm import tqdm
import datetime
import os
import inspect
# The following are the modules defined in the same directory
from MassLifetime import MassLifetime
from IMF import IMF
import SupernovaeIa
import primordial_gas as pg
import constants
from InterpolateYields import InterpolateYields

def ChemEvo(SFH, SFE, yield_files, imf_evolve=None, imf_dict=None, SNIaOn=True, p_preset=None, mass_lifetime_file=None,
            interp_kind="linear-log", solar_set='Default', Z_0=0, input_primodiral_gas=None,
            ElemNotice=["H", "He", "C", "N", "O", "Ne", "Si", "Mg", "Fe", "Other"], 
            output_dir="./outputs", out_file=None, comments=None):
    """
    Calculate how the abundance evolves with time with variable IMF

    
    Parameters
    ----------
    SFH: dict
        The star formation history.
        It should have the following keys:
        'Age':
            The age of the Universe in yr.
            Note that the unit should be yr!
        'SFR':
            The star formation rate in Msun/yr.
            Note that the unit should be Msun/yr!
            The star formation rate should be non-negative. 
    SFE: float
        The star formation efficiency.
        It should be a float.
        The total mass of gas is determined by the total mass of formed stars / SFE.
        The typical value is 0.1, and 0.3 is an enhanced value.
        It is unlikely that SFE > 0.5, and hence the code will EXIT if SFE > 0.5.
    yield_files: dict
        The dictionary that contains the yield file.
        It should have the following keys:
        'AGB+SNcc': str
            The yield file of AGB and SNcc.
        'SNIa': str
            The yield file of SNIa.
        For example, it could be {"AGB+SNcc": "./inputs/NuPyCEE/agb_and_massive_stars_C15_N13_0_5_HNe/yields1.h5",
                      "SNIa": "./inputs/NuPyCEE/sn1a_i99_W7/yields1.h5"}
    imf_evolve: function
        The default is None.
        It describes the IMF function evolves with ZGas.
        For example, it could be
        def imf_evolve(Z_gas, age):
            if Z_gas < 0.001:
                return lambda m: m**(-1.6) if (m>= constants.Mstar_min and m<= constants.Mstar_max) else 0
            else:
                return lambda m: m**(-2.3) if (m>= constants.Mstar_min and m<= constants.Mstar_max) else 0
        If it is not None, the imf_dict will be ignored.
    imf_dict: dict
        The IMF dictionary.
        It only takes effect when imf_evolve is None.
        It is {'IMF_type': parameter}, it should be one of the following:
        {'Salpeter': None}:
            The Salpeter IMF. No parameter is needed, and set the value to None.
        {'Kroupa': None}:
            The Kroupa IMF. No parameter is needed, and set the value to None.
        {'Chabrier': None}:
            The Chabrier IMF. No parameter is needed, and set the value to None.
        {'DietSalpeter': None}:
            The diet Salpeter IMF. No parameter is needed, and set the value to None.
        {'Custom': IMF_arrs}
            The custom IMF. IMF_arrs is a list of numpy arrays.
            IMF_arrs[i] is an array of shape (n, 2), where n is the number of mass bins.
            IMF_arrs[i][:, 0] is the mass of the stars in Msun.
            IMF_arrs[i][:, 1] is the sampled IMF, which can be non-normalized.
            The IMF class will normalize it.
            The length of IMF_arrs should be 1 or the SAME as the number of ages with non-zero SFR.
            If the length of IMF_arrs is 1, the code will use the same IMF for all the ages with non-zero SFR.
        {'PowerLaw': alpha}
            The power-law IMF. alpha is the power-law index.
            alpha should be a numpy array.
            If the length of alpha is 1, the code will use the same alpha for all the ages with non-zero SFR.
            If the length of alpha is not 1, the length should be the SAME as the number of ages with non-zero SFR.
            You can input an identical alpha array but it will be slower than inputting a single alpha.

    SNIaOn: bool, optional
        The default is True.
        If True, the code will include the SNIa events and the SNIa ejecta.
        If False, the code will not include the SNIa events and the SNIa ejecta.
    p_preset: float, optional
        The default is None.
        The calibration parameter of the SNIa rate.
        For more details, see the SupernovaeIa.py.
        I do not recommend you to set it unless you know what you are doing.
    mass_lifetime_file: str, optional
        The default is None, which means the code will use the default mass lifetime relation.
        For more details, see MassLifetime.py.
        If not None, the code will use the input hdf5 file that contains the mass-lifetime relation (see MassLifetime.py).
    interp_kind: str, optional
        The default is "linear-log", which the same as that of TNG,
        while the default of galIMF is "linear-linear".
        The kind of interpolation. 
        It should be in ['linear-linear', 'linear-log', 'log-linear', 'log-log', 'nearest', 'TNG-like'].
        The yields of the ith element from AGB and SNcc are a function of the initial metallicity and the mass of the star,
        namely,
        Y_i = Y_i(m, Z).
        The yields of the ith element from SNIa are a function of the initial metallicity,
        namely,
        Y_i = Y_i(Z).
        As is often the case, the SNIa yields are just one metallicity. However,
        there are still some yield tables of SNIa with multiple metallicities.
        Since the yields are tabulated, we need to do interpolation.
        For linear-log interpolation:
                Y = Y_low + (Y_high - Y_low)/(log(Z_high) - log(Z_low))*(log(ZGas) - log(Z_low))
                Y and log Z are linearly interpolated.
                In case of Z_low is zero, we set Z_low = 1e-10.
                This is the same as what is done in TNG simulation.
        linear means Y is linearly interpolated and log means log Z is linearly interpolated.
        For more details, please see the InterpolateYields.py.
        It is also used to interpolate the mass lifetime relation.
        Although the mass lifetime relations have weak dependence on Z.
    solar_set: str, optional
        The default is 'Default'.
        It should be in ['Default', 'AG', 'Allen', 'RA', 'Grevesse', 'GS', 'Lodders', 'solar']
        For more details, see the constants.py.
    Z_0: float, optional
        The default is 0.
        The primordial metallicity but without lithium.
        Together with the solar_set, it is used to add the abundances of the primordial gas except for lithium, 
        which obeys the ratio of the corresponding solar set.
    input_primodiral_gas: numpy array or None
        If it is None, use the default primordial gas from the table II in Cyburt et al. (2016).
        If it is a numpy array, it shuld have a shape of (len(constants.elem_names), ).
        For more details, see the primordial_gas.py.
    ElemNotice: list, optional
        The default is ["H", "He", "C", "N", "O", "Ne", "Si", "Mg", "Fe", "Other"].
        The list of elements that you want to track.
        Please note that
            1. You should include "H", "He", and "Other" in the list.
               Although the code will check it and add them automatically, it is better to include them.
            2. We only consider the first 30 elements in constants.py, if you need to track the elements 
               beyond the first 30 elements, please modify the constants.py.
            3. If some elements are just available in part of the AGB and SNcc yield tables, you cannot
               notice them.
            4. If some elements are just available in the SNIa yield tables, you cannot notice them.
               Although the code will remove them automatically, it is better to exclude them.
            5. The code will add up the elements you do not notice and name it as "Other".
    output_dir: str, optional
        The default is "./outputs".
        The directory to save the output files.
        The code will save the output files in the output_dir.
    out_file: str, optional
        The default is None.
        The name of the output file.
        If None, the code will use the default name, which is the current time ("YYYY-MM-DD-HH-MM.h5")
        If not None, the code will use the input name.
        The output file is in the HDF5 format and in the output_dir.
    comments: str, optional
        The default is None.
        You can add some comments to the output file.
        You may as well add something to distinguish the output file from others.
    """

    ###### Check the SFE ######
    if SFE > 0.5:
        print("The star formation efficiency is too high, please check the input value.")
        print("The typical value is 0.1, the enhanced value is 0.3.")
        print("EXIT!")
        sys.exit()
    ###### Check the SFE ######

    ###### Set the parameters of quad ######
    epsrel=1e-8
    limit=80
    ###### Set the parameters of quad ######

    ###### Initialize the IMF and SNIa class ######
    if imf_evolve is not None:
        # ZGas[0] has not been calculated yet, so we do not initialize the imf here.
        imf_key = "Evolve with ZGas"
        imf_num = 999  # Just set imf_num to a large number (larger than 1).
        if imf_dict is not None:
            print("The imf_dict will be ignored.")
            print("Only the imf_evolve will be used.")
        pass
    else:
        if imf_dict is None:
            print("You must specify one of the imf_evolve and imf_dict.")
            print("The imf_dict should be input because imf_evolve is None!")
            print("EXIT!")
            sys.exit()
        imf_key = list(imf_dict.keys())[0]
        imf_value = imf_dict[imf_key]
        imf_num = 1  # the number of IMF functions, 1 for invariant IMF, >1 for time-varying IMF
        imf_idx = 0  # the index of the IMF function, at first it is 0, and then it will be added by 1
                    # after a star formation epoch (SFR>0).
        if imf_key in ['Salpeter', 'Chabrier', 'Kroupa', 'DietSalpeter']:
            # The invariant IMF
            imf = IMF(IMF_type=imf_key).imf
        elif imf_key == "PowerLaw" and len(imf_value) == 1:
            # The invariant power-law IMF
            imf = IMF(IMF_type="PowerLaw", power_index=imf_value).imf
        elif imf_key == "PowerLaw" and len(imf_value) != 1:
            # The time-varying power-law IMF
            imf_num = len(imf_value)
            imf = IMF(IMF_type="PowerLaw", power_index=imf_value[imf_idx]).imf
            # Initialize the IMF for the time-varying power-law IMF
        elif imf_key == "Custom" and len(imf_value) == 1:
            # The invariant custom IMF
            imf = IMF(IMF_type="Custom", IMF_arr=imf_value).imf
        elif imf_key == "Custom" and len(imf_value) != 1:
            # The time-varying custom IMF
            imf_num = len(imf_value)
            imf = IMF(IMF_type="Custom", IMF_arr=imf_value[imf_idx]).imf
            # Initialize the IMF for the time-varying custom IMF
        else:
            print("Not supported IMF type %s" % imf_key)
            print("EXIT!")
            sys.exit()
        SNIa = SupernovaeIa.SNIa(imf, IMF_type=imf_key, p_preset=p_preset)  # SNIa class is always associated with the imf.
    ###### Initialize the IMF class ######

    ###### Check the SFH and add more keys ######
    if len(SFH['Age']) != len(SFH['SFR']):
        print("The length of the Age and SFR should be the same.")
        print("EXIT!")
        sys.exit()
    if np.diff(SFH['Age'])[0] < 0:
        print("The age should be in ascending order!")
        print("EXIT!")
        sys.exit()
    if SFH['SFR'][-1]>0:
        print("The SFR at the last time is non-zero.")
        print("It will be regarded as zero. Set it to zero.")
        SFH['SFR'][-1] = 0
    if (SFH['SFR'] < 0).any():
        print("The SFR should be non-negative.")
        print("EXIT!")
        sys.exit()
    if imf_key in ['PowerLaw', 'Custom'] and (SFH['SFR']>0).sum( ) > imf_num and imf_num !=1:
        print("The number of IMF functions is less than the number of ages with non-zero SFR.")
        print("EXIT!")
        sys.exit()
    elif imf_key in ['PowerLaw', 'Custom'] and (SFH['SFR']>0).sum( ) < imf_num:
        print("The number of IMF functions is larger than the number of ages with non-zero SFR.")
        print("The code will move on, but the extra IMF functions will not be used.")
    SFH['Mstar'] = np.zeros(len(SFH['Age']))
    SFH['TimeBin'] = np.zeros(len(SFH['Age']))
    SFH['TimeBin'][:-1] = np.diff(SFH['Age'])  # The time bin of each time, the time bin of the last time is set to ZERO.
    SFH['Mstar'][:-1] = SFH['SFR'][:-1] * np.diff(SFH['Age'])  # The mass of stars formed in each time bin
    # We do not consider the mass of stars formed in the last time!
    SFH['Mtot'] = np.sum(SFH['Mstar'])  # The total mass of stars formed.
    ###### Check the SFH and add more keys ######

    ###### Load the yields tables ######
    files = {key: h5py.File(value, 'r') for key, value in yield_files.items()}
    groups = {key: list(files[key].keys()) for key in files.keys()}
    if "Z_" in groups['AGB+SNcc'][0] or "Z_" in groups['SNIa'][0]:
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
    groups = {key: [groups[key][a] for a in Zyield[key].argsort()] for key in groups.keys()}
    groups_Z = {key: [groups_Z[key][a] for a in Zyield[key].argsort()] for key in groups_Z.keys()}
    # Check the ElemNotice
    # It should contain H, He and Other.
    # The finally selected elements are the intersection of ElemNotice and the elements avaliable in the yield table.
    # You may as well not to notice too many elements, which will speed down the calculation.
    if "H" not in ElemNotice:
        print("H should be in ElemNotice!")
        print("Add H to ElemNotice!")
        ElemNotice.append("H")
    if "He" not in ElemNotice:
        print("He should be in ElemNotice!")
        print("Add He to ElemNotice!")
        ElemNotice.append("He")
    if "Other" not in ElemNotice:
        print("Other should be in ElemNotice!")
        print("Add Other to ElemNotice!")
        ElemNotice.append("Other")
    dfs = { }
    ElemIntersect = { }
    if "Z_" not in groups["AGB+SNcc"][0]:
        for key in files.keys():
            dfs[key] = { }
            ElemIntersect[key] = { }
            for group in groups[key]:
                if key == "SNIa":
                    dfs[key][group] = pd.DataFrame(files[key][group]['Original'][...])
                else:
                    dfs[key][group] = pd.DataFrame(files[key][group]['Interpolated'][...])
                dfs[key][group].loc[:, 'M'] = dfs[key][group].loc[:, 'M'].astype(str)
                dfs[key][group].set_index('M', inplace=True)
                index = dfs[key][group].index
                index = [False if a in ElemNotice + ['Mrem'] else True for a in index]
                # Add the mass of the elements not in ElemNotice to the "Other" element
                if "Other" in dfs[key][group].index:
                    dfs[key][group].loc["Other"] += dfs[key][group].loc[index].sum(axis=0)
                else:
                    dfs[key][group].loc["Other"] = dfs[key][group].loc[index].sum(axis=0)
                ElemIntersect[key][group] = list(set(ElemNotice).intersection(set(dfs[key][group].index)))
    else:
        for key in files.keys():
            dfs[key] = { }
            ElemIntersect[key] = { }
            for group in groups[key]:
                if key == "SNIa":
                    dfs[key][group] = pd.read_hdf(yield_files["SNIa"], key="%s/Original"%group)
                else:
                    dfs[key][group] = pd.read_hdf(yield_files["AGB+SNcc"], key="%s/Interpolated"%group)
                index = dfs[key][group].index
                index = [False if a in ElemNotice + ['Mrem'] else True for a in index]
                # Add the mass of the elements not in ElemNotice to the "Other" element
                if "Other" in dfs[key][group].index:
                    dfs[key][group].loc["Other"] += dfs[key][group].loc[index].sum(axis=0)
                else:
                    dfs[key][group].loc["Other"] = dfs[key][group].loc[index].sum(axis=0)
                ElemIntersect[key][group] = list(set(ElemNotice).intersection(set(dfs[key][group].index)))
    for key in files.keys():
        for group in groups[key]:
            if len(set(ElemNotice)-set(ElemIntersect[key][group]))>0:
                print("The elements you notice but not available in the yield table of %s, %s are\n" %\
                    (key, group), set(ElemNotice)-set(ElemIntersect[key][group]))
    # If one element is only available in part of the metallicity values, you should not notice it.
    # For example, if "C" is only available in Z=0.001 of AGB+SNcc table, 
    # you should not notice "C", namely, you should not include "C" in ElemNotice.
    # The yields of "C" should be added to "Other".
    # It is the same for SNIa table.
    # It is OK if one element is just available in AGB+SNcc table but not in SNIa table,
    # which is the case for "H" and "He".
    # The following code is to check if you notice some elements that are just available in part of AGB+SNcc yield table 
    # or SNIa yield table.
    Elems = {key: [value for value in ElemIntersect[key].values()] for key in files.keys()}
    flag = False
    for key in files.keys():
        n = len(list(set.intersection(*map(set, Elems[key]))))
        for group in groups[key]:
            if len(ElemIntersect[key][group]) != n:
                flag = True
                print("Warning: the elements not available in the yield table of %s, %s are" %\
                    (key, group), set(ElemIntersect[key][group])-set.intersection(*map(set, Elems[key])))
                print("Although they are available in the %s yield table of the other metallicities,"%key)
                print("they should be and will be removed from the ElemNotice.")
                for elem in list(set(ElemIntersect[key][group])-set.intersection(*map(set, Elems[key]))):
                    ElemNotice.remove(elem)
    # Re-generate the dfs and ElemIntersect
    if flag and "Z_" not in groups["AGB+SNcc"][0]:
        dfs = { }
        ElemIntersect = { }
        for key in files.keys():
            dfs[key] = { }
            ElemIntersect[key] = { }
            for group in groups[key]:
                if key == "SNIa":
                    dfs[key][group] = pd.DataFrame(files[key][group]['Original'][...])
                else:
                    dfs[key][group] = pd.DataFrame(files[key][group]['Interpolated'][...])
                dfs[key][group].loc[:, 'M'] = dfs[key][group].loc[:, 'M'].astype(str)
                dfs[key][group].set_index('M', inplace=True)
                index = dfs[key][group].index
                index = [False if a in ElemNotice + ['Mrem'] else True for a in index]
                # Add the mass of the elements not in ElemNotice to the "Other" element
                if "Other" in dfs[key][group].index:
                    dfs[key][group].loc["Other"] += dfs[key][group].loc[index].sum(axis=0)
                else:
                    dfs[key][group].loc["Other"] = dfs[key][group].loc[index].sum(axis=0)
                ElemIntersect[key][group] = list(set(ElemNotice).intersection(set(dfs[key][group].index)))
        for key in files.keys():
            for group in groups[key]:
                if len(set(ElemNotice)-set(ElemIntersect[key][group]))>0:
                    print("The elements you notice but not available in the yield table of %s, %s are\n" %\
                        (key, group), set(ElemNotice)-set(ElemIntersect[key][group]))
    if flag and "Z_" in groups["AGB+SNcc"][0]:
        dfs = { }
        ElemIntersect = { }
        for key in files.keys():
            dfs[key] = { }
            ElemIntersect[key] = { }
            for group in groups[key]:
                if key == "SNIa":
                    dfs[key][group] = pd.read_hdf(yield_files["SNIa"], key="%s/Original"%group)
                else:
                    dfs[key][group] = pd.read_hdf(yield_files["AGB+SNcc"], key="%s/Interpolated"%group)
                index = dfs[key][group].index
                index = [False if a in ElemNotice + ['Mrem'] else True for a in index]
                # Add the mass of the elements not in ElemNotice to the "Other" element
                if "Other" in dfs[key][group].index:
                    dfs[key][group].loc["Other"] += dfs[key][group].loc[index].sum(axis=0)
                else:
                    dfs[key][group].loc["Other"] = dfs[key][group].loc[index].sum(axis=0)
                ElemIntersect[key][group] = list(set(ElemNotice).intersection(set(dfs[key][group].index)))
                dfs[key][group] = dfs[key][group].loc[ElemIntersect[key][group]]

        for key in files.keys():
            for group in groups[key]:
                if len(set(ElemNotice)-set(ElemIntersect[key][group]))>0:
                    print("The elements you notice but not available in the yield table of %s, %s are\n" %\
                        (key, group), set(ElemNotice)-set(ElemIntersect[key][group]))
    for key in dfs.keys():
        for group in dfs[key].keys():
            if (dfs[key][group]<0).any().any():
                # Create a boolean mask for negative values
                mask = dfs[key][group]<0
                # Find indices and columns where the elements are negative
                negative_positions = [(index, col) for col in dfs[key][group].columns \
                                      for index in dfs[key][group].index if mask.at[index, col]]
                print("Please check the yield file from %s, %s!"%(key, group))
                print("Indices and Columns of Negative Yields are")
                print(negative_positions)
                if interp_kind in ["log-linear", "log-log"]:
                    print("Since there are some negative yields,")
                    print("you cannot conduct interpolation in the logarithmic space of the yields.")
                    print("EXIT!")
                    sys.exit()

    ###### Load the yields tables ######

    ###### Load the mass-lifetime relation ######
    if mass_lifetime_file is None:
        lifetime = None
    else:
        file = h5py.File(mass_lifetime_file, "r")
        lifetime = { }
        for key in file.keys():
            key1 = key.replace("Z_","Z=").replace("_",".")  # in case it is Z_0_0001, as is the case for yields2.h5
            lifetime[key1] = { }
            lifetime[key1]['Mini'] = file[key]['MassLifetime'][:,0].astype(float)
            lifetime[key1]['lifetime'] = file[key]['MassLifetime'][:,1].astype(float)
        file.close()
    mass_lifetime = MassLifetime(lifetime=lifetime)
    ###### Load the mass-lifetime relation ######

    ###### Determine the primordial gas ######
    mass = SFH['Mtot'] / SFE  # The total mass of the primordial gas in the unit of Msun
    solar_set = "Default"  # The solar abundance set
    pr_gas = pg.Primordial_gas(mass=mass, Z_0=Z_0, 
                               input_array=input_primodiral_gas
                               ).add_metals(abund_table=constants.abund_tables[solar_set])
    ###### Determine the primordial gas ######

    ###### Initialize other essential variables ######
    YieldsTable = np.empty((len(SFH['Age']), 2), dtype=object)
    GasElement = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
    ZGas = np.zeros(len(SFH['Age']), dtype=np.float64)
    Nstar = np.zeros(len(SFH['Age']), dtype=np.float64)
    StarInitElement = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
    EjectElement = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
    SNIaElement = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
    SNccElement = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
    AGB_Element = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
    SNIaNum = np.zeros(len(SFH['Age']), dtype=np.float64)
    SNccNum = np.zeros(len(SFH['Age']), dtype=np.float64)
    AGB_Num = np.zeros(len(SFH['Age']), dtype=np.float64)
    interp_Y = InterpolateYields(dfs, kind=interp_kind)
    ###### Initialize other essential variables ######

    ###### Step 0: Initialize the first age, the primordial gas ######
    # The definition of the following variables is can be found in the excel file Variables.xlsx
    GasElement[0] = pr_gas['Gas']
    ZGas[0] = GasElement[0, 3:].sum() / GasElement.sum()  # The metallicity of the primordial gas
    StarInitElement[0] = SFH['Mstar'][0] * (GasElement[0]/GasElement[0].sum())
    if imf_evolve is not None:
        imf1 = imf_evolve(ZGas[0], SFH['Age'][0])
        # Normalize the IMF
        imf_norm = quad(imf1, constants.Mstar_min, constants.Mstar_max, 
                        epsrel=epsrel, limit=limit, full_output=1)[0]
        imf = lambda m: imf1(m)/imf_norm
        SNIa = SupernovaeIa.SNIa(imf, IMF_type=imf_key, p_preset=p_preset)
    Nstar[0] = SFH['Mstar'][0] / quad(lambda m: m*imf(m), constants.Mstar_min, constants.Mstar_max, 
                                      epsrel=epsrel, limit=limit, full_output=1)[0]
    StellarMass = SFH['Mstar'].cumsum()
    print("Step 0: Initialize the first age, the primordial gas")
    dEjectElement = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
    dSNIaElement = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
    dSNccElement = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
    dAGB_Element = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
    dSNIaNum = np.zeros(len(SFH['Age']), dtype=np.float64)
    dSNccNum = np.zeros(len(SFH['Age']), dtype=np.float64)
    dAGB_Num = np.zeros(len(SFH['Age']), dtype=np.float64)
    # mass boundary of each time bin
    dStellarMass = np.zeros(len(SFH['Age']), dtype=np.float64)
    mass_bounds = np.array([constants.Mstar_max]+
                                [mass_lifetime.lifetime_to_mass(SFH['TimeBin'][0:k+1].sum(), ZGas[0], kind=interp_kind) \
                                for k in range(0, len(SFH['Age'])-1)])
    for j in tqdm(range(0, len(SFH['Age'])-1)):
        # Calculate the mass of ejceta from AGB and SNcc happened in this time bin,
        # whose progenitor stars formed at Age[0].
        for elem in ElemIntersect["AGB+SNcc"][groups["AGB+SNcc"][0]]:
            ElemIndex = constants.elem_names.index(elem)
            interp = interp_Y.interps['AGB+SNcc'](ZGas[0], elem)
            dEjectElement[j, ElemIndex] = Nstar[0]*(quad(lambda m:
                                            imf(m)*interp(m),
                                            mass_bounds[j+1], mass_bounds[j],
                                            epsrel=epsrel, limit=limit, full_output=1)[0])
            ###### interp = interp_Y.interps['AGB+SNcc'](ZGas[0], elem) is a MUST! ######
            # You SHOULD NOT use the following code:
            # dEjectElement[j, ElemIndex] = Nstar[0]*(quad(lambda m:
            #                                 imf(m)*interp_Y.interps['AGB+SNcc'](ZGas[0], elem)(m),
            #                                 mass_bounds[j+1], mass_bounds[j],
            #                                 epsrel=epsrel, limit=limit, full_output=1)[0])
            # It will slow down the calculation.
            # quad will call the interp_Y.interps['AGB+SNcc'](ZGas[0], elem) again and again, which is not necessary
            # and slow down the calculation.
            ###### interp = interp_Y.interps['AGB+SNcc'](ZGas[0], elem) is a MUST! ######
        
        # Calculate the mass of ejecta from SNcc and the number of SNcc events happened in this time bin,
        # whose progenitor stars formed at Age[0].
        # Both mass_bounds[j+1] and mass_bounds[j] could be in three ranges:
        # >constants.SNcc_max, < constants.SNcc_max, [constants.SNcc_min, constants.SNcc_max]
        # 3*3=9 cases.
        # Actually there are only 6 cases because mass_bounds[j+1] < mass_bounds[j].
        # However, it is still complicated.
        # Here, I provide a simple way to calculate the mass of ejecta from SNcc and the number of SNcc events 
        # happened in this time bin.
        # 1. Calculate max(constants.SNcc_min, mass_bounds[j+1]) and min(constants.SNcc_max, mass_bounds[j])
        # 2. If max(constants.SNcc_min, mass_bounds[j+1])>=min(constants.SNcc_max, mass_bounds[j]), no ejecta come from SNcc.
        # 3. If max(constants.SNcc_min, mass_bounds[j+1])<min(constants.SNcc_max, mass_bounds[j]), 
        # all ejecta of SNcc come from the mass range 
        # [max(constants.SNcc_min, mass_bounds[j+1]), min(constants.SNcc_max, mass_bounds[j])].
        SNcc_low = np.maximum(constants.SNcc_min, mass_bounds[j+1])
        SNcc_high = np.minimum(constants.SNcc_max, mass_bounds[j])
        if SNcc_high > SNcc_low:
            dSNccNum[j] = Nstar[0]*(quad(lambda m: imf(m), SNcc_low, SNcc_high,
                                        epsrel=epsrel, limit=limit, full_output=1)[0])
            for elem in ElemIntersect["AGB+SNcc"][groups["AGB+SNcc"][0]]:
                interp = interp_Y.interps['AGB+SNcc'](ZGas[0], elem)
                ElemIndex = constants.elem_names.index(elem)
                dSNccElement[j, ElemIndex] = Nstar[0]*(quad(lambda m:
                                            imf(m)*interp(m),
                                            SNcc_low, SNcc_high,
                                            epsrel=epsrel, limit=limit, full_output=1)[0])
        else:
            pass
        # Calculate the mass of "ejecta" from AGB stars and the number of AGB events happened in this time bin,
        # whose progenitor stars formed at Age[0].
        # Again, both mass_bounds[j+1] and mass_bounds[j] could be in three ranges, ...
        # We adopt the same strategy as we did for SNcc.
        AGB_low = np.maximum(constants.AGB_min, mass_bounds[j+1])
        AGB_high = np.minimum(constants.AGB_max, mass_bounds[j])
        if AGB_high > AGB_low:
            dAGB_Num[j] = Nstar[0]*(quad(lambda m: imf(m), AGB_low, AGB_high,
                                        epsrel=epsrel, limit=limit, full_output=1)[0])
            for elem in ElemIntersect["AGB+SNcc"][groups["AGB+SNcc"][0]]:
                interp = interp_Y.interps['AGB+SNcc'](ZGas[0], elem)
                ElemIndex = constants.elem_names.index(elem)
                dAGB_Element[j, ElemIndex] = Nstar[0]*(quad(lambda m:
                                            imf(m)*interp(m),
                                            AGB_low, AGB_high,
                                            epsrel=epsrel, limit=limit, full_output=1)[0])
        else:
            pass
        # Calculate the mass of stars that explode or evolved as AGB stars in this time bin,
        # whose progenitor stars formed at Age[0].
        dStellarMass[j+1] = Nstar[0]*(quad(lambda m: m*imf(m), mass_bounds[j+1], mass_bounds[j],
                                            epsrel=epsrel, limit=limit, full_output=1)[0])
        dSNIaNum[j] = SNIa.number(SFH['Mstar'][0], SFH['Age'][j]-SFH['Age'][0], SFH['Age'][j+1]-SFH['Age'][0])

    StellarMass -= dStellarMass.cumsum()
    # The ejected mass of elements from AGB and SNcc, whose progenitor stars formed at Age[0]
    EjectElement += dEjectElement
    # The ejected mass of elements from SNcc, whose progenitor stars formed at Age[0]
    SNccElement += dSNccElement
    # The ejected mass of elements from AGB, whose progenitor stars formed at Age[0]
    AGB_Element += dAGB_Element
    # The number of SNcc events happened in each time bin, whose progenitor stars formed at Age[0]
    SNccNum += dSNccNum
    # The number of AGB events happened in each time bin, whose progenitor stars formed at Age[0]
    AGB_Num += dAGB_Num
    Zindex = interp_Y.Zindex
    if SNIaOn:
        # Calculate the number of SNIa events happened in this time bin, progenitor stars formed at Age[0]
        for elem in ElemIntersect["SNIa"][groups["SNIa"][0]]:
            y = interp_Y.interps['SNIa'](ZGas[0], elem)
            ElemIndex = constants.elem_names.index(elem)
            dSNIaElement[:, ElemIndex] += dSNIaNum * y
        SNIaNum += dSNIaNum
        dEjectElement += dSNIaElement
        EjectElement += dSNIaElement
        SNIaElement += dSNIaElement
        YieldsTable[0,0] = "AGB+SNcc: %s, SNIa: %s" % (groups_Z["AGB+SNcc"][Zindex["AGB+SNcc"]['low']],
                                                       groups_Z["SNIa"][Zindex["SNIa"]['low']])
        YieldsTable[0,1] = "AGB+SNcc: %s, SNIa: %s" % (groups_Z["AGB+SNcc"][Zindex["AGB+SNcc"]['high']],
                                                       groups_Z["SNIa"][Zindex["SNIa"]['high']])
    else:
        pass
        YieldsTable[0,0] = "AGB+SNcc: %s" % groups_Z["AGB+SNcc"][Zindex["AGB+SNcc"]['low']]
        YieldsTable[0,1] = "AGB+SNcc: %s" % groups_Z["AGB+SNcc"][Zindex["AGB+SNcc"]['high']]
    if imf_evolve is None:
        if imf_num>1 and SFH['Mstar'][0]>0:
            # It is time-varying IMF
            # and star formation rate at the first age is non-zero.
            # We need to update the IMF and SNIa class.
            imf_idx += 1
            if imf_key == "PowerLaw":
                imf = IMF(IMF_type=imf_key, power_index=imf_value[imf_idx]).imf
            elif imf_key == "Custom":
                imf = IMF(IMF_type=imf_key, IMF_arr=imf_value[imf_idx]).imf
            else:
                pass
            SNIa = SupernovaeIa.SNIa(imf, IMF_type=imf_key, p_preset=p_preset)
        else:
            pass
    else:
        pass
    ###### Step 0: Initialize the first age, the primordial gas ######

    ###### Step 1: Calculate the remaining ages ######
    print("Step 1: Calculate the remaining ages")
    for i in tqdm(range(1, len(SFH['Age'])-1)):
        GasElement[i] = GasElement[i-1] + EjectElement[i-1] - StarInitElement[i-1]
        ZGas[i] = GasElement[i, 3:].sum() / GasElement[i].sum()
        dEjectElement = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
        dEjectElement = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
        dSNIaElement = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
        dSNccElement = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
        dAGB_Element = np.zeros((len(SFH['Age']), len(constants.abund_tables[solar_set])), dtype=np.float64)
        dSNIaNum = np.zeros(len(SFH['Age']), dtype=np.float64)
        dSNccNum = np.zeros(len(SFH['Age']), dtype=np.float64)
        dAGB_Num = np.zeros(len(SFH['Age']), dtype=np.float64)
        dStellarMass = np.zeros(len(SFH['Age']), dtype=np.float64)
        if SFH['Mstar'][i] > 0:
            if imf_evolve is not None:
                imf1 = imf_evolve(ZGas[i], SFH['Age'][i])
                # Normalize the IMF
                imf_norm = quad(imf1, constants.Mstar_min, constants.Mstar_max, 
                                epsrel=epsrel, limit=limit, full_output=1)[0]
                imf = lambda m: imf1(m)/imf_norm
                SNIa = SupernovaeIa.SNIa(imf, IMF_type=imf_key, p_preset=p_preset)
            else:
                pass
            StarInitElement[i] = SFH['Mstar'][i] * (GasElement[i]/GasElement[i].sum())
            Nstar[i] = SFH['Mstar'][i] / quad(lambda m: m*imf(m), constants.Mstar_min, constants.Mstar_max,
                                              epsrel=epsrel, limit=limit, full_output=1)[0]
            mass_bounds = np.array([constants.Mstar_max]+
                                    [mass_lifetime.lifetime_to_mass(SFH['TimeBin'][i:k+1].sum(), ZGas[i], kind=interp_kind)\
                                    for k in range(i, len(SFH['Age'])-1)])
            # If the star formation rate at Age[i] is non-zero, 
            # the stars formed at Age[i] will explode and enrich the gas in the following ages.
            for j in range(i, len(SFH['Age'])-1):
                # Calculate the mass of ejceta from AGB and SNcc happened in this time bin
                for elem in ElemIntersect["AGB+SNcc"][groups["AGB+SNcc"][0]]:
                    # Interpolate the stellar yields of different initial metallicities
                    interp = interp_Y.interps['AGB+SNcc'](ZGas[i], elem)
                    ElemIndex = constants.elem_names.index(elem)
                    dEjectElement[j, ElemIndex] = Nstar[i]*(quad(lambda m:
                                                    imf(m)*interp(m),
                                                    mass_bounds[j+1-i], mass_bounds[j-i],
                                                    epsrel=epsrel, limit=limit, full_output=1)[0])
                # Calculate the mass of ejecta from SNcc and the number of SNcc events happened in this time bin,
                # whose progenitor stars formed at Age[i].
                SNcc_low = np.maximum(constants.SNcc_min, mass_bounds[j+1-i])
                SNcc_high = np.minimum(constants.SNcc_max, mass_bounds[j-i])
                if SNcc_high > SNcc_low:
                    dSNccNum[j] = Nstar[i]*(quad(lambda m: imf(m), SNcc_low, SNcc_high,
                                                epsrel=epsrel, limit=limit, full_output=1)[0])
                    for elem in ElemIntersect["AGB+SNcc"][groups["AGB+SNcc"][0]]:
                        interp = interp_Y.interps['AGB+SNcc'](ZGas[i], elem)
                        ElemIndex = constants.elem_names.index(elem)
                        dSNccElement[j, ElemIndex] = Nstar[i]*(quad(lambda m:
                                                    imf(m)*interp(m),
                                                    SNcc_low, SNcc_high,
                                                    epsrel=epsrel, limit=limit, full_output=1)[0])
                else:
                    pass
                # Calculate the mass of "ejecta" from AGB stars and the number of AGB events happened in this time bin,
                # whose progenitor stars formed at Age[i].
                AGB_low = np.maximum(constants.AGB_min, mass_bounds[j+1-i])
                AGB_high = np.minimum(constants.AGB_max, mass_bounds[j-i])
                if AGB_high > AGB_low:
                    dAGB_Num[j] = Nstar[i]*(quad(lambda m: imf(m), AGB_low, AGB_high,
                                                epsrel=epsrel, limit=limit, full_output=1)[0])
                    for elem in ElemIntersect["AGB+SNcc"][groups["AGB+SNcc"][0]]:
                        interp = interp_Y.interps['AGB+SNcc'](ZGas[i], elem)
                        ElemIndex = constants.elem_names.index(elem)
                        dAGB_Element[j, ElemIndex] = Nstar[i]*(quad(lambda m:
                                                    imf(m)*interp(m),
                                                    AGB_low, AGB_high,
                                                    epsrel=epsrel, limit=limit, full_output=1)[0])
                else:
                    pass
                # Calculate the mass of stars that explode or evolved as AGB stars in this time bin,
                # whose progenitor stars formed at Age[i].
                dStellarMass[j+1] = Nstar[i]*(quad(lambda m: m*imf(m), mass_bounds[j+1-i], mass_bounds[j-i],
                                                   epsrel=epsrel, limit=limit, full_output=1)[0])
                dSNIaNum[j] = SNIa.number(SFH['Mstar'][i], SFH['Age'][j]-SFH['Age'][i], SFH['Age'][j+1]-SFH['Age'][i])

            StellarMass -= dStellarMass.cumsum()
            # The ejected mass of elements from AGB and SNcc, whose progenitor stars formed at Age[i].
            EjectElement += dEjectElement
            # The ejected mass of elements from SNcc, whose progenitor stars formed at Age[i].
            SNccElement += dSNccElement
            # The ejected mass of elements
            AGB_Element += dAGB_Element
            # The number of SNcc events happened in each time bin, whose progenitor stars formed at Age[i].
            SNccNum += dSNccNum
            # The number of AGB events happened in each time bin, whose progenitor stars formed at Age[i].
            AGB_Num += dAGB_Num
            Zindex = interp_Y.Zindex
            if SNIaOn:
                # Calculate the number of SNIa events happened in this time bin, progenitor stars formed at Age[i]
                for elem in ElemIntersect["SNIa"][groups["SNIa"][0]]:
                    ElemIndex = constants.elem_names.index(elem)
                    y = interp_Y.interps['SNIa'](ZGas[i], elem)
                    dSNIaElement[:, ElemIndex] += dSNIaNum * y
                SNIaNum += dSNIaNum
                dEjectElement += dSNIaElement
                EjectElement += dSNIaElement
                SNIaElement += dSNIaElement
                YieldsTable[i,0] = "AGB+SNcc: %s, SNIa: %s" % (groups_Z["AGB+SNcc"][Zindex["AGB+SNcc"]['low']],
                                                               groups_Z["SNIa"][Zindex["SNIa"]['low']])
                YieldsTable[i,1] = "AGB+SNcc: %s, SNIa: %s" % (groups_Z["AGB+SNcc"][Zindex["AGB+SNcc"]['high']],
                                                               groups_Z["SNIa"][Zindex["SNIa"]['high']])
            else:
                YieldsTable[i,0] = "AGB+SNcc: %s" % groups_Z["AGB+SNcc"][Zindex["AGB+SNcc"]['low']]
                YieldsTable[i,1] = "AGB+SNcc: %s" % groups_Z["AGB+SNcc"][Zindex["AGB+SNcc"]['high']]
            if imf_evolve is None:
                if imf_num>1:
                    # It is time-varying IMF
                    # and star formation rate at the first age is non-zero.
                    # We need to update the IMF and SNIa class.
                    imf_idx += 1
                    if imf_key == "PowerLaw":
                        imf = IMF(IMF_type=imf_key, power_index=imf_value[imf_idx]).imf
                    elif imf_key == "Custom":
                        imf = IMF(IMF_type=imf_key, IMF_arr=imf_value[imf_idx]).imf
                    else:
                        pass
                    SNIa = SupernovaeIa.SNIa(imf, IMF_type=imf_key, p_preset=p_preset)
                else:
                    pass
            else:
                pass
        else:
            # If no star formation between this age and the next age, 
            # there is no ejecta from the stars formed in this time interval.
            YieldsTable[i,0] = "No Star Formation"
            YieldsTable[i,1] = "No Star Formation"
            continue

    # The last age
    GasElement[len(SFH['Age'])-1] = GasElement[len(SFH['Age'])-2] + \
        EjectElement[len(SFH['Age'])-2] - StarInitElement[len(SFH['Age'])-2]
    ZGas[len(SFH['Age'])-1] = GasElement[len(SFH['Age'])-1, 3:].sum() / GasElement[len(SFH['Age'])-1].sum()
    YieldsTable[len(SFH['Age'])-1,0] = "No Star Formation"
    YieldsTable[len(SFH['Age'])-1,1] = "No Star Formation"
    ###### Step 1: Calculate the remaining ages ######

    ###### Close the opened hdf5 files ######
    for file in files.values():
        file.close()
    ###### Close the opened hdf5 files ######

    ####### Step 2: Save the results #######
    # Get the current date and time
    current_time = datetime.datetime.now()
    # Format the date and time as "YYYY-MM-DD-HH-MM"
    formatted_time = current_time.strftime('%Y-%m-%d-%H-%M')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if out_file is None:
        out_file = f"{output_dir}/{formatted_time}.h5"
    else:
        out_file = f"{output_dir}/{out_file}"
    
    with h5py.File(out_file, 'w') as f:
        if comments is None:
            comments = "No comments."
        f.attrs["Comments"] = comments
        if imf_num>1:
            f.attrs["IMF"] = "Variant IMF"
        else:
            f.attrs["IMF"] = "Invariant IMF"
        f.attrs["Creation time"] = formatted_time
        f.attrs["AGB+SNcc yield file"] = yield_files["AGB+SNcc"]
        f.attrs["SNIa yield file"] = yield_files["SNIa"]
        f.attrs["SFE"] = SFE
        f.attrs["Z_0"] = Z_0
        f.attrs["Interpolation kind"] = interp_kind
        f.attrs["Solar set"] = solar_set
        if mass_lifetime_file is not None:
            f.attrs["Mass lifetime file"] = mass_lifetime_file
        else:
            f.attrs["Mass lifetime file"] = "None"
        f.attrs["SFH file"] = SFH['File']
        f.attrs['Concerned Elements'] = ElemNotice
        if SNIaOn:
            f.attrs["SNIaOn"] = "Yes"
        else:
            f.attrs["SNIaOn"] = "No"
        
        f.create_group("Star")
        f["Star"].attrs["Comments"] = "Mass in Msun, Time in Gyr, SFR in Msun/yr, Mtot in Msun"
        f["Star"].attrs["SFE"] = SFE
        f["Star"].attrs["Mtot"] = SFH["Mtot"]
        f.create_dataset("Star/StarInitElement", data=StarInitElement)
        f.create_dataset("Star/SNIaNum", data=SNIaNum)
        f.create_dataset("Star/SNccNum", data=SNccNum)
        f.create_dataset("Star/AGB_Num", data=AGB_Num)
        f.create_dataset("Star/Age", data=SFH['Age'])
        f.create_dataset("Star/SFR", data=SFH['SFR'])
        f.create_dataset("Star/Mstar", data=SFH['Mstar'])
        f.create_dataset("Star/TimeBin", data=SFH['TimeBin'])
        f.create_dataset("Star/Nstar", data=Nstar)
        f.create_dataset("Star/StellarMass", data=StellarMass)

        f.create_group("Gas")
        f.create_dataset("Gas/GasElement", data=GasElement)
        f.create_dataset("Gas/ZGas", data=ZGas)
        f.create_dataset("Gas/EjectElement", data=EjectElement)
        f.create_dataset("Gas/SNIaElement", data=SNIaElement)
        f.create_dataset("Gas/SNccElement", data=SNccElement)
        f.create_dataset("Gas/AGB_Element", data=AGB_Element)
        f.create_dataset("Gas/YieldsTable", data=np.array(YieldsTable, dtype='S'))

        f.create_group("IMF")
        if imf_evolve is None and imf_num>1:
            f["IMF"].attrs["IMF Type"] = imf_key
            if imf_key == "PowerLaw":
                PowerIndex_data = np.zeros((imf_idx, 2), dtype=np.float64)
                imf_idx = 0
                for i in range(len(SFH['Age'])-1):
                    if SFH['Mstar'][i]>0:
                        PowerIndex_data[imf_idx, 0] = SFH['Age'][i]
                        PowerIndex_data[imf_idx, 1] = imf_value[imf_idx]
                        imf_idx += 1
                    else:
                        continue
                f.create_dataset("IMF/PowerLawIndex", data=PowerIndex_data)
            elif imf_key == "Custom":
                imf_idx = 0
                for i in range(len(SFH['Age'])-1):
                    if SFH['Mstar'][i]>0:
                        f.create_dataset("IMF/%0.6eyr"%SFH['Age'][i], data=imf_value[imf_idx])
                        imf_idx += 1
                    else:
                        continue
            else:
                pass
        elif imf_evolve is None and imf_num<=1:
            f["IMF"].attrs["IMF Type"] = imf_key
            if imf_key == "PowerLaw":
                f.create_dataset("IMF/PowerLawIndex", data=imf_value[0])
            elif imf_key == "Custom":
                f.create_dataset("IMF/IMF", data=imf_value[0])
            else:
                pass
        else:
            # imf_evolve is not None
            f["IMF"].attrs["IMF Type"] = "IMF dependent on metallicity or age of the gas"
            f["IMF"].attrs["IMF"] = inspect.getsource(imf_evolve)
            m = np.logspace(np.log10(constants.Mstar_min), np.log10(constants.Mstar_max), 1000)
            for i in range(len(SFH['Age'])-1):
                if SFH['Mstar'][i]>0:
                    imf1 = imf_evolve(ZGas[i], SFH['Age'][i])
                    # Normalize the IMF
                    imf_norm = quad(imf1, constants.Mstar_min, constants.Mstar_max,
                                    epsrel=epsrel, limit=limit, full_output=1)[0]
                    imf = lambda m: imf1(m)/imf_norm
                    imf_data = np.zeros((len(m), 2), dtype=np.float64)
                    imf_data[:, 0] = m
                    imf_data[:, 1] = np.array([imf(a) for a in m])
                    f.create_dataset("IMF/%0.6eyr"%SFH['Age'][i], data=imf_data)
                else:
                    continue
    ####### Step 2: Save the results #######

    return GasElement, ZGas, EjectElement, SNIaElement, SNccElement, AGB_Element, StarInitElement, YieldsTable,\
           SNIaNum, SNccNum, AGB_Num, Nstar, StellarMass
