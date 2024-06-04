# -*- coding: utf-8 -*-
"""
Created on May 27, 2024
Last modified on May 27, 2024
@Author: Guan-Fu Liu

Define basic functions to read the results from the output files.
"""
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, './')
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy import interpolate
import h5py
from tqdm import tqdm
import datetime
import os
# The following are the modules defined in the same directory
from MassLifetime import MassLifetime
from IMF import IMF
import SupernovaeIa
import primordial_gas as pg
import constants
from InterpolateYields import InterpolateYields

def read_results(file_path):
    """
    Read the results from the output file.

    Parameters
    ----------
    file_path : str
        The name of the output file.
    """
    f =  h5py.File(file_path, 'r')
    yield_files = {a.split()[0]: f.attrs[a] for a in f.attrs.keys() if 'yield file' in a}
    files = {key: h5py.File(value, 'r') for key, value in yield_files.items()}
    groups = {key: list(files[key].keys()) for key in files.keys()}  # The names of the groups in the yield table
    ElemNotice = list(f.attrs['Concerned Elements'])
    ###### Load the yield tables ######
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
    
    if f.attrs['SNIaOn'] == 'Yes':
        print("The yields from SNIa are taken into consideration in the result to be analysed.")
        SNIaOn = True
    elif f.attrs['SNIaOn'] == 'No':
        print("The yields from SNIa are ignored in the result to be analysed.")
        SNIaOn = False
    else:
        raise ValueError("The SNIaOFF attribute should be either 'Yes' or 'No'.")
    ###### Load the yield tables ######

    ##### Load the mass lifetime relation ######
    mass_lifetime_file = f.attrs['Mass lifetime file']
    if mass_lifetime_file == "None":
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
    ##### Load the mass lifetime relation ######

    ###### Load the Star group ######
    SFH = { }
    for key in ['Age', 'SFR', 'TimeBin', 'Mstar']:
        SFH[key] = f['Star/'+key][()]

    StarInitElement = f['Star/StarInitElement'][()]
    StellarMass = f['Star/StellarMass'][()]
    SNccNum = f['Star/SNccNum'][()]
    SNIaNum = f['Star/SNIaNum'][()]
    AGB_Num = f['Star/AGB_Num'][()]
    ###### Load the Star group ######

    ###### Load the Gas group ######
    GasElement = f['Gas/GasElement'][()]
    ZGas = f['Gas/ZGas'][()]
    EjectElement = f['Gas/EjectElement'][()]
    SNccElement = f['Gas/SNccElement'][()]
    SNIaElement = f['Gas/SNIaElement'][()]
    AGB_Element = f['Gas/AGB_Element'][()]
    YieldsTable = f['Gas/YieldsTable'][()]
    ###### Load the Gas group ######

    ###### Load the IMF group ######
    IMF_type = f["IMF"].attrs["IMF Type"]
    IMF_data = { }
    if IMF_type == "PowerLaw":
        IMF_data[IMF_type] = f['IMF/PowerLawIndex'][()]
    elif IMF_type == "Custom":
        IMF_data[IMF_type] = [ ]
        for key in f['IMF'].keys():
            IMF_data[IMF_type].append(f['IMF/%s'+key][()])
    elif IMF_type in ['Salpeter', 'Chabrier', 'Kroupa', 'DietSalpeter']:
        IMF_data[IMF_type] = "invariant"
    else:
        IMF_data["Function"] =  [ ]
        for key in f['IMF'].keys():
            IMF_data["Function"].append(f['IMF/%s'%key][()])
    ###### Load the IMF group ######

    ###### Load the other attributes ######
    creation_time = f.attrs['Creation time']
    SFE = float(f.attrs['SFE'])
    interp_kind = f.attrs['Interpolation kind']
    comments = f.attrs['Comments']
    solar_set = f.attrs['Solar set']
    ###### Load the other attributes ######
    f.close()

    return yield_files, groups, ElemNotice, dfs, ElemIntersect, SNIaOn, mass_lifetime_file, mass_lifetime, SFH, StarInitElement,\
            StellarMass, SNccNum, SNIaNum, AGB_Num, GasElement, ZGas, EjectElement, SNccElement, SNIaElement,\
            AGB_Element, YieldsTable, IMF_type, IMF_data, creation_time, SFE, interp_kind, comments, solar_set


def GetZi2Zj(GasElement, Zi, Zj, solar_set, fill_value=None):
    """
    To get the ratio of the element j to the element i in the solar unit.

    If the mass of the element j at some ages are zero, the ratio of the element j to the element i is set to fill_value.

    Parameters
    ----------
    GasElement : numpy array
        The mass of elements. The shape of GasElement should be (N, 32). Here N is them number of ages,
          32 is the number of elements
        (the first and the last elemen are empty and "others", should not be used).
    Zi : str
        The element i. It should be in the list of constants.elem_names, namely, from H to Zn (Z=1, 2, 3, ..., 30).
    Zj : str
        The element j. It should be in the list of constants.elem_names, namely, from H to Zn (Z=1, 2, 3, ..., 30).
    solar_set : str
        The solar abundance set. It should be in the list of constants.abund_tables.keys().
    fill_value : float or None, optional
        The value to fill the ratio of the element j to the element i if the mass of the element j or the element i is zero.
        If it is None, mask the ages where the mass of the element j is zero, which will change the shape of the output.
        If it is a float, it should be a positive number.

    Returns
    -------
    Zi2Zj: dict
        The ratio of the element j to the element i in the solar unit.
        Zi2Zj['%s/%s'%(Zi, Zj)] is the ratio of the element j to the element i, in the solar unit.
        Zi2Zj['[%s/%s]'%(Zi, Zj)] is log10(Zi2Zj['%s/%s'%(Zi, Zj)]).
        Zi2Zj['%s-mask'%Zj] is the mask of the ages where BOTH of the mass of the element j and i is POSITIVE.
    """
    if Zi not in constants.elem_names or Zj not in constants.elem_names:
        print("The element should be in the list of constants.elem_names, which is ", constants.elem_names)
        return None
    if solar_set not in constants.abund_tables.keys():
        print("The solar abundance set should be in the list of constants.abund_tables.keys(), which is ",
               constants.abund_tables.keys())
        return None
    Zi2Zj_solar = constants.abund_tables[solar_set][constants.elem_names.index(Zi)] / \
                  constants.abund_tables[solar_set][constants.elem_names.index(Zj)]
    mask = GasElement[:, constants.elem_names.index(Zj)] > 0
    mask = mask & (GasElement[:, constants.elem_names.index(Zi)] > 0)
    Zi2Zj = { }
    if fill_value is None:
        Zi2Zj['%s/%s'%(Zi, Zj)] = GasElement[mask, constants.elem_names.index(Zi)] / \
                                  GasElement[mask, constants.elem_names.index(Zj)] / Zi2Zj_solar
        Zi2Zj['[%s/%s]'%(Zi, Zj)] = np.log10(Zi2Zj['%s/%s'%(Zi, Zj)])
    else:
        Zi2Zj['%s/%s'%(Zi, Zj)] = np.zeros(len(GasElement))
        Zi2Zj['%s/%s'%(Zi, Zj)][mask] = GasElement[mask, constants.elem_names.index(Zi)] / \
                                        GasElement[mask, constants.elem_names.index(Zj)] / Zi2Zj_solar
        Zi2Zj['%s/%s'%(Zi, Zj)][~mask] = fill_value
        Zi2Zj['[%s/%s]'%(Zi, Zj)] = np.log10(Zi2Zj['%s/%s'%(Zi, Zj)])
    Zi2Zj['%s/%s-mask'%(Zi,Zj)] = mask
    return Zi2Zj