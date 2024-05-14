# -*- coding: utf-8 -*-
"""
Created on May 11, 2024
Last modified on May 14, 2024
@Author: Guan-Fu Liu

Define the IMF class
"""
import constants
import numpy as np
from scipy import interpolate
from scipy.integrate import quad

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
        
        
        Returns
        -------
        self.imf: function
            The IMF function corresponding to the IMF type.
        self.Salpeter: function
            The Salpeter IMF function.
        self.Chabrier: function
            The Chabrier IMF function.
        self.Kroupa: function
            The Kroupa IMF function.
        self.Custom: function
            The custom IMF function, which is defined by the input array.

        
        Examples
        --------
        >>> from IMF import IMF
        >>> m = np.linspace(0.08, 150, 1000)
        >>> xi = np.array([IMF(IMF_type='Kroupa').imf(a) for a in m])
        # It is equivalent to
        >>> xi = np.array([IMF(IMF_type='Kroupa').Kroupa(a) for a in m])
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
            print("Error: The IMF type %s is not supported." % self.IMF_type)
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
        if m<constants.Mstar_min or m>constants.Mstar_max:
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
    
### This part is tested on May 14, 2024.
###### The IMF class ######