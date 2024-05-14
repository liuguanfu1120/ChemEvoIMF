# -*- coding: utf-8 -*-
"""
Created on May 12, 2024
Last modified on May 14, 2024
@Author: Guan-Fu Liu

Define the SNIa class
"""


###### Get the rate of SNIa per stellar mass ######
# The rate of SNIa per stellar mass can be prescribed as 
# the Equation 13 in Maoz & Mannucci (2012), https://ui.adsabs.harvard.edu/abs/2012PASA...29..447M/abstract
# which is a constraint from observation.
# Note that it assumes a diet Salpeter IMF, which should be calibrated to the adotped IMF.
# The calibration factor is the Equation 4 in Yan et al. (2019), https://www.aanda.org/10.1051/0004-6361/201936029

from scipy.integrate import quad
from IMF import IMF

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
            It is used to conduct the calibration
        
        Returns
        -------
        self.rate: function
            The rate of SNIa per stellar mass.
            The unit is 1/yr/solar mass.
        self.number: function
            The number of SNIa in the time interval [tmin, tmax] for the stellar population with the total mass Mstar.
        
        Examples
        --------
        >>> from IMF import IMF
        >>> import SupernovaeIa
        >>> SNIa = SupernovaeIa.SNIa(IMF(IMF_type='Kroupa').imf)
        >>> print(SNIa.number(1e9, 0, 4e7))
        0        
        # The total mass of the formed stars is 1e9 solar mass.
        # The number of SNIa in the time interval [0, 40 Myr] is 0.
        >>> print(SNIa.number(1e9, 4e7, 6e7))
        211325.61803158553       
        # The total mass of the formed stars is 1e9 solar mass.
        # The number of SNIa in the time interval [40 Myr, 60 Myr] is 211325.61803158553.
        # Here we do not truncate the number of SNIa to be an integer.
        """
        self.imf = imf
        imf_diet = IMF(IMF_type='DietSalpeter').imf  # The IMF function of the diet Salpeter IMF.
        p_denominator = (quad(imf_diet, 1.5, 8)[0])**2/\
                        (quad(imf_diet, 0.08, 150)[0]*quad(lambda m: m*imf_diet(m), 0.08, 150)[0])
        p_numerator = (quad(imf, 1.5, 8)[0])**2/\
                        (quad(imf, 0.08, 150)[0]*quad(lambda m: m*imf(m), 0.08, 150)[0])
        self.p = p_numerator/p_denominator
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

        Returns
        -------
        rate: float
            The rate of SNIa per stellar mass.
            The unit is 1/yr/solar mass.
        """
        y = self._rate(t)*self.p
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