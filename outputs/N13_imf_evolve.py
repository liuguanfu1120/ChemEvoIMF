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
