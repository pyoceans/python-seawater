# -*- coding: utf-8 -*-
import numpy as np

"""
Constants used
==============
"""
#TODO: make a commom constant file

db2Pascal  = 1e4
"""
Decibar to pascal
"""

gamma = 2.26e-7
"""
TODO: define Gamma
"""

def _specvol_SSO_0_CT25(p):
    """
    This function calcualtes specifc volume at the Standard Ocean Salinty, SSO, and at a Conservative Temperature of zero degrees C, as a function of pressure, p, in dbar, using a streamlined version of the 25-term CT version of specific volume, that is, a streamlined version of the code_rho_alpha_beta_CT25(SA,CT,p).

    Parameters
    ----------
    p : array_like
        pressure [db]

    Returns
    -------
    specvol_SSO_0_CT25 : array_like
                         Specific volume at (SSO, CT=0, p), 25-term equation.

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    TODO

    References
    ----------
    TODO

    Modifications
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    p = np.asarray(p)

    SSO = 35.16504 * np.ones( p.shape )
    specvol_SSO_0_CT25 = (1.00000000e+00 + SSO * ( 2.0777716085618458e-003 +np.sqrt(SSO) * 3.4688210757917340e-006) + p * 6.8314629554123324e-006) / (9.9984380290708214e+002 + SSO * ( 2.8925731541277653e+000 + SSO * 1.9457531751183059e-003) + p * ( 1.1930681818531748e-002 + SSO * 5.9355685925035653e-006 + p * -2.5943389807429039e-008) )

    return specvol_SSO_0_CT25

def  _enthalpy_SSO_0_CT25(p):
    """
     This function calcualtes enthalpy at the Standard Ocean Salinty, SSO, and at a Conservative Temperature of zero degrees C, as a function of pressure, p, in dbar, using a streamlined version of the 25-term CT version of the Gibbs function, that is, a streamlined version of the code _enthalpy_CT25(SA,CT,p)

     Parameters
     ----------
     p : array_like
        pressure [db]

    Returns
    -------
    enthalpy_CT25 : array_like
                    enthalpy_CT25 at (SSO, CT = 0, p), 25-term equation.

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    TODO

    References
    ----------
    TODO

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    p = np.asarray(p)

    SSO = 35.16504 * np.ones( p.shape )

    a0 = 1 + SSO * (2.0777716085618458e-3 + np.sqrt(SSO) * 3.4688210757917340e-6)
    a1 = 6.8314629554123324e-6
    b0 = 9.9984380290708214e2 + SSO * (2.8925731541277653e0 + SSO * 1.9457531751183059e-3)
    b1 = 0.5 * (1.1930681818531748e-2 + SSO * 5.9355685925035653e-6)
    b2 = -2.5943389807429039e-8
    A = b1 - np.sqrt(b1**2 - b0 * b2)
    B = b1 + np.sqrt(b1**2 - b0 * b2)

    part = ( a0 * b2 - a1 * b1) / (b2 * (B - A) )

    enthalpy_SSO_0_CT25 = db2Pascal * ( ( a1 / (2*b2) ) * np.log( 1 + p * ( 2 * b1 + b2 * p ) / b0 ) + part * np.log( 1 + ( b2 * p * (B - A) ) / (A * (B + b2 * p ) ) ) )

    return enthalpy_SSO_0_CT25

def  z_from_p(p, lat):
    """
    Calculates height from sea pressure using the computationally-efficient 25-term expression for density in terms of SA, CT and p.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90].
    p : array_like
        sea pressure [db]. Sea presure is absolute pressure - 10.1325 dbar (or minus atmospheric pressure)

    Returns
    -------
    z : array_like
        height [m]

    See Also
    --------
    TODO


    Examples
    --------
    TODO

    Notes
    -----
    At sea level z = 0, and since z (HEIGHT) is defined to be positive upwards, it follows that while z is positive in the atmosphere, it is NEGATIVE in the ocean.

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. Available from http://www.TEOS-10.org
    ,, [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater. FIXME: To be submitted to Ocean Science Discussions.
    .. [3] Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.

    Modifications:
    2010-08-26. Trevor McDougall, Claire Roberts-Thomson & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    p, lat = np.asarray(p), np.asarray(lat)

    X     = np.sin( np.deg2rad(lat) )
    sin2  = X**2
    B     = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2 )
    A     = -0.5 * gamma * B
    C     = _enthalpy_SSO_0_CT25(p)
    z     = -2 * C / ( B + np.sqrt( B**2 - 4 * A * C ) )

    return z

def grav(lat, p=0):
    """
    Calculates acceleration due to gravity as a function of latitude and as a function of pressure in the ocean.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90].
    p : number or array_like. Default p = 0
        pressure [db]

    Returns
    -------
    g : array_like
        gravity [m s :sup:`2`]

    See Also
    --------
    TODO

    Notes
    -----
    In the ocean z is negative.

    Examples
    --------
    TODO

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. Available from http://www.TEOS-10.org. Appendix D of this TEOS-10 Manual.
    .. [2] Moritz (2000) Goedetic reference system 1980. J. Geodesy,74,128-133.
    .. [3] Saunders, P.M., and N.P. Fofonoff (1976) Conversion of pressure to depth in the ocean. Deep-Sea Res.,pp. 109 - 111.

    Modifications:
    2010-07-23. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    lat, p, = np.asarray(lat), np.asarray(p)

    X = np.sin( np.deg2rad(lat) )
    sin2 = X**2
    gs = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2)
    z = z_from_p(p, lat)
    grav = gs * (1 - gamma * z) # z is the height corresponding to p
    return grav