# -*- coding: utf-8 -*-

"""
Absolute Salinity (SA) and Preformed Salinity (Sstar)

Functions:
----------
  SA_from_SP(SP, p, lon, lat)
        Absolute Salinity from Practical Salinity
  Sstar_from_SP(SP, p, lon, lat)
        Preformed Salinity from Practical Salinity
  SA_Sstar_from_SP(SP, p, lon, lat)
        Absolute Salinity & Preformed Salinity from Practical Salinity

This is part of the python Gibbs Sea Water library
http://code.google.com/p/python-seawater/

"""

from __future__ import division

import numpy as np
import seawater.constants as cte
import library as lib
from library import match_args_return, _delta_SA

# -----------------------

__all__ = ['SA_from_SP',
           'Sstar_from_SP',
           'SA_Sstar_from_SP']

# ----------------------

@match_args_return
def SA_from_SP(SP, p, lon, lat):
    r"""
    Calculates Absolute Salinity from Practical Salinity.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SA : masked array
         Absolute salinity [g kg :sup:`-1`]

    See Also
    --------
    _delta_SA, _SA_from_SP_Baltic

    Notes
    -----
    Since SP is non-negative by definition, this function changes
    any negative input values of SP to be zero.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SP = [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat = 188, 4
    >>> gsw.SA_from_SP(SP, p, lon, lat)
    array([ 34.71177971,  34.89152372,  35.02554774,  34.84723008,
            34.7366296 ,  34.73236186])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm
    for estimating Absolute Salinity in the global ocean. Submitted to Ocean
    Science. A preliminary version is available at Ocean Sci. Discuss.,
    6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    #SP[SP < 0] = 0
    SP = np.maximum(SP, 0)

    dSA = lib._delta_SA( p, lon, lat )

    SA = ( cte.SSO / 35 ) * SP + dSA
    SA_baltic = lib._SA_from_SP_Baltic( SP, lon, lat )

    if SA_baltic is not None:
        SA[~SA_baltic.mask] = SA_baltic[~SA_baltic.mask]

    return SA

# ------------------

@match_args_return
def Sstar_from_SP(SP, p, lon, lat):
    r"""
    Calculates Preformed Salinity from Absolute Salinity.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    Sstar : masked array
            Preformed Salinity [g kg :sup:`-1`]

    See Also
    --------
    _delta_SA, _SA_from_SP_Baltic

    Notes
    -----
    The mask is only set when the observation is well and truly on dry
    land; often the warning flag is not set until one is several hundred
    kilometers inland from the coast.

    Since SP is non-negative by definition, this function changes any negative
    input values of SP to be zero.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SP = [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat =  188, 4
    >>> gsw.Sstar_from_SP(SP, p, lon, lat)
    array([ 34.7115532 ,  34.89116101,  35.02464926,  34.84359277,
            34.7290336 ,  34.71967638])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm
    for estimating Absolute Salinity in the global ocean. Submitted to Ocean
    Science. A preliminary version is available at Ocean Sci. Discuss.,
    6, 215-242.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    lon, lat, p, SP = np.broadcast_arrays(lon, lat, p, SP)

    #SP[SP < 0] = 0
    SP.clip(0, np.inf)

    dSA = lib._delta_SA( p, lon, lat )
    Sstar = (cte.SSO/35.) * SP - cte.r1 * dSA

    # In the Baltic Sea, SA == Sstar.
    Sstar_baltic = lib._SA_from_SP_Baltic( SP, lon, lat )

    if Sstar_baltic is not None:
        Sstar[~Sstar_baltic.mask] = Sstar_baltic[~Sstar_baltic.mask]

    return Sstar

# -----------------

def SA_Sstar_from_SP(SP, p, lon, lat):
    r"""
    Calculates Absolute Salinity and Preformed Salinity from Practical Salinity.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SA : masked array
         Absolute salinity [g kg :sup:`-1`]
    Sstar : masked array
            Preformed Salinity [g kg :sup:`-1`]

    See Also
    --------
    _delta_SA, _SA_from_SP_Baltic

    Notes
    -----
    In the Baltic Sea, Sstar == SA.

    The mask is only set when the observation is well and truly on dry
    land; often the warning flag is not set until one is several hundred
    kilometers inland from the coast.

    Since SP is non-negative by definition, this function changes any negative
    input values of SP to be zero.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SP = [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat =  188, 4
    >>> gsw.SA_Sstar_from_SP(SP, p, lon, lat)
    array([[ 34.71177971,  34.89152372,  35.02554774,  34.84723008,
             34.7366296 ,  34.73236186],
           [ 34.7115532 ,  34.89116101,  35.02464926,  34.84359277,
             34.7290336 ,  34.71967638]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm
    for estimating Absolute Salinity in the global ocean. Submitted to Ocean
    Science. A preliminary version is available at Ocean Sci. Discuss.,
    6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    return SA_from_SP(SP, p, lon, lat), Sstar_from_SP(SP, p, lon, lat)

# -------------------

if __name__=='__main__':
    import doctest
    doctest.testmod()
