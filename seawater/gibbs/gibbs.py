# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np
import seawater.constants as cte
import library as lib
from library import match_args_return
from conversions import *




@match_args_return
def entropy_derivative_SA(SA, CT):
    r"""
    Calculates the derivatives of specific entropy (eta_SA) with respect to
    Absolute Salinity at constant Conservative Temperature.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_SA : array_like
             The derivative of specific entropy with respect to SA at constant
             CT [J g :sup:`-1` K :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_derivative_SA(SA, CT)
    array([ -0.2632868 ,  -0.26397728,  -0.2553675 ,  -0.23806659,
            -0.23443826,  -0.23282068])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.8) and (P.14a,c).

    Modifications:
    2010-08-21. Trevor McDougall.
    2011-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1

    pt = pt_from_CT(SA, CT)

    eta_SA = -(lib._gibbs(n1, n0, n0, SA, pt, 0) ) / (cte.Kelvin + pt)

    return eta_SA


@match_args_return
def entropy_derivative_CT(SA, CT):
    r"""
    Calculates the derivatives of specific entropy (eta_CT) with respect to
    Conservative Temperature at constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_CT : array_like
             The derivative of specific entropy with respect to CT at constant
             SA [ J (kg K :sup:`-2`) :sup:`-1` ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_derivative_CT(SA, CT)
    array([13.22103121,  13.23691119,  13.48900463,  14.08659902,
           14.25772958,  14.38642995])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.8) and (P.14a,c).

    Modifications:
    2010-08-21. Trevor McDougall.
    2011-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    pt = pt_from_CT(SA, CT)

    eta_CT = cte.cp0 / (cte.Kelvin + pt)

    return eta_CT


def entropy_first_derivatives(SA, CT):
    r"""
    Calculates the following two partial derivatives of specific entropy (eta)
    (1) eta_SA, the derivative with respect to Absolute Salinity at constant
        Conservative Temperature, and
    (2) eta_CT, the derivative with respect to Conservative Temperature at
        constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_SA : array_like
             The derivative of specific entropy with respect to SA at constant
             CT [J g :sup:`-1` K :sup:`-1`]
    eta_CT : array_like
             The derivative of specific entropy with respect to CT at constant
             SA [ J (kg K :sup:`-2`) :sup:`-1` ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_first_derivatives(SA, CT)
    array([[ -0.2632868 ,  -0.26397728,  -0.2553675 ,  -0.23806659,
             -0.23443826,  -0.23282068],
           [ 13.22103121,  13.23691119,  13.48900463,  14.08659902,
             14.25772958,  14.38642995]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.8) and (P.14a,c).

    Modifications:
    2010-08-21. Trevor McDougall.
    2011-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    return ( entropy_derivative_SA(SA, CT),
             entropy_derivative_CT(SA, CT) )


@match_args_return
def entropy_derivative_CT_CT(SA, CT):
    r"""
    Calculates the the second-order partial derivatives of specific entropy
    with respect to Conservative Temperature at constant Absolute
    Salinity (eta_CT_CT).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_CT_CT : array_like
                The second derivative of specific entropy with respect to CT at
                constant SA [J (kg K :sup:`3`) :sup:`-1` ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_derivative_CT_CT(SA, CT)
    array([-0.04366502, -0.04378134, -0.04550611, -0.04970894, -0.05093869,
           -0.05187502])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (P.14b) and (P.15a,b).

    Modifications:
    2010-08-23. Trevor McDougall and Paul Barker.
    2011-02-28. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n2 = 0, 2

    pt = pt_from_CT(SA, CT)
    abs_pt = 273.15 + pt

    CT_pt = -(abs_pt * lib._gibbs(n0, n2, n0, SA, pt, 0) ) / cte.cp0

    eta_CT_CT = -cte.cp0 / ( CT_pt * abs_pt**2 )

    return eta_CT_CT


@match_args_return
def entropy_derivative_SA_CT(SA, CT):
    r"""
    Calculates the the second-order partial derivatives of specific entropy
    with to Absolute Salinity and Conservative Temperature (eta_SA_CT).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_SA_CT : array_like
                The second derivative of specific entropy with respect to
                SA and CT [J (kg (g kg :sup:`-1` ) K :sup:`2`) :sup:`-1` ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_derivative_SA_CT(SA, CT)
    array([-0.0018331 , -0.00181947, -0.00158084, -0.00093011, -0.00071701,
           -0.00054841])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (P.14b) and (P.15a,b).

    Modifications:
    2010-08-23. Trevor McDougall and Paul Barker.
    2011-02-28. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1

    pt = pt_from_CT(SA, CT)
    abs_pt = 273.15 + pt

    CT_SA = (lib._gibbs(n1, n0, n0, SA, pt, 0) -
                    ( abs_pt * lib._gibbs(n1, n1, n0, SA, pt, 0 ) ) ) / cte.cp0

    eta_SA_CT = - CT_SA * entropy_derivative_CT_CT(SA, CT)

    return eta_SA_CT


@match_args_return
def entropy_derivative_SA_SA(SA, CT):
    r"""
    Calculates the the second-order partial derivatives of specific entropy
    with respect to Absolute Salinity at constant Conservative
    Temperature (eta_SA_SA).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_SA_SA : array_like
                The second derivative of specific entropy with respect to SA at
                constant CT [J (kg K (g kg :sup:`-1` ) :sup:`2`) :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_derivative_SA_SA(SA, CT)
    TODO

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (P.14b) and (P.15a,b).

    Modifications:
    2010-08-23. Trevor McDougall and Paul Barker.
    2011-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1, n2 = 0, 1, 2

    pt = pt_from_CT(SA, CT)
    abs_pt = 273.15 + pt

    CT_SA = (lib._gibbs(n1, n0, n0, SA, pt, 0) -
                    ( abs_pt * lib._gibbs(n1, n1, n0, SA, pt, 0 ) ) ) / cte.cp0

    eta_SA_CT = entropy_derivative_SA_CT(SA, CT)

    eta_SA_SA = -lib._gibbs(n2, n0, n0, SA, pt, 0) / abs_pt - CT_SA * eta_SA_CT

    return eta_SA_SA


def entropy_second_derivatives(SA, CT):
    r"""
    Calculates the following three second-order partial derivatives of specific
    entropy (eta)
    (1) eta_SA_SA, the second derivative with respect to Absolute Salinity at
        constant Conservative Temperature, and
    (2) eta_SA_CT, the derivative with respect to Absolute Salinity and
        Conservative Temperature.
    (3) eta_CT_CT, the second derivative with respect to Conservative
        Temperature at constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_SA_SA : array_like
                The second derivative of specific entropy with respect to SA at
                constant CT [J (kg K (g kg :sup:`-1` ) :sup:`2`) :sup:`-1`]
    eta_SA_CT : array_like
                The second derivative of specific entropy with respect to
                SA and CT [J (kg (g kg :sup:`-1` ) K :sup:`2`) :sup:`-1` ]
    eta_CT_CT : array_like
                The second derivative of specific entropy with respect to CT at
                constant SA [J (kg K :sup:`3`) :sup:`-1` ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_second_derivatives(SA, CT)
    TODO

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (P.14b) and (P.15a,b).

    Modifications:
    2010-08-23. Trevor McDougall and Paul Barker.
    2011-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    return ( entropy_derivative_SA_SA(SA, CT),
             entropy_derivative_SA_CT(SA, CT),
             entropy_derivative_CT_CT(SA, CT),
           )


@match_args_return
def CT_derivative_SA(SA, pt):
    r"""
    Calculates the derivatives of Conservative Temperature (CT_SA) with respect
    to Absolute Salinity at constant potential temperature (with pr = 0 dbar).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    CT_SA : array_like
            The derivative of CT with respect to SA at constant potential
            temperature reference sea pressure of 0 dbar.
            [K (g kg :sup:`-1`) :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> pt = [28.7832, 28.4209, 22.7850, 10.2305, 6.8292, 4.3245]
    >>> gsw.CT_derivative_SA(SA, pt)
    array([-0.04198109, -0.04155814, -0.03473921, -0.0187111 , -0.01407594,
           -0.01057172])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.3) and (A.12.9a,b).

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
    Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    Modifications:
    2010-08-05. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1
    abs_pt = cte.Kelvin + pt

    g100 = lib._gibbs(n1, n0, n0, SA, pt, 0)
    g110 = lib._gibbs(n1, n1, n0, SA, pt, 0)
    CT_SA = ( g100 - abs_pt * g110 ) / cte.cp0

    return CT_SA

@match_args_return
def CT_derivative_pt(SA, pt):
    r"""
    Calculates the derivatives of Conservative Temperature (CT_pt) with respect
    to potential temperature (the regular potential temperature which is
    referenced to 0 dbar) at constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    CT_pt : array_like
            The derivative of CT with respect to pt at constant SA.
            [ unitless ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> pt = [28.7832, 28.4209, 22.7850, 10.2305, 6.8292, 4.3245]
    >>> gsw.CT_derivative_pt(SA, pt)
    array([1.00281494,  1.00255482,  1.00164514,  1.00000377,  0.99971636,
           0.99947433])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.3) and (A.12.9a,b).

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
    Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    Modifications:
    2010-08-05. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n2 = 0, 2
    abs_pt = cte.Kelvin + pt

    g020 = lib._gibbs(n0, n2, n0, SA, pt, 0)
    CT_pt = - (abs_pt * g020 ) / cte.cp0

    return CT_pt

def CT_first_derivatives(SA, pt):
    r"""
    Calculates the following two derivatives of Conservative Temperature
    (1) CT_SA, the derivative with respect to Absolute Salinity at constant
        potential temperature (with pr = 0 dbar), and
    (2) CT_pt, the derivative with respect to potential temperature (the
        regular potential temperature which is referenced to 0 dbar) at
        constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    CT_SA : array_like
            The derivative of CT with respect to SA at constant potential
            temperature reference sea pressure of 0 dbar.
            [K (g kg :sup:`-1`) :sup:`-1`]

    CT_pt : array_like
            The derivative of CT with respect to pt at constant SA.
            [ unitless ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> pt = [28.7832, 28.4209, 22.7850, 10.2305, 6.8292, 4.3245]
    >>> gsw.CT_first_derivatives(SA, pt)
    array([[-0.04198109, -0.04155814, -0.03473921, -0.0187111 , -0.01407594,
            -0.01057172],
           [ 1.00281494,  1.00255482,  1.00164514,  1.00000377,  0.99971636,
             0.99947433]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.3) and (A.12.9a,b).

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
    Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    Modifications:
    2010-08-05. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    return CT_derivative_SA(SA, pt), CT_derivative_pt(SA, pt)


"""
Section D: extra functions for Depth, Pressure and Distance
"""



@match_args_return
def grav(lat, p=0):
    r"""
    Calculates acceleration due to gravity as a function of latitude and as a
    function of pressure in the ocean.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    p : number or array_like. Default p = 0
        pressure [dbar]

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
    >>> import seawater.gibbs as gsw
    >>> lat = [-90, -60, -30, 0]
    >>> p = 0
    >>> gsw.grav(lat, p)
    array([ 9.83218621,  9.81917886,  9.79324926,  9.780327  ])
    >>> gsw.grav(45)
    9.8061998770458008

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    .. [2] Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.

    .. [3] Saunders, P.M., and N.P. Fofonoff (1976) Conversion of pressure to
    depth in the ocean. Deep-Sea Res.,pp. 109 - 111.

    Modifications:
    2010-07-23. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    X = np.sin( np.deg2rad(lat) )
    sin2 = X**2
    gs = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2)
    z = z_from_p(p, lat)
    grav = gs * (1 - cte.gamma * z) # z is the height corresponding to p
    return grav

@match_args_return
def distance(lon, lat, p=0):
    r"""
    Calculates the distance in met res between successive points in the vectors
    lon and lat, computed using the Haversine formula on a spherical earth of
    radius 6,371 km, being the radius of a sphere having the same volume as
    Earth. For a spherical Earth of radius 6,371,000 m, one nautical mile is
    1,853.2488 m, thus one degree of latitude is 111,194.93 m.

    Haversine formula:
        R = earth's radius (mean radius = 6,371 km)

    .. math::
        a = \sin^2(\delta \text{lat}/2) + \cos(\text{lat}_1) \cos(\text{lat}_2) \sin^2(\delta \text{lon}/2)

        c = 2 \times \text{atan2}(\sqrt{a}, \sqrt{(1-a)})

        d = R \times c

    Parameters
    ----------
    lon : array_like
          decimal degrees east [0..+360] or [-180 ... +180]
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    p : number or array_like. Default p = 0
        pressure [dbar]

    Returns
    -------
    dist: array_like
          distance between points on a spherical Earth at pressure (p) [m]

    See Also
    --------
    TODO

    Notes
    -----
    z is height and is negative in the oceanographic.

    Distances are probably good to better than 1\% of the "true" distance on the
    ellipsoidal earth.

    The check value below differ from the original online docs at
    "http://www.teos-10.org/pubs/gsw/html/gsw_distance.html" but agree with the
    result.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> lon = [159, 220]
    >>> lat = [-35, 35]
    >>> gsw.distance(lon, lat)
    array([[ 10030974.652916]])
    >>> p = [200, 1000]
    >>> gsw.distance(lon, lat, p)
    array([[ 10030661.63878009]])
    >>> p = [[200], [1000]]
    >>> gsw.distance(lon, lat, p)
    array([[ 10030661.63878009],
           [ 10029412.58776001]])

    References
    ----------
    .. [1] http://www.eos.ubc.ca/~rich/map.html

    Modifications:
    2000-11-06. Rich Pawlowicz
    2010-07-28. Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """
    #FIXME? The argument handling seems much too complicated.
    # Maybe we can come up with some simple specifications of
    # what argument combinations are permitted, and handle everything
    # with broadcasting. - EF

    #FIXME: Eric what do you think? This assume p(stations, depth)
    lon, lat, = np.atleast_2d(lon), np.atleast_2d(lat)

    if (lon.size == 1) & (lat.size == 1):
        raise ValueError('more than one point is needed to compute distance')
    elif lon.ndim != lat.ndim:
        raise ValueError('lon, lat must have the same dimension')

    lon, lat, p = np.broadcast_arrays(lon, lat, p)

    dlon = np.deg2rad( np.diff(lon) )
    dlat = np.deg2rad( np.diff(lat) )

    a = ( ( np.sin(dlat/2.) )**2 + np.cos( np.deg2rad( lat[:,:-1] ) ) *
    np.cos( np.deg2rad( lat[:,1:] ) ) * ( np.sin(dlon/2.) )**2 )

    angles = 2. * np.arctan2( np.ma.sqrt(a), np.ma.sqrt(1-a) )

    p_mid = 0.5 * (   p[:,0:-1] +   p[:,0:-1] )
    lat_mid = 0.5 * ( lat[:,:-1] + lat[:,1:] )

    z = z_from_p(p_mid, lat_mid)

    distance = (cte.a + z) * angles

    return distance

"""
Section E: extra functions for Salinity
"""

"""cndr_from_SP == sw.cndr
Calculates Practical Salinity from conductivity ratio (R), using the PSS-78
algorithm. Note that the PSS-78 algorithm for Practical Salinity
is only valid in the range 2 < SP < 42.  The output, SP, of this function is
constrained to be non-negative.
"""
from  seawater.csiro import cndr as cndr_from_SP

"""SP_from_cndr == sw.salt
Calculates conductivity ratio (R) from (SP,t,p) using PSS-78. Note that the
PSS-78 algorithm for Practical Salinity is only valid in the range 2 < SP < 42.
"""
from  seawater.csiro import salt as SP_from_cndr




@match_args_return
def SA_from_rho(rho, t, p):
    r"""
    Calculates the Absolute Salinity of a seawater sample, for given values of
    its density, in situ temperature and sea pressure (in dbar).

    One use for this function is in the laboratory where a measured value of
    the in situ density :math:`\rho` of a seawater sample may have been made at
    the laboratory temperature :math:`t` and at atmospheric pressure :math:`p`.
    The present function will return the Absolute Salinity SA of this seawater
    sample.

    Parameters
    ----------
    rho : array_like
          in situ density [kg m :sup:`-3`]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    This is expressed on the Reference-Composition Salinity Scale of
    Millero et al. (2008).

    After two iterations of a modified Newton-Raphson iteration, the error in SA
    is typically no larger than 2 :math:`^\times` 10 :sup:`-13` [g kg :sup:`-1`]

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> rho = [1021.839, 1022.262, 1024.426, 1027.792, 1029.839, 1032.002]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.SA_from_rho(rho, t, p)
    array([ 34.71022966,  34.89057683,  35.02332421,  34.84952096,
            34.73824809,  34.73188384])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.5.

    .. [2] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008:
    The composition of Standard Seawater and the definition of the
    Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.

    Modifications:
    2010-08-23. Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """
    n0, n1 = 0, 1
    v_lab = np.ones( rho.shape ) / rho
    v_0 = lib._gibbs(n0, n0, n1, 0, t, p)
    v_120 = lib._gibbs(n0, n0, n1, 120, t, p)
    SA = 120 * ( v_lab - v_0 ) / ( v_120 - v_0 ) # initial estimate of SA
    Ior = (SA < 0) | (SA > 120)
    v_SA = ( v_120 - v_0 ) / 120 # initial estimate of v_SA, SA derivative of v

    for iter in range(0,2):
        SA_old = SA
        delta_v = lib._gibbs(n0, n0, n1, SA_old, t, p) - v_lab
        SA = SA_old - delta_v / v_SA # half way through the modified N-R method
        SA_mean = 0.5 * ( SA + SA_old )
        v_SA = lib._gibbs(n1, n0, n1, SA_mean, t, p)
        SA = SA_old - delta_v / v_SA

    SA[Ior] = np.ma.masked

    return SA

# ----------------------------


if __name__=='__main__':
    import doctest
    doctest.testmod()
