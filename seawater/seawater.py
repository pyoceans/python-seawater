# -*- coding: utf-8 -*-

import numpy as np
import sys

"""
The International Practical Temperature Scale of 1968 (IPTS-68)
#:math:`T68  = 1.00024 * T90`
this linear transformation is accurate within 0.5 m C for conversion between IPTS-68 and ITS-90 over the oceanographic temperature range (Saunders,et al 1991).
"""
T68conv = 1.00024

"""
0.017453292519943295
"""
DEG2RAD = np.pi/180.

"""
A.E.Gill p.597
..:math:
  \Omega = \frac{2*\\pi}{\\textrm{sidereal day}}

1 sidereal day = 23.9344696 hours
units : radians/sec
"""
OMEGA = 7.292e-5

"""
Conductivity at S=35 psu , T=15 C [ITPS 68] and P=0 db)
units : mmho cm :sup:`-1` == mS cm :sup:`-1`
Reference: R.C. Millard and K. Yang 1992. "CTD Calibration and Processing Methods used by Woods Hole Oceanographic Institution"  Draft April 14, 1992 (Personal communication)
"""
C3515 = 42.914

"""
acceleration of gravity in m s :sup:`2`
"""
g = 9.8

""" Used by sw.dist
nautical miles to kilometers is
defined in Pond & Pickard p303
"""
NM2KM   = 1.8520
DEG2MIN = 60
DEG2NM  = 60

""" Used by solubility functions
"""
Kelvin = 273.15

""" Decibar to pascal
"""
db2Pascal  = 1e4

def adtg(s, t, p):
    """
    Calculates adiabatic temperature gradient as per UNESCO 1983 routines.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"

    Returns
    -------
    adtg : array_like
           adiabatic temperature gradient [:math:`^\\circ` C db :sup:`-1`]

    See Also
    --------
    ptmp

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested

    Examples
    --------
    Data from UNESCO 1983 p45

    >>> import numpy as np
    >>> import seawater as sw
    >>> from seawater import T68conv
    >>> t = np.array([[ 0,  0,  0,  0,  0,  0], [10, 10, 10, 10, 10, 10], [20, 20, 20, 20, 20, 20], [30, 30, 30, 30, 30, 30], [40, 40, 40, 40, 40, 40]])
    >>> t = t / T68conv
    >>> s = np.array([[25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35]])
    >>> p = np.array([[0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000]])
    >>> sw.adtg(s, t, p)
    array([[  1.68710000e-05,   1.04700000e-04,   1.69426000e-04,
              3.58030000e-05,   1.17956500e-04,   1.77007000e-04],
           [  1.00194580e-04,   1.60959050e-04,   2.06874170e-04,
              1.14887280e-04,   1.71364200e-04,   2.12991770e-04],
           [  1.73819840e-04,   2.13534000e-04,   2.44483760e-04,
              1.84273240e-04,   2.21087800e-04,   2.49137960e-04],
           [  2.41720460e-04,   2.64764100e-04,   2.82959590e-04,
              2.47934560e-04,   2.69466550e-04,   2.86150390e-04],
           [  3.07870120e-04,   3.16988600e-04,   3.23006480e-04,
              3.09844920e-04,   3.18839700e-04,   3.24733880e-04]])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater. UNESCO Tech. Pap. in Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39

    Bryden, H. 1973.
    "New Polynomials for thermal expansion, adiabatic temperature gradient
    and potential temperature of sea water."
    DEEP-SEA RES., 1973, Vol20,401-408.

    Authors
    -------
    Phil Morgan, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-16. Filipe Fernandes, Reformulated docstring.
    """

    T68 = T68conv * t

    a0 =  3.5803E-5
    a1 =  8.5258E-6
    a2 = -6.836E-8
    a3 =  6.6228E-10

    b0 =  1.8932E-6
    b1 = -4.2393E-8

    c0 =  1.8741E-8
    c1 = -6.7795E-10
    c2 =  8.733E-12
    c3 = -5.4481E-14

    d0 = -1.1351E-10
    d1 =  2.7759E-12

    e0 = -4.6206E-13
    e1 =  1.8676E-14
    e2 = -2.1687E-16

    adtg = a0 + ( a1 + ( a2 + a3 * T68 ) * T68) * T68 \
            + ( b0 + b1 * T68 ) * ( s-35 ) \
            + ( ( c0 + ( c1 + ( c2 + c3 * T68 ) * T68 ) * T68 ) \
            + ( d0 + d1 * T68 ) * ( s-35 ) ) * p \
            + (  e0 + (e1 + e2 * T68) * T68 )*p*p

    return adtg

def alpha(s, t, p, pt=False):
    """
    Calculate the thermal expansion coefficient.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature or potential temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"
    pt : bool
         True if temperature is potential, default is False

    Returns
    -------
    alpha : array_like
            thermal expansion coeff (:math:`\\alpha`) [:math:`^\\circ` C :sup:`-1`]

    See Also
    --------
    beta, aonb

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested

    Examples
    --------
    Data from McDougall 1987

    >>> import seawater as sw
    >>> s, t, p = 40., 10., 4000.
    >>> sw.alpha(s, t, p, pt=True)
    0.00025061316481624323

    Reference
    ---------
    McDougall, T.J. 1987.  "Neutral Surfaces"
    Journal of Physical Oceanography vol 17 pages 1950-1964,

    Authors
    -------
    N.L. Bindoff  1993, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    93-04-22. Phil Morgan,  Help display modified to suit library
    93-04-23. Phil Morgan,  Input argument checking
    94-10-15. Phil Morgan,  Pass S,T,P and keyword for 'ptmp'
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-16. Filipe Fernandes, Reformulated docstring.
    """

    alpha = aonb(s, t, p, pt) * beta(s, t, p, pt)
    return alpha

def aonb(s, t, p, pt=False):
    """
    Calculate :math:`\\alpha/\\beta`.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature or potential temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"
    pt : bool
         True if temperature is potential, default is False

    Returns
    -------
    aonb : array_like
           :math:`\\alpha/\\beta` [psu :math:`^\\circ` C :sup:`-1`]

    See Also
    --------
    alpha, beta

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested
    TODO: Test pt=False

    Examples
    --------
    Data from McDouogall 1987

    >>> import seawater as sw
    >>> s, t, p = 40.0, 10.0, 4000
    >>> sw.aonb(s, t, p, pt=True)
    0.347650567047807

    Reference
    ---------
    McDougall, T.J. 1987. "Neutral Surfaces"
    Journal of Physical Oceanography vol 17 pages 1950-1964,

    Authors
    -------
    N.L. Bindoff  1993, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    93-04-22. Phil Morgan,  Help display modified to suit library
    93-04-23. Phil Morgan,  Input argument checking
    94-10-15. Phil Morgan,  Pass S,T,P and keyword for 'ptmp'
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-16. Filipe Fernandes, Reformulated docstring.
    """

    # Ensure we use ptmp in calculations
    if pt==False:
        t = ptmp(s, t, p, 0) # now we have ptmp

    p = np.float32(p)
    t = t * T68conv

    c1  = np.array([-0.255019e-7, 0.298357e-5, -0.203814e-3, \
                    0.170907e-1, 0.665157e-1])
    c2  = np.array([-0.846960e-4, 0.378110e-2])
    c2a = np.array([-0.251520e-11, -0.164759e-6, 0.0])
    c3  = -0.678662e-5
    c4  = np.array([0.791325e-8, -0.933746e-6, 0.380374e-4])
    c5  =  0.512857e-12
    c6  = -0.302285e-13
    # Now calculate the thermal expansion saline contraction ratio adb
    sm35  = s - 35.0
    aonb  = np.polyval(c1, t) + sm35 * ( np.polyval(c2, t) \
            + np.polyval(c2a, p) ) \
            + sm35**2 * c3 + p * np.polyval(c4, t) \
            + c5 * (p**2) * (t**2) + c6 * p**3

    return aonb

def beta(s, t, p, pt=False):
    """
    Calculate the saline contraction coefficient :math:`\\beta` as defined by T.J. McDougall.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature or potential temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"
    pt : bool
         True if temperature is potential, default is False

    Returns
    -------
    beta : array_like
           saline Contraction Coefficient  [psu :sup:`-1`]

    Examples
    --------
    Data from McDouogall 1987

    >>> import seawater as sw
    >>> s, t, p = 40.0, 10.0, 4000
    >>> sw.beta(s, t, p, pt=True)
    0.00072087661741618932

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested
    TODO: Test pt=False for alpha, beta and aonb

    Authors
    -------
    N.L. Bindoff  1993, Lindsay Pender (Lindsay.pender@csiro.au)

    Reference
    ---------
    McDougall, T.J. 1987. "Neutral Surfaces"
    Journal of Physical Oceanography vol 17 pages 1950-1964,

    Modifications
    -------------
    93-04-22. Phil Morgan,  Help display modified to suit library
    93-04-23. Phil Morgan,  Input argument checking
    94-10-15. Phil Morgan,  Pass S,T,P and keyword for 'ptmp'
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-16. Filipe Fernandes, Reformulated docstring.
    """

    # Ensure we use ptmp in calculations
    if pt==False:
        t = ptmp(s, t, p, 0) # now we have ptmp

    p = np.float32(p)
    t = t * T68conv

    c1 = np.array([-0.415613e-9, 0.555579e-7, -0.301985e-5, 0.785567e-3])
    c2 = np.array([0.788212e-8, -0.356603e-6])
    c3 = np.array([-0.602281e-15, 0.408195e-10, 0.0])
    c4 = 0.515032e-8
    c5 = np.array([-0.213127e-11, 0.192867e-9, -0.121555e-7])
    c6 = np.array([-0.175379e-14, 0.176621e-12])
    c7 = 0.121551e-17

    # Now calculate the thermal expansion saline contraction ratio adb
    sm35  = s - 35
    beta  = np.polyval(c1, t) + sm35 * (np.polyval(c2, t) + \
            np.polyval(c3, p) ) + c4 * (sm35**2) + \
            p * np.polyval(c5, t) + (p**2) * np.polyval(c6, t) \
            + c7 * (p**3)

    return beta

def bfrq(s, t, p, lat=None):
    """
    Calculates Brunt-Vaisala Frequency squared (N :sup:`2`) at the mid depths from the equation:

    .. math::
        N^{2} = \\frac{-g}{\\sigma_{\\theta}} \\frac{d\\sigma_{\\theta}}{dz}

    Also calculates Potential Vorticity from:

    .. math::
        q=f \\frac{N^2}{g}

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature or potential temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"

    lat : number or array_like, optional
          latitude in decimal degrees north [-90..+90].
          Will grav instead of the default g = 9.8 m :sup:`2` s :sup:`-1`) and d(z) instead of d(p)

    Returns
    -------
    n2 : array_like
           Brünt-Väisälä Frequency squared (M-1xN)  [s :sup:`-2`]
    q : array_like
           planetary potential vorticity (M-1xN)  [ m s :sup:`-1`]
    p_ave : array_like
            mid pressure between P grid (M-1xN) [db]

    See Also
    --------
    pden, dens

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested
    The value of gravity is a global constant

    Examples
    --------
    TODO: add a test here
    >>> import numpy as np
    >>> import seawater as sw
    >>> s = np.array([[0, 0, 0], [15, 15, 15], [30, 30, 30],[35,35,35]])
    >>> t = np.repeat(15, s.size).reshape(s.shape)
    >>> p = np.array([0, 250, 500, 1000])
    >>> lat = np.array([30,32,35])
    >>> sw.bfrq(s, t, p, lat)[0]
    array([[  4.51543648e-04,   4.51690708e-04,   4.51920753e-04],
           [  4.45598092e-04,   4.45743207e-04,   4.45970207e-04],
           [  7.40996788e-05,   7.41238078e-05,   7.41615525e-05]])

    References
    ----------
    A.E. Gill 1982. p.54  eqn 3.7.15
    "Atmosphere-Ocean Dynamics"
    Academic Press: New York.  ISBN: 0-12-283522-0

    Jackett, D.R. and McDougall, T.J. 1994.
    Minimal adjustment of hydrographic properties to achieve static
    stability. Aubmitted J.Atmos.Ocean.Tech.

    Authors
    -------
    Phil Morgan 93-06-24, Lindsay Pender (Lindsay.Pender@csiro.au)
    Greg Johnson (gjohnson@pmel.noaa.gov)
    added potential vorticity calculation

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    06-04-19. Lindsay Pender, Corrected sign of PV.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-17. Filipe Fernandes, Reformulated docstring.
    """

    # if pressure is a vector make it a array of the same size as t/s
    if p.ndim == 1:
        p = np.repeat(p[np.newaxis,:], t.shape[1], axis=1).reshape(t.shape)

    if lat is not None:
        z = depth(p, lat)
        g = grav(lat, -z) # note that grav expects height as argument
        f = cor(lat)
    else: # TODO: test this if logic
        z = p
        f = NaN

    m,n   = p.shape # TODO: check where depth increases to automagically find which dimension to operate
    iup   = np.arange(0, m-1)
    ilo   = np.arange(1, m)

    p_ave    = ( p[iup,:] + p[ilo,:] )/2.
    pden_up  = pden( s[iup,:], t[iup,:], p[iup,:], p_ave )
    pden_lo  = pden( s[ilo,:], t[ilo,:], p[ilo,:], p_ave )
    mid_pden = ( pden_up + pden_lo )/2
    dif_pden = pden_up - pden_lo
    mid_g    = ( g[iup,:] + g[ilo,:] )/2
    dif_z    = np.diff(z, axis=0) # TODO: check where depth increases to automagically find which dimension to operate
    n2       = -mid_g * dif_pden / ( dif_z * mid_pden )
    q        = -f * dif_pden / ( dif_z * mid_pden )

    return n2, q, p_ave

def depth(p, lat):
    """
    Calculates depth in metres from pressure in dbars.

    Parameters
    ----------
    p : array_like
        pressure [db].
    lat : number or array_like
          latitude in decimal degress north [-90..+90]. The shape can be "broadcasted"

    Returns
    -------
    z : array_like
        depth [metres]

    Examples
    --------
    UNESCO 1983 data p30

    >>> import seawater as sw
    >>> lat = np.array([0, 30, 45, 90])
    >>> p   = np.array([[  500,   500,   500,  500], [ 5000,  5000,  5000, 5000], [10000, 10000, 10000, 10000]])
    >>> sw.depth(p, lat)
    array([[  496.65299239,   495.99772917,   495.3427354 ,   494.03357499],
           [ 4915.04099112,  4908.55954332,  4902.08075214,  4889.13132561],
           [ 9725.47087508,  9712.6530721 ,  9699.84050403,  9674.23144056]])

    Notes
    -----
    original seawater name is dpth

    References
    ----------
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Authors
    -------
    Phil Morgan 92-04-06  (morgan@ml.csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    # Eqn 25, p26.  UNESCO 1983.
    c1 =  9.72659
    c2 = -2.2512E-5
    c3 =  2.279E-10
    c4 = -1.82E-15

    gam_dash = 2.184e-6

    lat = abs(lat)
    X   = np.sin( lat * DEG2RAD ) # convert to radians
    X   = X * X

    bot_line = 9.780318 * ( 1.0 + ( 5.2788E-3 + 2.36E-5 * X ) * X ) + \
               gam_dash * 0.5 * p
    top_line = ( ( ( c4 * p + c3 ) * p + c2 ) * p + c1 ) * p
    depthm   = top_line / bot_line
    return depthm

def grav(lat, z=0):
    """
    Calculates acceleration due to gravity as function of latitude.

    Parameters
    ----------
    lat : array_like
         latitude in decimal degrees north [-90..+90].

    z : number or array_like. Default z = 0
        height in metres (+ve above sea surface, -ve below).

    Returns
    -------
    g : array_like
        gravity [m s :sup:`2`]

    See Also
    --------
    bfrq

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater as sw
    >>> lat = 45.
    >>> sw.grav(lat, z=0)
    9.8061898752053995

    References
    ----------
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    A.E. Gill 1982. p.597
    "Atmosphere-Ocean Dynamics"
    Academic Press: New York.  ISBN: 0-12-283522-0

    Authors
    -------
    Phil Morgan 93-04-20  (morgan@ml.csiro.au)

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    # Eqn p27.  UNESCO 1983.
    a       = 6371000. # mean radius of earth  A.E.Gill
    lat     = abs(lat)
    X       = np.sin( lat * DEG2RAD )  # convert to radians
    sin2    = X * X
    grav    = 9.780318 * ( 1.0 + ( 5.2788E-3 + 2.36E-5 * sin2 ) * sin2 )
    grav    = grav / ( ( 1 + z/a )**2 )    # from A.E.Gill p.597
    return  grav

def cor(lat):
    """
    Calculates the Coriolis factor :math:`f` defined by:


    .. math::
        f = 2 \\Omega \\sin(lat)

    where:


    .. math::
        \\Omega = \\frac{2 \\pi}{\\textrm{sidereal day}} = 7.292e^{-5} \\textrm{ radians sec}^{-1}


    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90].

    Returns
    -------
    f : array_like
        Coriolis factor [s :sup:`-1`]

    See Also
    --------
    TODO: inertial period

    Notes
    -----
    1 sidereal day = 23.9344696 hours

    Examples
    --------
    >>> import seawater as sw
    >>> sw.cor(45)
    0.00010312445296824608

    Notes
    -----
    The value of Omega is a global constant

    References
    ----------
    S. Pond & G.Pickard  2nd Edition 1986
    Introductory Dynamical Oceanogrpahy
    Pergamon Press Sydney.  ISBN 0-08-028728-X

    A.E. Gill 1982. p.597
    "Atmosphere-Ocean Dynamics"
    Academic Press: New York.  ISBN: 0-12-283522-0

    Authors
    -------
    Phil Morgan 93-04-20  (morgan@ml.csiro.au)

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    # Eqn p27.  UNESCO 1983.
    f = 2 * OMEGA * np.sin( lat * DEG2RAD )
    return f

def cndr(s, t, p):
    """
    Calculates conductivity ratio from S, T, P.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"

    Returns
    -------
    cndr : array_like
           conductivity ratio. R = C(s,t,p)/C(35,15(IPTS-68),0) [no units]

    See Also
    --------
    salds, sals, salrt

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested

    Examples
    --------
    Data from UNESCO 1983 p9

    >>> import seawater as sw
    >>> t    = np.array([0, 10, 0, 10, 10, 30]) / T68conv
    >>> p    = np.array([0, 0, 1000, 1000, 0, 0])
    >>> s    = np.array([25, 25, 25, 25, 40, 40])
    >>> sw.cndr(s, t, p)
    array([ 0.49800825,  0.65499015,  0.50624434,  0.66297496,  1.00007311,
            1.52996697])

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Authors
    -------
    Phil Morgan 93-04-21, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    T68 = t * T68conv
    DT  = T68 - 15

    Rx  = (s/35.0)**0.5 # first guess at Rx = sqrt(Rt)
    SInc  = sals(Rx*Rx, t)

    # DO A NEWTON-RAPHSON ITERATION FOR INVERSE INTERPOLATION OF Rt FROM S.
    for n in range(100): # TODO: is a 100 important? changed to 10?
        Rx   = Rx + (s - SInc) / salds(Rx, DT)
        SInc = sals(Rx*Rx, t)
        DELS = abs(SInc - s)
        DELS = np.array(DELS) # TODO: "any" is an array method
        if (DELS.any() < 1.0E-4):
            break
    """ original matlab TODO: implement this above and check the difference
    for i = 1:ms
    for j = 1:ns
        ---------------------------------------------------------------------
        DO A NEWTON-RAPHSON ITERATION FOR INVERSE INTERPOLATION OF Rt FROM S.
        ---------------------------------------------------------------------
        S_loop   = S(i,j) # S in the loop
        T_loop   = T(i,j) # T in the loop
        Rx_loop  = sqrt(S_loop/35.0) #first guess at Rx = sqrt(Rt)
        SInc     = sals(Rx_loop*Rx_loop, T_loop) # S INCrement (guess) from Rx
        iloop    = 0
        end_loop = 0
        while ~end_loop:
                Rx_loop = Rx_loop + (S_loop - SInc) / salds(Rx_loop, T_loop - 15)
        SInc    = sals(Rx_loop * Rx_loop, T_loop)
        iloop   = iloop + 1
        dels    = abs(SInc-S_loop)
        if (dels>1.0e-4 & iloop<10) :
            end_loop = 0
        else:
            end_loop = 1

        Rx(i,j) = Rx_loop
    """
    # ONCE Rt FOUND, CORRESPONDING TO EACH (S,T) EVALUATE R
    # eqn(4) p.8 UNESCO 1983
    d1 =  3.426e-2
    d2 =  4.464e-4
    d3 =  4.215e-1
    d4 = -3.107e-3

    e1 =  2.070e-5
    e2 = -6.370e-10
    e3 =  3.989e-15

    A  = ( d3 + d4 * T68 )
    B  = 1 + d1 * T68 + d2 * T68**2
    C  = p * ( e1 + e2 * p + e3 * p**2 )

    # eqn(6) p.9 UNESCO 1983.
    Rt    = Rx * Rx
    rt    = salrt(t)
    Rtrt  = rt * Rt
    D     = B - A * rt * Rt
    E     = rt * Rt * A * ( B + C )
    r     = ( abs( D**2 + 4 * E ) )**0.5 - D
    r     = 0.5 * r/A
    return r

def sals(rt, t):
    """
    Salinity of sea water as a function of Rt and T.
    UNESCO 1983 polynomial.

    Parameters
    ----------
    rt : array_like
         :math:`rt(s,t) = \\frac{C(s,t,0)}{C(35, t(\\textrm{IPTS-68}), 0)}`
    t : array_like
        temperature [:math:`^\\circ` C (ITS-90)]

    Returns
    -------
    s : array_like
        salinity [psu (PSS-78)]

    See Also
    --------
    salt

    Notes
    -----
    TODO

    Examples
    --------
    Data from UNESCO 1983 p9

    >>> import seawater as sw
    >>> t    = np.array([15, 20, 5]) / T68conv
    >>> rt   = np.array([  1, 1.0568875, 0.81705885])
    >>> sw.sals(rt,t)
    array([ 35.        ,  37.24562718,  27.99534701])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Authors
    -------
    Phil Morgan 93-04-17, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    # eqn (1) & (2) p6,7 unesco
    del_T68 = t * T68conv - 15

    a0 =  0.0080
    a1 = -0.1692
    a2 = 25.3851
    a3 = 14.0941
    a4 = -7.0261
    a5 =  2.7081

    b0 =  0.0005
    b1 = -0.0056
    b2 = -0.0066
    b3 = -0.0375
    b4 =  0.0636
    b5 = -0.0144

    k  =  0.0162

    Rtx   = (rt)**0.5
    del_S = ( del_T68 / ( 1 + k * del_T68 ) ) * \
            ( b0 + ( b1 + ( b2 + ( b3 + ( b4 + b5 * Rtx ) \
            * Rtx ) * Rtx ) * Rtx ) * Rtx)

    s = a0 + ( a1 + ( a2 + ( a3 + ( a4 + a5 * Rtx) * \
               Rtx) * Rtx ) * Rtx ) * Rtx

    s = s + del_S

    return s

def salds(rtx, delt):
    """
    Calculates Salinity differential dS/d(sqrt(Rt)) at constant T.
    UNESCO 1983 polynomial.

    Parameters
    ----------
    rtx : array_like
          :math:`\\sqrt{rt}`
    delt : array_like
           t-15 [:math:`^\\circ` C (IPTS-68)]

    Returns
    -------
    ds : array_like
         S differential :math:`\\frac{dS}{d(\\sqrt{(Rt)})} at constant T.

    See Also
    --------
    cndr, salt

    Notes
    -----
    TODO

    Examples
    --------
    Data from UNESCO 1983 p9

    >>> import seawater as sw
    >>> delt = np.array([15, 20, 5]) / T68conv  - 15
    >>> rtx  = np.array([  1, 1.0568875, 0.81705885])**0.5
    >>> sw.salds(rtx, delt)
    array([ 78.31921607,  81.5689307 ,  68.19023687])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Authors
    -------
    Phil Morgan 93-04-21, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    a0 =  0.0080
    a1 = -0.1692
    a2 = 25.3851
    a3 = 14.0941
    a4 = -7.0261
    a5 =  2.7081

    b0 =  0.0005
    b1 = -0.0056
    b2 = -0.0066
    b3 = -0.0375
    b4 =  0.0636
    b5 = -0.0144

    k  =  0.0162

    ds =  a1 + ( 2 * a2 + ( 3 * a3 + ( 4 * a4 + 5 * a5 * rtx ) * rtx ) \
          * rtx ) * rtx + ( delt / ( 1 + k * delt ) ) * \
          ( b1 + ( 2 * b2 + ( 3 * b3 + ( 4 * b4 + 5 * b5 * rtx ) \
          * rtx ) * rtx) * rtx)

    return ds

def salrt(t):
    """
    Equation for rt used in calculating salinity. UNESCO 1983 polynomial.

    .. math::
        rt(t) = \\frac{C(35,t,0)}{C(35,15(\textrm{IPTS-68}), 0)}


    Parameters
    ----------
      t : array_like
          temperature [:math:`^\\circ` C (ITS-90)]

    Returns
    -------
    rt : array_like
    conductivity ratio  [no units]

    See Also
    --------
    salt

    Notes
    -----
    TODO

    Examples
    --------
    Data from UNESCO 1983 p9

    >>> import seawater as sw
    >>> t = np.array([15, 20, 5]) / T68conv
    >>> sw.salrt(t)
    array([ 1.        ,  1.11649272,  0.77956585])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Authors
    -------
    Phil Morgan 93-04-17, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    #Eqn (3) p.7 UNESCO.
    T68 = t * T68conv

    c0 =  0.6766097
    c1 =  2.00564e-2
    c2 =  1.104259e-4
    c3 = -6.9698e-7
    c4 =  1.0031e-9

    rt = c0 + ( c1 + ( c2 + ( c3 + c4 * T68) * T68) * T68 ) * T68
    return rt

def salt(r, t, p):
    """
    Calculates Salinity from conductivity ratio. UNESCO 1983 polynomial.

    Parameters
    ----------
    r : array_like
            conductivity ratio :math:`R = \\frac{C(S,T,P)}{C(35,15(IPTS-68),0)}` [no units]
    t : array_like
        temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]

    Returns
    -------
    s : array_like
        salinity [psu (PSS-78)]

    See Also
    --------
    sals, salrt, salrp

    Notes
    -----
    TODO

    Examples
    --------
    Data from UNESCO 1983 p9

    >>> import seawater as sw
    >>> r = np.array([1, 1.2, 0.65])
    >>> t = np.array([15, 20, 5]) / T68conv
    >>> p = np.array([0, 2000, 1500])
    >>> sw.salt(r, t, p)
    array([ 34.99999992,  37.24562765,  27.99534693])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Authors
    -------
    Phil Morgan 93-04-17, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    rt = salrt(t)
    rp = salrp(r, t, p )
    rt = r / ( rp * rt )
    s  = sals(rt, t)

    return s

def salrp(r, t, p):
    """
    Equation for Rp used in calculating salinity. UNESCO 1983 polynomial.

    .. math::
        Rp(S,T,P) = \\frac{C(S,T,P)}{C(S,T,0)}


    Parameters
    ----------
    r : array_like
        conductivity ratio :math:`R = \\frac{C(S,T,P)}{C(35,15(IPTS-68),0)}` [no units]
    t : array_like
        temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]

    Returns
    -------
    rp : array_like
        conductivity ratio :math:`Rp(S,T,P) = \\frac{C(S,T,P)}{C(S,T,0)}` [no units]

    See Also
    --------
    salt

    Notes
    -----
    TODO

    Examples
    --------

    >>> import seawater as sw
    >>> r = np.array([1, 1.2, 0.65])
    >>> t = np.array([15, 20, 5]) / T68conv
    >>> p = np.array([0, 2000, 1500])
    >>> sw.salrp(r, t, p)
    array([ 1.        ,  1.01694294,  1.02048638])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Authors
    -------
    Phil Morgan 93-04-17, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    # eqn (4) p.8 unesco.
    T68 = t * T68conv

    d1 =  3.426e-2
    d2 =  4.464e-4
    d3 =  4.215e-1
    d4 = -3.107e-3

    e1 =  2.070e-5
    e2 = -6.370e-10
    e3 =  3.989e-15

    rp = 1 + ( p * ( e1 + e2 * p + e3 * p**2 ) ) \
         / ( 1 + d1 * T68 + d2 * T68**2 + ( d3 + d4 * T68 ) * r)

    return rp

def fp(s, p):
    """
    Freezing point of Sea Water using UNESCO 1983 polynomial.

    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    p : array_like
        pressure [db]

    Returns
    -------
    fp : array_like
        freezing point temperature [:math:`^\\circ` C (ITS-90)]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    UNESCO DATA p.30

    >>> import seawater as sw
    >>> s = np.array([[5, 10, 15, 20, 25, 30, 35, 40], [5, 10, 15, 20, 25, 30, 35, 40]])
    >>> p = np.array([[ 0, 0, 0, 0, 0, 0, 0, 0], [500, 500, 500, 500, 500, 500, 500, 500]])
    >>> sw.fp(s, p)
    array([[-0.27369757, -0.54232831, -0.81142026, -1.0829461 , -1.35804594,
            -1.63748903, -1.9218401 , -2.2115367 ],
           [-0.65010724, -0.91873798, -1.18782992, -1.45935577, -1.73445561,
            -2.01389869, -2.29824976, -2.58794636]])

    References
    ----------
    Fofonff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Authors
    -------
    Phil Morgan 93-04-20, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """
    #TODO: P = P/10; # to convert db to Bar as used in UNESCO routines (was commented in the original)

    # eqn  p.29
    a0 = -0.0575
    a1 =  1.710523e-3
    a2 = -2.154996e-4
    b  = -7.53e-4

    fp = ( a0 * s + a1 * s * (s)**0.5 + a2 * s**2 + b * p ) / T68conv

    return fp

def svel(s, t, p):
    """
    Sound Velocity in sea water using UNESCO 1983 polynomial.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"

    Returns
    -------
    svel : array_like
           sound velocity  [m/s]


    See Also
    --------
    TODO

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested
    TODO: Add equation to docstring

    Examples
    --------
    Data from Pond and Pickard Intro. Dynamical Oceanography 2nd ed. 1986

    >>> import numpy as np
    >>> import seawater as sw
    >>> from seawater import T68conv
    >>> t = np.array([[  0,  0,  0,  0,  0,  0], [ 10, 10, 10, 10, 10, 10], [ 20, 20, 20, 20, 20, 20], [ 30, 30, 30, 30, 30, 30], [ 40, 40, 40, 40, 40, 40]]) / T68conv
    >>> s = np.array([[ 25, 25, 25, 35, 35, 35], [ 25, 25, 25, 35, 35, 35], [ 25, 25, 25, 35, 35, 35], [ 25, 25, 25, 35, 35, 35], [ 25, 25, 25, 35, 35, 35]])
    >>> p = np.array([[ 0, 5000, 10000, 0, 5000, 10000], [ 0, 5000, 10000, 0, 5000, 10000], [ 0, 5000, 10000, 0, 5000, 10000], [ 0, 5000, 10000, 0, 5000, 10000], [ 0, 5000, 10000, 0, 5000, 10000]])
    >>> sw.svel(s, t, p)
    array([[ 1435.789875  ,  1520.358725  ,  1610.4074    ,  1449.13882813,
             1533.96863705,  1623.15007097],
           [ 1477.68316464,  1561.30635914,  1647.39267114,  1489.82233602,
             1573.40946928,  1658.99115504],
           [ 1510.31388348,  1593.59671798,  1676.80967748,  1521.4619731 ,
             1604.4762822 ,  1687.18305631],
           [ 1535.21434752,  1618.95631952,  1700.60547902,  1545.59485539,
             1628.97322783,  1710.06294277],
           [ 1553.44506636,  1638.02522336,  1719.15088536,  1563.20925247,
             1647.29949576,  1727.83176404]])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Authors
    -------
    Phil Morgan 93-04-20, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    # UNESCO 1983. eqn.33  p.46
    p = p/10  # convert db to bars as used in UNESCO routines
    T68 = t * T68conv

    # eqn 34 p.46
    c00 = 1402.388
    c01 =    5.03711
    c02 =   -5.80852e-2
    c03 =    3.3420e-4
    c04 =   -1.47800e-6
    c05 =    3.1464e-9

    c10 =  0.153563
    c11 =  6.8982e-4
    c12 = -8.1788e-6
    c13 =  1.3621e-7
    c14 = -6.1185e-10

    c20 =  3.1260e-5
    c21 = -1.7107e-6
    c22 =  2.5974e-8
    c23 = -2.5335e-10
    c24 =  1.0405e-12

    c30 = -9.7729e-9
    c31 =  3.8504e-10
    c32 = -2.3643e-12

    Cw  = ( ( ( ( c32 * T68 + c31 ) * T68 + c30 ) * p + \
          ( ( ( ( c24 * T68 + c23 ) * T68 + c22 ) * T68 + c21 ) * T68 + c20 ) ) * p + \
          ( ( ( ( c14 * T68 + c13 ) * T68 + c12 ) * T68 + c11 ) * T68 + c10 ) ) * p + \
          ( ( ( ( c05 * T68 + c04 ) * T68 + c03 ) * T68 + c02 ) * T68 + c01 ) * T68 + c00

    # eqn 35. p.47
    a00 =  1.389
    a01 = -1.262e-2
    a02 =  7.164e-5
    a03 =  2.006e-6
    a04 = -3.21e-8

    a10 =  9.4742e-5
    a11 = -1.2580e-5
    a12 = -6.4885e-8
    a13 =  1.0507e-8
    a14 = -2.0122e-10

    a20 = -3.9064e-7
    a21 =  9.1041e-9
    a22 = -1.6002e-10
    a23 =  7.988e-12

    a30 =  1.100e-10
    a31 =  6.649e-12
    a32 = -3.389e-13

    A = ( ( ( ( a32 * T68 + a31 ) * T68 + a30 ) * p + \
        ( ( ( a23 * T68 + a22 ) * T68 + a21 ) * T68 + a20 ) ) * p + \
        ( ( ( ( a14 * T68 + a13 ) * T68 + a12 ) * T68 + a11 ) * T68 + a10 ) ) * p + \
        ( ( ( a04 * T68 + a03 ) * T68 + a02 ) * T68 + a01 ) * T68 + a00

    # eqn 36 p.47
    b00 = -1.922e-2
    b01 = -4.42e-5
    b10 =  7.3637e-5
    b11 =  1.7945e-7

    B = b00 + b01 * T68 + ( b10 + b11 * T68) * p

    # eqn 37 p.47
    d00 =  1.727e-3
    d10 = -7.9836e-6

    D = d00 + d10 * p

    # eqn 33 p.46
    svel = Cw + A * s + B * s * (s)**0.5 + D * s**2

    return svel

def pres(depth, lat):
    """
    Calculates pressure in dbars from depth in meters.

    Parameters
    ----------
    depth : array_like
            depth [metres]
    lat : array_like
          latitude in decimal degress north [-90..+90]
          The shape can be "broadcasted"

    Returns
    -------
    p : array_like
           pressure [db]

    See Also
    --------
    TODO: pressure from depth


    Notes
    -----
    TODO: lat broadcast feature need to be tested

    Examples
    --------
    >>> import seawater as sw
    >>> depth, lat = 7321.45, 30
    >>> sw.pres(depth,lat)
    7500.0065130118019

    References
    ----------
    Saunders, P.M. 1981
    "Practical conversion of Pressure to Depth"
    Journal of Physical Oceanography, 11, 573-574

    Authors
    -------
    Phil Morgan 93-06-25  (morgan@ml.csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    X       = np.sin( abs(lat) * DEG2RAD )  # convert to radians
    C1      = 5.92E-3 + X**2 * 5.25E-3
    pres    = ( ( 1 - C1 ) - ( ( ( 1 - C1 )**2 ) - ( 8.84E-6 * depth ) )**0.5 ) / 4.42E-6
    return pres

def dist(lon, lat, units='km'):
    """
    Calculate distance between two positions on globe using the "Plane
    Sailing" method.  Also uses simple geometry to calculate the bearing of
    the path between position pairs.

    Parameters
    ----------
    lon : array_like
          decimal degrees (+ve E, -ve W) [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [- 90.. +90]
    units : string, optional
            default kilometers

    Returns
    -------
    dist : array_like
           distance between positions in units
    phaseangle : array_like
                 angle of line between stations with x axis (East). Range of values are -180..+180. (E=0, N=90, S=-90)

    See Also
    --------
    TODO

    Notes
    -----
    Usually used to creat a distance vector to plot hydrographic data. However, pay attention to the phaseangle to aviod apples and oranges!

    Also not that the input order for the matlab version is lat,lon (alphabetic order), while this version is lon,lat (geometric order).

    Examples
    --------
    >>> import numpy as np
    >>> import seawater as sw
    >>> lon = np.array([35, 35])
    >>> lat = np.array([41, 40])
    >>> sw.dist(lon, lat)
    (array([ 111.12]), array([-90.]))

    Create a distance vector

    >>> lon = np.arange(30,40,1)
    >>> lat = np.array(35)
    >>> np.cumsum(np.append(0, sw.dist(lon, lat, units='km')[0]))
    array([   0.        ,   91.02417516,  182.04835032,  273.07252548,
            364.09670065,  455.12087581,  546.14505097,  637.16922613,
            728.19340129,  819.21757645])

    References
    ----------
    The PLANE SAILING method as described in "CELESTIAL NAVIGATION" 1989 by Dr. P. Gormley. The Australian Antartic Division.

    Authors
    -------
    Phil Morgan and Steve Rintoul 92-02-10

    Modifications
    -------------
    99-06-25. Lindsay Pender, Function name change from distance to sw_dist.
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    if lat.size == 1:
        lat = np.repeat(lat, lon.size)
    elif lon.size == 1:
        lon = np.repeat(lon, lat.size)

    npositions = max(lon.shape)

    ind = np.arange( 0, npositions-1, 1) # index to first of position pairs

    dlon = np.diff(lon, axis=0)
    if any( abs(dlon) > 180 ):
        flag = abs(dlon) > 180
        dlon[flag] = -np.sign( dlon[flag] ) * ( 360 - abs( dlon[flag] ) )

    latrad = abs( lat * DEG2RAD )
    dep    = np.cos( ( latrad [ind+1] + latrad[ind] ) / 2 ) * dlon
    dlat   = np.diff( lat, axis=0 )
    dist   = DEG2NM * ( dlat**2 + dep**2 )**0.5

    if units == 'km':
        dist = dist * NM2KM

    # Calcualte angle to x axis
    RAD2DEG     = 1/DEG2RAD
    phaseangle  = np.angle( dep + dlat * 1j ) * RAD2DEG
    return dist, phaseangle

def satAr(s, t):
    """
    Solubility (satuaration) of Argon (Ar) in sea water.

    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [:math:`^\\circ` C (ITS-90)]

    Returns
    -------
    satAr : array_like
            solubility of Ar [ml l :sup:`-1`]

    See Also
    --------
    satO2, satN2

    Notes
    -----
    TODO

    Examples
    --------
    Data from Weiss 1970

    >>> import numpy
    >>> import seawater as sw
    >>> from seawater import T68conv
    >>> t = np.array([[ -1, -1], [ 10, 10], [ 20, 20], [ 40, 40]]) / T68conv
    >>> s = np.array([[ 20, 40], [ 20, 40], [ 20, 40], [ 20, 40]])
    >>> sw.satAr(s, t)
    array([[ 0.4455784 ,  0.38766011],
           [ 0.33970659,  0.29887756],
           [ 0.27660227,  0.24566428],
           [ 0.19861429,  0.17937698]])

    References
    ----------
    Weiss, R. F. 1970
    "The solubility of nitrogen, oxygen and argon in water and seawater."
    Deap-Sea Research., 1970, Vol 17, pp721-735.

    Authors
    -------
    Phil Morgan 97-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    # convert T to Kelvin
    t = Kelvin + t * T68conv

    # constants for Eqn (4) of Weiss 1970
    a1 = -173.5146
    a2 =  245.4510
    a3 =  141.8222
    a4 =  -21.8020
    b1 =   -0.034474
    b2 =    0.014934
    b3 =   -0.0017729

    # Eqn (4) of Weiss 1970
    lnC = a1 + a2 * ( 100/t ) + a3 * np.log( t/100 ) + a4 * ( t/100 ) + \
          s * ( b1 + b2 * ( t/100 ) + b3 * ( ( t/100 )**2) )

    c = np.exp(lnC)

    return c

def satN2(s, t):
    """
    Solubility (satuaration) of Nitrogen (N2) in sea water.

    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [:math:`^\\circ` C (ITS-90)]

    Returns
    -------
    satN2 : array_like
            solubility of N2  [ml l :sup:`-1`]

    See Also
    --------
    satO2, satAr

    Notes
    -----
    TODO

    Examples
    --------
    Data from Weiss 1970

    >>> import numpy
    >>> import seawater as sw
    >>> from seawater import T68conv
    >>> t = np.array([[ -1, -1], [ 10, 10], [ 20, 20], [ 40, 40]]) / T68conv
    >>> s = np.array([[ 20, 40], [ 20, 40], [ 20, 40], [ 20, 40]])
    >>> sw.satN2(s, t)
    array([[ 16.27952432,  14.00784526],
           [ 12.64036196,  11.01277257],
           [ 10.46892822,   9.21126859],
           [  7.78163876,   6.95395099]])


    References
    ----------
    Weiss, R. F. 1970
    "The solubility of nitrogen, oxygen and argon in water and seawater."
    Deap-Sea Research., 1970, Vol 17, pp721-735.

    Authors
    -------
    Phil Morgan 97-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    # convert T to Kelvin
    t = Kelvin + t * T68conv

    # constants for Eqn (4) of Weiss 1970
    a1 = -172.4965
    a2 =  248.4262
    a3 =  143.0738
    a4 =  -21.7120
    b1 =   -0.049781
    b2 =    0.025018
    b3 =   -0.0034861

    # Eqn (4) of Weiss 1970
    lnC = a1 + a2 * ( 100/t ) + a3 * np.log( t/100 ) + a4 * ( t/100 ) + \
          s * ( b1 + b2 * ( t/100 ) + b3 * ( ( t/100 )**2 ) )

    c = np.exp(lnC)
    return c

def satO2(s, t):
    """
    Solubility (satuaration) of Oxygen (O2) in sea water.

    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [:math:`^\\circ` C (ITS-68)]

    Returns
    -------
    satO2 : array_like
            solubility of O2  [ml l :sup:`-1` ]

    Notes
    -----
    TODO

    Examples
    --------
    Data from Weiss 1970

    >>> import numpy
    >>> import seawater as sw
    >>> from seawater import T68conv
    >>> t = np.array([[ -1, -1], [ 10, 10], [ 20, 20], [ 40, 40]]) / T68conv
    >>> s = np.array([[ 20, 40], [ 20, 40], [ 20, 40], [ 20, 40]])
    >>> sw.satO2(s, t)
    array([[ 9.162056  ,  7.98404249],
           [ 6.95007741,  6.12101928],
           [ 5.64401453,  5.01531004],
           [ 4.0495115 ,  3.65575811]])


    References
    ----------
    Weiss, R. F. 1970
    "The solubility of nitrogen, oxygen and argon in water and seawater."
    Deap-Sea Research., 1970, Vol 17, pp721-735.

    Authors
    -------
    Phil Morgan 97-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    # convert T to Kelvin
    t = Kelvin + t * T68conv

    # constants for Eqn (4) of Weiss 1970
    a1 = -173.4292
    a2 =  249.6339
    a3 =  143.3483
    a4 =  -21.8492
    b1 =   -0.033096
    b2 =    0.014259
    b3 =   -0.0017000

    # Eqn (4) of Weiss 1970
    lnC = a1 + a2 * ( 100/t ) + a3 * np.log( t/100 ) + a4 * ( t/100 ) + \
          s * ( b1 + b2 * ( t/100 ) + b3 * ( ( t/100 )**2 ) )

    c = np.exp(lnC)
    return c

def dens0(s, t):
    """
    Density of Sea Water at atmospheric pressure using
    UNESCO 1983 (EOS 1980) polynomial.

    Parameters
    ----------
    s(p=0) : array_like
             salinity [psu (PSS-78)]
    t(p=0) : array_like
             temperature [:math:`^\\circ` C (ITS-90)]

    Returns
    -------
    dens0(s, t) : array_like
                  density  [kg m :sup:`3`] of salt water with properties (s, t, p=0) 0 db gauge pressure

    See Also
    --------
    svan, dens, smow, seck, pden

    Notes
    -----
    Dens0 is the density as a function of salinity and temperature only used to TODO

    Examples
    --------
    Data from UNESCO Tech. Paper in Marine Sci. No. 44, p22

    >>> import numpy as np
    >>> import seawater as sw
    >>> s = np.array([0, 0, 0, 0, 35, 35, 35, 35])
    >>> t = np.array([0, 0, 30, 30, 0, 0, 30, 30]) / T68conv
    >>> sw.dens0(s, t)
    array([  999.842594  ,   999.842594  ,   995.65113374,   995.65113374,
            1028.10633141,  1028.10633141,  1021.72863949,  1021.72863949])

    References
    ----------
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Millero, F.J. and  Poisson, A.
    International one-atmosphere equation of state of seawater.
    Deep-Sea Res. 1981. Vol28A(6) pp625-629.

    Authors
    -------
    Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    T68 = t * T68conv

    # UNESCO 1983 eqn(13) p17
    b0 =  8.24493e-1
    b1 = -4.0899e-3
    b2 =  7.6438e-5
    b3 = -8.2467e-7
    b4 =  5.3875e-9

    c0 = -5.72466e-3
    c1 =  1.0227e-4
    c2 = -1.6546e-6

    d0 = 4.8314e-4
    dens0 = smow(t) + ( b0 + ( b1 + ( b2 + ( b3 + b4 * T68 ) * T68 ) * T68 ) * T68 ) * \
    s + ( c0 + ( c1 + c2 * T68 ) * T68 ) * s * (s)**0.5 + d0 * s**2
    return dens0

def smow(t):
    """
    Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980.

    Parameters
    ----------
    t : array_like
        temperature [:math:`^\\circ` C (ITS-90)]

    Returns
    -------
    dens(t) : array_like
              density  [kg m :sup:`3`]

    See Also
    --------
    svan, dens, dens0, seck, pden

    Notes
    -----
    Standard Mean Ocean Water (SMOW) is the water collected in the deep ocean used as a reference.

    Examples
    --------
    Data from UNESCO Tech. Paper in Marine Sci. No. 44, p22

    >>> import numpy as np
    >>> import seawater as sw
    >>> t = np.array([0, 0, 30, 30, 0, 0, 30, 30]) / T68conv
    >>> sw.smow(t)
    array([ 999.842594  ,  999.842594  ,  995.65113374,  995.65113374,
            999.842594  ,  999.842594  ,  995.65113374,  995.65113374])

    References
    ----------
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
    UNESCO 1983 p17  Eqn(14)

    Millero, F.J & Poisson, A.
    INternational one-atmosphere equation of state for seawater.
    Deep-Sea Research Vol28A No.6. 1981 625-629. Eqn (6)

    Authors
    -------
    Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    a0 = 999.842594
    a1 =   6.793952e-2
    a2 =  -9.095290e-3
    a3 =   1.001685e-4
    a4 =  -1.120083e-6
    a5 =   6.536332e-9

    T68  = t * T68conv
    dens = a0 + ( a1 + ( a2 + ( a3 + ( a4 + a5 * T68 ) * T68 ) * T68 ) * T68 ) * T68
    return dens

def seck(s, t, p=0):
    """
    Secant Bulk Modulus (K) of Sea Water using Equation of state 1980. UNESCO polynomial implementation.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"

    Returns
    -------
    k : array_like
        secant bulk modulus [bars]

    See Also
    --------
    dens

    Notes
    -----
    TODO

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22

    >>> import numpy as np
    >>> import seawater as sw
    >>> s = np.array([0, 0, 0, 0, 35, 35, 35, 35])
    >>> t = np.array([0, 0, 30, 30, 0, 0, 30, 30]) / T68conv
    >>> p = np.array([0, 10000, 0, 10000, 0, 10000, 0, 10000])
    >>> sw.seck(s, t, p)
    array([ 19652.21      ,  22977.2115    ,  22336.0044572 ,  25656.8196222 ,
            21582.27006823,  24991.99729129,  23924.21823158,  27318.32472464])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
    Eqn.(15) p.18

    Millero, F.J. and  Poisson, A.
    International one-atmosphere equation of state of seawater.
    Deep-Sea Res. 1981. Vol28A(6) pp625-629.

    Authors
    -------
    Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    # Compute compression terms
    p   = p/10.0 # convert from db to atmospheric pressure units
    T68 = t * T68conv

    # Pure water terms of the secant bulk modulus at atmos pressure.
    # UNESCO eqn 19 p 18
    h3 = -5.77905E-7
    h2 =  1.16092E-4
    h1 =  1.43713E-3
    h0 =  3.239908   #[-0.1194975]

    AW = h0 + ( h1 + ( h2 + h3 * T68 ) * T68 ) * T68

    k2 =  5.2787E-8
    k1 = -6.12293E-6
    k0 =  8.50935E-5  #[+3.47718E-5]

    BW = k0 + ( k1 + k2 * T68 ) * T68

    e4 =    -5.155288E-5
    e3 =     1.360477E-2
    e2 =    -2.327105
    e1 =   148.4206
    e0 = 19652.21    #[-1930.06]

    KW  = e0 + ( e1 + ( e2 + ( e3 + e4 * T68 ) * T68 ) * T68 ) * T68 # eqn 19

    # Sea water terms of secant bulk modulus at atmos. pressure
    j0 = 1.91075E-4

    i2 = -1.6078E-6
    i1 = -1.0981E-5
    i0 =  2.2838E-3

    SR = (s)**0.5

    A  = AW + ( i0 + ( i1 + i2 * T68 ) * T68 + j0 * SR ) * s

    m2 =  9.1697E-10
    m1 =  2.0816E-8
    m0 = -9.9348E-7

    B = BW + ( m0 + ( m1 + m2 * T68 ) * T68 ) * s # eqn 18

    f3 =  -6.1670E-5
    f2 =   1.09987E-2
    f1 =  -0.603459
    f0 =  54.6746

    g2 = -5.3009E-4
    g1 =  1.6483E-2
    g0 =  7.944E-2

    K0 = KW + ( f0 + ( f1 + ( f2 + f3 * T68 ) * T68 ) * T68 \
            + ( g0 + ( g1 + g2 * T68 ) * T68 ) * SR ) * s # eqn 16

    K = K0 + ( A + B * p ) * p # eqn 15
    return K

def dens(s, t, p):
    """
    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"

    Returns
    -------
    dens : array_like
           density  [kg m :sup:`3`]

    See Also
    --------
    dens0, seck

    Notes
    -----
    TODO

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22

    >>> import numpy as np
    >>> import seawater as sw
    >>> s = np.array([0, 0, 0, 0, 35, 35, 35, 35])
    >>> t = np.array([0, 0, 30, 30, 0, 0, 30, 30]) / T68conv
    >>> p = np.array([0, 10000, 0, 10000, 0, 10000, 0, 10000])
    >>> sw.dens(s, t, p)
    array([  999.842594  ,  1045.33710972,   995.65113374,  1036.03148891,
            1028.10633141,  1070.95838408,  1021.72863949,  1060.55058771])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Millero, F.J., Chen, C.T., Bradshaw, A., and Schleicher, K.
    " A new high pressure equation of state for seawater"
    Deap-Sea Research., 1980, Vol27A, pp255-264.

    Authors
    -------
    Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    # UNESCO 1983. eqn.7  p.15
    densP0 = dens0(s, t)
    K      = seck(s, t, p)
    p      = p / 10.0 # convert from db to atm pressure units
    dens   = densP0 / ( 1-p / K )
    return dens

def pden(s, t, p, pr=0):
    """
    Calculates potential density of water mass relative to the specified reference pressure by pden = dens(S, ptmp, PR).

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"
    pr : number
         reference pressure [db], default = 0

    Returns
    -------
    pden : array_like
           potential denisty relative to the ref. pressure [kg m :sup:3]

    See Also
    --------
    ptmp, dens

    Notes
    -----
    The reference pressure is in "oceanographic" standards, so 0 db means at surface or 1 atm.

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22

    >>> import numpy as np
    >>> import seawater as sw
    >>> s = np.array([0, 0, 0, 0, 35, 35, 35, 35])
    >>> t = np.array([0, 0, 30, 30, 0, 0, 30, 30]) / T68conv
    >>> p = np.array([0, 10000, 0, 10000, 0, 10000, 0, 10000])
    >>> sw.pden(s, t, p)
    array([  999.842594  ,   999.79523994,   995.65113374,   996.36115932,
            1028.10633141,  1028.15738545,  1021.72863949,  1022.59634627])

    :math:`\\sigma_{4}` (at 4000 db)

    >>> sw.pden(s, t, p, 4000) - 1000
    array([ 19.2895493 ,  19.33422519,  12.43271053,  13.27563816,
            46.30976432,  46.48818851,  37.76150878,  38.74500757])

    References
    ----------
    A.E. Gill 1982. p.54
    "Atmosphere-Ocean Dynamics"
    Academic Press: New York.  ISBN: 0-12-283522-0

    Authors
    -------
    Phil Morgan 1992/04/06, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    pt   = ptmp(s, t, p, pr)
    pden = dens(s, pt, pr)
    return pden

def svan(s, t, p=0):
    """
    Specific Volume Anomaly calculated as svan = 1/dens(s, t, p) - 1/dens(35, 0, p).

    Note that it is often quoted in literature as 1e8*units.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"

    Returns
    -------
    svan : array_like
           specific volume anomaly  [m :sup:`3` kg :sup:`-1`]

    See Also
    --------
    dens

    Notes
    -----
    TODO

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22

    >>> import numpy as np
    >>> import seawater as sw
    >>> s = np.array([0, 0, 0, 0, 35, 35, 35, 35])
    >>> t = np.array([0, 0, 30, 30, 0, 0, 30, 30]) / T68conv
    >>> p = np.array([0, 10000, 0, 10000, 0, 10000, 0, 10000])
    >>> sw.svan(s, t, p)
    array([  2.74953924e-05,   2.28860986e-05,   3.17058231e-05,
             3.14785290e-05,   0.00000000e+00,   0.00000000e+00,
             6.07141523e-06,   9.16336113e-06])

    References
    ----------
    Fofonoff, N.P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
    #Eqn (9) p.15.

    S. Pond & G.Pickard  2nd Edition 1986. Introductory Dynamical Oceanography Pergamon Press Sydney. ISBN 0-08-028728-X

    Authors
    -------
    Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    svan = 1/dens( s, t, p ) - 1/dens( 35, 0.0, p )
    return svan

def gpan(s, t, p):
    """
    Geopotential Anomaly calculated as the integral of svan from the
    the sea surface to the bottom. THUS RELATIVE TO SEA SURFACE.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"
    axis : int, optional
           Axis along which the means are computed. The default is to compute using axis=1 or stations.
           Usually arrays are 2D MxN where M is depth and M stations.

    Returns
    -------
    gpan : array_like
           geopotential anomaly [m :sup:`3` kg :sup:`-1` Pa == m  :sup:`2` s :sup:`-2` == J kg :sup:`-1`]

    See Also
    --------
    svan, gvel

    Notes
    -----
    TODO: example with values that make some sense
    TODO: pass axis as argument

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22

    >>> import numpy as np
    >>> import seawater as sw
    >>> s = np.array([[0, 0, 0], [15, 15, 15], [30, 30, 30],[35,35,35]])
    >>> t = np.repeat(15, s.size).reshape(s.shape)
    >>> p = np.array([0, 250, 500, 1000])
    >>> sw.gpan(s, t, p)
    array([[   0.        ,    0.        ,    0.        ],
           [  56.35465209,   56.35465209,   56.35465209],
           [  84.67266947,   84.67266947,   84.67266947],
           [ 104.95799186,  104.95799186,  104.95799186]])

    References
    ----------
    S. Pond & G.Pickard  2nd Edition 1986. Introductory Dynamical Oceanogrpahy Pergamon Press Sydney. ISBN 0-08-028728-X

    Note that older literature may use units of "dynamic decimeter" for above.

    Adapted method from Pond and Pickard (p76) to calc gpan relative to sea surface whereas P&P calculated relative to the deepest common depth.

    Authors
    -------
    Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    if t.ndim == 1:
        print "need at least to profiles to compute geopotential anomaly"
        return

    if p.ndim == 1:
        p = np.repeat(p[np.newaxis,:], t.shape[1], axis=1).reshape(t.shape)

    svn        = svan(s, t, p)
    mean_svan  = 0.5*(svn[1:,:] + svn[0:-1,:])

    top = svn[0,:] * p[0,:] * db2Pascal

    delta_ga   = ( mean_svan * np.diff(p, axis=0) ) * db2Pascal
    ga         = np.cumsum( np.vstack( ( top, delta_ga ) ), axis=0 ) # TODO: test different integration methods

    return ga

def gvel(ga, distm, lat):
    """
    Calculates geostrophic velocity given the geopotential anomaly
    and position of each station.

    Parameters
    ----------
    ga : array_like
         geopotential anomaly relative to the sea surface.
    dist : array_like
           distance between stations [meters]
    lat : array_like
          lat to use for constant f determination

    Returns
    -------
    vel : array_like
           geostrophic velocity relative to the sea surface.
           dimension will be MxN-1 (N: stations)

    See Also
    --------
    svan, gpan

    Notes
    -----
    The original matlab version had gvel and gvel2, here the logic is "gvel2", where one must compute the distance first.
    TODO: dim(m, nstations-1) or pass axis?
    TODO: add example with a referece level
    TODO: example with values that make some sense
    TODO: pass axis as argument

    Examples
    --------
    >>> import numpy as np
    >>> import seawater as sw
    >>> lon = np.array([30,30,30])
    >>> lat = np.array([30,32,35])
    >>> s = np.array([[0, 1, 2], [15, 16, 17], [30, 31, 32],[35,35,35]])
    >>> t = np.repeat(15, s.size).reshape(s.shape)
    >>> p = np.array([0, 250, 500, 1000])
    >>> ga = sw.gpan(s,t,p)
    >>> distm = 1000.0 * sw.dist(lon, lat, units='km')[0]
    >>> sw.gvel(ga, distm, lat)
    array([[-0.        , -0.        ],
           [ 0.11385857,  0.07154328],
           [ 0.22436908,  0.14112984],
           [ 0.33366938,  0.20996603]])


    References
    ----------
    S. Pond & G.Pickard  2nd Edition 1986 Introductory Dynamical Oceanogrpahy Pergamon Press Sydney. ISBN 0-08-028728-X
    Equation 8.9A p73  Pond & Pickard

    Authors
    -------
    Phil Morgan 1992/03/26 (morgan@ml.csiro.au)

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    f     = cor( ( lat[0:-1] + lat[1:] )/2 )
    lf    = f * distm
    vel   = -np.diff(ga, axis=1) / lf

    return vel

def cp(s, t, p):
    """
    Heat Capacity of Sea Water using UNESCO 1983 polynomial.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"

    Returns
    -------
    cp : array_like
         specific heat capacity [J kg :sup:`-1` C :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    Data from Pond and Pickard Intro. Dynamical Oceanography 2nd ed. 1986
    >>> import numpy as np
    >>> import seawater as sw
    >>> from seawater import T68conv
    >>> t = np.array([[0, 0, 0, 0, 0, 0], [10, 10, 10, 10, 10, 10], [20, 20, 20, 20, 20, 20], [30, 30, 30, 30, 30, 30], [40, 40, 40, 40, 40, 40]]) / T68conv
    >>> s = np.array([[25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35]])
    >>> p = np.array([[0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000]])
    >>> sw.cp(s, t, p)
    array([[ 4048.4405375 ,  3896.25585   ,  3807.7330375 ,  3986.53309476,
             3849.26094605,  3769.11791286],
           [ 4041.8276691 ,  3919.5550066 ,  3842.3111366 ,  3986.34061786,
             3874.72665865,  3804.415624  ],
           [ 4044.8438591 ,  3938.5978466 ,  3866.7400391 ,  3993.85441786,
             3894.99294519,  3828.29059113],
           [ 4049.0984351 ,  3952.0375476 ,  3882.9855526 ,  4000.68382238,
             3909.24271128,  3844.32151784],
           [ 4051.2244911 ,  3966.1132036 ,  3905.9162711 ,  4003.46192541,
             3923.89463092,  3868.28959814]])

    References
    ----------
    Fofonff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Authors
    -------
    Phil Morgan, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    p = p/10. # to convert [db] to [bar] as used in UNESCO routines
    T68 = t * T68conv

    # eqn 26 p.32
    c0 = 4217.4
    c1 =   -3.720283
    c2 =    0.1412855
    c3 =   -2.654387e-3
    c4 =    2.093236e-5

    a0 = -7.64357
    a1 =  0.1072763
    a2 = -1.38385e-3

    b0 =  0.1770383
    b1 = -4.07718e-3
    b2 =  5.148e-5

    Cpst0 = ( ( ( c4 * T68 + c3 ) * T68 + c2 ) * T68 + c1 ) * T68 + c0 + \
            ( a0 + a1 * T68 + a2 * T68**2 ) * s + \
            ( b0 + b1 * T68 + b2 * T68**2 ) * s *(s)**0.5

    # eqn 28 p.33
    a0 = -4.9592e-1
    a1 =  1.45747e-2
    a2 = -3.13885e-4
    a3 =  2.0357e-6
    a4 =  1.7168e-8

    b0 =  2.4931e-4
    b1 = -1.08645e-5
    b2 =  2.87533e-7
    b3 = -4.0027e-9
    b4 =  2.2956e-11

    c0 = -5.422e-8
    c1 =  2.6380e-9
    c2 = -6.5637e-11
    c3 =  6.136e-13

    del_Cp0t0 = ( ( ( ( ( c3 * T68 + c2 ) * T68 + c1 ) * T68 + c0 ) * p + \
                ( ( ( ( b4 * T68 + b3 ) * T68 + b2 ) * T68 + b1 ) * T68 + b0 ) ) * p + \
                ( ( ( ( a4 * T68 + a3 ) * T68 + a2 ) * T68 + a1 ) * T68 + a0 ) ) * p

    # eqn 29 p.34
    d0 =  4.9247e-3
    d1 = -1.28315e-4
    d2 =  9.802e-7
    d3 =  2.5941e-8
    d4 = -2.9179e-10

    e0 = -1.2331e-4
    e1 = -1.517e-6
    e2 =  3.122e-8

    f0 = -2.9558e-6
    f1 =  1.17054e-7
    f2 = -2.3905e-9
    f3 =  1.8448e-11

    g0 =  9.971e-8

    h0 =  5.540e-10
    h1 = -1.7682e-11
    h2 =  3.513e-13

    j1 = -1.4300e-12

    S3_2 = s * (s)**0.5

    """
    del_Cpstp = [ ( ( ( ( d4 * T68 + d3 ) * T68 + d2 ) * T68 + d1 ) * T68 + d0 ) * S + \
                ( ( e2 * T68 + e1 ) * T68 + e0 ) * S3_2 ] * P + \
                [ ( ( ( f3 * T68 + f2 ) * T68 + f1 ) * T68 + f0 ) * S + \
                g0 * S3_2 ] * P**2 + \
                [ ( ( h2 * T68 + h1 ) * T68 + h0 ) * S + \
                j1 * T68 * S3_2 ] * P**3
    """

    del_Cpstp = ( ( ( ( ( d4 * T68 + d3 ) * T68 + d2 ) * T68 + d1 ) * T68 + d0 ) * s + \
                ( ( e2 * T68 + e1 ) * T68 + e0 ) * S3_2 ) * p + \
                ( ( ( ( f3 * T68 + f2 ) * T68 + f1 ) * T68 + f0 ) * s + \
                 g0 * S3_2 ) * p**2 + \
                ( ( ( h2 * T68 + h1 ) * T68 + h0 ) * s + \
                j1 * T68 * S3_2 ) * p**3

    cp = Cpst0 + del_Cp0t0 + del_Cpstp

    return cp

def ptmp(s, t, p, pr=0):
    """
    Calculates potential temperature as per UNESCO 1983 report.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"
    pr : array_like
        reference pressure [db], default = 0

    Returns
    -------
    pt : array_like
         potential temperature relative to PR [:math:`^\\circ` C (ITS-90)]

    See Also
    --------
    adtg, pden

    Notes
    -----
    TODO

    Examples
    --------
    >>> import numpy as np
    >>> import seawater as sw
    >>> from seawater import T68conv
    >>> t = np.array([[0, 0, 0, 0, 0, 0], [10, 10, 10, 10, 10, 10], [20, 20, 20, 20, 20, 20], [30, 30, 30, 30, 30, 30], [40, 40, 40, 40, 40, 40]]) / T68conv
    >>> s = np.array([[25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35],  [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35]])
    >>> p = np.array([[0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000]])
    >>> sw.ptmp(s, t, p, pr=0) * T68conv
    array([[  0.        ,  -0.30614418,  -0.96669485,   0.        ,
             -0.3855565 ,  -1.09741136],
           [ 10.        ,   9.35306331,   8.46840949,  10.        ,
              9.29063461,   8.36425752],
           [ 20.        ,  19.04376281,  17.94265   ,  20.        ,
             18.99845171,  17.86536441],
           [ 30.        ,  28.75124632,  27.43529911,  30.        ,
             28.72313484,  27.38506197],
           [ 40.        ,  38.46068173,  36.92544552,  40.        ,
             38.44979906,  36.90231661]])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
    Eqn.(31) p.39

    Bryden, H. 1973.
    "New Polynomials for thermal expansion, adiabatic temperature gradient and potential temperature of sea water."
    DEEP-SEA RES., 1973, Vol20,401-408.

    Authors
    -------
    Phil Morgan 92-04-06, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    # theta1
    del_P  = pr - p
    del_th = del_P * adtg(s, t, p)
    th     = t * T68conv + 0.5 * del_th
    q      = del_th

    # theta2
    del_th = del_P * adtg(s, th/T68conv, p + 0.5 * del_P )
    th     = th + ( 1 - 1/(2)**00.5 ) * ( del_th - q )
    q      = ( 2 - (2)**0.5 ) * del_th + ( -2 + 3/(2)**0.5 ) * q

    # theta3
    del_th = del_P * adtg( s, th/T68conv, p + 0.5 * del_P )
    th     = th + ( 1 + 1/(2)**0.5 ) * ( del_th - q )
    q      = ( 2 + (2)**0.5 ) * del_th + ( -2 -3/(2)**0.5 ) * q

    # theta4
    del_th = del_P * adtg( s, th/T68conv, p + del_P )
    pt     = ( th + ( del_th - 2 * q ) / 6 ) / T68conv
    return pt

def temp(s, pt, p, pr=0):
    """
    Calculates temperature from potential temperature at the reference pressure PR and in-situ pressure P.

    Parameters
    ----------
    s(p) : array_like
        salinity [psu (PSS-78)]
    pt(p) : array_like
         potential temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"
    pr : array_like
         reference pressure [db]

    Returns
    -------
    temp : array_like
           temperature [:math:`^\\circ` C (ITS-90)]

    See Also
    --------
    ptmp

    Notes
    -----
    TODO

    Examples
    --------
    >>> import numpy as np
    >>> import seawater as sw
    >>> from seawater import T68conv
    >>> s, t, p = 35, 15, 100
    >>> sw.temp(s, sw.ptmp(s, t, p), p)
    15.0

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    UNESCO 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _UNESCO Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
    Eqn.(31) p.39

    Bryden, H. 1973.
    "New Polynomials for thermal expansion, adiabatic temperature gradient and potential temperature of sea water."
    DEEP-SEA RES., 1973, Vol20,401-408.

    Authors
    -------
    Phil Morgan 92-04-06, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    # Carry out inverse calculation by swapping p0 & pr
    t = ptmp(s, pt, pr, p);

    return t

def swvel(lenth, depth):
    """
    Calculates surface wave velocity.

    lenth : array_like
            wave length
    depth : array_like
            water depth [meters]

    Returns
    -------
    speed : array_like
            surface wave speed [m s :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO: add my wave function

    Examples
    --------
    >>> import seawater as sw
    >>> sw.swvel(10,100)
    3.9493270848342941

    Authors
    ------
    Lindsay Pender 2005

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
    10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    k = 2.0 * np.pi / lenth
    speed = (g * np.tanh(k * depth) / k)**0.5
    return speed