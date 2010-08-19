# -*- coding: utf-8 -*-
"""
seawater.py
Translated from matlab CSIRO seawater toolbox Version 3.2

Filipe P. A. Fernandes
e-mail:   ocefpaf@gmail.com
web:      http://ocefpaf.tiddlyspot.com/
date:     14-Jan-2010
modified: 17-Aug-2010
obs:      flag: TODO
          some keywords and default values are hardcoded!!!
          create docstrings
          shear/richnumber are in "extras"
"""

# IMPORTS: (TODO: change to import numpy as np)
from numpy import array, ones, polyval, float32, diff, pi, sin, cos, \
                  angle, arange, log, exp, cumsum, vstack, hstack, tanh, sign

#(TODO: move these to test routine)
from os    import uname
from time  import asctime, localtime
from sys   import version

# COMMON:
DEG2RAD = pi/180

# FUNCTIONS:
def adtg(s, t, p):
    """
    Calculates adiabatic temperature gradient as per UNESCO 1983 routines.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [ :math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"

    Returns
    -------
    adtg : array_like
           adiabatic temperature gradient [ :math:`^\\circ` C db :sup:`-1` ]

    See Also
    --------
    ptmp : calls adtg

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested

    Examples
    --------
    Data from Unesco 1983 p45

    >>> t = array([[ 0,  0,  0,  0,  0,  0], [10, 10, 10, 10, 10, 10], [20, 20, 20, 20, 20, 20], [30, 30, 30, 30, 30, 30], [40, 40, 40, 40, 40, 40]])
    >>> t = t / 1.00024
    >>> s = array([[25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35]])
    >>> p = array([[0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000]])
    >>> adtg(s, t, p)
     array([[ 1.68710000e-05, 1.04700000e-04, 1.69426000e-04,
              3.58030000e-05, 1.17956500e-04, 1.77007000e-04],
            [ 1.00213255e-04, 1.60972016e-04, 2.06883149e-04,
              1.14904937e-04, 1.71376482e-04, 2.13000398e-04],
            [ 1.73853488e-04, 2.13558726e-04, 2.44501964e-04,
              1.84304853e-04, 2.21111157e-04, 2.49155462e-04],
            [ 2.41768241e-04, 2.64801063e-04, 2.82987774e-04,
              2.47979289e-04, 2.69501460e-04, 2.86177520e-04],
            [ 3.07934056e-04, 3.17039963e-04, 3.23045906e-04,
              3.09904786e-04, 3.18888326e-04, 3.24771901e-04]])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater. Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39

    Bryden, H. 1973.
    "New Polynomials for thermal expansion, adiabatic temperature gradient
    and potential temperature of sea water."
    DEEP-SEA RES., 1973, Vol20,401-408.

    Author
    ------
    Phil Morgan, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-16. Filipe Fernandes, Reformulated docstring.
    """

    T68 = 1.00024 * t

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
           temperature or potential temperature [ :math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"
    pt : bool
         True if temperature is potential, default is False

    Returns
    -------
    alpha : array_like
            thermal expansion coeff ( :math:`\\alpha` ) [ :math:`^\\circ` C :sup:`-1` ]

    See Also
    --------
    beta and aonb

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested

    Examples
    --------
    Data from McDougall 1987

    >>> s, pt, p = 40., 10., 4000.
    >>> alpha(s, pt, p)
    0.00025061316481624323

    Reference
    ---------
    McDougall, T.J. 1987.  "Neutral Surfaces"
    Journal of Physical Oceanography vol 17 pages 1950-1964,

    Author
    ------
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
           temperature or potential temperature [:math:`^\\circ`  C (ITS-90)]
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
    alpha and beta

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested
    TODO: Test pt=False

    Examples
    --------
    >>> s, pt, p = 40.0, 10.0, 4000
    >>> aonb(s, pt, p)
    0.347650567047807

    Reference
    ---------
    McDougall, T.J. 1987. "Neutral Surfaces"
    Journal of Physical Oceanography vol 17 pages 1950-1964,

    Author
    ------
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

    p = float32(p)
    t = t * 1.00024

    c1  = array([-0.255019e-7, 0.298357e-5, -0.203814e-3, \
                    0.170907e-1, 0.665157e-1])
    c2  = array([-0.846960e-4, 0.378110e-2])
    c2a = array([-0.251520e-11, -0.164759e-6, 0.0])
    c3  = -0.678662e-5
    c4  = array([0.791325e-8, -0.933746e-6, 0.380374e-4])
    c5  =  0.512857e-12
    c6  = -0.302285e-13
    # Now calculate the thermal expansion saline contraction ratio adb
    sm35  = s - 35.0
    aonb  = polyval(c1, t) + sm35 * ( polyval(c2, t) \
            + polyval(c2a, p) ) \
            + sm35**2 * c3 + p * polyval(c4, t) \
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
           temperature or potential temperature [:math:`^\\circ`  C (ITS-90)]
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
    >>> s, pt, p = 40.0, 10.0, 4000
    >>> sw.beta(s,pt,p)
    0.00072087661741618932

    Notes
    -----
    Pressure broadcast feature need to be tested
    TODO: Test pt=False for alpha, beta and aonb

    Author
    ------
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

    p = float32(p)
    t = t * 1.00024

    c1 = array([-0.415613e-9, 0.555579e-7, -0.301985e-5, 0.785567e-3])
    c2 = array([0.788212e-8, -0.356603e-6])
    c3 = array([-0.602281e-15, 0.408195e-10, 0.0])
    c4 = 0.515032e-8
    c5 = array([-0.213127e-11, 0.192867e-9, -0.121555e-7])
    c6 = array([-0.175379e-14, 0.176621e-12])
    c7 = 0.121551e-17

    # Now calaculate the thermal expansion saline contraction ratio adb
    sm35  = S - 35
    beta  = polyval(c1, pt) + sm35 * (polyval(c2, pt) + \
            polyval(c3, P) ) + c4 * (sm35**2) + \
            P * polyval(c5, pt) + (P**2) * polyval(c6, pt) \
            + c7 * (P**3)

    return beta

def bfrq(s, t, p, lat=None):
    """
    Calculates Brunt-Vaisala Frequency squared (N^2) at the mid depths
    from the equation, :math:`N^{2} = \\frac{-g}{\\sigma_{\\theta}} \\frac{d\\sigma_{\\theta}}{dz}`

    Also calculates Potential Vorticity from, :math:`q=f \\frac{N^2}{g}`

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature or potential temperature [:math:`^\\circ`  C (ITS-90)]
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

    Examples
    --------
    TODO: add a test here
    >>> n2, q, p_ave = bfrq(s, t, p, lat)

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested

    See Also
    --------
    pden, dens

    Author
    ------
    Phil Morgan 93-06-24, Lindsay Pender (Lindsay.Pender@csiro.au)
    Greg Johnson (gjohnson@pmel.noaa.gov)
    added potential vorticity calculation

    References
    ----------
    A.E. Gill 1982. p.54  eqn 3.7.15
    "Atmosphere-Ocean Dynamics"
    Academic Press: New York.  ISBN: 0-12-283522-0

    Jackett, D.R. and McDougall, T.J. 1994.
    Minimal adjustment of hydrographic properties to achieve static
    stability.  submitted J.Atmos.Ocean.Tech.

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    06-04-19. Lindsay Pender, Corrected sign of PV.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-17. Filipe Fernandes, Reformulated docstring.
    """

    if lat is not None:
        z = depth(p, lat)
        g = grav(lat, -z) # note that grav expects height as argument
        f = cor(lat)
    else: # TODO: test this if logic
        z = p
        g = 9.8
        f = NaN

    m,n   = p.shape # TODO: check where depth increases to automagically find which dimension to operate
    iup   = arange(0, m-1)
    ilo   = arange(1, m)

    p_ave    = ( p[iup,:] + p[ilo,:] )/2.
    pden_up  = pden( s[iup,:], t[iup,:], p[iup,:], p_ave )
    pden_lo  = pden( s[ilo,:], t[ilo,:], p[ilo,:], p_ave )
    mid_pden = ( pden_up + pden_lo )/2
    dif_pden = pden_up - pden_lo
    mid_g    = ( g[iup,:] + g[ilo,:] )/2
    dif_z    = diff(z, axis=0) # TODO: check where depth increases to automagically find which dimension to operate
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
    Unesco 1983 data p30

    >>> lat = array([0, 30, 45, 90])
    >>> p   = array([[  500,   500,   500,  500], [ 5000,  5000,  5000, 5000], [10000, 10000, 10000, 10000]])
    >>> z = depth(p, lat)
    array([[  496.65299239,   495.99772917,   495.3427354 ,   494.03357499],
       [ 4915.04099112,  4908.55954332,  4902.08075214,  4889.13132561],
       [ 9725.47087508,  9712.6530721 ,  9699.84050403,  9674.23144056]])

    Notes
    -----
    original seawater name is dpth

    Author
    ------
    Phil Morgan 92-04-06  (morgan@ml.csiro.au)

    References
    ----------
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    # Eqn 25, p26.  Unesco 1983.
    c1 =  9.72659
    c2 = -2.2512E-5
    c3 =  2.279E-10
    c4 = -1.82E-15

    gam_dash = 2.184e-6

    lat = abs(lat)
    X   = sin( lat * DEG2RAD ) # convert to radians
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

    Examples
    --------
    >>> g = grav(lat, z=0)
    9.8061898752053995

    See Also
    --------
    bfrq

    Author
    ------
    Phil Morgan 93-04-20  (morgan@ml.csiro.au)

    References
    ----------
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    A.E. Gill 1982. p.597
    "Atmosphere-Ocean Dynamics"
    Academic Press: New York.  ISBN: 0-12-283522-0

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    # Eqn p27.  Unesco 1983.
    a       = 6371000. # mean radius of earth  A.E.Gill
    lat     = abs(lat)
    X       = sin( lat * DEG2RAD )  # convert to radians
    sin2    = X * X
    grav    = 9.780318 * ( 1.0 + ( 5.2788E-3 + 2.36E-5 * sin2 ) * sin2 )
    grav    = grav / ( ( 1 + z/a )**2 )    # from A.E.Gill p.597
    return  grav