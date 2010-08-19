# -*- coding: utf-8 -*-

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
    None

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested

    Examples
    --------
    Data from Unesco 1983 p45

    >>> t = np.array([[ 0,  0,  0,  0,  0,  0], [10, 10, 10, 10, 10, 10], [20, 20, 20, 20, 20, 20], [30, 30, 30, 30, 30, 30], [40, 40, 40, 40, 40, 40]])
    >>> t = t / T68conv
    >>> s = np.array([[25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35]])
    >>> p = np.array([[0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000], [0, 5000, 10000, 0, 5000, 10000]])
    >>> sw.adtg(s, t, p)
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

    >>> s, pt, p = 40., 10., 4000.
    >>> sw.alpha(s, pt, p)
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
    >>> s, pt, p = 40.0, 10.0, 4000
    >>> sw.aonb(s, pt, p)
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
    >>> s, pt, p = 40.0, 10.0, 4000
    >>> sw.beta(s,pt,p)
    0.00072087661741618932

    Notes
    -----
    Pressure broadcast feature need to be tested
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

    p = np.float32(p)
    t = t * T68conv

    c1 = np.array([-0.415613e-9, 0.555579e-7, -0.301985e-5, 0.785567e-3])
    c2 = np.array([0.788212e-8, -0.356603e-6])
    c3 = np.array([-0.602281e-15, 0.408195e-10, 0.0])
    c4 = 0.515032e-8
    c5 = np.array([-0.213127e-11, 0.192867e-9, -0.121555e-7])
    c6 = np.array([-0.175379e-14, 0.176621e-12])
    c7 = 0.121551e-17

    # Now calaculate the thermal expansion saline contraction ratio adb
    sm35  = S - 35
    beta  = np.polyval(c1, pt) + sm35 * (np.polyval(c2, pt) + \
            np.polyval(c3, P) ) + c4 * (sm35**2) + \
            P * np.polyval(c5, pt) + (P**2) * np.polyval(c6, pt) \
            + c7 * (P**3)

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
    >>> n2, q, p_ave = sw.bfrq(s, t, p, lat)

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
    Unesco 1983 data p30

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
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Authors
    -------
    Phil Morgan 92-04-06  (morgan@ml.csiro.au)

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
    None

    Examples
    --------
    >>> sw.grav(lat, z=0)
    9.8061898752053995

    References
    ----------
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

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

    # Eqn p27.  Unesco 1983.
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

    # Eqn p27.  Unesco 1983.
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
    Data from Unesco 1983 p9

    >>> t    = np.array([0, 10, 0, 10, 10, 30]) / T68conv
    >>> p    = np.array([0, 0, 1000, 1000, 0, 0])
    >>> s    = np.array([25, 25, 25, 25, 40, 40])
    >>> sw.cndr(s, t, p)
    array([0.49800825, 0.65499015, 0.50624434, 0.66297496, 1.00007311, 1.52996697])

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

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
    # eqn(4) p.8 Unesco 1983
    d1 =  3.426e-2
    d2 =  4.464e-4
    d3 =  4.215e-1
    d4 = -3.107e-3

    e1 =  2.070e-5
    e2 = -6.370e-10
    e3 =  3.989e-15

    A  = ( d3 + d4 * T68 )
    B  = 1 + d1 * T68 + d2 * T68**2
    C  = P * ( e1 + e2 * p + e3 * P**2 )

    # eqn(6) p.9 UNESCO 1983.
    Rt    = Rx * Rx
    rt    = salrt(t)
    Rtrt  = rt * Rt
    D     = B - A * rt * Rt
    E     = rt * Rt * A * ( B + C )
    r     = ( abs( D**2 + 4 * E ) )**0.5 - D
    r     = 0.5 * R/A
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
    None

    Examples
    --------
    Data from Unesco 1983 p9

    >>> t    = np.array([15, 20, 5]) / T68conv
    >>> rt   = np.array([  1, 1.0568875, 0.81705885])
    >>> sw.sals(rt,t)
    array([ 35.        ,  37.24562718,  27.99534701])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

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
    None

    Examples
    --------
    Sata from Unesco 1983 p9
    >>> delt = np.array([15, 20, 5]) / T68conv  - 15
    >>> rtx  = np.array([  1, 1.0568875, 0.81705885])**0.5
    >>> sw.salds(rtx, delt)
    array([ 78.31921607,  81.5689307 ,  68.19023687])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

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
    None

    Examples
    --------
    Data from Unesco 1983 p9

    >>> t = np.array([15, 20, 5]) / T68conv
    >>> sw.saltrt(t)
    array([ 1.,  1.11649272,  0.77956585])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Authors
    -------
    Phil Morgan 93-04-17, Lindsay Pender (Lindsay.Pender@csiro.au)

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    #Eqn (3) p.7 Unesco.
    T68 = T * T68conv

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
    None

    Examples
    --------
    Data from Unesco 1983 p9

    >>> r = np.array([1, 1.2, 0.65])
    >>> t = np.array([15, 20, 5]) / T68conv
    >>> p = np.array([0, 2000, 1500])
    >>> sw.salt(r, t, p)
    array([ 34.99999992,  37.24562765,  27.99534693])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

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
    None

    Examples
    --------
    >>> r = np.array([1, 1.2, 0.65])
    >>> t = np.array([15, 20, 5]) / T68conv
    >>> p = np.array([0, 2000, 1500])
    >>> sw.salrp(r, t, p)
    array([1., 1.01694294, 1.02048638])

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

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
    None

    Notes
    -----
    None

    Examples
    --------
    UNESCO DATA p.30

    >>> s = np.array([[5, 10, 15, 20, 25, 30, 35, 40], [5, 10, 15, 20, 25, 30, 35, 40]])
    >>> p = np.array([[ 0, 0, 0, 0, 0, 0, 0, 0], [500, 500, 500, 500, 500, 500, 500, 500]])
    >>> sw.fp(s, p)
    array([[-0.27369757, -0.54232831, -0.81142026, -1.0829461, -1.35804594, -1.63748903, -1.9218401, -2.2115367], [-0.65010724, -0.91873798, -1.18782992, -1.45935577, -1.73445561, -2.01389869, -2.29824976, -2.58794636]])

    References
    ----------
    Fofonff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

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
    #TODO: P = P/10; # to convert db to Bar as used in Unesco routines (was commented in the original)

    # eqn  p.29
    a0 = -0.0575
    a1 =  1.710523e-3
    a2 = -2.154996e-4
    b  = -7.53e-4

    fp = ( a0 * s + a1 * s * (s)**0.5 + a2 * s**2 + b * p ) / T68conv

    return fp