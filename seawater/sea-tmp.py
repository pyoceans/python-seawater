def depth(P, LAT):
    """
    USAGE:  z = depth(P, LAT)

    DESCRIPTION:
    Calculates depth in metres from pressure in dbars.

    INPUT:  (all must have same dimensions)
    P   = Pressure    [db]
    LAT = Latitude in decimal degress north [-90..+90]
            (lat may have dimensions 1x1 or 1xn where P(mxn).

    OUTPUT:
      z = depth [metres]

    Author
    ------
    Phil Morgan 92-04-06  (morgan@ml.csiro.au)

    References
    ----------
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    CALLER:  general purpose
    CALLEE:  none

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    10-01-14. Filipe Fernandes, Python translation.
    """
    # Eqn 25, p26.  Unesco 1983.
    c1 =  9.72659
    c2 = -2.2512E-5
    c3 =  2.279E-10
    c4 = -1.82E-15

    gam_dash = 2.184e-6

    LAT = abs(LAT)
    X   = sin( LAT * DEG2RAD );  # convert to radians
    X   = X * X

    bot_line = 9.780318 * ( 1.0 + ( 5.2788E-3 + 2.36E-5 * X ) * X ) + \
               gam_dash * 0.5 * P
    top_line = ( ( ( c4 * P + c3 ) * P + c2 ) * P + c1 ) * P
    DEPTHM   = top_line / bot_line
    return DEPTHM

def grav(LAT, z=0):
    """
    USAGE:  g = grav(LAT, z=0)

    DESCRIPTION:
    Calculates acceleration due to gravity as function of latitude.

    INPUT:  (all must have same dimensions)
    LAT = Latitude in decimal degress north [-90..+90]
    z   = height in metres (+ve above sea surface, -ve below)
          default z = 0

    OUTPUT:
    g    = gravity [m/s^2]

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

    CALLER:  general purpose
    CALLEE:  none

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
    """
    # Eqn p27.  Unesco 1983.
    a       = 6371000.    # mean radius of earth  A.E.Gill
    LAT     = abs(LAT)
    X       = sin( LAT * DEG2RAD )  # convert to radians
    sin2    = X * X
    grav    = 9.780318 * ( 1.0 + ( 5.2788E-3 + 2.36E-5 * sin2 ) * sin2 )
    grav    = grav / ( ( 1 + z/a )**2 )    # from A.E.Gill p.597
    return  grav

def cor(LAT):
    """
    USAGE:  f = cor(LAT)

    DESCRIPTION:
    Calculates the Coriolis factor "f" defined by
        f = 2*Omega*Sin(lat)  where Omega = 7.292e-5 radians/sec.

    INPUT:
    LAT = Latitude in decimal degress north [-90..+90]

    OUTPUT:
    f    = Coriolis Factor "f" [s-1]

    Author
    ------
    Phil Morgan 93-04-20  (morgan@ml.csiro.au)

    References
    ----------
    S. Pond & G.Pickard  2nd Edition 1986
    Introductory Dynamical Oceanogrpahy
    Pergamon Press Sydney.  ISBN 0-08-028728-X

    A.E. Gill 1982. p.597
    "Atmosphere-Ocean Dynamics"
    Academic Press: New York.  ISBN: 0-12-283522-0

    CALLER:  general purpose
    CALLEE:  none

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
    """
    # Eqn p27.  Unesco 1983.
    OMEGA   = 7.292e-5   # s**-1   A.E.Gill p.597
    f       = 2 * OMEGA * sin( LAT * DEG2RAD )
    return f

def c3515():
    """
    USAGE:  c3515 = c3515

    DESCRIPTION:
    Returns conductivity at S=35 psu , T=15 C [ITPS 68] and P=0 db).

    INPUT: (none)

    OUTPUT:
    c3515  = Conductivity   [mmho/cm == mS/cm]

    Author
    ------
    Phil Morgan 93-04-17  (morgan@ml.csiro.au)

    References
    ----------
    R.C. Millard and K. Yang 1992.
    "CTD Calibration and Processing Methods used by Woods Hole
    Oceanographic Institution"  Draft April 14, 1992
    (Personal communication)

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
    """
    return 42.914

def cndr(S, T, P):
    """
    USAGE:  cndr = cndr(S, T, P)

    DESCRIPTION:
    Calculates conductivity ratio from S, T, P.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78) ]
    T = temperature [:math:`^\\circ`  C (ITS-90)]
    P = pressure    [db]

    OUTPUT:
    cndr = Conductivity ratio     R =  C(S,T,P)/C(35,15(IPTS-68),0) [no units]

    Author
    ------
    Phil Morgan 93-04-21, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    CALLER: general purpose
    CALLEE: salds sals salrt

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    T68 = T * 1.00024
    DT  = T68 - 15

    Rx  = (S/35.0)**0.5 # first guess at Rx = sqrt(Rt)
    SInc  = sals(Rx*Rx, T)

    # DO A NEWTON-RAPHSON ITERATION FOR INVERSE INTERPOLATION OF Rt FROM S.
    for n in range(10): # TODO: is a 100 important? changed to 10
        Rx   = Rx + (S - SInc) / salds(Rx, DT)
        SInc = sals(Rx*Rx, T)
        #try:
            #DELS = abs(SInc - S)
        #except TypeError:
        DELS = abs(SInc - S)
        DELS = array(DELS) # TODO: "any" is an array method
        if (DELS.any() < 1.0E-4):
            break

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
    C  = P * ( e1 + e2 * P + e3 * P**2 )

    # eqn(6) p.9 UNESCO 1983.
    Rt    = Rx * Rx
    rt    = salrt(T)
    Rtrt  = rt * Rt
    D     = B - A * rt * Rt
    E     = rt * Rt * A * ( B + C )
    R     = ( abs( D**2 + 4 * E ) )**0.5 - D
    R     = 0.5 * R/A
    return R

def sals(Rt, T):
    """
    USAGE:  S = sals(Rt, T)

    DESCRIPTION:
    Salinity of sea water as a function of Rt and T.
    UNESCO 1983 polynomial.

    INPUT:
    Rt = Rt(S,T) = C(S,T,0) / C(35, T(IPTS-68), 0)
    T  = temperature [:math:`^\\circ`  C (ITS-90)]

    OUTPUT:
    S  = salinity    [psu      (PSS-78)]

    Author
    ------
    Phil Morgan 93-04-17, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    CALLER: salt
    CALLEE: none

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """
    # eqn (1) & (2) p6,7 unesco
    del_T68 = T * 1.00024 - 15

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

    Rtx   = (Rt)**0.5
    del_S = ( del_T68 / ( 1 + k * del_T68 ) ) * \
            ( b0 + ( b1 + ( b2 + ( b3 + ( b4 + b5 * Rtx ) \
            * Rtx ) * Rtx ) * Rtx ) * Rtx)

    S = a0 + ( a1 + ( a2 + ( a3 + ( a4 + a5 * Rtx) * \
               Rtx) * Rtx ) * Rtx ) * Rtx

    S = S + del_S

    return S

def salds(Rtx, delT):
    """
    USAGE:  dS = salds(Rtx,delT)

    DESCRIPTION:
    Calculates Salinity differential dS/d(sqrt(Rt)) at constant T.
    UNESCO 1983 polynomial.

    INPUT: (all must have same dimensions)
    Rtx   = sqrt(Rt) where Rt defined in salt
    delT  = T-15     [:math:`^\\circ`  C (IPTS-68)]

    OUTPUT:
    dS = S differential dS/d(sqrt(Rt)) at constant T.

    Author
    ------
    Phil Morgan 93-04-21, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    CALLER: cndr
    CALLEE: none

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
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

    dS =  a1 + ( 2 * a2 + ( 3 * a3 + ( 4 * a4 + 5 * a5 * Rtx ) * Rtx ) \
          * Rtx ) * Rtx + ( delT / ( 1 + k * delT ) ) * \
          ( b1 + ( 2 * b2 + ( 3 * b3 + ( 4 * b4 + 5 * b5 * Rtx ) \
          * Rtx ) * Rtx) * Rtx)

    return dS

def salrt(T):
    """
    USAGE:  rt = salrt(T)

    DESCRIPTION:
       Equation rt(T) = C(35,T,0) / C(35,15(IPTS-68), 0)
       used in calculating salinity. UNESCO 1983 polynomial..

    INPUT:
      T = temperature [:math:`^\\circ`  C (ITS-90)]

    OUTPUT:
      rt = conductivity ratio  [no units]

    Author
    ------
    Phil Morgan 93-04-17, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    CALLER: salt
    CALLEE: none

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """
    #rt = rt(T) = C(35,T,0)/C(35,15,0)
    #Eqn (3) p.7 Unesco.
    T68 = T * 1.00024

    c0 =  0.6766097
    c1 =  2.00564e-2
    c2 =  1.104259e-4
    c3 = -6.9698e-7
    c4 =  1.0031e-9

    rt = c0 + ( c1 + ( c2 + ( c3 + c4 * T68) * T68) * T68 ) * T68
    return rt

def salt(cndr, T, P):
    """
    USAGE: S = salt(cndr, T, P)

    DESCRIPTION:
    Calculates Salinity from conductivity ratio. UNESCO 1983 polynomial.

    INPUT:
    cndr = Conductivity ratio     R =  C(S,T,P)/C(35,15(IPTS-68),0) [no units]
    T    = temperature [:math:`^\\circ`  C (ITS-90)]
    P    = pressure    [db]

    OUTPUT:
    S    = salinity    [psu      (PSS-78)]

    Author
    ------
    Phil Morgan 93-04-17, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    CALLER: general purpose
    CALLEE: sals salrt salrp

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """
    R  = cndr
    rt = salrt(T)
    Rp = salrp( R, T, P )
    Rt = R / ( Rp * rt )
    S  = sals( Rt, T )

    return S

def salrp(R, T, P):
    """
    USAGE:  Rp = salrp(R, T, P)

    DESCRIPTION:
    Equation Rp(S,T,P) = C(S,T,P)/C(S,T,0) used in calculating salinity.
    UNESCO 1983 polynomial.

    INPUT: (All must have same shape)
    R = Conductivity ratio  R =  C(S,T,P)/C(35,15(IPTS-68),0) [no units]
    T = temperature [:math:`^\\circ`  C (ITS-90)]
    P = pressure    [db]

    OUTPUT:
    Rp = conductivity ratio  Rp(S,T,P) = C(S,T,P)/C(S,T,0)  [no units]

    Author
    ------
    Phil Morgan 93-04-17, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    CALLER: salt
    CALLEE: none

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """
    # eqn (4) p.8 unesco.
    T68 = T * 1.00024

    d1 =  3.426e-2
    d2 =  4.464e-4
    d3 =  4.215e-1
    d4 = -3.107e-3

    e1 =  2.070e-5
    e2 = -6.370e-10
    e3 =  3.989e-15

    Rp = 1 + ( P * ( e1 + e2 * P + e3 * P**2 ) ) \
         / ( 1 + d1 * T68 + d2 * T68**2 + ( d3 + d4 * T68 ) * R)

    return Rp

def fp(S, P):
    """
    USAGE:  fp = fp(S, P)

    DESCRIPTION:
    Freezing point of Sea Water using UNESCO 1983 polynomial.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78)]
    P = pressure    [db]

    OUTPUT:
    fp = Freezing Point temperature [:math:`^\\circ`  C (ITS-90)]

    Author
    ------
    Phil Morgan 93-04-20, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    CALLER: general purpose
    CALLEE: none

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    #P = P/10; # to convert db to Bar as used in Unesco routines
    # eqn  p.29
    a0 = -0.0575
    a1 =  1.710523e-3
    a2 = -2.154996e-4
    b  = -7.53e-4

    fp = ( a0 * S + a1 * S * (S)**0.5 + a2 * S**2 + b * P ) / 1.00024

    return fp

def svel(S, T, P):
    """
    USAGE:  svel = svel(S, T, P)

    DESCRIPTION:
    Sound Velocity in sea water using UNESCO 1983 polynomial.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78)]
    T = temperature [:math:`^\\circ`  C (ITS-90)]
    P = pressure    [db]

    OUTPUT:
    svel = sound velocity  [m/s]

    Author
    ------
    Phil Morgan 93-04-20, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    CALLER: general purpose
    CALLEE: none

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """
    #UNESCO 1983. eqn.33  p.46
    P = P/10  # convert db to bars as used in UNESCO routines
    T68 = T * 1.00024

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

    Cw  = ( ( ( ( c32 * T68 + c31 ) * T68 + c30 ) * P + \
          ( ( ( ( c24 * T68 + c23 ) * T68 + c22 ) * T68 + c21 ) * T68 + c20 ) ) * P + \
          ( ( ( ( c14 * T68 + c13 ) * T68 + c12 ) * T68 + c11 ) * T68 + c10 ) ) * P + \
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

    A = ( ( ( ( a32 * T68 + a31 ) * T68 + a30 ) * P + \
        ( ( ( a23 * T68 + a22 ) * T68 + a21 ) * T68 + a20 ) ) * P + \
        ( ( ( ( a14 * T68 + a13 ) * T68 + a12 ) * T68 + a11 ) * T68 + a10 ) ) * P + \
        ( ( ( a04 * T68 + a03 ) * T68 + a02 ) * T68 + a01 ) * T68 + a00

    # eqn 36 p.47
    b00 = -1.922e-2
    b01 = -4.42e-5
    b10 =  7.3637e-5
    b11 =  1.7945e-7

    B = b00 + b01 * T68 + ( b10 + b11 * T68) * P

    # eqn 37 p.47
    d00 =  1.727e-3
    d10 = -7.9836e-6

    D = d00 + d10 * P

    # eqn 33 p.46
    svel = Cw + A * S + B * S * (S)**0.5 + D * S**2

    return svel

def pres(DEPTH, LAT):
    """
    USAGE:  pres = pres(DEPTH, LAT)

    DESCRIPTION:
    Calculates pressure in dbars from depth in meters.

    INPUT:  (all must have same dimensions)
    depth = depth [metres]
    lat   = Latitude in decimal degress north [-90..+90]

    OUTPUT:
    pres   = Pressure    [db]

    Author
    ------
    Phil Morgan 93-06-25  (morgan@ml.csiro.au)

    References
    ----------
    Saunders, P.M. 1981
    "Practical conversion of Pressure to Depth"
    Journal of Physical Oceanography, 11, 573-574

    CHECK VALUE:
    P=7500.00 db for LAT=30 deg, depth=7321.45 meters

    CALLER:  general purpose
    CALLEE:  none

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    10-01-14. Filipe Fernandes, Python translation.
    """
    X       = sin( abs(LAT) * DEG2RAD )  # convert to radians
    C1      = 5.92E-3 + X**2 * 5.25E-3
    pres    = ( ( 1 - C1 ) - ( ( ( 1 - C1 )**2 ) - ( 8.84E-6 * DEPTH ) )**0.5 ) / 4.42E-6
    return pres

def dist(lon, lat): # TODO: add keywords options for units
    """
    USAGE:  dist, phaseangle = dist(lon, lat)

    DESCRIPTION:
    Calculate distance between two positions on globe using the "Plane
    Sailing" method.  Also uses simple geometry to calculate the bearing of
    the path between position pairs.

    INPUT:
    lon      = decimal degrees (+ve E, -ve W) [-180..+180]
    lat      = decimal degrees (+ve N, -ve S) [- 90.. +90]

    OUTPUT:
    dist        = distance between positions in units
    phaseangle  = angle of line between stations with x axis (East).
                    Range of values are -180..+180. (E=0, N=90, S=-90)

    Author
    ------
    Phil Morgan and Steve Rintoul 92-02-10

    References
    ----------
    The PLANE SAILING method as descriibed in "CELESTIAL NAVIGATION" 1989 by
    Dr. P. Gormley. The Australian Antartic Division.

    CALLER:   general purpose
    CALLEE:   none

    Modifications
    -------------
    99-06-25. Lindsay Pender, Function name change from distance to sw_dist.
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    10-01-14. Filipe Fernandes, Python translation.
    """
    DEG2MIN = 60
    DEG2NM  = 60
    NM2KM   = 1.8520    # Defined in Pond & Pickard p303.

    npositions = max(lat.shape)
    ind = arange( 0, npositions-1, 1) # index to first of position pairs

    dlon = diff(lon, axis=0)
    if any( abs(dlon) > 180 ):
        flag = abs(dlon) > 180
        dlon[flag] = -sign( dlon[flag] ) * ( 360 - abs( dlon[flag] ) )

    latrad = abs( lat * DEG2RAD )
    dep    = cos( ( latrad [ind+1] + latrad[ind] ) / 2 ) * dlon
    dlat   = diff( lat, axis=0 )
    dist   = DEG2NM * ( dlat**2 + dep**2 )**0.5
    # TODO: add keyword options defaults is in miles
    #if strcmp(units,'km') # TODO: add keyword options defaults to n.miles
    dist = dist * NM2KM
    #end

    # CALCUALTE ANGLE TO X AXIS
    RAD2DEG     = 1/DEG2RAD
    phaseangle  = angle( dep + dlat * 1j ) * RAD2DEG
    return dist, phaseangle

def satAr(S, T):
    """
    USAGE:  satAr = satAr(S, T)

    DESCRIPTION:
    Solubility (satuaration) of Argon (Ar) in sea water.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78)]
    T = temperature [:math:`^\\circ`  C (ITS-90)]

    OUTPUT:
    satAr = solubility of Ar  [ml/l]

    Author
    ------
    Phil Morgan 97-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Weiss, R. F. 1970
    "The solubility of nitrogen, oxygen and argon in water and seawater."
    Deap-Sea Research., 1970, Vol 17, pp721-735.

    CALLER: general purpose
    CALLEE:

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    # convert T to Kelvin
    T = 273.15 + T * 1.00024

    # constants for Eqn (4) of Weiss 1970
    a1 = -173.5146
    a2 =  245.4510
    a3 =  141.8222
    a4 =  -21.8020
    b1 =   -0.034474
    b2 =    0.014934
    b3 =   -0.0017729

    # Eqn (4) of Weiss 1970
    lnC = a1 + a2 * ( 100/T ) + a3 * log( T/100 ) + a4 * ( T/100 ) + \
          S * ( b1 + b2 * ( T/100 ) + b3 * ( ( T/100 )**2) )

    c = exp(lnC)

    return c

def satN2(S, T):
    """
    USAGE:  satN2 = satN2(S, T)

    DESCRIPTION:
    Solubility (satuaration) of Nitrogen (N2) in sea water.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78)]
    T = temperature [:math:`^\\circ`  C (ITS-90)]

    OUTPUT:
    satN2 = solubility of N2  [ml/l]

    Author
    ------
    Phil Morgan 97-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Weiss, R. F. 1970
    "The solubility of nitrogen, oxygen and argon in water and seawater."
    Deap-Sea Research., 1970, Vol 17, pp721-735.

    CALLER: general purpose
    CALLEE:

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    # convert T to Kelvin
    T = 273.15 + T * 1.00024

    # constants for Eqn (4) of Weiss 1970
    a1 = -172.4965
    a2 =  248.4262
    a3 =  143.0738
    a4 =  -21.7120
    b1 =   -0.049781
    b2 =    0.025018
    b3 =   -0.0034861

    # Eqn (4) of Weiss 1970
    lnC = a1 + a2 * ( 100/T ) + a3 * log( T/100 ) + a4 * ( T/100 ) + \
          S * ( b1 + b2 * ( T/100 ) + b3 * ( ( T/100 )**2 ) )

    c = exp(lnC)
    return c

def satO2(S,T):
    """
    USAGE:  satO2 = _satO2(S, T)

    DESCRIPTION:
    Solubility (satuaration) of Oxygen (O2) in sea water.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78)]
    T = temperature [:math:`^\\circ`  C (ITS-68)]

    OUTPUT:
    satO2 = solubility of O2  [ml/l]

    Author
    ------
    Phil Morgan 97-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Weiss, R. F. 1970
    "The solubility of nitrogen, oxygen and argon in water and seawater."
    Deap-Sea Research., 1970, Vol 17, pp721-735.

    CALLER: general purpose
    CALLEE:

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    # convert T to Kelvin
    T = 273.15 + T * 1.00024

    # constants for Eqn (4) of Weiss 1970
    a1 = -173.4292
    a2 =  249.6339
    a3 =  143.3483
    a4 =  -21.8492
    b1 =   -0.033096
    b2 =    0.014259
    b3 =   -0.0017000

    # Eqn (4) of Weiss 1970
    lnC = a1 + a2 * ( 100/T ) + a3 * log( T/100 ) + a4 * ( T/100 ) + \
          S * ( b1 + b2 * ( T/100 ) + b3 * ( ( T/100 )**2 ) )

    c = exp(lnC)
    return c

def dens0(S,T):
    """
    USAGE:  dens0 = dens0(S, T)

    DESCRIPTION:
    Density of Sea Water at atmospheric pressure using
    UNESCO 1983 (EOS 1980) polynomial.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78)]
    T = temperature [:math:`^\\circ`  C (ITS-90)]

    OUTPUT:
    dens0 = density  [kg/m^3] of salt water with properties S,T,
            P=0 (0 db gauge pressure)

    Author
    ------
    Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Millero, F.J. and  Poisson, A.
    International one-atmosphere equation of state of seawater.
    Deep-Sea Res. 1981. Vol28A(6) pp625-629.

    CALLER: general purpose, dens
    CALLEE: smow

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    T68 = T * 1.00024

    #     UNESCO 1983 eqn(13) p17
    b0 =  8.24493e-1
    b1 = -4.0899e-3
    b2 =  7.6438e-5
    b3 = -8.2467e-7
    b4 =  5.3875e-9

    c0 = -5.72466e-3
    c1 =  1.0227e-4
    c2 = -1.6546e-6

    d0 = 4.8314e-4
    dens = smow(T) + ( b0 + ( b1 + ( b2 + ( b3 + b4 * T68 ) * T68 ) * T68 ) * T68 ) * \
    S + ( c0 + ( c1 + c2 * T68 ) * T68 ) * S * (S)**0.5 + d0 * S**2
    return dens

def smow(T):
    """
    USAGE:  dens = smow(T)

    DESCRIPTION:
    Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980.

    INPUT:
    T = temperature [:math:`^\\circ`  C (ITS-90)]

    OUTPUT:
    dens = density  [kg/m^3]

    Author
    ------
    Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
    UNESCO 1983 p17  Eqn(14)

    Millero, F.J & Poisson, A.
    INternational one-atmosphere equation of state for seawater.
    Deep-Sea Research Vol28A No.6. 1981 625-629.    Eqn (6)

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    a0 = 999.842594
    a1 =   6.793952e-2
    a2 =  -9.095290e-3
    a3 =   1.001685e-4
    a4 =  -1.120083e-6
    a5 =   6.536332e-9

    T68  = T * 1.00024
    dens = a0 + ( a1 + ( a2 + ( a3 + ( a4 + a5 * T68 ) * T68 ) * T68 ) * T68 ) * T68
    return dens

def seck(S, T, P=0):
    """
    USAGE:  dens = seck(S, T, P)

    DESCRIPTION:
    Secant Bulk Modulus (K) of Sea Water using Equation of state 1980.
    UNESCO polynomial implementation.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78) ]
    T = temperature [:math:`^\\circ`  C (ITS-90)]
    P = pressure    [db]

    OUTPUT:
    K = Secant Bulk Modulus  [bars]

    Author
    ------
    Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
    Eqn.(15) p.18

    Millero, F.J. and  Poisson, A.
    International one-atmosphere equation of state of seawater.
    Deep-Sea Res. 1981. Vol28A(6) pp625-629.

    CALLER: dens
    CALLEE: none

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    # COMPUTE COMPRESSION TERMS
    P   = P/10.0 # convert from db to atmospheric pressure units
    T68 = T * 1.00024

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

    # SEA WATER TERMS OF SECANT BULK MODULUS AT ATMOS PRESSURE.
    j0 = 1.91075E-4

    i2 = -1.6078E-6
    i1 = -1.0981E-5
    i0 =  2.2838E-3

    SR = (S)**0.5

    A  = AW + ( i0 + ( i1 + i2 * T68 ) * T68 + j0 * SR ) * S

    m2 =  9.1697E-10
    m1 =  2.0816E-8
    m0 = -9.9348E-7

    B = BW + ( m0 + ( m1 + m2 * T68 ) * T68 ) * S # eqn 18

    f3 =  -6.1670E-5
    f2 =   1.09987E-2
    f1 =  -0.603459
    f0 =  54.6746

    g2 = -5.3009E-4
    g1 =  1.6483E-2
    g0 =  7.944E-2

    K0 = KW + ( f0 + ( f1 + ( f2 + f3 * T68 ) * T68 ) * T68 \
            + ( g0 + ( g1 + g2 * T68 ) * T68 ) * SR ) * S # eqn 16

    K = K0 + ( A + B * P ) * P # eqn 15
    return K

def dens(S, T, P):
    """
    USAGE:  dens = dens(S, T, P)

    DESCRIPTION:
    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78)]
    T = temperature [:math:`^\\circ`  C (ITS-90)]`
    P = pressure    [db]

    OUTPUT:
    dens = density  [kg/m^3]

    Author
    ------
    Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    Millero, F.J., Chen, C.T., Bradshaw, A., and Schleicher, K.
    " A new high pressure equation of state for seawater"
    Deap-Sea Research., 1980, Vol27A, pp255-264.

    CALLER: general purpose
    CALLEE: dens0 seck

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    #UNESCO 1983. eqn.7  p.15
    densP0 = dens0(S, T)
    K      = seck(S, T, P)
    P      = P / 10.0  #convert from db to atm pressure units
    dens   = densP0 / ( 1-P / K )
    return dens

def pden(S, T, P, PR=0):
    """
    USAGE:  pden = pden(S, T, P, PR)

    DESCRIPTION:
    Calculates potential density of water mass relative to the specified
    reference pressure by pden = dens(S, ptmp, PR).

    INPUT:  (all must have same dimensions)
    S  = salinity    [psu      (PSS-78) ]
    T  = temperature [:math:`^\\circ`  C (ITS-90)]
    P  = pressure    [db]
    PR = Reference pressure  [db]
         default = 0

    OUTPUT:
    pden = Potential denisty relative to the ref. pressure [kg/m^3]

    Author
    ------
    Phil Morgan 1992/04/06, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    A.E. Gill 1982. p.54
    "Atmosphere-Ocean Dynamics"
    Academic Press: New York.  ISBN: 0-12-283522-0

    CALLER:  general purpose
    CALLEE:  ptmp dens

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    PT   = ptmp(S, T, P, PR)
    pden = dens(S, PT, PR)
    return pden

def svan(S, T, P=0):
    """
    USAGE:  svan = svan(S, T, P)

    DESCRIPTION:
    Specific Volume Anomaly calculated as
    svan = 1/dens(s, t, p) - 1/dens(35, 0, p).
    Note that it is often quoted in literature as 1e8*units.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78) ]
    T = temperature [:math:`^\\circ`  C (ITS-90)]
    P = Pressure    [db]

    OUTPUT:
    svan = Specific Volume Anomaly  [m^3 kg^-1]

    Author
    ------
    Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonoff, N.P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
    Eqn (9) p.15.

    S. Pond & G.Pickard  2nd Edition 1986
    Introductory Dynamical Oceanogrpahy
    Pergamon Press Sydney.  ISBN 0-08-028728-X

    CALLER: general purpose
    CALLEE: dens

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    svan = 1/dens( S, T, P ) - 1/dens( 35, 0.0, P )

    return svan

def gpan(S, T, P):
    """
    USAGE:  gpan = gpan(S, T, P)

    DESCRIPTION:
    Geopotential Anomaly calculated as the integral of svan from the
    the sea surface to the bottom. Thus RELATIVE TO SEA SURFACE.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78)]
    T = temperature [:math:`^\\circ`  C (ITS-90)]
    P = Pressure    [db]

    OUTPUT:
    gpan = Geopotential Anomaly  [m^3 kg^-1 Pa == m^2 s^-2 == J kg^-1]

    Author
    ------
    Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    S. Pond & G.Pickard  2nd Edition 1986
    Introductory Dynamical Oceanogrpahy
    Pergamon Press Sydney.  ISBN 0-08-028728-X

    Note that older literature may use units of "dynamic decimeter" for above.

    Adapted method from Pond and Pickard (p76) to calc gpan rel to sea
    surface whereas P&P calculated relative to the deepest common depth.

    CALLER: general purpose
    CALLEE: svan

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    db2Pascal  = 1e4

    if P.ndim == 1:
        m   = P.size
    else:
        m,n = P.shape

    svn        = svan(S, T, P)
    mean_svan  = 0.5 * ( svn[1:m,:] + svn[0:m-1,:] )

    if n == 1:
        top = svn[0,0] * P[0,0] * db2Pascal
    else:
        top = svn[0,:] * P[0,:] * db2Pascal

    delta_ga   = ( mean_svan * diff(P,axis=0) ) * db2Pascal
    ga         = cumsum( vstack( ( top, delta_ga ) ),  axis=0 ) # TODO: I do not remember why?

    return ga

def gvel(ga, lon, lat):
    """
    USAGE:  vel = gvel(ga, lat, lon)

    DESCRIPTION:
    Calculates geostrophic velocity given the geopotential anomaly
    and position of each station.

    INPUT:
    ga   = geopotential anomaly relative to the sea surface.
            dim(mxnstations)
    lon  = longitude of each station (+ve = E, -ve = W) [-180..+180]
    lat  = latitude  of each station (+ve = N, -ve = S) [ -90.. +90]

    OUTPUT:
    vel  = geostrophic velocity RELATIVE to the sea surface.
            dim(m,nstations-1)

    Author
    -----
    Phil Morgan   1992/03/26  (morgan@ml.csiro.au)

    References
    ----------
    S. Pond & G.Pickard  2nd Edition 1986
    Introductory Dynamical Oceanography
    Pergamon Press Sydney.  ISBN 0-08-028728-X
    Equation 8.9A p73  Pond & Pickard

    CALLER:   general purpose
    CALLEE:   dist

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
    """

    # You may replace the call to dist if you have
    # a more appropriate distance routine.
    dista, phase = dist(lon,lat) # TODO: if I implement keyword for NM/KM this must be changed!
    distm = 1000.0 * dista # meters to 'km'
    m,n   = ga.shape
    f     = cor( ( lat[0:n-1] + lat[1:n] )/2 )
    lf    = f * distm
    #LF    = lf[(ones((m,1)),:] # TODO: not sure why?
    vel   = -( ga[:,1:n] - ga[:,0:n-1] ) / LF

    return vel

def gvel2(ga, dist, lat):
    """
    USAGE:  vel = gvel(ga, dist, lat)

    DESCRIPTION:
    Calculates geostrophic velocity given the geopotential anomaly
    and position of each station.

    INPUT:
    ga   = geopotential anomaly relative to the sea surface.
            dim(mxnstations)
    dist = distance between stations, in kilometers
    lat  = lat to use for constant f determination

    OUTPUT:
    vel  = geostrophic velocity RELATIVE to the sea surface.
            dim(m,nstations-1)

    Author
    ------
    Phil Morgan   1992/03/26  (morgan@ml.csiro.au)

    References
    ----------
    S. Pond & G.Pickard  2nd Edition 1986
    Introductory Dynamical Oceanogrpahy
    Pergamon Press Sydney.  ISBN 0-08-028728-X
    Equation 8.9A p73  Pond & Pickard

    CALLER:   general purpose

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
    """
    distm = 1000.0 * dista # meters to 'km'
    m,n   = ga.shape
    f     = cor( ( lat[0:n-1] + lat[1:n] )/2 )
    lf    = f * distm
    vel   = -( ga[:,1:n] - ga[:,0:n-1] ) / LF

    return vel

def cp(S, T, P):
    """
    USAGE: cp = cp(S, T, P)

    DESCRIPTION:
    Heat Capacity of Sea Water using UNESCO 1983 polynomial.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78)]
    T = temperature [:math:`^\\circ`  C (ITS-90)]
    P = pressure    [db]

    OUTPUT:
    cp = Specific Heat Capacity  [J kg^-1 C^-1]

    Author
    ------
    Phil Morgan, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

    CALLER: general purpose
    CALLEE: none

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """

    P = P/10 # to convert db to Bar as used in Unesco routines
    T68 = T * 1.00024

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
            ( a0 + a1 * T68 + a2 * T68**2 ) * S + \
            ( b0 + b1 * T68 + b2 * T68**2 ) * S *(S)**0.5

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

    del_Cp0t0 = ( ( ( ( ( c3 * T68 + c2 ) * T68 + c1 ) * T68 + c0 ) * P + \
                ( ( ( ( b4 * T68 + b3 ) * T68 + b2 ) * T68 + b1 ) * T68 + b0 ) ) * P + \
                ( ( ( ( a4 * T68 + a3 ) * T68 + a2 ) * T68 + a1 ) * T68 + a0 ) ) * P

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

    S3_2 = S * (S)**0.5

    """
    del_Cpstp = [ ( ( ( ( d4 * T68 + d3 ) * T68 + d2 ) * T68 + d1 ) * T68 + d0 ) * S + \
                ( ( e2 * T68 + e1 ) * T68 + e0 ) * S3_2 ] * P + \
                [ ( ( ( f3 * T68 + f2 ) * T68 + f1 ) * T68 + f0 ) * S + \
                g0 * S3_2 ] * P**2 + \
                [ ( ( h2 * T68 + h1 ) * T68 + h0 ) * S + \
                j1 * T68 * S3_2 ] * P**3
    """

    del_Cpstp = ( ( ( ( ( d4 * T68 + d3 ) * T68 + d2 ) * T68 + d1 ) * T68 + d0 ) * S + \
                ( ( e2 * T68 + e1 ) * T68 + e0 ) * S3_2 ) * P + \
                ( ( ( ( f3 * T68 + f2 ) * T68 + f1 ) * T68 + f0 ) * S + \
                 g0 * S3_2 ) * P**2 + \
                ( ( ( h2 * T68 + h1 ) * T68 + h0 ) * S + \
                j1 * T68 * S3_2 ) * P**3

    cp = Cpst0 + del_Cp0t0 + del_Cpstp

    return cp

def ptmp(S, T, P, PR=0):
    """
    USAGE:  ptmp = sw_ptmp(S,T,P,PR)

    DESCRIPTION:
    Calculates potential temperature as per UNESCO 1983 report.

    INPUT:  (all must have same dimensions)
    S  = salinity    [psu      (PSS-78) ]
    T  = temperature [:math:`^\\circ`  C (ITS-90)]
    P  = pressure    [db]
    PR = Reference pressure  [db]
         default = 0

    OUTPUT:
    ptmp = Potential temperature relative to PR [:math:`^\\circ`  C (ITS-90)]

    Author
    ------
    Phil Morgan 92-04-06, Lindsay Pender (Lindsay.Pender@csiro.au)

    DISCLAIMER:
    This software is provided "as is" without warranty of any kind.
    See the file sw_copy.m for conditions of use and licence.

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
    Eqn.(31) p.39

    Bryden, H. 1973.
    "New Polynomials for thermal expansion, adiabatic temperature gradient
    and potential temperature of sea water."
    DEEP-SEA RES., 1973, Vol20,401-408.

    CALLER:  general purpose
    CALLEE:  adtg

    Modifications
    -------------
    99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """
    # theta1
    del_P  = PR - P
    del_th = del_P * adtg(S, T, P)
    th     = T * 1.00024 + 0.5 * del_th
    q      = del_th

    # theta2
    del_th = del_P * adtg(S, th/1.00024, P + 0.5 * del_P )
    th     = th + ( 1 - 1/(2)**00.5 ) * ( del_th - q )
    q      = ( 2 - (2)**0.5 ) * del_th + ( -2 + 3/(2)**0.5 ) * q

    # theta3
    del_th = del_P * adtg( S, th/1.00024, P + 0.5 * del_P )
    th     = th + ( 1 + 1/(2)**0.5 ) * ( del_th - q )
    q      = ( 2 + (2)**0.5 ) * del_th + ( -2 -3/(2)**0.5 ) * q

    # theta4
    del_th = del_P * adtg( S, th/1.00024, P + del_P )
    PT     = ( th + ( del_th - 2 * q ) / 6 ) / 1.00024
    return PT

def temp(S, PTMP, P, PR):
    """
    USAGE:  temp = temp(S, PTMP, P, PR)

    DESCRIPTION:
    Calculates temperature from potential temperature at the reference
    pressure PR and in-situ pressure P.

    INPUT:  (all must have same dimensions)
    S     = salinity              [psu      (PSS-78) ]
    PTMP  = potential temperature [:math:`^\\circ`  C (ITS-90)]
    P     = pressure              [db]
    PR    = Reference pressure    [db]

    OUTPUT:
    temp = temperature [:math:`^\\circ`  C (ITS-90)]

    Author
    ------
    Phil Morgan 92-04-06, Lindsay Pender (Lindsay.Pender@csiro.au)

    References
    ----------
    Fofonoff, P. and Millard, R.C. Jr
    Unesco 1983. Algorithms for computation of fundamental properties of
    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
    Eqn.(31) p.39

    Bryden, H. 1973.
    "New Polynomials for thermal expansion, adiabatic temperature gradient
    and potential temperature of sea water."
    DEEP-SEA RES., 1973, Vol20,401-408.

    CALLER:  general purpose
    CALLEE:  ptmp

    Modifications
    -------------
    03-12-12. Lindsay Pender, Converted to ITS-90.
    10-01-14. Filipe Fernandes, Python translation.
    """
    # CARRY OUT INVERSE CALCULATION BY SWAPPING P0 & PR
    T = ptmp(S, PTMP, PR, P);

    return T

def swvel(lenth, depth):
    """
    USAGE:  speed = swvel(lenth, depth)

    DESCRIPTION:
    Calculates surface wave velocity.

    INPUT:  (all must have same dimensions)
    lenth = wave length
    depth = water depth [metres]

    OUTPUT:
    speed   = Surface wave speed (m/s)

    Author
    ------
    Lindsay Pender 2005

    CALLER:  general purpose
    CALLEE:  none

    Modifications
    -------------
    10-01-14. Filipe Fernandes, Python translation.
    """
    g = 9.8 # TODO: use grav function
    k = 2.0 * pi / lenth
    speed = (g * tanh(k * depth) / k)**0.5
    return speed