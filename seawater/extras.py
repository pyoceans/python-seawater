# -*- coding: utf-8 -*-

from numpy import sqrt, abs, arange, sign, diff, pi # fixme
from seawater import dens, pden, cor  # fixme

def sigmat(S, T, P):
    """
    USAGE:  sgmt = sigmat(S, T, P)

    DESCRIPTION:
    Sigma-T is the remainder of subtracting 1000 kg m**-3 from the density of a sea water sample at atmospheric pressure

    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78)]
    T = temperature [degree C (ITS-90)]
    P = pressure    [db]

    OUTPUT:
    sgmt = density  [kg / m**3]

    AUTHOR:  Filipe Fernandes, 2010

    MODIFICATIONS:
    10-01-26. Filipe Fernandes, first version.
    """
    sgmt = dens(S, T, P) - 1000.0
    return sgmt

def sigmatheta(S, T, P):
    """
    USAGE:  sgmte = sigmatheta(S, T, P)

    DESCRIPTION:
    Sigma-Theta is a measure of the density of ocean water where the quantity sigma-t is calculated using the potential temperature (theta) rather than the in situ temperature

    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78)]
    T = temperature [degree C (ITS-90)]
    P = pressure    [db]

    OUTPUT:
    sgmte = density  [kg / m**3]

    AUTHOR:  Filipe Fernandes, 2010

    MODIFICATIONS:
    10-01-26. Filipe Fernandes, first version.
    """
    sgmte = pden(S, T, P) - 1000.0
    return sgmte

def N(bvfr2):
    """
    USAGE:  bvfr = N(Nsqrd)

    DESCRIPTION:
    Buoyancy frequency is the frequency with which a parcel or particle of fluid displaced a small vertical distance from its equilibrium position in a stable environment will oscillate. It will oscillate in simple harmonic motion with an angular frequency defined by

                -g      d(pdens)
         N =  ( ----- x -------- )**1/2
                pdens     d(z)

    INPUT:
    Nsqrd = frequency [s**-2]

    OUTPUT:
    Nsqrd = frequency [s**-1]

    AUTHOR:  Filipe Fernandes, 2010

    MODIFICATIONS:
    10-01-28. Filipe Fernandes, first version.
    """
    bvfr  = sqrt( abs( bvfr2 ) ) * sign( bvfr2 )
    return bvfr

def shear(P, U, V=0):
    """
    USAGE:  shr = shear(P, U, V)

    DESCRIPTION:
    Calculates the vertical shear from U, V velocity section

    INPUT:
    P     = pressure    [db]
    U, V  = [m / s]

    OUTPUT:
    shr   = frequency [s**-1]
    p_ave = Mid pressure between P grid (M-1xN)  [db]

    AUTHOR:  Filipe Fernandes, 2010

    MODIFICATIONS:
    10-01-28. Filipe Fernandes, first version.
    """

    m,n      = P.shape #fixme *check where depth increases to find dimension
    iup      = arange(0,m-1)
    ilo      = arange(1,m)
    p_ave    = ( P[iup,:] + P[ilo,:] )/2. # fixme
    vel      = sqrt( U**2 + V**2 )
    diff_vel = diff( vel, axis=0 )
    diff_z   = diff(   P, axis=0 ) # fixme to Z ?
    shr      = diff_vel / diff_z

    return shr, p_ave

def richnumb(N, S):
    """
    USAGE:  ri = richnumb( N, S)

    DESCRIPTION:
    Calculates  the ratio of buoyancy to inertial forces which measures the stability of a fluid layer.
    this functions computes the gradient Richardson number in the form of:

    Ri =  N**2 / S**2

    Representing a dimensionless number that expresses the ratio of the energy extracted by buoyancy forces to the energy gained from the shear of the large-scale velocity field.

    INPUT:
    N   = frequency [s**-1]
    shr = frequency [s**-1]

    OUTPUT:
    ri  = non-dimensional

    AUTHOR:  Filipe Fernandes, 2010

    MODIFICATIONS:
    10-01-28. Filipe Fernandes, first version.
    """

    N2 = N**2 * sign( N ) # fixme
    S2 = S**2
    ri   =  N2 / S2

    return ri

def intertialfrq(lat):
    """
    USAGE:  Ti = intertialfrq(lat)

    DESCRIPTION:
    The frequency "f" of rotation of inertial waves.

    Calculates local inertial period as:

    2*pi/cor(lat)/60/60


    INPUT:
    LAT = Latitude in decimal degress north [-90..+90]

    OUTPUT:
    Period in hours

    AUTHOR:  Filipe Fernandes, 2010

    MODIFICATIONS:
    10-04-29. Filipe Fernandes, first version.
    """

    Ti = 2*pi / cor(lat)/3600

    return Ti