# -*- coding: utf-8 -*-

def sigma():
    """
    USAGE:  sgm = sigma(S, T, P)

    DESCRIPTION:
    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial in
    Sigma Units or dens(S,T,P) - 1000.0

    INPUT:  (all must have same dimensions)
    S = salinity    [psu      (PSS-78)]
    T = temperature [degree C (ITS-90)]
    P = pressure    [db]

    OUTPUT:
    sgm = density  [kg / m**3]

    AUTHOR:  Filipe Fernandes, 2010

    MODIFICATIONS:
    10-01-26. Filipe Fernandes, first version.
    """
    sgm = dens(S,T,P) - 1000.0
    return sgm

def N(bvfr):
    """
    USAGE:  bvfr = N(Nsqrd)

    DESCRIPTION:
    Calculates Brunt-Vaisala Frequency (N)
    obs: Not the N square!!!

    INPUT:
    Nsqrd = frequency [s**-2]

    OUTPUT:
    Nsqrd = frequency [s**-1]

    AUTHOR:  Filipe Fernandes, 2010

    MODIFICATIONS:
    10-01-28. Filipe Fernandes, first version.
    """
    bvfr  = np.sqrt( np.abs( bvfr2 ) ) * np.sign( bvfr2 )
    return bvfr

def shear(U, V=0, P):
    """
    USAGE:  shr = shear(U, V, P)

    DESCRIPTION:
    Calculates Shear from velocity data.

    INPUT:
    U, V  = [m / s]
    P   = pressure    [db]

    OUTPUT:
    shr = frequency [s**-1]

    AUTHOR:  Filipe Fernandes, 2010

    MODIFICATIONS:
    10-01-28. Filipe Fernandes, first version.
    """

    m,n      = P.shape #fixme *check where depth increases to find dimension
    iup      = arange(0,m-1)
    ilo      = arange(1,m)
    p_ave    = ( P[iup,:] + P[ilo,:] )/2.
    vel      = np.sqrt( U**2 + V**2 )
    diff_vel = diff( vel, axis=0 )
    diff_z   = diff(   Z, axis=0 )
    shr      = diff_vel / diff_z

    return shr

def richnumb(shr, N):
    """
    USAGE:  ri = richnumb(shr, N)

    DESCRIPTION:
    Calculates gradient Richardson number in the form of:

    Ri =  shr**2 / N**2


    INPUT:
    shr = frequency [s**-1]
    N   = frequency [s**-1]

    OUTPUT:
    ri  = non-dimensional

    AUTHOR:  Filipe Fernandes, 2010

    MODIFICATIONS:
    10-01-28. Filipe Fernandes, first version.
    """
    shr2 = shr**2
    N2   =   N**2
    ri   = shr2 / N2
    
    return ri
