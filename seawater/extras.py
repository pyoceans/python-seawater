# -*- coding: utf-8 -*-

# extras
def sigma():
    """
    USAGE:  sgm = sigma(S, T, P)

    DESCRIPTION:
    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.

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

#The following functions are from:
#http://www.imr.no/~bjorn/python/seawater/index.html

def drhodt(S,T,P=0):
    """
    USAGE: drhodt(S, T, [P])

    DESCRIPTION:
    Compute temperature derivative of density

    INPUT:
        S = Salinity,     [PSS-78]
        T = Temperature,  [degC]
        P = Pressure,     [dbar]
    P is optional, with a default value = zero

    OUTPUT:
        Temperature derivative of density  [kg /(K m**3)]
    """

    a1 =  6.793952e-2
    a2 = -1.819058e-2
    a3 =  3.005055e-4
    a4 = -4.480332e-6
    a5 =  3.268166e-8

    b1 = -4.0899e-3
    b2 =  1.52876e-4
    b3 = -2.47401e-6
    b4 =  2.155e-8 

    c1 =  1.0227e-4
    c2 = -3.3092e-6

    e1 =  148.4206
    e2 = -4.65421
    e3 =  4.081431e-2
    e4 = -2.0621152e-4

    f1 = -0.603459
    f2 =  2.19974e-2
    f3 = -1.8501e-4
  
    g1 =  1.6483e-2
    g2 = -1.06018e-3

    h1 =  1.43713e-3
    h2 =  2.32184e-4
    h3 = -1.733715e-6
 
    i1 = -1.0981e-5
    i2 = -3.2156e-6

    k1 = -6.12293e-6
    k2 =  1.05574e-7
  
    m1 = 2.0816e-8
    m2 = 1.83394e-9

    P = P/10.0

    DSMOV = a1 + ( a2 + ( a3 + ( a4 + a5 * T ) * T ) * T ) * T
    DRHO0 = DSMOV + ( b1 + ( b2 + ( b3 + b4 * T ) * T ) * T ) * S \
                  + ( c1 + c2 * T ) * S**1.5

    DAW = h1 + ( h2 + h3 * T )
    DA  = DAW + ( i1 + i2 * T ) * S

    DBW = k1 + k2 * T
    DB  = DBW + ( m1 + m2 * T ) * S

    DKW = e1 + ( e2 + ( e3 + e4 * T ) * T ) * T
    DK0 = DKW + ( f1 + ( f2 + f3 * T ) * T ) * S + ( g1 + g2 * T ) * S**1.5
    DK  = DK0 + ( DA + DB * P ) * P

    K      = seck(S, T, P)
    RHO0   = dens0(S, T) 
    denom  = 1.0 - P/K
    DRHODT = ( DRHO0 * denom - RHO0 * P * DK/( K * K ) ) / ( denom * denom )
    return DRHODT

def drhods(S, T, P=0):
    """
    USAGE: drhodt(S, T, [P])
    
    DESCRIPTION: Compute salinity derivative of density

    INPUT:
        S = Salinity,     [PSS-78]
        T = Temperature,  [degC]
        P = Pressure,     [dbar]
    P is optional, with a default value = zero

    OUTPUT:
        Salinity derivative of density [kg/m**3]
    """
    
    b0 =  8.24493e-1
    b1 = -4.0899e-3
    b2 =  7.6438e-5
    b3 = -8.2467e-7
    b4 =  5.3875e-9

    c0 = -5.72466e-3
    c1 =  1.0227e-4
    c2 = -1.6546e-6

    d0 =  9.6628e-4

    f0 = 54.6746
    f1 = -0.603459
    f2 =  1.09987e-2
    f3 = -6.1670e-5

    g0 =  7.944e-2
    g1 =  1.6483e-2
    g2 = -5.3009e-4

    i0 =  2.2838e-3
    i1 = -1.0981e-5
    i2 = -1.6078e-6

    j0 =  2.866125e-4

    m0 = -9.9348e-7
    m1 =  2.0816e-8
    m2 =  9.1697e-10

    P = 0.1*P # Convert to bar

    DRHO0 = b0 + T * ( b1 + T * ( b2 + T * ( b3 + T * b4 ) ) ) +   \
            1.5 * S**0.5 * ( c0 + T * ( c1 + T * c2 ) ) + S * d0
    DK0   = f0 + T * ( f1 + T * ( f2 + T * f3 ) ) + \
            1.5 * S**0.5 * ( g0 + T * ( g1 + T * g2 ) )
    DA    = i0 + T * ( i1 + T * i2 ) + j0 * S**0.5
    DB    = m0 + T * ( m1 + T * m2 )
    DK    = DK0 + P * ( DA + P * DB )
    RHO0  = dens0(S, T)
    K     = seck(S, T, P)
    denom = 1.0 - P/K
    DRHO  = ( DRHO0 * denom - RHO0 * P * DK / ( K * K ) ) / ( denom * denom )
    return DRHODS