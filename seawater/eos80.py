# -*- coding: utf-8 -*-


from __future__ import division, absolute_import

import numpy as np

from .constants import deg2rad, earth_radius
from .library import T90conv, T68conv, salrt, salrp, sals, seck, smow


__all__ = ['adtg',
           'alpha',
           'aonb',
           'beta',
           'dpth',
           'g',
           'salt',
           'fp',
           'svel',
           'pres',
           'dens0',
           'dens',
           'pden',
           'cp',
           'ptmp',
           'temp']


def adtg(s, t, p):
    """
    Calculates adiabatic temperature gradient as per UNESCO 1983 routines.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [℃ (ITS-90)]
    p : array_like
        pressure [db]

    Returns
    -------
    adtg : array_like
           adiabatic temperature gradient [℃ db :sup:`-1`]

    Examples
    --------
    >>> # Data from UNESCO 1983 p45.
    >>> import seawater as sw
    >>> from seawater.library import T90conv
    >>> t = T90conv([[ 0,  0,  0,  0,  0,  0],
    ...                 [10, 10, 10, 10, 10, 10],
    ...                 [20, 20, 20, 20, 20, 20],
    ...                 [30, 30, 30, 30, 30, 30],
    ...                 [40, 40, 40, 40, 40, 40]])
    >>> s = [[25, 25, 25, 35, 35, 35],
    ...      [25, 25, 25, 35, 35, 35],
    ...      [25, 25, 25, 35, 35, 35],
    ...      [25, 25, 25, 35, 35, 35],
    ...      [25, 25, 25, 35, 35, 35]]
    >>> p = [0, 5000, 10000, 0, 5000, 10000]
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
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Bryden, H. 1973. New Polynomials for thermal expansion, adiabatic
       temperature gradient and potential temperature of sea water. Deep-Sea
       Res. Vol20,401-408. doi:10.1016/0011-7471(73)90063-6

    """
    s, t, p = list(map(np.asanyarray, (s, t, p)))

    T68 = T68conv(t)

    a = [3.5803e-5, 8.5258e-6, -6.836e-8, 6.6228e-10]
    b = [1.8932e-6, -4.2393e-8]
    c = [1.8741e-8, -6.7795e-10, 8.733e-12, -5.4481e-14]
    d = [-1.1351e-10, 2.7759e-12]
    e = [-4.6206e-13, 1.8676e-14, -2.1687e-16]
    return (a[0] + (a[1] + (a[2] + a[3] * T68) * T68) * T68 +
            (b[0] + b[1] * T68) * (s - 35) +
            ((c[0] + (c[1] + (c[2] + c[3] * T68) * T68) * T68) +
             (d[0] + d[1] * T68) * (s - 35)) * p +
            (e[0] + (e[1] + e[2] * T68) * T68) * p * p)


def alpha(s, t, p, pt=False):
    """
    Calculate the thermal expansion coefficient.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature or potential temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].
    pt : bool
         True if temperature is potential, default is False

    Returns
    -------
    alpha : array_like
            thermal expansion coeff :math:`\alpha` [℃ :sup:`-1`]

    Examples
    --------
    >>> # Data from McDougall 1987
    >>> import seawater as sw
    >>> s, t, p = 40, 10, 4000
    >>> sw.alpha(s, t, p, pt=True)
    0.00025061316481624323

    References
    ----------
    .. [1] McDougall, Trevor J., 1987: Neutral Surfaces. J. Phys. Oceanogr.,
       17, 1950-1964. doi: 10.1175/1520-0485(1987)017<1950:NS>2.0.CO;2

    """
    s, t, p, pt = list(map(np.asanyarray, (s, t, p, pt)))
    return aonb(s, t, p, pt) * beta(s, t, p, pt)


def aonb(s, t, p, pt=False):
    """
    Calculate :math:`\alpha/\beta`.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature or potential temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].
    pt : bool
         True if temperature is potential, default is False

    Returns
    -------
    aonb : array_like
           :math:`\alpha/\beta` [psu ℃ :sup:`-1`]

    Examples
    --------
    >>> # Data from McDougall 1987.
    >>> import seawater as sw
    >>> s, t, p = 40, 10, 4000
    >>> sw.aonb(s, t, p, pt=True)
    0.347650567047807

    References
    ----------
    .. [1] McDougall, Trevor J., 1987: Neutral Surfaces. J. Phys. Oceanogr.,
       17, 1950-1964. doi: 10.1175/1520-0485(1987)017<1950:NS>2.0.CO;2

    """

    # Ensure we use ptmp in calculations.
    s, t, p, pt = list(map(np.asanyarray, (s, t, p, pt)))

    if not pt:
        t = ptmp(s, t, p, 0)  # Now we have ptmp.

    p = np.float_(p)
    t = T68conv(t)

    c1 = np.array([-0.255019e-7, 0.298357e-5, -0.203814e-3,
                   0.170907e-1, 0.665157e-1])
    c2 = np.array([-0.846960e-4, 0.378110e-2])
    c2a = np.array([-0.251520e-11, -0.164759e-6, 0.0])
    c3 = -0.678662e-5
    c4 = np.array([0.791325e-8, -0.933746e-6, 0.380374e-4])
    c5 = 0.512857e-12
    c6 = -0.302285e-13

    # Now calculate the thermal expansion saline contraction ratio aonb.
    sm35 = s - 35.0
    return (np.polyval(c1, t) + sm35 *
            (np.polyval(c2, t) + np.polyval(c2a, p)) +
            sm35 ** 2 * c3 + p * np.polyval(c4, t) +
            c5 * (p ** 2) * (t ** 2) + c6 * p ** 3)


def beta(s, t, p, pt=False):
    """
    Calculate the saline contraction coefficient :math:`\beta` as defined
    by T.J. McDougall.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature or potential temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].
    pt : bool
         True if temperature is potential, default is False

    Returns
    -------
    beta : array_like
           saline Contraction Coefficient [psu :sup:`-1`]

    Examples
    --------
    >>> # Data from McDougall 1987
    >>> import seawater as sw
    >>> s, t, p = 40, 10, 4000
    >>> sw.beta(s, t, p, pt=True)
    0.00072087661741618932

    References
    ----------
    .. [1] McDougall, Trevor J., 1987: Neutral Surfaces. J. Phys. Oceanogr.,
       17, 1950-1964. doi: 10.1175/1520-0485(1987)017<1950:NS>2.0.CO;2

    """

    s, t, p, pt = list(map(np.asanyarray, (s, t, p, pt)))

    # Ensure we use ptmp in calculations
    if not pt:
        t = ptmp(s, t, p, 0)  # Now we have ptmp.

    t = T68conv(t)

    c1 = np.array([-0.415613e-9, 0.555579e-7, -0.301985e-5, 0.785567e-3])
    c2 = np.array([0.788212e-8, -0.356603e-6])
    c3 = np.array([-0.602281e-15, 0.408195e-10, 0.0])
    c4 = 0.515032e-8
    c5 = np.array([-0.213127e-11, 0.192867e-9, -0.121555e-7])
    c6 = np.array([-0.175379e-14, 0.176621e-12])
    c7 = 0.121551e-17

    # Now calculate the thermal expansion saline contraction ratio adb
    sm35 = s - 35
    return (np.polyval(c1, t) + sm35 *
            (np.polyval(c2, t) + np.polyval(c3, p)) +
            c4 * (sm35 ** 2) + p * np.polyval(c5, t) +
            (p ** 2) * np.polyval(c6, t) + c7 * (p ** 3))


def cp(s, t, p):
    """
    Heat Capacity of Sea Water using UNESCO 1983 polynomial.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    cp : array_like
         specific heat capacity [J kg :sup:`-1` C :sup:`-1`]

    Examples
    --------
    >>> # Data from Pond and Pickard Intro. Dyn. Oceanography 2nd ed. 1986.
    >>> import seawater as sw
    >>> from seawater.library import T90conv
    >>> t = T90conv([[0, 0, 0, 0, 0, 0],
    ...              [10, 10, 10, 10, 10, 10],
    ...              [20, 20, 20, 20, 20, 20],
    ...              [30, 30, 30, 30, 30, 30],
    ...             [40, 40, 40, 40, 40, 40]])
    >>> s = [[25, 25, 25, 35, 35, 35],
    ...      [25, 25, 25, 35, 35, 35],
    ...      [25, 25, 25, 35, 35, 35],
    ...      [25, 25, 25, 35, 35, 35],
    ...      [25, 25, 25, 35, 35, 35]]
    >>> p = [0, 5000, 10000, 0, 5000, 10000]
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
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp. Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    """

    s, t, p = list(map(np.asanyarray, (s, t, p)))

    p = p / 10.  # To convert [db] to [bar] as used in UNESCO routines.
    T68 = T68conv(t)

    # Eqn. 26 p.32.
    a = (-7.64357, 0.1072763, -1.38385e-3)
    b = (0.1770383, -4.07718e-3, 5.148e-5)
    c = (4217.4, -3.720283, 0.1412855, -2.654387e-3, 2.093236e-5)

    Cpst0 = ((((c[4] * T68 + c[3]) * T68 + c[2]) * T68 + c[1]) * T68 + c[0] +
             (a[0] + a[1] * T68 + a[2] * T68 ** 2) * s +
             (b[0] + b[1] * T68 + b[2] * T68 ** 2) * s * s ** 0.5)

    # Eqn. 28 p.33.
    a = (-4.9592e-1, 1.45747e-2, -3.13885e-4, 2.0357e-6, 1.7168e-8)
    b = (2.4931e-4, -1.08645e-5, 2.87533e-7, -4.0027e-9, 2.2956e-11)
    c = (-5.422e-8, 2.6380e-9, -6.5637e-11, 6.136e-13)

    del_Cp0t0 = ((((((c[3] * T68 + c[2]) * T68 + c[1]) * T68 + c[0]) * p +
                   ((((b[4] * T68 + b[3]) * T68 + b[2]) * T68 + b[1]) *
                    T68 + b[0])) * p + ((((a[4] * T68 + a[3]) * T68 + a[2]) *
                                         T68 + a[1]) * T68 + a[0])) * p)

    # Eqn 29 p.34.
    d = (4.9247e-3, -1.28315e-4, 9.802e-7, 2.5941e-8, -2.9179e-10)
    e = (-1.2331e-4, -1.517e-6, 3.122e-8)
    f = (-2.9558e-6, 1.17054e-7, -2.3905e-9, 1.8448e-11)
    g0 = 9.971e-8
    h = (5.540e-10, -1.7682e-11, 3.513e-13)
    j1 = -1.4300e-12

    S3_2 = s * s ** 0.5

    del_Cpstp = ((((((d[4] * T68 + d[3]) * T68 + d[2]) * T68 + d[1]) *
                   T68 + d[0]) * s + ((e[2] * T68 + e[1]) * T68 + e[0]) *
                  S3_2) * p +
                 ((((f[3] * T68 + f[2]) * T68 + f[1]) * T68 + f[0]) * s +
                  g0 * S3_2) * p ** 2 + (((h[2] * T68 + h[1]) * T68 + h[0]) *
                                         s + j1 * T68 * S3_2) * p ** 3)

    return Cpst0 + del_Cp0t0 + del_Cpstp


def dens0(s, t):
    """
    Density of Sea Water at atmospheric pressure.

    Parameters
    ----------
    s(p=0) : array_like
             salinity [psu (PSS-78)]
    t(p=0) : array_like
             temperature [℃ (ITS-90)]

    Returns
    -------
    dens0(s, t) : array_like
                  density  [kg m :sup:`3`] of salt water with properties
                  (s, t, p=0) 0 db gauge pressure

    Examples
    --------
    >>> # Data from UNESCO Tech. Paper in Marine Sci. No. 44, p22
    >>> import seawater as sw
    >>> from seawater.library import T90conv
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> sw.dens0(s, t)
    array([  999.842594  ,   999.842594  ,   995.65113374,   995.65113374,
            1028.10633141,  1028.10633141,  1021.72863949,  1021.72863949])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere
       equation of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
       doi:10.1016/0198-0149(81)90122-9

    """

    s, t = list(map(np.asanyarray, (s, t)))

    T68 = T68conv(t)

    # UNESCO 1983 Eqn.(13) p17.
    b = (8.24493e-1, -4.0899e-3, 7.6438e-5, -8.2467e-7, 5.3875e-9)
    c = (-5.72466e-3, 1.0227e-4, -1.6546e-6)
    d = 4.8314e-4
    return (smow(t) + (b[0] + (b[1] + (b[2] + (b[3] + b[4] * T68) * T68) *
            T68) * T68) * s + (c[0] + (c[1] + c[2] * T68) * T68) * s *
            s ** 0.5 + d * s ** 2)


def dens(s, t, p):
    """
    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    dens : array_like
           density  [kg m :sup:`3`]

    Examples
    --------
    >>> # Data from Unesco Tech. Paper in Marine Sci. No. 44, p22.
    >>> import seawater as sw
    >>> from seawater.library import T90conv
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> p = [0, 10000, 0, 10000, 0, 10000, 0, 10000]
    >>> sw.dens(s, t, p)
    array([  999.842594  ,  1045.33710972,   995.65113374,  1036.03148891,
            1028.10633141,  1070.95838408,  1021.72863949,  1060.55058771])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J., Chen, C.T., Bradshaw, A., and Schleicher, K. A new
       high pressure equation of state for seawater. Deap-Sea Research., 1980,
       Vol27A, pp255-264. doi:10.1016/0198-0149(80)90016-3

    """

    s, t, p = list(map(np.asanyarray, (s, t, p)))

    # UNESCO 1983. Eqn..7  p.15.
    densP0 = dens0(s, t)
    K = seck(s, t, p)
    p = p / 10.  # Convert from db to atm pressure units.
    return densP0 / (1 - p / K)


def dpth(p, lat):
    """
    Calculates depth in meters from pressure in dbars.

    Parameters
    ----------
    p : array_like
        pressure [db].
    lat : number or array_like
          latitude in decimal degrees north [-90..+90].

    Returns
    -------
    z : array_like
        depth [meters]

    Examples
    --------
    >>> # UNESCO 1983 data p30.
    >>> import seawater as sw
    >>> lat = [0, 30, 45, 90]
    >>> p = [[  500,   500,   500,  500],
    ...      [ 5000,  5000,  5000, 5000],
    ...      [10000, 10000, 10000, 10000]]
    >>> sw.dpth(p, lat)
    array([[  496.65299239,   495.99772917,   495.3427354 ,   494.03357499],
           [ 4915.04099112,  4908.55954332,  4902.08075214,  4889.13132561],
           [ 9725.47087508,  9712.6530721 ,  9699.84050403,  9674.23144056]])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    """
    p, lat = list(map(np.asanyarray, (p, lat)))

    # Eqn 25, p26.  UNESCO 1983.
    c = [9.72659, -2.2512e-5, 2.279e-10, -1.82e-15]
    gam_dash = 2.184e-6

    lat = abs(lat)
    X = np.sin(lat * deg2rad)
    X = X * X

    bot_line = (9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * X) * X) +
                gam_dash * 0.5 * p)
    top_line = (((c[3] * p + c[2]) * p + c[1]) * p + c[0]) * p
    return top_line / bot_line


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
        freezing point temperature [℃ (ITS-90)]

    Examples
    --------
    >>> # UNESCO DATA p.30.
    >>> import seawater as sw
    >>> s = [[5, 10, 15, 20, 25, 30, 35, 40],
    ...      [5, 10, 15, 20, 25, 30, 35, 40]]
    >>> p = [[ 0, 0, 0, 0, 0, 0, 0, 0],
    ...      [500, 500, 500, 500, 500, 500, 500, 500]]
    >>> sw.fp(s, p)
    array([[-0.27369757, -0.54232831, -0.81142026, -1.0829461 , -1.35804594,
            -1.63748903, -1.9218401 , -2.2115367 ],
           [-0.65010724, -0.91873798, -1.18782992, -1.45935577, -1.73445561,
            -2.01389869, -2.29824976, -2.58794636]])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    """

    s, p = list(map(np.asanyarray, (s, p)))

    # NOTE: P = P/10 # to convert db to Bar as used in UNESCO routines.
    # Eqn  p.29.
    a = [-0.0575, 1.710523e-3, -2.154996e-4]
    b = -7.53e-4
    return T90conv(a[0] * s + a[1] * s * s ** 0.5 + a[2] * s ** 2 + b * p)


def g(lat, z=0):
    """
    Calculates acceleration due to gravity as function of latitude.

    Parameters
    ----------
    lat : array_like
         latitude in decimal degrees north [-90..+90].

    z : number or array_like. Default z = 0
        height in meters (+ve above sea surface, -ve below).

    Returns
    -------
    g : array_like
        gravity [m s :sup:`-2`]

    Examples
    --------
    >>> import seawater as sw
    >>> sw.g(45, z=0)
    9.8061898752053995

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] A.E. Gill 1982. p.54  Eqn. 3.7.15 "Atmosphere-Ocean Dynamics"
       Academic Press: New York. ISBN: 0-12-283522-0

    """

    lat, z = list(map(np.asanyarray, (lat, z)))

    # Eqn p27.  UNESCO 1983.
    lat = np.abs(lat)
    X = np.sin(lat * deg2rad)
    sin2 = X * X
    grav = 9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * sin2) * sin2)
    return grav / ((1 + z / earth_radius) ** 2)  # From A.E.Gill p.597.


def pden(s, t, p, pr=0):
    """
    Calculates potential density of water mass relative to the specified
    reference pressure by pden = dens(S, ptmp, PR).

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].
    pr : number
         reference pressure [db], default = 0

    Returns
    -------
    pden : array_like
           potential density relative to the ref. pressure [kg m :sup:3]

    Examples
    --------
    >>> # Data from Unesco Tech. Paper in Marine Sci. No. 44, p22.
    >>> import seawater as sw
    >>> from seawater.library import T90conv
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> p = [0, 10000, 0, 10000, 0, 10000, 0, 10000]
    >>> sw.pden(s, t, p)
    array([  999.842594  ,   999.79523994,   995.65113374,   996.36115932,
            1028.10633141,  1028.15738545,  1021.72863949,  1022.59634627])

    :math:`\sigma_{4}` (at 4000 db)

    >>> sw.pden(s, t, p, 4000) - 1000
    array([ 19.2895493 ,  19.33422519,  12.43271053,  13.27563816,
            46.30976432,  46.48818851,  37.76150878,  38.74500757])

    References
    ----------
    .. [1] A.E. Gill 1982. p.54  Eqn. 3.7.15 "Atmosphere-Ocean Dynamics"
       Academic Press: New York. ISBN: 0-12-283522-0

    """

    s, t, p, pr = list(map(np.asanyarray, (s, t, p, pr)))

    pt = ptmp(s, t, p, pr)
    return dens(s, pt, pr)


def pres(depth, lat):
    """
    Calculates pressure in dbars from depth in meters.

    Parameters
    ----------
    depth : array_like
            depth [meters]
    lat : array_like
          latitude in decimal degrees north [-90..+90]

    Returns
    -------
    p : array_like
           pressure [db]

    Examples
    --------
    >>> import seawater as sw
    >>> depth, lat = 7321.45, 30
    >>> sw.pres(depth,lat)
    7500.0065130118019

    References
    ----------
    .. [1] Saunders, Peter M., 1981: Practical Conversion of Pressure to Depth.
       J. Phys. Oceanogr., 11, 573-574.
       doi: 10.1175/1520-0485(1981)011<0573:PCOPTD>2.0.CO;2

    """
    depth, lat = list(map(np.asanyarray, (depth, lat)))

    X = np.sin(np.abs(lat * deg2rad))
    C1 = 5.92e-3 + X ** 2 * 5.25e-3
    return ((1 - C1) - (((1 - C1) ** 2) - (8.84e-6 * depth)) ** 0.5) / 4.42e-6


def ptmp(s, t, p, pr=0):
    """
    Calculates potential temperature as per UNESCO 1983 report.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].
    pr : array_like
        reference pressure [db], default = 0

    Returns
    -------
    pt : array_like
         potential temperature relative to PR [℃ (ITS-90)]

    Examples
    --------
    >>> import seawater as sw
    >>> from seawater.library import T90conv, T68conv
    >>> t = T90conv([[0, 0, 0, 0, 0, 0],
    ...              [10, 10, 10, 10, 10, 10],
    ...              [20, 20, 20, 20, 20, 20],
    ...              [30, 30, 30, 30, 30, 30],
    ...              [40, 40, 40, 40, 40, 40]])
    >>> s = [[25, 25, 25, 35, 35, 35],
    ...      [25, 25, 25, 35, 35, 35],
    ...      [25, 25, 25, 35, 35, 35],
    ...      [25, 25, 25, 35, 35, 35],
    ...      [25, 25, 25, 35, 35, 35]]
    >>> p = [0, 5000, 10000, 0, 5000, 10000]
    >>> T68conv(sw.ptmp(s, t, p, pr=0))
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
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Bryden, H. 1973. New Polynomials for thermal expansion, adiabatic
       temperature gradient and potential temperature of sea water. Deep-Sea
       Res. Vol20,401-408. doi:10.1016/0011-7471(73)90063-6

    """

    s, t, p, pr = list(map(np.asanyarray, (s, t, p, pr)))

    # Theta1.
    del_P = pr - p
    del_th = del_P * adtg(s, t, p)
    th = T68conv(t) + 0.5 * del_th
    q = del_th

    # Theta2.
    del_th = del_P * adtg(s, T90conv(th), p + 0.5 * del_P)
    th = th + (1 - 1 / 2 ** 0.5) * (del_th - q)
    q = (2 - 2 ** 0.5) * del_th + (-2 + 3 / 2 ** 0.5) * q

    # Theta3.
    del_th = del_P * adtg(s, T90conv(th), p + 0.5 * del_P)
    th = th + (1 + 1 / 2 ** 0.5) * (del_th - q)
    q = (2 + 2 ** 0.5) * del_th + (-2 - 3 / 2 ** 0.5) * q

    # Theta4.
    del_th = del_P * adtg(s, T90conv(th), p + del_P)
    return T90conv(th + (del_th - 2 * q) / 6)


def salt(r, t, p):
    """
    Calculates Salinity from conductivity ratio. UNESCO 1983 polynomial.

    Parameters
    ----------
    r : array_like
        conductivity ratio :math:`R = \frac{C(S,T,P)}{C(35,15(IPTS-68),0)}`
    t : array_like
        temperature [℃ (ITS-90)]
    p : array_like
        pressure [db]

    Returns
    -------
    s : array_like
        salinity [psu (PSS-78)]

    Examples
    --------
    Data from UNESCO 1983 p9.
    >>> import seawater as sw
    >>> from seawater.library import T90conv
    >>> r = [1, 1.2, 0.65]
    >>> t = T90conv([15, 20, 5])
    >>> p = [0, 2000, 1500]
    >>> sw.salt(r, t, p)
    array([ 34.99999992,  37.24562765,  27.99534693])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp. Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    """
    r, t, p = list(map(np.asanyarray, (r, t, p)))

    rt = salrt(t)
    rp = salrp(r, t, p)
    rt = r / (rp * rt)
    return sals(rt, t)


def svel(s, t, p):
    """
    Sound Velocity in sea water using UNESCO 1983 polynomial.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    svel : array_like
           sound velocity  [m/s]

    Examples
    --------
    Data from Pond and Pickard Intro. Dynamical Oceanography 2nd ed. 1986

    >>> import seawater as sw
    >>> from seawater.library import T90conv
    >>> t = T90conv([[  0,  0,  0,  0,  0,  0],
    ...              [ 10, 10, 10, 10, 10, 10],
    ...              [ 20, 20, 20, 20, 20, 20],
    ...              [ 30, 30, 30, 30, 30, 30],
    ...              [ 40, 40, 40, 40, 40, 40]])
    >>> s = [[ 25, 25, 25, 35, 35, 35],
    ...      [ 25, 25, 25, 35, 35, 35],
    ...      [ 25, 25, 25, 35, 35, 35],
    ...      [ 25, 25, 25, 35, 35, 35],
    ...      [ 25, 25, 25, 35, 35, 35]]
    >>> p = [ 0, 5000, 10000, 0, 5000, 10000]
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
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    """
    s, t, p = list(map(np.asanyarray, (s, t, p)))

    # UNESCO 1983. Eqn..33  p.46.
    p = p / 10  # Convert db to bars as used in UNESCO routines.
    T68 = T68conv(t)

    # Eqn 34 p.46.
    c00, c01, c02, c03, c04, c05 = (1402.388, 5.03711, -5.80852e-2, 3.3420e-4,
                                    -1.47800e-6, 3.1464e-9)
    c10, c11, c12, c13, c14 = (0.153563, 6.8982e-4, -8.1788e-6, 1.3621e-7,
                               -6.1185e-10)
    c20, c21, c22, c23, c24 = (3.1260e-5, -1.7107e-6, 2.5974e-8, -2.5335e-10,
                               1.0405e-12)
    c30, c31, c32 = (-9.7729e-9, 3.8504e-10, -2.3643e-12)

    Cw = (((((c32 * T68 + c31) * T68 + c30) * p +
            ((((c24 * T68 + c23) * T68 + c22) * T68 + c21) * T68 + c20)) * p +
           ((((c14 * T68 + c13) * T68 + c12) * T68 + c11) * T68 + c10)) *
          p + ((((c05 * T68 + c04) * T68 + c03) * T68 + c02) * T68 + c01) *
          T68 + c00)

    # Eqn. 35. p.47
    a00, a01, a02, a03, a04 = (1.389, -1.262e-2, 7.164e-5, 2.006e-6, -3.21e-8)
    a10, a11, a12, a13, a14 = (9.4742e-5, -1.2580e-5, -6.4885e-8, 1.0507e-8,
                               -2.0122e-10)
    a20, a21, a22, a23 = (-3.9064e-7, 9.1041e-9, -1.6002e-10, 7.988e-12)
    a30, a31, a32 = (1.100e-10, 6.649e-12, -3.389e-13)

    A = (((((a32 * T68 + a31) * T68 + a30) * p +
           (((a23 * T68 + a22) * T68 + a21) * T68 + a20)) * p +
          ((((a14 * T68 + a13) * T68 + a12) * T68 + a11) * T68 + a10)) * p +
         (((a04 * T68 + a03) * T68 + a02) * T68 + a01) * T68 + a00)

    # Eqn 36 p.47.
    b00, b01, b10, b11 = -1.922e-2, -4.42e-5, 7.3637e-5, 1.7945e-7
    B = b00 + b01 * T68 + (b10 + b11 * T68) * p

    # Eqn 37 p.47.
    d00, d10 = 1.727e-3, -7.9836e-6
    D = d00 + d10 * p

    # Eqn 33 p.46.
    return Cw + A * s + B * s * s ** 0.5 + D * s ** 2


def temp(s, pt, p, pr=0):
    """
    Calculates temperature from potential temperature at the reference
    pressure PR and in situ pressure P.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    pt(p) : array_like
            potential temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].
    pr : array_like
         reference pressure [db]

    Returns
    -------
    temp : array_like
           temperature [℃ (ITS-90)]

    Examples
    --------
    >>> import seawater as sw
    >>> s, t, p = 35, 15, 100
    >>> sw.temp(s, sw.ptmp(s, t, p), p)
    15.0

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Bryden, H. 1973. New Polynomials for thermal expansion, adiabatic
       temperature gradient and potential temperature of sea water. Deep-Sea
       Res.  Vol20,401-408. doi:10.1016/0011-7471(73)90063-6

    """
    s, pt, p, pr = list(map(np.asanyarray, (s, pt, p, pr)))
    # Carry out inverse calculation by swapping p0 & pr.
    return ptmp(s, pt, pr, p)
