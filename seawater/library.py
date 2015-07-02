# -*- coding: utf-8 -*-


from __future__ import division, absolute_import

import numpy as np


__all__ = ['cndr',
           'salds',
           'salrp',
           'salrt',
           'seck',
           'sals',
           'smow',
           'T68conv',
           'T90conv']


# Constants.
a = [0.0080, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081]
b = [0.0005, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144]
c = [0.6766097, 2.00564e-2, 1.104259e-4, -6.9698e-7, 1.0031e-9]
d = [3.426e-2, 4.464e-4, 4.215e-1, -3.107e-3]
e = [2.070e-5, -6.370e-10, 3.989e-15]
k = 0.0162


def cndr(s, t, p):
    """
    Calculates conductivity ratio.

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
    cndr : array_like
           conductivity ratio. R = C(s,t,p) / C(35,15(IPTS-68),0) [no units]

    Examples
    --------
    >>> # Data from UNESCO 1983 p9.
    >>> import seawater as sw
    >>> t = T90conv([0, 10, 0, 10, 10, 30])
    >>> p = [0, 0, 1000, 1000, 0, 0]
    >>> s = [25, 25, 25, 25, 40, 40]
    >>> sw.cndr(s, t, p)
    array([ 0.49800825,  0.65499015,  0.50624434,  0.66297496,  1.00007311,
            1.52996697])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.UNESCO.org/images/0005/000598/059832eb.pdf

    """
    s, t, p = list(map(np.asanyarray, (s, t, p)))

    T68 = T68conv(t)

    shape = s.shape
    s, t = list(map(np.ravel, (s, t)))
    # Do a Newton-Raphson iteration for inverse interpolation of Rt from s.
    Rx = []
    for S, T in zip(s, t):
        Rx_loop = np.sqrt(S / 35.0)  # first guess at Rx = sqrt(Rt).
        SInc = sals(Rx_loop * Rx_loop, T)  # S Increment (guess) from Rx.
        iloop = 0
        while True:
            # FIXME: I believe that T / 1.00024 isn't correct here.  But I'm
            # reproducing seawater up to its bugs!
            Rx_loop = Rx_loop + (S - SInc) / salds(Rx_loop, T / 1.00024 - 15)
            SInc = sals(Rx_loop * Rx_loop, T)
            iloop += 1
            dels = abs(SInc - S)
            if (dels > 1.0e-10) and (iloop < 100):
                pass
            else:
                break

        Rx.append(Rx_loop)

    Rx = np.array(Rx).reshape(shape)

    # Once Rt found, corresponding to each (s,t) evaluate r.
    # Eqn(4) p.8 UNESCO 1983.
    A = (d[2] + d[3] * T68)
    B = 1 + d[0] * T68 + d[1] * T68 ** 2
    C = p * (e[0] + e[1] * p + e[2] * p ** 2)

    # Eqn(6) p.9 UNESCO 1983.
    Rt = Rx ** 2
    rt = salrt(t)
    # Rtrt  = rt * Rt # NOTE: unused in the code, but present in the original.
    D = B - A * rt * Rt
    E = rt * Rt * A * (B + C)
    r = np.sqrt(np.abs(D ** 2 + 4 * E)) - D
    r = 0.5 * r / A
    return r


def salds(rtx, delt):
    """
    Calculates Salinity differential (:math:`\\frac{dS}{d(\\sqrt{Rt})}`) at
    constant temperature.

    Parameters
    ----------
    rtx : array_like
          :math:`\\sqrt{rt}`
    delt : array_like
           t-15 [℃ (IPTS-68)]

    Returns
    -------
    ds : array_like
         :math:`\\frac{dS}{d rtx}`

    Examples
    --------
    >>> # Data from UNESCO 1983 p9.
    >>> import numpy as np
    >>> import seawater as sw
    >>> delt = T90conv([15, 20, 5]) - 15
    >>> rtx  = np.array([ 1, 1.0568875, 0.81705885]) ** 0.5
    >>> sw.salds(rtx, delt)
    array([ 78.31921607,  81.5689307 ,  68.19023687])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.UNESCO.org/images/0005/000598/059832eb.pdf

    """
    rtx, delt = list(map(np.asanyarray, (rtx, delt)))

    ds = (a[1] +
          (2 * a[2] + (3 * a[3] + (4 * a[4] + 5 * a[5] * rtx) * rtx) * rtx) *
          rtx + (delt / (1 + k * delt)) *
          (b[1] +
           (2 * b[2] + (3 * b[3] + (4 * b[4] + 5 * b[5] * rtx) * rtx) * rtx) *
           rtx))

    return ds


def salrp(r, t, p):
    """
    Equation for Rp used in calculating salinity. UNESCO 1983 polynomial.

    .. math::
        Rp(S,T,P) = \\frac{C(S,T,P)}{C(S,T,0)}


    Parameters
    ----------
    r : array_like
        conductivity ratio :math:`R = \\frac{C(S,T,P)}{C(35,15(IPTS-68),0)}`
    t : array_like
        temperature [℃ (ITS-90)]
    p : array_like
        pressure [db]

    Returns
    -------
    rp : array_like
        conductivity ratio :math:`Rp(S,T,P) = \\frac{C(S,T,P)}{C(S,T,0)}`

    Examples
    --------
    >>> import seawater as sw
    >>> r = [1, 1.2, 0.65]
    >>> t = T90conv([15, 20, 5])
    >>> p = [0, 2000, 1500]
    >>> sw.salrp(r, t, p)
    array([ 1.        ,  1.01694294,  1.02048638])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    """
    r, t, p = list(map(np.asanyarray, (r, t, p)))

    # Eqn(4) p.8 UNESCO.
    T68 = T68conv(t)

    rp = (1 + (p * (e[0] + e[1] * p + e[2] * p ** 2)) /
          (1 + d[0] * T68 + d[1] * T68 ** 2 + (d[2] + d[3] * T68) * r))

    return rp


def salrt(t):
    """
    Equation for rt used in calculating salinity. UNESCO 1983 polynomial.

    .. math::
        rt(t) = \\frac{C(35,t,0)}{C(35,15(\\textrm{IPTS-68}), 0)}


    Parameters
    ----------
      t : array_like
          temperature [℃ (ITS-90)]

    Returns
    -------
    rt : array_like
    conductivity ratio  [no units]

    Examples
    --------
    >>> # Data from UNESCO 1983 p9.
    >>> import seawater as sw
    >>> t = T90conv([15, 20, 5])
    >>> sw.salrt(t)
    array([ 1.        ,  1.11649272,  0.77956585])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    """
    t = np.asanyarray(t)
    T68 = T68conv(t)
    # Eqn (3) p.7 UNESCO.
    return c[0] + (c[1] + (c[2] + (c[3] + c[4] * T68) * T68) * T68) * T68


def seck(s, t, p=0):
    """
    Secant Bulk Modulus (K) of Sea Water using Equation of state 1980.
    UNESCO polynomial implementation.

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
    k : array_like
        secant bulk modulus [bars]

    Examples
    --------
    >>> # Data from Unesco Tech. Paper in Marine Sci. No. 44, p22.
    >>> import seawater as sw
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> p = [0, 10000, 0, 10000, 0, 10000, 0, 10000]
    >>> sw.seck(s, t, p)
    array([ 19652.21      ,  22977.2115    ,  22336.0044572 ,  25656.8196222 ,
            21582.27006823,  24991.99729129,  23924.21823158,  27318.32472464])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere equation
       of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
       doi:10.1016/0198-0149(81)90122-9

    """
    s, t, p = list(map(np.asanyarray, (s, t, p)))

    # Compute compression terms.
    p = p / 10.0  # Convert from db to atmospheric pressure units.
    T68 = T68conv(t)

    # Pure water terms of the secant bulk modulus at atmos pressure.
    # UNESCO Eqn 19 p 18.
    # h0 = -0.1194975
    h = [3.239908, 1.43713e-3, 1.16092e-4, -5.77905e-7]
    AW = h[0] + (h[1] + (h[2] + h[3] * T68) * T68) * T68

    # k0 = 3.47718e-5
    k = [8.50935e-5, -6.12293e-6, 5.2787e-8]
    BW = k[0] + (k[1] + k[2] * T68) * T68

    # e0 = -1930.06
    e = [19652.21, 148.4206, -2.327105, 1.360477e-2, -5.155288e-5]
    KW = e[0] + (e[1] + (e[2] + (e[3] + e[4] * T68) * T68) * T68) * T68

    # Sea water terms of secant bulk modulus at atmos. pressure.
    j0 = 1.91075e-4
    i = [2.2838e-3, -1.0981e-5, -1.6078e-6]
    A = AW + (i[0] + (i[1] + i[2] * T68) * T68 + j0 * s ** 0.5) * s

    m = [-9.9348e-7, 2.0816e-8, 9.1697e-10]
    B = BW + (m[0] + (m[1] + m[2] * T68) * T68) * s  # Eqn 18.

    f = [54.6746, -0.603459, 1.09987e-2, -6.1670e-5]
    g = [7.944e-2, 1.6483e-2, -5.3009e-4]
    K0 = (KW + (f[0] + (f[1] + (f[2] + f[3] * T68) * T68) * T68 +
                (g[0] + (g[1] + g[2] * T68) * T68) * s ** 0.5) * s)  # Eqn 16.
    return K0 + (A + B * p) * p  # Eqn 15.


def sals(rt, t):
    """
    Salinity of sea water as a function of Rt and T.  UNESCO 1983 polynomial.

    Parameters
    ----------
    rt : array_like
         :math:`rt(s,t) = \\frac{C(s,t,0)}{C(35, t(\\textrm{IPTS-68}), 0)}`
    t : array_like
        temperature [℃ (ITS-90)]

    Returns
    -------
    s : array_like
        salinity [psu (PSS-78)]

    Examples
    --------
    >>> # Data from UNESCO 1983 p9.
    >>> import seawater as sw
    >>> t = T90conv([15, 20, 5])
    >>> rt = [1, 1.0568875, 0.81705885]
    >>> sw.sals(rt, t)
    array([ 35.        ,  37.24562718,  27.99534701])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.UNESCO.org/images/0005/000598/059832eb.pdf

    """
    rt, t = list(map(np.asanyarray, (rt, t)))

    # Eqn (1) & (2) p6,7 UNESCO.
    del_T68 = T68conv(t) - 15

    Rtx = (rt) ** 0.5
    del_S = ((del_T68 / (1 + k * del_T68)) *
             (b[0] + (b[1] + (b[2] + (b[3] + (b[4] + b[5] * Rtx) *
                                      Rtx) * Rtx) * Rtx) * Rtx))
    s = a[0] + (a[1] + (a[2] + (a[3] + (a[4] + a[5] * Rtx) *
                                Rtx) * Rtx) * Rtx) * Rtx
    s += del_S

    return s


def smow(t):
    """
    Density of Standard Mean Ocean Water (Pure Water) using EOS 1980.

    Parameters
    ----------
    t : array_like
        temperature [℃ (ITS-90)]

    Returns
    -------
    dens(t) : array_like
              density  [kg m :sup:`3`]

    Examples
    --------
    >>> # Data from UNESCO Tech. Paper in Marine Sci. No. 44, p22.
    >>> import seawater as sw
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> sw.smow(t)
    array([ 999.842594  ,  999.842594  ,  995.65113374,  995.65113374,
            999.842594  ,  999.842594  ,  995.65113374,  995.65113374])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere equation
       of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
       doi:10.1016/0198-0149(81)90122-9

    """
    t = np.asanyarray(t)

    a = (999.842594, 6.793952e-2, -9.095290e-3, 1.001685e-4, -1.120083e-6,
         6.536332e-9)

    T68 = T68conv(t)
    return (a[0] + (a[1] + (a[2] + (a[3] + (a[4] + a[5] * T68) * T68) * T68) *
            T68) * T68)


def T68conv(T90):
    """
    Convert ITS-90 temperature to IPTS-68

    :math:`T68  = T90 * 1.00024`

    Parameters
    ----------
    t : array_like
           temperature [℃ (ITS-90)]

    Returns
    -------
    t : array_like
           temperature [℃ (IPTS-68)]

    Notes
    -----
    The International Practical Temperature Scale of 1968 (IPTS-68) need to be
    correct to the ITS-90. This linear transformation is accurate within
    0.5 ℃ for conversion between IPTS-68 and ITS-90 over the
    oceanographic temperature range.

    Examples
    --------
    >>> import seawater as sw
    >>> T68conv(19.995201151723585)
    20.0

    References
    ----------
    .. [1] Saunders, P. M., 1991: The International Temperature Scale of 1990,
       ITS-90. WOCE Newsletter, No. 10, WOCE International Project Office,
       Southampton, United Kingdom, 10.

    """
    T90 = np.asanyarray(T90)
    return T90 * 1.00024


def T90conv(t, t_type='T68'):
    """
    Convert IPTS-68 or IPTS-48 to temperature to ITS-90.

    T48 apply to all data collected prior to 31/12/1967.
    T68 apply to all data collected between 01/10/1968 and 31/12/1989.

    .. math::
        T90 = T68 / 1.00024

        T90 = T48 - (4.4e-6) * T48 * (100-T48) ) / 1.00024

    Parameters
    ----------
    t : array_like
           temperature [℃ (IPTS-68) or (IPTS-48)]
    t_type : string, optional
            'T68' (default) or 'T48'

    Returns
    -------
    T90 : array_like
           temperature [℃ (ITS-90)]

    Notes
    -----
    The International Practical Temperature Scale of 1968 (IPTS-68) need to be
    correct to the ITS-90. This linear transformation is accurate within
    0.5 ℃ for conversion between IPTS-68 and ITS-90 over the
    oceanographic temperature range.

    Examples
    --------
    >>> T90conv(20.004799999999999)
    20.0
    >>> T90conv(20., t_type='T48')
    19.988162840918179

    References
    ----------
    .. [1] Saunders, P. M., 1991: The International Temperature Scale of 1990,
       ITS-90. WOCE Newsletter, No. 10, WOCE International Project Office,
       Southampton, United Kingdom, 10.

    .. [2] International Temperature Scales of 1948, 1968 and 1990, an ICES
       note, available from http://www.ices.dk/ocean/procedures/its.htm

    """
    t = np.asanyarray(t)

    if t_type == 'T68':
        T90 = t / 1.00024
    elif t_type == 'T48':
        T90 = (t - 4.4e-6 * t * (100 - t)) / 1.00024
    else:
        raise NameError("Unrecognized temperature type.  Try 'T68'' or 'T48'")

    return T90


def atleast_2d(*arys):
    """
    Same as numpy atleast_2d, but with the single dimension last, instead of
    first.

    """
    res = []
    for ary in arys:
        ary = np.asanyarray(ary)
        if len(ary.shape) == 0:
            result = ary.reshape(1, 1)
        elif len(ary.shape) == 1:
            result = ary[:, np.newaxis]
        else:
            result = ary
        res.append(result)
    if len(res) == 1:
        return res[0]
    else:
        return res


if __name__ == '__main__':
    import doctest
    doctest.testmod()
