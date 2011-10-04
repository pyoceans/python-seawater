# -*- coding: utf-8 -*-

"""
Practical Salinity (SP), PSS-78

Functions:
----------

  gsw_SP_from_C
    Practical Salinity from conductivity, C (inc. for SP < 2)
  gsw_C_from_SP
    conductivity, C, from Practical Salinity (inc. for SP < 2)
  gsw_SP_from_R
    Practical Salinity from conductivity ratio, R (inc. for SP < 2)
  gsw_R_from_SP
    conductivity ratio, R, from Practical Salinity (inc. for SP < 2)
  gsw_SP_salinometer
    Practical Salinity from a laboratory salinometer (inc. for SP < 2)


This is part of the python Gibbs Sea Water library
http://code.google.com/p/python-seawater/

"""

from __future__ import division

import numpy as np

from library import gsw_Hill_ratio_at_SP2
from utilities import match_args_return

__all__ = ['SP_from_C',
           #'gsw_C_from_SP',
           #'gsw_SP_from_R',
           #'gsw_R_from_SP',
           'SP_salinometer']


@match_args_return
def SP_from_C(C, t, p):
    """
    Calculates Practical Salinity, SP, from conductivity, C, primarily using
    the PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical
    Salinity is only valid in the range 2 < SP < 42.  If the PSS-78 algorithm
    produces a Practical Salinity that is less than 2 then the Practical
    Salinity is recalculated with a modified form of the Hill et al. (1986)
    formula. The modification of the Hill et al. (1986) expression is to ensure
    that it is exactly consistent with PSS-78 at SP = 2.  Note that the input
    values of conductivity need to be in units of mS/cm (not S/m).

    Parameters
    ----------
    C : array
        conductivity [ mS/cm ]
    t : array
        in-situ temperature [:math:`^\circ` C (ITS-90)]
    p : array
        sea pressure [dbar]
        (i.e. absolute pressure - 10.1325 dbar)

    Returns
    -------
    SP : array
         Practical Salinity [psu (PSS-78), unitless]

    Examples
    --------
    TODO

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    References
    ----------
    .. [1] Culkin and Smith, 1980:  Determination of the Concentration of Potassium
    Chloride Solution Having the Same Electrical Conductivity, at 15C and
    Infinite Frequency, as Standard Seawater of Salinity 35.0000 (Chlorinity
    19.37394), IEEE J. Oceanic Eng, 5, 22-23.

    .. [2] Hill, K.D., T.M. Dauphinee & D.J. Woods, 1986: The extension of the
    Practical Salinity Scale 1978 to low salinities. IEEE J. Oceanic Eng., 11,
    109 - 112.

    .. [3] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
    seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.  Appendix E.

    .. [4] Unesco, 1983: Algorithms for computation of fundamental properties of
    seawater. Unesco Technical Papers in Marine Science, 44, 53 pp.

    Modifications:
    2011-04-01. Paul Barker, Trevor McDougall and Rich Pawlowicz.
    """

    C, t, p = np.broadcast_arrays(C, t, p)

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

    c0 =  0.6766097
    c1 =  2.00564e-2
    c2 =  1.104259e-4
    c3 = -6.9698e-7
    c4 =  1.0031e-9

    d1 =  3.426e-2
    d2 =  4.464e-4
    d3 =  4.215e-1
    d4 = -3.107e-3

    e1 =  2.070e-5
    e2 = -6.370e-10
    e3 =  3.989e-15

    k  =  0.0162

    t68 = t * 1.00024
    ft68 = (t68 - 15) / (1 + k*(t68 - 15))

    # The dimensionless conductivity ratio, R, is the conductivity input, C,
    # divided by the present estimate of C(SP=35, t_68=15, p=0) which is
    # 42.9140 mS/cm (=4.29140 S/m), (Culkin and Smith, 1980).

    R = 0.023302418791070513 * C  # 0.023302418791070513 = 1./42.9140

    # rt_lc corresponds to rt as defined in the UNESCO 44 (1983) routines.
    rt_lc = c0 + (c1 + (c2 + (c3 + c4*t68)*t68)*t68)*t68
    Rp = ( 1 + (p*(e1 + e2*p + e3*p*p)) /
          (1 + d1*t68 + d2*t68*t68 + (d3 + d4*t68)*R) )
    Rt = R / (Rp*rt_lc)

    Rt[Rt < 0] = np.nan
    Rtx = np.sqrt(Rt)

    SP = ( a0 + (a1 + (a2 + (a3 + (a4 + a5*Rtx)*Rtx)*Rtx)*Rtx)*Rtx +
        ft68*(b0 + (b1 + (b2 + (b3 + (b4 + b5*Rtx)*Rtx)*Rtx)*Rtx)*Rtx) )

    # The following section of the code is designed for SP < 2 based on the
    # Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
    # exactly equal to the PSS-78 algorithm at SP = 2.

    I2, = np.nonzero(np.ravel(SP) < 2)  # find
    if len(I2) > 0:
        Hill_ratio = Hill_ratio_at_SP2(t[I2])
        x = 400 * Rt[I2]
        sqrty = 10 * Rtx[I2]
        part1 = 1 + x * (1.5 + x)
        part2 = 1 + sqrty * (1 + sqrty * (1 + sqrty))
        SP_Hill_raw = SP[I2] - a0/part1 - b0*ft68[I2]/part2
        SP[I2] = Hill_ratio * SP_Hill_raw

    # Ensure that SP is non-negative.
    SP[SP < 0] = 0.0

    return SP


@match_args_return
def SP_salinometer(Rt, t):
    r"""
    Calculates Practical Salinity SP from a salinometer, primarily using the
    PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical Salinity is
    only valid in the range 2 < SP < 42.  If the PSS-78 algorithm produces a
    Practical Salinity that is less than 2 then the Practical Salinity is
    recalculated with a modified form of the Hill et al. (1986) formula. The
    modification of the Hill et al. (1986) expression is to ensure that it is
    exactly consistent with PSS-78 at SP = 2.

    A laboratory salinometer has the ratio of conductivities, Rt, as an output,
    and the present function uses this conductivity ratio and the temperature t
    of the salinometer bath as the two input variables.

    Parameters
    ----------
    Rt : array
         C(SP,t_68,0)/C(SP=35,t_68,0) [unitless]
         conductivity ratio :math:`R = \frac{C(S, t_68, 0)}{C(35, 15(IPTS-68),0)} [unitless]
    t : array
        Temperature of the bath of the salinometer [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    SP : array
         Practical Salinity [psu (PSS-78), unitless]

    See Also
    --------
    TODO: sw.sals

    Notes
    -----
    TODO

    Examples
    --------
    FIXME

    Data from UNESCO 1983 p9

    >>> import seawater.csiro as sw
    >>> t = T90conv([15, 20, 5])
    >>> rt   = [  1, 1.0568875, 0.81705885]
    >>> sw.sals(rt, t)
    array([ 35.        ,  37.24562718,  27.99534701])


    References
    -----------
    ..[1] Fofonoff, P. and R.C. Millard Jr. 1983: Algorithms for computation of
    fundamental properties of seawater. Unesco Tech. Pap. in Mar. Sci., 44,
    53 pp.

    ..[2] Hill, K.D., T.M. Dauphinee & D.J. Woods, 1986: The extension of the
    Practical Salinity Scale 1978 to low salinities. IEEE J. Oceanic Eng., 11,
    109 - 112.

    .. [3] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
    seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See appendix E of this TEOS-10 Manual, and in
    particular, Eqns. (E.2.1) and (E.2.6).

    Modifications:
    2011-04-30. Paul Barker, Trevor McDougall and Rich Pawlowicz. Version 3.0
    """

    C, t = np.broadcast_arrays(C, t)

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

    t68 = t * 1.00024
    ft68 = (t68 - 15) / (1 + k * (t68 - 15))

    #NOTE: Should we use numpy.NA or masked_array?
    Rt[Rt < 0] = np.NaN
    Rtx = np.sqrt(Rt)

    SP = (a0 + (a1 + (a2 + (a3 + (a4 + a5 * Rtx) * Rtx) * Rtx) * Rtx) * Rtx +
          ft68 * (b0 + (b1 + (b2 + (b3 + (b4 + b5 * Rtx) * Rtx) * Rtx) * Rtx)
                                                                        * Rtx))

    # The following section of the code is designed for SP < 2 based on the
    # Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
    # exactly equal to the PSS-78 algorithm at SP = 2.

    I2 = SP < 2
    if I2.any():
        Hill_ratio = Hill_ratio_at_SP2(t[I2])
        x = 400 * Rt[I2]
        sqrty = 10 * Rtx[I2]
        part1 = 1 + x * (1.5 + x)
        part2 = 1 + sqrty * (1 + sqrty * (1 + sqrty))
        SP_Hill_raw = SP[I2] - a0 / part1 - b0 * ft68[I2] / part2
        SP[I2] = Hill_ratio * SP_Hill_raw

    # Ensure that SP is non-negative.
    SP[SP < 0] = 0

    return SP


if __name__ == '__main__':
    import doctest
    doctest.testmod()
