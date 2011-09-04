# -*- coding: utf-8 -*-

"""
Other conversions between temperatures, salinities, pressure and height

Functions:
----------

  t_from_CT(SA, CT, p)
      in-situ temperature from Conservative Temperature
  pt_from_t(SA, t, p, pr=0)
      potential temperature
  CT_from_pt(SA, pt)
      Conservative Temperature from potential temperature
  pot_enthalpy_from_pt(SA, pt)
      potential enthalpy from potential temperature
  pt0_from_t(SA, t, p)
      potential temperature with a reference pressure of zero dbar
  pt_from_CT(SA, CT)
      potential temperature from Conservative Temperature
  SP_from_SA(SA, p, lon, lat)
      Practical Salinity from Absolute Salinity
  Sstar_from_SA(SA, p, lon, lat)
      Preformed Salinity from Absolute Salinity
  SA_from_Sstar(Sstar, p, lon, lat)
      Absolute Salinity from Preformed Salinity
  SP_from_Sstar(Sstar, p, lon, lat)
      Practical Salinity from Preformed Salinity
  z_from_p(p, lat)
      height from pressure
  p_from_z(z, lat)
      pressure from height
  t90_from_t48(t48)
      ITS-90 temperature from IPTS-48 temperature
  t90_from_t68(t68)
      ITS-90 temperature from IPTS-68 temperature

This is part of the python Gibbs Sea Water library
http://code.google.com/p/python-seawater/

"""

from __future__ import division

import numpy as np

from seawater.constants import SSO, cp0, Kelvin, sfac
from library import entropy_part, gibbs
from utilities import match_args_return, strip_mask

__all__ = ['t_from_CT',
           'pt_from_t',
           'CT_from_pt',
           'pot_enthalpy_from_pt',
           #'pt0_from_t',
           #'pt_from_CT',
           #'SP_from_SA',
           #'Sstar_from_SA',
           #'SA_from_Sstar',
           #'SP_from_Sstar',
           #'z_from_p',
           #'p_from_z',
           #'t90_from_t48',
           #'t90_from_t68'
          ]

rad = np.pi / 180.0


@match_args_return
def t_from_CT(SA, CT, p):
    r"""
    Calculates *in-situ* temperature from Conservative Temperature of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.t_from_CT(SA, CT, p)
    array([ 28.78558023,  28.43287225,  22.81032309,  10.26001075,
             6.8862863 ,   4.40362445])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See sections 3.1 and 3.3.

    Modifications:
    2010-08-26. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    pt0 = pt_from_CT(SA, CT)

    t = pt_from_t(SA, pt0, 0, p)

    return t


@match_args_return
def pt_from_t(SA, t, p, p_ref=0):
    r"""
    Calculates potential temperature with the general reference pressure, pr,
    from in situ temperature.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]
    p_ref : int, float, optional
         reference pressure, default = 0

    Returns
    -------
    pt : array_like
         potential temperature [:math:`^\circ` C (ITS-90)]

    See Also
    --------
    _entropy_part

    Notes
    -----
    This function calls `entropy_part` which evaluates entropy except for the
    parts which are a function of Absolute Salinity alone. A faster routine
    exists pt0_from_t(SA,t,p) if p_ref is indeed zero dbar.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.pt_from_t(SA, t, p)
    array([ 28.78319682,  28.42098334,  22.7849304 ,  10.23052366,
             6.82923022,   4.32451057])
    >>> gsw.pt_from_t(SA, t, p, pr = 1000)
    array([ 29.02665528,  28.662375  ,  22.99149634,  10.35341725,
             6.92732954,   4.4036    ])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 3.1.

    .. [2] McDougall T.J., P.M. Barker, R. Feistel and D.R. Jackett, 2011:  A
    computationally efficient 48-term expression for the density of seawater
    in terms of Conservative Temperature, and related properties of seawater.

    Modifications:
    2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and
    Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    p_ref = np.asanyarray(p_ref)

    n0, n2 = 0, 2

    SA[SA < 0] = 0

    s1 = SA * 35. / SSO

    pt = (t + (p - p_ref) * (8.65483913395442e-6 -
                       s1 * 1.41636299744881e-6 -
              (p + p_ref) * 7.38286467135737e-9 +
                        t * (-8.38241357039698e-6 +
                       s1 * 2.83933368585534e-8 +
                        t * 1.77803965218656e-8 +
              (p + p_ref) * 1.71155619208233e-10)))

    dentropy_dt = cp0 / ((Kelvin + pt) * (1 - 0.05 *
                                       (1 - SA / SSO)))

    true_entropy_part = entropy_part(SA, t, p)

    for Number_of_iterations in range(0, 2, 1):
        pt_old = pt
        dentropy = entropy_part(SA, pt_old, p_ref) - true_entropy_part
        pt = pt_old - dentropy / dentropy_dt  # half way through the method
        ptm = 0.5 * (pt + pt_old)
        dentropy_dt = -gibbs(n0, n2, n0, SA, ptm, p_ref)
        pt = pt_old - dentropy / dentropy_dt

    # maximum error of 6.3x10^-9 degrees C for one iteration.
    # maximum error is 1.8x10^-14 degrees C for two iterations (two iterations
    # is the default, "for Number_of_iterations = 1:2).
    # These errors are over the full "oceanographic funnel" of
    # McDougall et al. (2010), which reaches down to p = 8000 dbar.

    return pt


@match_args_return
def CT_from_pt(SA, pt):
    r"""
    Calculates Conservative Temperature of seawater from potential temperature
    (whose reference sea pressure is zero dbar).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> pt = [28.7832, 28.4209, 22.7850, 10.2305, 6.8292, 4.3245]
    >>> gsw.CT_from_pt(SA, pt)
    array([ 28.80992302,  28.43914426,  22.78624661,  10.22616561,
             6.82718342,   4.32356518])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 3.3.

    Modifications:
    2010-08-05. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, pt, mask = strip_mask(SA, pt)

    pot_enthalpy = pot_enthalpy_from_pt(SA, pt)

    CT = pot_enthalpy / cp0

    return np.ma.array(CT, mask=mask, copy=False)


@match_args_return
def pot_enthalpy_from_pt(SA, pt):
    r"""
    Calculates the potential enthalpy of seawater from potential temperature
    (whose reference sea pressure is zero dbar).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    pot_enthalpy : array_like
                   potential enthalpy [J kg :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> pt = [28.7832, 28.4209, 22.7850, 10.2305, 6.8292, 4.3245]
    >>> gsw.pot_enthalpy_from_pt(SA, pt)
    array([ 115005.40853458,  113525.30870246,   90959.68769935,
             40821.50280454,   27253.21472227,   17259.10131183])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 3.2.

    Modifications:
    2010-08-26. David Jackett, Trevor McDougall and Paul Barker
    2011-02-16. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, pt, mask = strip_mask(SA, pt)

    x2 = sfac * SA
    x = np.sqrt(x2)
    y = pt * 0.025  # Normalize for F03 and F08

    pot_enthalpy = (61.01362420681071 + y * (168776.46138048015 +
    y * (-2735.2785605119625 + y * (2574.2164453821433 +
    y * (-1536.6644434977543 + y * (545.7340497931629 +
    (-50.91091728474331 - 18.30489878927802 * y) * y))))) +
    x2 * (268.5520265845071 + y * (-12019.028203559312 +
    y * (3734.858026725145 + y * (-2046.7671145057618 +
    y * (465.28655623826234 + (-0.6370820302376359 -
    10.650848542359153 * y) * y)))) +
    x * (937.2099110620707 + y * (588.1802812170108 +
    y * (248.39476522971285 + (-3.871557904936333 -
    2.6268019854268356 * y) * y)) +
    x * (-1687.914374187449 + x * (246.9598888781377 +
    x * (123.59576582457964 - 48.5891069025409 * x)) +
    y * (936.3206544460336 +
    y * (-942.7827304544439 + y * (369.4389437509002 +
    (-33.83664947895248 - 9.987880382780322 * y) * y)))))))

    """
    The above polynomial for pot_enthalpy is the full expression for potential
    enthalpy in terms of SA and pt, obtained from the Gibbs function as below.

    It has simply collected like powers of x and y so that it is
    computationally faster than calling the Gibbs function twice as is done in
    the commented code below. When this code below is run, the results are
    identical to calculating pot_enthalpy as above, to machine precision.

    n0, n1 = 0, 1
    g000 = gibbs(n0, n0, n0, SA, pt, 0)
    g010 = gibbs(n0, n1, n0, SA, pt, 0)
    pot_enthalpy = g000 - (Kelvin + pt) * g010

    This is the end of the alternative code
    %timeit gsw.CT_from_pt(SA, pt)
    1000 loops, best of 3: 1.34 ms per loop <- calling gibbs
    1000 loops, best of 3: 254 us per loop <- standard
    """

    return np.ma.array(pot_enthalpy, mask=mask, copy=False)
