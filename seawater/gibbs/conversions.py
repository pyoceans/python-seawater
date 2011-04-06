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
import seawater.constants as cte
from library import match_args_return
import library as lib

# -----------

__all__ = ['t_from_CT',
           'pt_from_t',
           'CT_from_pt',
           'pot_enthalpy_from_pt',
           'pt0_from_t',
           'pt_from_CT',
           'SP_from_SA',
           'Sstar_from_SA',
           'SA_from_Sstar',
           'SP_from_Sstar',
           'z_from_p',
           'p_from_z',
           't90_from_t48',
           't90_from_t68']

# -----------

rad = np.pi / 180.0

@match_args_return
def t_from_CT(SA, CT, p):
    r"""
    Calculates in situ temperature from Conservative Temperature of seawater.

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

# -------------------
@match_args_return
def pt_from_t(SA, t, p, pr=0):
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
    pr : int, float, optional
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
    This function calls "entropy_part" which evaluates entropy except for the
    parts which are a function of Absolute Salinity alone. A faster routine
    exists pt0_from_t(SA,t,p) if pr is indeed zero dbar.

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

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson,
    R. Feistel and R. W. Hallberg, 2010: A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative
    Temperature, and related properties of seawater.

    Modifications:
    2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    pr = np.asanyarray(pr)

    n0, n2 = 0, 2

    SA[SA < 0] = 0

    s1 = SA * 35. / cte.SSO

    pt = t + ( p - pr ) * ( 8.65483913395442e-6  -
    s1 * 1.41636299744881e-6 -
    ( p + pr ) * 7.38286467135737e-9 +
    t * ( -8.38241357039698e-6 +
    s1 * 2.83933368585534e-8 +
    t * 1.77803965218656e-8 +
    ( p + pr ) * 1.71155619208233e-10 ) )

    dentropy_dt = cte.cp0 / ( (cte.Kelvin + pt) * ( 1 - 0.05 *
                                       ( 1 - SA / cte.SSO ) ) )

    true_entropy_part = lib._entropy_part(SA, t, p)

    for Number_of_iterations in range(0,2,1):
        pt_old = pt
        dentropy = lib._entropy_part(SA, pt_old, pr) - true_entropy_part
        pt = pt_old - dentropy / dentropy_dt # half way through the mod. method
        ptm = 0.5 * (pt + pt_old)
        dentropy_dt = -lib._gibbs(n0, n2, n0, SA, ptm, pr)
        pt = pt_old - dentropy / dentropy_dt

    # maximum error of 6.3x10^-9 degrees C for one iteration.
    # maximum error is 1.8x10^-14 degrees C for two iterations (two iterations
    # is the default, "for Number_of_iterations = 1:2).
    # These errors are over the full "oceanographic funnel" of
    # McDougall et al. (2010), which reaches down to p = 8000 dbar.

    return pt


# ------------------------------

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

    SA, pt, mask = lib.strip_mask(SA, pt)

    pot_enthalpy =  pot_enthalpy_from_pt(SA, pt)

    CT = pot_enthalpy / cte.cp0

    return np.ma.array(CT, mask=mask, copy=False)

# ---------------------------


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

    SA, pt, mask = lib.strip_mask(SA, pt)

    x2 = cte.sfac * SA
    x = np.ma.sqrt(x2)
    y = pt * 0.025 # normalize for F03 and F08

    pot_enthalpy =  ( 61.01362420681071 + y * ( 168776.46138048015 +
    y * ( -2735.2785605119625 + y * ( 2574.2164453821433 +
    y * ( -1536.6644434977543 + y * ( 545.7340497931629 +
    ( -50.91091728474331 - 18.30489878927802 * y ) * y ) ) ) ) ) +
    x2 * ( 268.5520265845071 + y * ( -12019.028203559312 +
    y * ( 3734.858026725145 + y * ( -2046.7671145057618 +
    y * ( 465.28655623826234 + ( -0.6370820302376359 -
    10.650848542359153 * y ) * y ) ) ) ) +
    x * ( 937.2099110620707 + y * ( 588.1802812170108 +
    y * ( 248.39476522971285 + ( -3.871557904936333 -
    2.6268019854268356 * y ) * y ) ) +
    x * ( -1687.914374187449 + x * ( 246.9598888781377 +
    x * ( 123.59576582457964 - 48.5891069025409 * x ) ) +
    y * ( 936.3206544460336 +
    y * ( -942.7827304544439 + y * ( 369.4389437509002 +
    ( -33.83664947895248 - 9.987880382780322 * y ) * y ) ) ) ) ) ) )

    """
    The above polynomial for pot_enthalpy is the full expression for potential
    enthalpy in terms of SA and pt, obtained from the Gibbs function as below.

    It has simply collected like powers of x and y so that it is
    computationally faster than calling the Gibbs function twice as is done in
    the commented code below. When this code below is run, the results are
    identical to calculating pot_enthalpy as above, to machine precision.

    n0, n1 = 0, 1
    g000 = lib._gibbs(n0, n0, n0, SA, pt, 0)
    g010 = lib._gibbs(n0, n1, n0, SA, pt, 0)
    pot_enthalpy = g000 - (cte.Kelvin + pt) * g010

    #----------------This is the end of the alternative code------------------
    #%timeit gsw.CT_from_pt(SA, pt)
    #1000 loops, best of 3: 1.34 ms per loop <- calling gibbs
    #1000 loops, best of 3: 254 us per loop <- standard
    """

    return np.ma.array(pot_enthalpy, mask=mask, copy=False)

# -------------------------

@match_args_return
def pt0_from_t(SA, t, p):
    r"""
    Calculates potential temperature with reference pressure, pr = 0 dbar.
    The present routine is computationally faster than the more general
    function "pt_from_t(SA, t, p, pr)" which can be used for any reference
    pressure value.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    pt0 : array_like
          potential temperature relative to 0 dbar [:math:`^\circ` C (ITS-90)]

    See Also
    --------
    _entropy_part, _gibbs_pt0_pt0, _entropy_part_zerop

    Notes
    -----
    pt_from_t  has the same result (only slower)

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.pt0_from_t(SA, t, p)
    array([ 28.78319682,  28.42098334,  22.7849304 ,  10.23052366,
             6.82923022,   4.32451057])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 3.1.

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson,
    R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    Modifications:
    2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    s1 = SA * (35. / cte.SSO)

    pt0 = t + p * ( 8.65483913395442e-6 -
             s1 *   1.41636299744881e-6 -
              p *   7.38286467135737e-9 +
              t * (-8.38241357039698e-6 +
             s1 *   2.83933368585534e-8 +
              t *   1.77803965218656e-8 +
              p *   1.71155619208233e-10 ) )

    dentropy_dt = cte.cp0 / ( (cte.Kelvin + pt0) * ( 1 - 0.05 *
                                        ( 1 - SA / cte.SSO ) ) )

    true_entropy_part = lib._entropy_part(SA, t, p)

    for Number_of_iterations in range(0,2,1):
        pt0_old = pt0
        dentropy = lib._entropy_part_zerop(SA, pt0_old) - true_entropy_part
        pt0 = pt0_old - dentropy / dentropy_dt # half way through mod. method
        pt0m = 0.5 * (pt0 + pt0_old)
        dentropy_dt = -lib._gibbs_pt0_pt0(SA, pt0m)
        pt0 = pt0_old - dentropy / dentropy_dt

    """
    maximum error of 6.3x10^-9 degrees C for one iteration.
    maximum error is 1.8x10^-14 degrees C for two iterations (two iterations is
    the default, "for Number_of_iterations = 1:2")

    These errors are over the full "oceanographic funnel" of McDougall et al.
    (2010), which reaches down to p = 8000 dbar.
    """

    return pt0




# -----------------------


@match_args_return
def pt_from_CT(SA, CT):
    r"""
    Calculates potential temperature (with a reference sea pressure of zero
    dbar) from Conservative Temperature.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]

    See Also
    --------
    specvol_anom

    Notes
    -----
    This function uses 1.5 iterations through a modified Newton-Raphson (N-R)
    iterative solution procedure, starting from a rational-function-based
    initial condition for both pt and dCT_dpt.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.pt_from_CT(SA, CT)
    array([ 28.78317705,  28.4209556 ,  22.78495347,  10.23053439,
             6.82921659,   4.32453484])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See sections 3.1 and 3.3.

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson,
    R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    Modifications:
    2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, CT, mask = lib.strip_mask(SA, CT)

    s1 = SA * 35. / cte.SSO

    a0 = -1.446013646344788e-2
    a1 = -3.305308995852924e-3
    a2 =  1.062415929128982e-4
    a3 =  9.477566673794488e-1
    a4 =  2.166591947736613e-3
    a5 =  3.828842955039902e-3

    b0 =  1.000000000000000e+0
    b1 =  6.506097115635800e-4
    b2 =  3.830289486850898e-3
    b3 =  1.247811760368034e-6

    a5CT = a5 * CT
    b3CT = b3 * CT
    CT_factor = ( a3 + a4 * s1 + a5CT )
    pt_num = a0 + s1 * ( a1 + a2 * s1 ) + CT * CT_factor
    pt_den = b0 + b1 * s1 + CT * ( b2 + b3CT )
    pt = pt_num / pt_den

    dCT_dpt = pt_den / (CT_factor + a5CT - (b2 + b3CT + b3CT) * pt)

    # 1.5 iterations through the modified Newton-Rapshon iterative method
    CT_diff = CT_from_pt(SA, pt) - CT
    pt_old = pt
    pt = pt_old - CT_diff /dCT_dpt # 1/2-way through the 1st modified N-R loop
    ptm = 0.5 * (pt + pt_old)

    # This routine calls gibbs_pt0_pt0(SA,pt0) to get the second derivative of
    # the Gibbs function with respect to temperature at zero sea pressure.

    dCT_dpt = -(ptm + cte.Kelvin) * lib._gibbs_pt0_pt0(SA, ptm) / cte.cp0
    pt = pt_old - CT_diff / dCT_dpt # end of 1st full modified N-R iteration
    CT_diff = CT_from_pt(SA, pt) - CT
    pt_old = pt
    pt = pt_old - CT_diff / dCT_dpt # 1.5 iterations of the modified N-R method

    return np.ma.array(pt, mask=mask, copy=False)

# -------------------------------

@match_args_return
def SP_from_SA(SA, p, lon, lat):
    r"""
    Calculates Practical Salinity from Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SP : masked array
         salinity [psu (PSS-78)], unitless

    See Also
    --------
    _delta_SA, _SP_from_SA_Baltic

    Notes
    -----
    The mask is only set when the observation is well and truly on dry
    land; often the warning flag is not set until one is several hundred
    kilometers inland from the coast.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat =  188, 4
    >>> gsw.SP_from_SA(SA, p, lon, lat)
    array([ 34.54872019,  34.72747639,  34.86055202,  34.68097006,
            34.56797054,  34.56003796])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    lon, lat, p, SA = np.broadcast_arrays(lon, lat, p, SA)

    dSA = lib._delta_SA( p, lon, lat )

    SP = (35./cte.SSO) * ( SA - dSA )

    SP_baltic = lib._SP_from_SA_Baltic( SA, lon, lat )

    if SP_baltic is not None:
        SP[~SP_baltic.mask] = SP_baltic[~SP_baltic.mask]

    return SP

# ------------------------------


@match_args_return
def Sstar_from_SA(SA, p, lon, lat):
    r"""
    Converts Preformed Salinity from Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    Sstar : masked array
            Preformed Salinity [g kg :sup:`-1`]

    See Also
    --------
    _delta_SA

    Notes
    -----
    In the Baltic Sea, SA = Sstar, and note that _delta_SA returns zero for dSA
    in the Baltic.

    The mask is only set when the observation is well and truly on dry
    land; often the warning flag is not set until one is several hundred
    kilometers inland from the coast.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat =  188, 4
    >>> gsw.Sstar_from_SA(SA, p, lon, lat)
    array([ 34.71157349,  34.89113729,  35.02470152,  34.84356269,
            34.729004  ,  34.71971452])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    lon, lat, p, SA = np.broadcast_arrays(lon, lat, p, SA)

    dSA = lib._delta_SA( p, lon, lat )

    Sstar =  SA - ( 1 + cte.r1 ) * dSA

    return Sstar

# -----------------------------

@match_args_return
def SA_from_Sstar(Sstar, p, lon, lat):
    r"""
    Calculates Absolute Salinity from Preformed Salinity.

    Parameters
    ----------
    Sstar : array_like
            Preformed Salinity [g kg :sup:`-1`]
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SA : masked array
         Absolute salinity [g kg :sup:`-1`]

    See Also
    --------
    _delta_SA

    Notes
    -----
    In the Baltic Sea, SA = Sstar, and note that _delta_SA returns zero for dSA
    in the Baltic.

    The mask is only set when the observation is well and truly on dry
    land; often the warning flag is not set until one is several hundred
    kilometers inland from the coast.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> Sstar = [34.7115, 34.8912, 35.0247, 34.8436, 34.7291, 34.7197]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat = 188, 4
    >>> gsw.SA_from_Sstar(Sstar, p, lon, lat)
    array([ 34.71172651,  34.89156271,  35.02559848,  34.84723731,
            34.736696  ,  34.73238548])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    lon, lat, p, Sstar = np.broadcast_arrays(lon, lat, p, Sstar)

    dSA = lib._delta_SA( p, lon, lat )

    SA = Sstar + ( 1 + cte.r1 ) * dSA

    return SA

# -----------------------

@match_args_return
def SP_from_Sstar(Sstar, p, lon, lat):
    r"""
    Calculates Practical Salinity from Preformed Salinity.

    Parameters
    ----------
    Sstar : array_like
            Preformed Salinity [g kg :sup:`-1`]
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SP : masked array
         salinity [psu (PSS-78)], unitless

    See Also
    --------
    _delta_SA, _SP_from_SA_Baltic

    Notes
    -----
    The mask is only set when the observation is well and truly on dry
    land; often the warning flag is not set until one is several hundred
    kilometers inland from the coast.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> Sstar = [34.7115, 34.8912, 35.0247, 34.8436, 34.7291, 34.7197]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat =  188, 4
    >>> gsw.SP_from_Sstar(Sstar, p, lon, lat)
    array([ 34.54864705,  34.72753881,  34.8605505 ,  34.68100719,
            34.56806609,  34.56002351])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    lon, lat, p, Sstar = np.broadcast_arrays(lon, lat, p, Sstar)

    dSA = lib._delta_SA( p, lon, lat )
    SP = (35/cte.SSO) * ( Sstar + cte.r1 * dSA )

    # In the Baltic Sea, SA = Sstar.
    SP_baltic = lib._SP_from_SA_Baltic( Sstar, lon, lat )

    if SP_baltic is not None:
        SP[~SP_baltic.mask] = SP_baltic[~SP_baltic.mask]

    return SP

# -----------------------

@match_args_return
def  z_from_p(p, lat):
    r"""
    Calculates height from sea pressure using the computationally-efficient
    25-term expression for density in terms of SA, CT and p.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    p : array_like
        pressure [dbar]

    Returns
    -------
    z : array_like
        height [m]

    See Also
    --------
    _enthalpy_SSO_0_CT25


    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lat = 4
    >>> gsw.z_from_p(p, lat)
    array([  -9.94460074,  -49.71817465, -124.2728275 , -248.47044828,
           -595.82618014, -992.0931748 ])

    Notes
    -----
    At sea level z = 0, and since z (HEIGHT) is defined to be positive upwards,
    it follows that while z is positive in the atmosphere, it is NEGATIVE in
    the ocean.

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson,
    R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    .. [3] Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.

    Modifications:
    2010-08-26. Trevor McDougall, Claire Roberts-Thomson & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    X     = np.sin(lat*rad)
    sin2  = X**2
    B     = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2 )
    A     = -0.5 * cte.gamma * B
    C     = lib._enthalpy_SSO_0_CT25(p)
    z     = -2 * C / ( B + np.sqrt( B**2 - 4 * A * C ) )

    return z

# -------------------------

@match_args_return
def  p_from_z(z, lat):
    r"""
    Calculates sea pressure from height using computationally-efficient 25-term
    expression for density, in terms of SA, CT and p.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    z : array_like
        height [m]

    Returns
    -------
    p : array_like
        pressure [dbar]

    See Also
    --------
    _specvol_SSO_0_CT25, _enthalpy_SSO_0_CT25

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> z = [10, 50, 125, 250, 600, 1000]
    >>> lat = 4.
    >>> gsw.p_from_z(z, lat)
    array([  -10.05521794,   -50.2711751 ,  -125.6548857 ,  -251.23284504,
            -602.44050752, -1003.07609807])
    >>> z = [-9.94460074, -49.71817465, -124.2728275, -248.47044828, -595.82618014, -992.0931748]
    >>> gsw.p_from_z(z, lat)
    array([   10.,    50.,   125.,   250.,   600.,  1000.])

    Notes
    -----
    Height (z) is NEGATIVE in the ocean. Depth is -z. Depth is not used in the
    gibbs library.

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson,
    R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    .. [3] Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.

    .. [4] Saunders, P. M., 1981: Practical conversion of pressure to depth.
    Journal of Physical Oceanography, 11, 573-574.

    Modifications:
    2010-08-26. Trevor McDougall, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    X     = np.sin(lat*rad)
    sin2  = X**2
    gs    = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2 )
    # get the first estimate of p from Saunders (1981)
    c1 =  5.25e-3 * sin2 + 5.92e-3
    p  = -2 * z / ( (1-c1) + np.ma.sqrt( (1-c1) * (1-c1) + 8.84e-6 * z ) )
    # end of the first estimate from Saunders (1981)
    df_dp = cte.db2Pascal * lib._specvol_SSO_0_CT25(p) # initial of df_dp
    f     = lib._enthalpy_SSO_0_CT25(p) + gs * ( z - 0.5 * cte.gamma * ( z**2 ) )
    p_old = p
    p     = p_old - f / df_dp
    pm    = 0.5 * (p + p_old)
    df_dp = cte.db2Pascal * lib._specvol_SSO_0_CT25(pm)
    p     = p_old - f / df_dp

    return p

# ------------------------

@match_args_return
def t90_from_t48(t48):
    """
    Converts IPTS-48 temperature to ITS-90

    Parameters
    ---------
    t48 : array-like
          in-situ temperature (ITPS-48)    [deg C]

    Returns
    -------
    t90 : array-like
          in-situ temperature  (ITS-90)    [deg C]

    """

    return (t48 - (4.4e-6)*t48*(100 - t48)) / 1.00024

# -------------------------------

@match_args_return
def t90_from_t68(t68):
    """
    Converts IPTS-68 temperature to ITS-90

    Parameters
    ---------
    t68 : array-like
          in-situ temperature (ITPS-68)    [deg C]

    Returns
    -------
    t90 : array-like
          in-situ temperature  (ITS-90)    [deg C]

    """

    # t90 = t68 / 1.00024
    return t68 * 0.999760057586179

# -------------------

if __name__=='__main__':
    import doctest
    doctest.testmod()
