# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np
import seawater.constants as cte
import seawater.library as lib
from seawater.library import match_args_return

"""
Section B: functions
"""

@match_args_return
def entropy(SA, t, p):
    r"""
    Calculates specific entropy of seawater.

    The specific entropy of seawater :math:`\eta` is given by:

    .. math::
        \eta(SA, t, p) = -g_T = \frac{\partial g}{\partial T}\Big|_{SA,p}

    When taking derivatives with respect to *in situ* temperature, the symbol
    :math:`T` will be used for temperature in order that these derivatives not
    be confused with time derivatives.

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
    entropy : array_like
              specific entropy [J kg :sup:`-1` K :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.entropy(SA, t, p)
    array([ 400.38942528,  395.43817843,  319.8664982 ,  146.79088159,
             98.64734087,   62.79150873])

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

    n0, n1 = 0, 1
    entropy = -lib._gibbs(n0, n1, n0, SA, t, p)

    return entropy

@match_args_return
def rho(SA, t, p):
    r"""
    Calculates in situ density of seawater from Absolute Salinity and in situ
    temperature.

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
    rho : array_like
          in situ density [kg m :sup:`-3`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.rho(SA, t, p)
    array([ 1021.84017319,  1022.26268993,  1024.42771594,  1027.79020181,
            1029.83771473,  1032.00240412])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.8.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1
    rho = 1. / lib._gibbs(n0, n0, n1, SA, t, p)

    return rho

@match_args_return
def cp(SA, t, p):
    r"""
    Calculates the isobaric heat capacity of seawater.

    SA : array_like
        Absolute salinity [g kg :sup:`-1`]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    cp : array_like
         heat capacity of seawater [J kg :sup:`-1` K :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.cp(SA, t, p)
    array([ 4002.88800396,  4000.98028393,  3995.54646889,  3985.07676902,
            3973.59384348,  3960.18408479])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n2 = 0, 2
    cp = -( t + cte.Kelvin ) * lib._gibbs(n0, n2, n0, SA, t, p)

    return cp

@match_args_return
def helmholtz_energy(SA, t, p):
    r"""
    Calculates the Helmholtz energy of seawater.

    The specific Helmholtz energy of seawater :math:`f` is given by:

    .. math::
        f(SA, t, p) = g - (p + P_0) \nu = g - (p + P_0) \frac{\partial g}{\partial P}\Big|_{SA,T}

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
    Helmholtz_energy : array_like
                       Helmholtz energy [J kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.helmholtz_energy(SA, t, p)
    array([-5985.58288209, -5830.81845224, -3806.96617841,  -877.66369421,
            -462.17033905,  -245.50407205])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.13.

    Modifications:
    2010-08-26. Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1
    g000 = lib._gibbs(n0, n0, n0, SA, t, p)
    g001 = lib._gibbs(n0, n0, n1, SA, t, p)
    helmholtz_energy =( g000 - ( cte.db2Pascal * p + 101325 ) * g001 )

    return helmholtz_energy

@match_args_return
def internal_energy(SA, t, p):
    r"""
    Calculates the Helmholtz energy of seawater.

    The specific internal energy of seawater :math:`u` is given by:

    .. math::
        u(SA, t, p) = g + (T_0 + t)\eta - (p + P_0)\nu = g - (T_0 + t)\frac{\partial g}{\partial T}\Big|_{SA,p} - (p + P_0)\frac{\partial g}{\partial P}\Big|_{SA,T}

    where :math:`T_0` is the Celsius zero point, 273.15 K and
    :math:`P_0` = 101 325 Pa is the standard atmosphere pressure.

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
    internal_energy (u) : array_like
                          specific internal energy [J kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.internal_energy(SA, t, p)
    array([ 114906.23847309,  113426.57417062,   90860.81858842,
             40724.34005719,   27162.66600185,   17182.50522667])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqn. (2.11.1)

    Modifications:
    2010-08-22. Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1
    g000 = lib._gibbs(n0, n0, n0, SA, t, p)
    g010 = lib._gibbs(n0, n1, n0, SA, t, p)
    g001 = lib._gibbs(n0, n0, n1, SA, t, p)
    internal_energy = g000 - (cte.Kelvin + t) * g010 - ( cte.db2Pascal * p +
                                                             101325 ) * g001

    return internal_energy

@match_args_return
def sound_speed(SA, t, p):
    r"""
    Calculates the speed of sound in seawater.

    The speed of sound in seawater :math:`c` is given by:

    .. math::
        c(SA, t, p) = \sqrt{ \partial P  / \partial \rho |_{SA,\eta}} = \sqrt{(\rho\kappa)^{-1}} = g_P \sqrt{g_{TT}/(g^2_{TP} - g_{TT}g_{PP})}

    Note that in these expressions, since sound speed is in m s :sup`-1` and
    density has units of kg m :sup:`-3` it follows that the pressure of the
    partial derivatives must be in Pa and the isentropic compressibility
    :math:`kappa` must have units of Pa :sup:`-1`. The sound speed c produced
    by both the SIA and the GSW software libraries (appendices M and N) has
    units of m s :sup:`-1`.

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
    sound_speed : array_like
                  speed of sound in seawater [m s :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.sound_speed(SA, t, p)
    array([ 1542.61580359,  1542.70353407,  1530.84497914,  1494.40999692,
            1487.37710252,  1483.93460908])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqn. (2.17.1)

    Modifications:
    2010-07-23. David Jackett, Paul Barker and Trevor McDougall.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1, n2 = 0, 1, 2
    g020 = lib._gibbs(n0, n2, n0, SA, t, p)
    g011 = lib._gibbs(n0, n1, n1, SA, t, p)
    g001 = lib._gibbs(n0, n0, n1, SA, t, p)
    g002 = lib._gibbs(n0, n0, n2, SA, t, p)
    sound_speed = g001 * np.ma.sqrt( g020 / ( g011**2  - g020 * g002 ) )

    return sound_speed

@match_args_return
def adiabatic_lapse_rate(SA, t, p):
    r"""
    Calculates the adiabatic lapse rate of sea water.

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
    adiabatic_lapse_rate : array_like
                           Adiabatic lapse rate [K Pa :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    The output is in unit of degrees Celsius per Pa, (or equivalently K/Pa) not
    in units of K/dbar

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.adiabatic_lapse_rate(SA, t, p)
    array([  2.40350282e-08,   2.38496700e-08,   2.03479880e-08,
             1.19586543e-08,   9.96170718e-09,   8.71747270e-09])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqn. (2.22.1).

    Modifications:
    2010-08-26. Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1, n2 = 0, 1, 2
    g011 = lib._gibbs(n0, n1, n1, SA, t, p)
    g020 = lib._gibbs(n0, n2, n0, SA, t, p)
    adiabatic_lapse_rate = - g011 / g020

    return adiabatic_lapse_rate

@match_args_return
def chem_potential_relative(SA, t, p):
    r"""
    Calculates the adiabatic lapse rate of sea water.

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
    chem_potential_relative : array_like
                              relative chemical potential [J kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.chem_potential_relative(SA, t, p)
    array([ 79.4254481 ,  79.25989214,  74.69154859,  65.64063719,
            61.22685656,  57.21298557])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    Modifications:
    2010-08-26. Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1
    chem_potential_relative = lib._gibbs(n1, n0, n0, SA, t, p)

    return chem_potential_relative

@match_args_return
def specvol(SA, t, p):
    r"""
    Calculates the specific volume of seawater.

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
    specvol : array_like
              specific volume [m :sup:`3` kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.specvol(SA, t, p)
    array([ 0.00097863,  0.00097822,  0.00097615,  0.00097296,  0.00097103,
            0.00096899])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.7.

    Modifications:
    2010-08-26. David Jackett & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1
    specvol = lib._gibbs(n0, n0, n1, SA, t, p)

    return specvol

@match_args_return
def conservative_t(SA, t, p):
    r"""
    Calculates Conservative Temperature of seawater from in situ temperature.

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.conservative_t(SA, t, p)
    array([ 28.80991983,  28.43922782,  22.78617689,  10.22618927,
             6.82721363,   4.32357575])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 3.3.

    Modifications:
    2010-08-26. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    #pt0 = potential_t(SA, t, p) #NOTE: original call a
                                 # faster version (pt0_from_t)
                                 # instead of potential_t
    pt0 = pt0_from_t(SA, t, p)
    CT = CT_from_pt(SA, pt0)

    return CT

@match_args_return
def potential_t(SA, t, p, pr=0):
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
    >>> gsw.potential_t(SA, t, p)
    array([ 28.78319682,  28.42098334,  22.7849304 ,  10.23052366,
             6.82923022,   4.32451057])
    >>> gsw.potential_t(SA, t, p, pr = 1000)
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

@match_args_return
def enthalpy(SA, t, p):
    r"""
    Calculates the specific enthalpy of seawater.

    The specific enthalpy of seawater :math:`h` is given by:

    .. math::
        h(SA, t, p) = g + (T_0 + t)\eta = g - (T_0 + t) \frac{\partial g}{\partial T}\Big|_{SA,p}

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
    enthalpy : array_like
               specific enthalpy [J kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.enthalpy(SA, t, p)
    array([ 115103.26047838,  114014.8036012 ,   92179.9209311 ,
             43255.32838089,   33087.21597002,   26970.5880448 ])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See appendix A.11.

    Modifications:
    2010-08-26. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1
    g000 = lib._gibbs(n0, n0, n0, SA, t, p)
    g010 = lib._gibbs(n0, n1, n0, SA, t, p)
    enthalpy = g000 - (t + cte.Kelvin) * g010

    return enthalpy

@match_args_return
def chem_potential_water(SA, t, p):
    r"""
    Calculates the chemical potential of water in seawater.

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
    chem_potential_water : array_like
                           chemical potential of water in seawater
                           [J kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.chem_potential_water(SA, t, p)
    array([-8545.56114628, -8008.08554834, -5103.98013987,  -634.06778275,
            3335.56680347,  7555.43444597])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    Modifications:
    2010-09-28. Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """
    SA, t, p, mask = lib.strip_mask(SA, t, p)

    x2 = cte.sfac * SA

    x = np.sqrt(x2)
    y = t * 0.025
    z = p * 1e-4 # pressure (p) is sea pressure in units of dbar

    g03_g = ( 101.342743139674 + z * ( 100015.695367145 +
        z * ( -2544.5765420363 + z * ( 284.517778446287 +
        z * ( -33.3146754253611 + ( 4.20263108803084 -
        0.546428511471039 * z ) * z ) ) ) ) +
        y * ( 5.90578347909402 + z * ( -270.983805184062 +
        z * ( 776.153611613101 + z * ( -196.51255088122 +
        ( 28.9796526294175 - 2.13290083518327 * z ) * z ) ) ) +
        y * ( -12357.785933039 + z * ( 1455.0364540468 +
        z * ( -756.558385769359 + z * ( 273.479662323528 +
        z * ( -55.5604063817218 + 4.34420671917197 * z ) ) ) ) +
        y * ( 736.741204151612 + z * ( -672.50778314507 +
        z * ( 499.360390819152 + z * ( -239.545330654412 +
        ( 48.8012518593872 - 1.66307106208905 * z ) * z ) ) ) +
        y * ( -148.185936433658 + z * ( 397.968445406972 +
        z * ( -301.815380621876 + ( 152.196371733841 -
        26.3748377232802 * z ) * z ) ) +
        y * ( 58.0259125842571 + z * ( -194.618310617595 +
        z * ( 120.520654902025 + z * ( -55.2723052340152 +
        6.48190668077221 * z ) ) ) +
        y * ( -18.9843846514172 + y * ( 3.05081646487967 -
        9.63108119393062 * z ) +
        z * ( 63.5113936641785 + z * ( -22.2897317140459 +
        8.17060541818112 * z ) ) ) ) ) ) ) ) )

    g08_g = x2 * ( 1416.27648484197 +
        x * ( -2432.14662381794 + x * ( 2025.80115603697 +
        y * ( 543.835333000098 + y * ( -68.5572509204491 +
        y * ( 49.3667694856254 + y * ( -17.1397577419788 +
        2.49697009569508 * y ) ) ) - 22.6683558512829 * z ) +
        x * ( -1091.66841042967 - 196.028306689776 * y +
        x * ( 374.60123787784 - 48.5891069025409 * x +
        36.7571622995805 * y ) + 36.0284195611086 * z ) +
        z * ( -54.7919133532887 + ( -4.08193978912261 -
        30.1755111971161 * z ) * z ) ) +
        z * ( 199.459603073901 + z * ( -52.2940909281335 +
        ( 68.0444942726459 - 3.41251932441282 * z ) * z ) ) +
        y * ( -493.407510141682 + z * ( -175.292041186547 +
        ( 83.1923927801819 - 29.483064349429 * z ) * z ) +
        y * ( -43.0664675978042 + z * ( 383.058066002476 +
        z * ( -54.1917262517112 + 25.6398487389914 * z ) ) +
        y * ( -10.0227370861875 - 460.319931801257 * z + y *
        ( 0.875600661808945 + 234.565187611355 * z ) ) ) ) ) +
        y * ( 168.072408311545 ) )

    g_SA_part = ( 8645.36753595126 +
        x * ( -7296.43987145382 + x * ( 8103.20462414788 +
        y * ( 2175.341332000392 + y * ( -274.2290036817964 +
        y * ( 197.4670779425016 + y * ( -68.5590309679152 +
        9.98788038278032 * y ) ) ) - 90.6734234051316 * z ) +
        x * ( -5458.34205214835 - 980.14153344888 * y +
        x * ( 2247.60742726704 - 340.1237483177863 * x +
        220.542973797483 * y ) + 180.142097805543 * z ) +
        z * ( -219.1676534131548 + ( -16.32775915649044 -
        120.7020447884644 * z ) * z ) ) +
        z * ( 598.378809221703 + z * ( -156.8822727844005 +
        ( 204.1334828179377 - 10.23755797323846 * z ) * z ) ) +
        y * ( -1480.222530425046 + z * ( -525.876123559641 +
        ( 249.57717834054571 - 88.449193048287 * z ) * z ) +
        y * ( -129.1994027934126 + z * ( 1149.174198007428 +
        z * ( -162.5751787551336 + 76.9195462169742 * z ) ) +
        y * ( -30.0682112585625 - 1380.9597954037708 * z + y *
        ( 2.626801985426835 + 703.695562834065 * z ) ) ) ) ) +
        y * ( 1187.3715515697959) )

    chem_potential_water =  g03_g + g08_g - 0.5 * cte.sfac * SA * g_SA_part

    return np.ma.array(chem_potential_water, mask=mask, copy=False)

@match_args_return
def chem_potential_salt(SA, t, p):
    r"""
    Calculates the chemical potential of salt in seawater.

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
    chem_potential_salt : array_like
                          chemical potential of salt in seawater [J kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.chem_potential_salt(SA, t, p)
    array([-8466.13569818, -7928.8256562 , -5029.28859129,  -568.42714556,
            3396.79366004,  7612.64743154])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.9.

    Modifications:
    2010-09-28. Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    chem_potential_salt = ( chem_potential_relative(SA, t, p) +
                            chem_potential_water(SA, t, p) )

    return chem_potential_salt

@match_args_return
def isochoric_heat_cap(SA, t, p):
    r"""
    Calculates the isochoric heat capacity of seawater.

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
    isochoric_heat_cap : array_like
                         isochoric heat capacity [J kg :sup:`-1` K :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.isochoric_heat_cap(SA, t, p)
    array([ 3928.13708702,  3927.27381633,  3941.36418525,  3966.26126146,
            3960.50903222,  3950.13901342])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.21.

    Modifications:
    2010-08-26. Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1, n2 = 0, 1, 2

    g020 = lib._gibbs(n0, n2, n0, SA, t, p)
    g011 = lib._gibbs(n0, n1, n1, SA, t, p)
    g002 = lib._gibbs(n0, n0, n2, SA, t, p)

    isochoric_heat_cap = -(cte.Kelvin + t) * (g020 - g011 * g011 / g002)

    return isochoric_heat_cap

@match_args_return
def kappa(SA, t, p):
    r"""
    Calculates the isentropic compressibility of seawater.

    When the entropy and Absolute Salinity are held constant while the pressure
    is changed, the isentropic and isohaline compressibility
    :math:`kappa` is obtained:

    .. math::
        \kappa(SA, t, p) = \rho^{-1}\frac{\partial \rho}{\partial P}\Big|_{SA,\eta} = -\nu^{-1}\frac{\partial \nu}{\partial P}\Big|_{SA,\eta} = \rho^{-1}\frac{\partial \rho}{\partial P}\Big|_{SA,\theta} = -\nu^{-1}\frac{\partial \nu}{\partial P}\Big|_{SA,\theta} =
        -\frac{ (g_{TP}^2 - g_{TT} g_{PP} ) }{g_P g_{TT}}

    The isentropic and isohaline compressibility is sometimes called simply the
    isentropic compressibility (or sometimes the "adiabatic compressibility"),
    on the unstated understanding that there is also no transfer of salt during
    the isentropic or adiabatic change in pressure.

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
    kappa : array_like
            Isentropic compressibility [Pa :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    The output is Pascal and not dbar.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.kappa(SA, t, p)
    array([  4.11245799e-10,   4.11029072e-10,   4.16539558e-10,
             4.35668338e-10,   4.38923693e-10,   4.40037576e-10])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (2.16.1) and the row for kappa in
    Table P.1 of appendix P

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1, n2 = 0, 1, 2

    g020 = lib._gibbs(n0, n2, n0, SA, t, p)
    g011 = lib._gibbs(n0, n1, n1, SA, t, p)
    g002 = lib._gibbs(n0, n0, n2, SA, t, p)
    g001 = lib._gibbs(n0, n0, n1, SA, t, p)

    kappa = ( g011 * g011 - g020 * g002 ) / ( g001 * g020 )

    return kappa

@match_args_return
def kappa_const_t(SA, t, p):
    r"""
    Calculates isothermal compressibility of seawater at constant in situ
    temperature.

    .. math::
        \kappa^t(SA, t, p) = \rho^{-1}\frac{\partial \rho}{\partial P}\Big|_{SA,T} = -\nu^{-1}\frac{\partial \nu}{\partial P}\Big|_{SA,T} = -\frac{g_{PP}}{g_P}

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
    kappa : array_like
            Isothermal compressibility [Pa :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    This is the compressibility of seawater at constant in situ temperature.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.kappa_const_t(SA, t, p)
    array([  4.19071646e-10,   4.18743202e-10,   4.22265764e-10,
             4.37735100e-10,   4.40373818e-10,   4.41156577e-10])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqn. (2.15.1)

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1, n2 = 0, 1, 2
    g002 = lib._gibbs(n0, n0, n2, SA, t, p)
    g001 = lib._gibbs(n0, n0, n1, SA, t, p)
    kappa = -g002 / g001

    return kappa

@match_args_return
def pot_rho(SA, t, p, pr=0):
    r"""
    Calculates potential density of seawater.

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
    pot_rho : array_like
              potential density  [kg m :sup:`-3`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.pot_rho(SA, t, p)
    array([ 1021.79814581,  1022.05248442,  1023.89358365,  1026.66762112,
            1027.10723087,  1027.40963126])
    >>> gsw.pot_rho(SA, t, p, pr=1000)
    array([ 1025.95554512,  1026.21306986,  1028.12563226,  1031.1204547 ,
            1031.63768355,  1032.00240412])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 3.4.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    pt = potential_t(SA, t, p, pr=pr)

    pot_rho = rho(SA, pt, pr)

    return pot_rho

@match_args_return
def specvol_anom(SA, t, p):
    r"""
    Calculates specific volume anomaly from Absolute Salinity, in situ
    temperature and pressure, using the full TEOS-10 Gibbs function.

    The reference value of Absolute Salinity is SSO and the reference value of
    Conservative Temperature is equal to 0 :math:`^\circ` C.

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
    specvol_anom : array_like
                   specific volume anomaly  [m :sup:`3` kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.specvol_anom(SA, t, p)
    array([  6.01044463e-06,   5.78602432e-06,   4.05564999e-06,
             1.42198662e-06,   1.04351837e-06,   7.63964850e-07])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqn. (3.7.3)

    Modifications:
    2010-08-26. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1

    SSO = cte.SSO

    CT0 = 0

    pt_zero = pt_from_CT(SSO, CT0)

    t_zero = potential_t(SSO, pt_zero, 0, p)

    specvol_anom = ( lib._gibbs(n0, n0, n1, SA, t, p) -
                     lib._gibbs(n0, n0, n1, SSO, t_zero, p) )

    return specvol_anom

@match_args_return
def alpha_wrt_t(SA, t, p):
    r"""
    Calculates the thermal expansion coefficient of seawater with respect to in
    situ temperature.

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
    alpha_wrt_t : array_like
                  thermal expansion coefficient [K :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.alpha_wrt_t(SA, t, p)
    array([ 0.0003256 ,  0.00032345,  0.00028141,  0.00017283,  0.00014557,
            0.00012836])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqn. (2.18.1)

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm
    for estimating Absolute Salinity in the global ocean. Submitted to Ocean
    Science. A preliminary version is available at Ocean Sci. Discuss.,
    6, 215-242.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1

    g011 = lib._gibbs(n0, n1, n1, SA, t, p)
    g001 = lib._gibbs(n0, n0, n1, SA, t, p)
    alpha_wrt_t = g011 / g001

    return alpha_wrt_t

@match_args_return
def alpha_wrt_CT(SA, t, p):
    r"""
    Calculates the thermal expansion coefficient of seawater with respect to
    Conservative Temperature.

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
    alpha_wrt_CT : array_like
                   thermal expansion coefficient [K :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.alpha_wrt_CT(SA, t, p)
    array([ 0.00032471,  0.00032272,  0.00028118,  0.00017314,  0.00014627,
            0.00012943])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqn. (2.18.3).

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1, n2 = 0, 1, 2

    # pt0 = potential_t(SA, t, p) #NOTE: original call a faster version
    #"pt0_from_t" instead of potential_t
    pt0 = pt0_from_t(SA, t, p)

    fac = -cte.cp0 / ( (cte.Kelvin + pt0) * lib._gibbs(n0, n2, n0, SA, t, p ) )

    g011 = lib._gibbs(n0, n1, n1, SA, t, p)
    g001 = lib._gibbs(n0, n0, n1, SA, t, p)
    alpha_wrt_CT = fac * ( g011 / g001 )

    return alpha_wrt_CT

@match_args_return
def alpha_wrt_pt(SA, t, p):
    r"""
    Calculates the thermal expansion coefficient of seawater with respect to
    potential temperature, with a reference pressure of zero.

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
    alpha_wrt_pt : array_like
                   thermal expansion coefficient [K :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.alpha_wrt_pt(SA, t, p)
    array([ 0.00032562,  0.00032355,  0.00028164,  0.00017314,  0.00014623,
            0.00012936])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqn. (2.18.2).

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1, n2 = 0, 1, 2

    #pt0 = potential_t(SA, t, p) #NOTE: original call a faster version
    # (pt0_from_t) instead of potential_t
    pt0 = pt0_from_t(SA, t, p)

    fac = lib._gibbs(n0, n2, n0, SA, pt0, 0) / lib._gibbs(n0, n2, n0, SA, t, p)
    g011 = lib._gibbs(n0, n1, n1, SA, t, p)
    g001 = lib._gibbs(n0, n0, n1, SA, t, p)
    alpha_wrt_pt = fac * ( g011 / g001 )

    return alpha_wrt_pt

@match_args_return
def beta_const_t(SA, t, p):
    r"""
    Calculates the saline (i.e. haline) contraction coefficient of seawater at
    constant in situ temperature.

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
    beta_const_t : array_like
                   saline contraction coefficient [kg g :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.beta_const_t(SA, t, p)
    array([ 0.00073112,  0.00073107,  0.00073602,  0.00075381,  0.00075726,
            0.00075865])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqn. (2.19.1)

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1

    g101 = lib._gibbs(n1, n0, n1, SA, t, p)
    g001 = lib._gibbs(n0, n0, n1, SA, t, p)
    beta_const_t = -g101 / g001

    return beta_const_t

@match_args_return
def beta_const_CT(SA, t, p):
    r"""
    Calculates the saline (i.e. haline) contraction coefficient of seawater at
    constant Conservative Temperature.

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
    beta_const_CT : array_like
                    saline contraction coefficient [kg g :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.beta_const_CT(SA, t, p)
    array([ 0.00071749,  0.00071765,  0.00072622,  0.00075051,  0.00075506,
            0.00075707])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqn. (2.19.3)

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1, n2 = 0, 1, 2

    # pt0 = potential_t(SA, t, p) #NOTE: original call a faster version
    # (pt0_from_t) instead of potential_t
    pt0 = pt0_from_t(SA, t, p)

    g001 = lib._gibbs(n0, n0, n1, SA, t, p)
    g110 = lib._gibbs(n1, n1, n0, SA, t, p)
    g100 = lib._gibbs(n1, n0, n0, SA, pt0, 0)

    factora = g110 - g100 / (cte.Kelvin + pt0)
    g020 = lib._gibbs(n0, n2, n0, SA, t, p)
    factor = factora / ( g001 * g020 )

    g011 = lib._gibbs(n0, n1, n1, SA, t, p)
    g101 = lib._gibbs(n1, n0, n1, SA, t, p)
    beta_const_CT = g011 * factor - g101 / g001

    return beta_const_CT

@match_args_return
def beta_const_pt(SA, t, p):
    r"""
    Calculates the saline (i.e. haline) contraction coefficient of seawater at
    constant potential temperature with a reference pressure of 0 dbar.

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
    beta_const_pt : array_like
                    saline contraction coefficient [kg g :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.beta_const_pt(SA, t, p)
    array([ 0.00073112,  0.00073106,  0.00073599,  0.00075375,  0.00075712,
            0.00075843])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqn. (2.19.2)

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1, n2 = 0, 1, 2

    # pt0 = potential_t(SA, t, p) #NOTE: original call a faster version
    # "pt0_from_t" instead of potential_t
    pt0 = pt0_from_t(SA, t, p)

    g001 = lib._gibbs(n0, n0, n1, SA, t, p)
    g110 = lib._gibbs(n1, n1, n0, SA, t, p)
    factora = g110 - lib._gibbs(n1, n1, n0, SA, pt0, 0)

    g020 = lib._gibbs(n0, n2, n0, SA, t, p)
    factor = factora / ( g001 * g020 )

    g011 = lib._gibbs(n0, n1, n1, SA, t, p)
    g101 = lib._gibbs(n1, n0, n1, SA, t, p)

    beta_const_pt = g011 * factor - g101 / g001

    return beta_const_pt

@match_args_return
def osmotic_coefficient(SA, t, p):
    r"""
    Calculates the osmotic coefficient of seawater.

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
    osmotic_coefficient : array_like
                          osmotic coefficient of seawater [unitless]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.osmotic_coefficient(SA,t , p)
    array([ 0.90284718,  0.90298624,  0.90238866,  0.89880927,  0.89801054,
            0.89767912])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    Modifications:
    2010-09-28. Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0 = 0

    molal = molality(SA) # molality of seawater in mol kg :sup:`-1`
    part = molal * cte.R * ( cte.Kelvin + t )

    SAzero = 0
    g000 = lib._gibbs(n0, n0, n0, SAzero, t, p)
    osmotic_coefficient = ( g000 - chem_potential_water(SA, t, p) ) / part

    return osmotic_coefficient

@match_args_return
def molality(SA):
    r"""
    Calculates the molality of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]

    Returns
    -------
    molal : array_like
            seawater molality [mol kg :sup:`-1`]

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
    >>> gsw.molality(SA)
    array([ 1.14508476,  1.15122708,  1.15581223,  1.14971265,  1.14593231,
            1.14578877])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    Modifications:
    2010-09-28. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # molality of seawater in mol kg :sup:`-1`
    molal = SA / (cte.M_S * ( 1000 - SA ) )
    molal[SA < 0] = np.ma.masked

    return molal

@match_args_return
def ionic_strength(SA):
    r"""
    Calculates the ionic strength of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]

    Returns
    -------
    ionic_strength : array_like
                     ionic strength of seawater [mol kg :sup:`-1`]

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
    >>> gsw.ionic_strength(SA)
    array([ 0.71298118,  0.71680567,  0.71966059,  0.71586272,  0.71350891,
            0.71341953])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Table L.1.

    .. [2] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008:
    The composition of Standard Seawater and the definition of the
    Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.
    See Eqns. 5.9 and 5.12.

    Modifications:
    2010-09-28. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    Z_2 = 1.2452898 # the valence factor of sea salt

    molal = molality(SA) # molality of seawater in mol kg :sup:`-1`

    ionic_strength = 0.5 * Z_2 * molal

    return ionic_strength

"""
Section C: extra functions for Temperature
NOTE: section B calls: CT_from_pt(SA, pt0), pt_from_CT(SSO, CT0), pt0_from_t(SA,t, p)
pt_from_t(SA, t, p, pr=0) was renamed to potential_t(SA, t, p, pr=0) in Section B
"""

@match_args_return
def CT_from_pt(SA, pt): #NOTE: used in conservative_t(SA, t, p)
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

    t = potential_t(SA, pt0, 0, p)

    return t

@match_args_return
def pt0_from_t(SA, t, p):
    r"""
    Calculates potential temperature with reference pressure, pr = 0 dbar.
    The present routine is computationally faster than the more general
    function "potential_t(SA, t, p, pr)" which can be used for any reference
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
    potential_t  has the same result (only slower)

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

@match_args_return
def pt_from_entropy(SA, entropy):
    r"""
    Calculates potential temperature with reference pressure pr = 0 dbar
    from entropy.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    entropy : array_like
              specific entropy [J kg :sup:`-1` K :sup:`-1`]

    Returns
    -------
    pt : array_like
         potential temperature [:math:`^\circ` C (ITS-90)] (Default) or

    See Also
    --------
    _gibbs_pt0_pt0

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> entropy = [400.3892, 395.4378, 319.8668, 146.7910, 98.6471, 62.7919]
    >>> gsw.pt_from_entropy(SA, entropy)
    array([ 28.78317983,  28.42095483,  22.78495274,  10.23053207,
             6.82921333,   4.32453778])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See appendix  A.10.

    Modifications:
    2010-10-13. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA.clip(0, np.inf)
    #SA[SA < 0] = 0

    n0, n1 = 0, 1

    part1 = 1 - SA / cte.SSO
    part2 = 1 - 0.05 * part1
    ent_SA = (cte.cp0 / cte.Kelvin) * part1 * ( 1 - 1.01 * part1)
    c = (entropy - ent_SA) * part2 / cte.cp0
    pt = cte.Kelvin * (np.exp(c) - 1)
    dentropy_dt = cte.cp0 / ( (cte.Kelvin + pt) * part2) # initial dentropy_dt

    for Number_of_iterations in range(0,3):
        pt_old = pt
        dentropy = entropy_from_pt(SA, pt_old) - entropy
        pt = pt_old - dentropy / dentropy_dt # half way through mod. method
        ptm = 0.5 * (pt + pt_old)
        dentropy_dt = -lib._gibbs_pt0_pt0(SA, ptm)
        pt = pt_old - dentropy / dentropy_dt

    return pt

@match_args_return
def CT_from_entropy(SA, entropy):
    r"""
    Calculates potential temperature with reference pressure pr = 0 dbar or
    Conservative temperature from entropy.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    entropy : array_like
              specific entropy [J kg :sup:`-1` K :sup:`-1`]

    Returns
    -------
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    See Also
    --------
    _gibbs_pt0_pt0

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> entropy = [400.3892, 395.4378, 319.8668, 146.7910, 98.6471, 62.7919]
    >>> gsw.CT_from_entropy(SA, entropy)
    array([ 28.80990279,  28.43919923,  22.78619927,  10.22619767,
             6.82719674,   4.32360295])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See appendix  A.10.

    Modifications:
    2010-10-13. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA.clip(0, np.inf)
    #SA[SA < 0] = 0

    n0, n1 = 0, 1

    pt = pt_from_entropy(SA, entropy)
    CT = CT_from_pt(SA, pt)

    return CT

@match_args_return
def entropy_from_CT(SA, CT):
    r"""
    Calculates specific entropy of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    entropy : array_like
              specific entropy [J kg :sup:`-1` K :sup:`-1`]

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
    >>> gsw.entropy_from_CT(SA, CT)
    array([ 400.38916315,  395.43781023,  319.86680989,  146.79103279,
             98.64714648,   62.79185763])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See appendix A.10.

    Modifications:
    2010-10-13. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    #SA[SA < 0] = 0
    SA.clip(0, np.inf)

    n0, n1 = 0, 1

    pt0 = pt_from_CT(SA, CT)
    entropy = entropy_from_pt(SA, pt0)

    return entropy

@match_args_return
def entropy_from_pt(SA, pt):
    r"""
    Calculates specific entropy of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    entropy : array_like
              specific entropy [J kg :sup:`-1` K :sup:`-1`]

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
    >>> pt = [28.7832, 28.4210, 22.7850, 10.2305, 6.8292, 4.3245]
    >>> gsw.entropy_from_pt(SA, pt)
    array([ 400.38946744,  395.43839949,  319.86743859,  146.79054828,
             98.64691006,   62.79135672])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See appendix A.10.

    Modifications:
    2010-10-13. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    #SA[SA < 0] = 0
    SA.clip(0, np.inf)

    n0, n1 = 0, 1

    entropy = -lib._gibbs(n0, n1, n0, SA, pt, 0)

    return entropy

@match_args_return
def temps_maxdensity(SA, p):
    r"""
    Calculates the temperatures of maximum density of seawater. This function
    returns the in-situ, potential, and Conservative temperatures at which the
    density of seawater is a maximum, at given Absolute Salinity, SA, and sea
    pressure, p (in dbar).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    p : array_like
        pressure [dbar]

    Returns
    -------
    t_maxden : array_like
               in situ temperature [:math:`^\circ` C (ITS-90)]
    pt_maxden : array_like
                potential temperature [:math:`^\circ` C (ITS-90)]
    CT_maxden : array_like
                Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    See Also
    --------
    TODO

    Notes
    -----
    Returning only in-situ t for now. I don't see the point in returning all 3.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.temps_maxdensity(SA, p)
    array([-3.72500891, -3.85371429, -4.05143318, -4.29554251, -5.07474662,
           -6.01170197])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 3.42

    Modifications:
    2010-08-26. Trevor McDougall & Paul Barker.
    2010-02-21. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, p, mask = lib.strip_mask(SA, p)

    n0, n1 = 0, 1

    dt = 0.001 # the temperature increment for the gibbs_PTT derivative
    t = 3.978 - 0.22072*SA # the initial guess of t_maxden
    gibbs_PTT = 1.1e-8 # the initial guess for g_PTT

    for Number_of_iterations in range(0,3,1):
        t_old = t
        gibbs_PT = lib._gibbs(n0, n1, n1, SA, t_old, p)
        t = t_old - gibbs_PT / gibbs_PTT # half way through the modified method
        t_mean = 0.5*(t + t_old)
        gibbs_PTT = ( (lib._gibbs(n0, n1, n1, SA, t_mean + dt,p) -
        lib._gibbs(n0, n1, n1, SA, t_mean - dt, p) ) / (dt + dt) )
        t = t_old - gibbs_PT / gibbs_PTT

    # After three iterations of this modified Newton-Raphson iteration, the
    #  error in t_maxden is typically no larger than 1x10^-15 degC
    #t_maxden  = t
    #pt_maxden = gsw_pt0_from_t(SA,t_maxden,p)
    #CT_maxden = gsw_CT_from_pt(SA,pt_maxden)

    #return t_maxden, pt_maxden, CT_maxden
    return np.ma.array(t, mask=mask, copy=False)


@match_args_return
def entropy_first_derivative_SA(SA, CT):
    r"""
    Calculates the derivatives of specific entropy (eta_SA) with respect to
    Absolute Salinity at constant Conservative Temperature.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_SA : array_like
             The derivative of specific entropy with respect to SA at constant
             CT [J g :sup:`-1` K :sup:`-1`]

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
    >>> gsw.entropy_first_derivative_SA(SA, CT)
    array([ -0.2632868 ,  -0.26397728,  -0.2553675 ,  -0.23806659,
            -0.23443826,  -0.23282068])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.8) and (P.14a,c).

    Modifications:
    2010-08-21. Trevor McDougall.
    2010-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1

    pt = pt_from_CT(SA, CT)

    eta_SA = -(lib._gibbs(n1, n0, n0, SA, pt, 0) ) / (cte.Kelvin + pt)

    return eta_SA


@match_args_return
def entropy_first_derivative_CT(SA, CT):
    r"""
    Calculates the derivatives of specific entropy (eta_CT) with respect to
    Conservative Temperature at constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_CT : array_like
             The derivative of specific entropy with respect to CT at constant
             SA [ J (kg K :sup:`-2`) :sup:`-1` ]

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
    >>> gsw.entropy_first_derivative_CT(SA, CT)
    array([13.22103121,  13.23691119,  13.48900463,  14.08659902,
           14.25772958,  14.38642995])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.8) and (P.14a,c).

    Modifications:
    2010-08-21. Trevor McDougall.
    2010-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    pt = pt_from_CT(SA, CT)

    eta_CT = cte.cp0 / (cte.Kelvin + pt)

    return eta_CT

def entropy_first_derivatives(SA, CT):
    r"""
    Calculates the following two partial derivatives of specific entropy (eta)
    (1) eta_SA, the derivative with respect to Absolute Salinity at constant
        Conservative Temperature, and
    (2) eta_CT, the derivative with respect to Conservative Temperature at
        constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_SA : array_like
             The derivative of specific entropy with respect to SA at constant
             CT [J g :sup:`-1` K :sup:`-1`]
    eta_CT : array_like
             The derivative of specific entropy with respect to CT at constant
             SA [ J (kg K :sup:`-2`) :sup:`-1` ]

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
    >>> gsw.entropy_first_derivatives(SA, CT)
    array([[ -0.2632868 ,  -0.26397728,  -0.2553675 ,  -0.23806659,
             -0.23443826,  -0.23282068],
           [ 13.22103121,  13.23691119,  13.48900463,  14.08659902,
             14.25772958,  14.38642995]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.8) and (P.14a,c).

    Modifications:
    2010-08-21. Trevor McDougall.
    2010-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    return ( entropy_first_derivative_SA(SA, CT),
             entropy_first_derivative_CT(SA, CT) )

@match_args_return
def CT_first_derivative_SA(SA, pt):
    r"""
    Calculates the derivatives of Conservative Temperature (CT_SA) with respect
    to Absolute Salinity at constant potential temperature (with pr = 0 dbar).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    CT_SA : array_like
            The derivative of CT with respect to SA at constant potential
            temperature reference sea pressure of 0 dbar.
            [K (g kg :sup:`-1`) :sup:`-1`]

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
    >>> gsw.CT_first_derivative_SA(SA, pt)
    array([-0.04198109, -0.04155814, -0.03473921, -0.0187111 , -0.01407594,
           -0.01057172])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.3) and (A.12.9a,b).

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
    Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    Modifications:
    2010-08-05. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1
    abs_pt = cte.Kelvin + pt

    g100 = lib._gibbs(n1, n0, n0, SA, pt, 0)
    g110 = lib._gibbs(n1, n1, n0, SA, pt, 0)
    CT_SA = ( g100 - abs_pt * g110 ) / cte.cp0

    return CT_SA

@match_args_return
def CT_first_derivative_pt(SA, pt):
    r"""
    Calculates the derivatives of Conservative Temperature (CT_pt) with respect
    to potential temperature (the regular potential temperature which is
    referenced to 0 dbar) at constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    CT_pt : array_like
            The derivative of CT with respect to pt at constant SA.
            [ unitless ]

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
    >>> gsw.CT_first_derivative_pt(SA, pt)
    array([1.00281494,  1.00255482,  1.00164514,  1.00000377,  0.99971636,
           0.99947433])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.3) and (A.12.9a,b).

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
    Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    Modifications:
    2010-08-05. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n2 = 0, 2
    abs_pt = cte.Kelvin + pt

    g020 = lib._gibbs(n0, n2, n0, SA, pt, 0)
    CT_pt = - (abs_pt * g020 ) / cte.cp0

    return CT_pt

def CT_first_derivatives(SA, pt):
    r"""
    Calculates the following two derivatives of Conservative Temperature
    (1) CT_SA, the derivative with respect to Absolute Salinity at constant
        potential temperature (with pr = 0 dbar), and
    (2) CT_pt, the derivative with respect to potential temperature (the
        regular potential temperature which is referenced to 0 dbar) at
        constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    CT_SA : array_like
            The derivative of CT with respect to SA at constant potential
            temperature reference sea pressure of 0 dbar.
            [K (g kg :sup:`-1`) :sup:`-1`]

    CT_pt : array_like
            The derivative of CT with respect to pt at constant SA.
            [ unitless ]

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
    >>> gsw.CT_first_derivatives(SA, pt)
    array([[-0.04198109, -0.04155814, -0.03473921, -0.0187111 , -0.01407594,
            -0.01057172],
           [ 1.00281494,  1.00255482,  1.00164514,  1.00000377,  0.99971636,
             0.99947433]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.3) and (A.12.9a,b).

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
    Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    Modifications:
    2010-08-05. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    return CT_first_derivative_SA(SA, pt), CT_first_derivative_pt(SA, pt)


"""
Section D: extra functions for Depth, Pressure and Distance
"""

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

    X     = np.sin( np.deg2rad(lat) )
    sin2  = X**2
    B     = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2 )
    A     = -0.5 * cte.gamma * B
    C     = lib._enthalpy_SSO_0_CT25(p)
    z     = -2 * C / ( B + np.sqrt( B**2 - 4 * A * C ) )

    return z

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

    X     = np.sin( np.deg2rad(lat) )
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

@match_args_return
def grav(lat, p=0):
    r"""
    Calculates acceleration due to gravity as a function of latitude and as a
    function of pressure in the ocean.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    p : number or array_like. Default p = 0
        pressure [dbar]

    Returns
    -------
    g : array_like
        gravity [m s :sup:`2`]

    See Also
    --------
    TODO

    Notes
    -----
    In the ocean z is negative.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> lat = [-90, -60, -30, 0]
    >>> p = 0
    >>> gsw.grav(lat, p)
    array([ 9.83218621,  9.81917886,  9.79324926,  9.780327  ])
    >>> gsw.grav(45)
    9.8061998770458008

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    .. [2] Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.

    .. [3] Saunders, P.M., and N.P. Fofonoff (1976) Conversion of pressure to
    depth in the ocean. Deep-Sea Res.,pp. 109 - 111.

    Modifications:
    2010-07-23. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    X = np.sin( np.deg2rad(lat) )
    sin2 = X**2
    gs = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2)
    z = z_from_p(p, lat)
    grav = gs * (1 - cte.gamma * z) # z is the height corresponding to p
    return grav

@match_args_return
def distance(lon, lat, p=0):
    r"""
    Calculates the distance in met res between successive points in the vectors
    lon and lat, computed using the Haversine formula on a spherical earth of
    radius 6,371 km, being the radius of a sphere having the same volume as
    Earth. For a spherical Earth of radius 6,371,000 m, one nautical mile is
    1,853.2488 m, thus one degree of latitude is 111,194.93 m.

    Haversine formula:
        R = earth's radius (mean radius = 6,371 km)

    .. math::
        a = \sin^2(\delta \text{lat}/2) + \cos(\text{lat}_1) \cos(\text{lat}_2) \sin^2(\delta \text{lon}/2)

        c = 2 \times \text{atan2}(\sqrt{a}, \sqrt{(1-a)})

        d = R \times c

    Parameters
    ----------
    lon : array_like
          decimal degrees east [0..+360] or [-180 ... +180]
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    p : number or array_like. Default p = 0
        pressure [dbar]

    Returns
    -------
    dist: array_like
          distance between points on a spherical Earth at pressure (p) [m]

    See Also
    --------
    TODO

    Notes
    -----
    z is height and is negative in the oceanographic.

    Distances are probably good to better than 1\% of the "true" distance on the
    ellipsoidal earth.

    The check value below differ from the original online docs at
    "http://www.teos-10.org/pubs/gsw/html/gsw_distance.html" but agree with the
    result.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> lon = [159, 220]
    >>> lat = [-35, 35]
    >>> gsw.distance(lon, lat)
    array([[ 10030974.652916]])
    >>> p = [200, 1000]
    >>> gsw.distance(lon, lat, p)
    array([[ 10030661.63878009]])
    >>> p = [[200], [1000]]
    >>> gsw.distance(lon, lat, p)
    array([[ 10030661.63878009],
           [ 10029412.58776001]])

    References
    ----------
    .. [1] http://www.eos.ubc.ca/~rich/map.html

    Modifications:
    2000-11-06. Rich Pawlowicz
    2010-07-28. Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """
    #FIXME? The argument handling seems much too complicated.
    # Maybe we can come up with some simple specifications of
    # what argument combinations are permitted, and handle everything
    # with broadcasting. - EF

    #FIXME: Eric what do you think? This assume p(stations, depth)
    lon, lat, = np.atleast_2d(lon), np.atleast_2d(lat)

    if (lon.size == 1) & (lat.size == 1):
        raise ValueError('more than one point is needed to compute distance')
    elif lon.ndim != lat.ndim:
        raise ValueError('lon, lat must have the same dimension')

    lon, lat, p = np.broadcast_arrays(lon, lat, p)

    dlon = np.deg2rad( np.diff(lon) )
    dlat = np.deg2rad( np.diff(lat) )

    a = ( ( np.sin(dlat/2.) )**2 + np.cos( np.deg2rad( lat[:,:-1] ) ) *
    np.cos( np.deg2rad( lat[:,1:] ) ) * ( np.sin(dlon/2.) )**2 )

    angles = 2. * np.arctan2( np.ma.sqrt(a), np.ma.sqrt(1-a) )

    p_mid = 0.5 * (   p[:,0:-1] +   p[:,0:-1] )
    lat_mid = 0.5 * ( lat[:,:-1] + lat[:,1:] )

    z = z_from_p(p_mid, lat_mid)

    distance = (cte.a + z) * angles

    return distance

"""
Section E: extra functions for Salinity
"""

"""cndr_from_SP == sw.cndr
Calculates Practical Salinity from conductivity ratio (R), using the PSS-78
algorithm. Note that the PSS-78 algorithm for Practical Salinity
is only valid in the range 2 < SP < 42.  The output, SP, of this function is
constrained to be non-negative.
"""
from  seawater.csiro import cndr as cndr_from_SP

"""SP_from_cndr == sw.salt
Calculates conductivity ratio (R) from (SP,t,p) using PSS-78. Note that the
PSS-78 algorithm for Practical Salinity is only valid in the range 2 < SP < 42.
"""
from  seawater.csiro import salt as SP_from_cndr

@match_args_return
def SA_from_SP(SP, p, lon, lat):
    r"""
    Calculates Absolute Salinity from Practical Salinity.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
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
    _delta_SA, _SA_from_SP_Baltic

    Notes
    -----
    Since SP is non-negative by definition, this function changes
    any negative input values of SP to be zero.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SP = [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat = 188, 4
    >>> gsw.SA_from_SP(SP, p, lon, lat)
    array([ 34.71177971,  34.89152372,  35.02554774,  34.84723008,
            34.7366296 ,  34.73236186])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm
    for estimating Absolute Salinity in the global ocean. Submitted to Ocean
    Science. A preliminary version is available at Ocean Sci. Discuss.,
    6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    #SP[SP < 0] = 0
    SP = np.maximum(SP, 0)

    dSA = lib._delta_SA( p, lon, lat )

    SA = ( cte.SSO / 35 ) * SP + dSA
    SA_baltic = lib._SA_from_SP_Baltic( SP, lon, lat )

    if SA_baltic is not None:
        SA[~SA_baltic.mask] = SA_baltic[~SA_baltic.mask]

    return SA


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
        SP[indsbaltic] = SP_baltic[indsbaltic]

    return SP

@match_args_return
def Sstar_from_SP(SP, p, lon, lat):
    r"""
    Calculates Preformed Salinity from Absolute Salinity.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
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
    _delta_SA, _SA_from_SP_Baltic

    Notes
    -----
    The mask is only set when the observation is well and truly on dry
    land; often the warning flag is not set until one is several hundred
    kilometers inland from the coast.

    Since SP is non-negative by definition, this function changes any negative
    input values of SP to be zero.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SP = [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat =  188, 4
    >>> gsw.Sstar_from_SP(SP, p, lon, lat)
    array([ 34.7115532 ,  34.89116101,  35.02464926,  34.84359277,
            34.7290336 ,  34.71967638])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm
    for estimating Absolute Salinity in the global ocean. Submitted to Ocean
    Science. A preliminary version is available at Ocean Sci. Discuss.,
    6, 215-242.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    lon, lat, p, SP = np.broadcast_arrays(lon, lat, p, SP)

    SP[SP < 0] = 0

    dSA = lib._delta_SA( p, lon, lat )
    Sstar = (cte.SSO/35.) * SP - cte.r1 * dSA

    # In the Baltic Sea, SA == Sstar.
    Sstar_baltic = lib._SA_from_SP_Baltic( SP, lon, lat )

    if Sstar_baltic is not None:
        Sstar[~Sstar_baltic.mask] = Sstar_baltic[~Sstar_baltic.mask]

    return Sstar

def SA_Sstar_from_SP(SP, p, lon, lat):
    r"""
    Calculates Absolute Salinity and Preformed Salinity from Practical Salinity.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
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
    Sstar : masked array
            Preformed Salinity [g kg :sup:`-1`]

    See Also
    --------
    _delta_SA, _SA_from_SP_Baltic

    Notes
    -----
    In the Baltic Sea, Sstar == SA.

    The mask is only set when the observation is well and truly on dry
    land; often the warning flag is not set until one is several hundred
    kilometers inland from the coast.

    Since SP is non-negative by definition, this function changes any negative
    input values of SP to be zero.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SP = [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat =  188, 4
    >>> gsw.SA_Sstar_from_SP(SP, p, lon, lat)
    array([[ 34.71177971,  34.89152372,  35.02554774,  34.84723008,
             34.7366296 ,  34.73236186],
           [ 34.7115532 ,  34.89116101,  35.02464926,  34.84359277,
             34.7290336 ,  34.71967638]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm
    for estimating Absolute Salinity in the global ocean. Submitted to Ocean
    Science. A preliminary version is available at Ocean Sci. Discuss.,
    6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    return SA_from_SP(SP, p, lon, lat), Sstar_from_SP(SP, p, lon, lat)


@match_args_return
def SA_from_rho(rho, t, p):
    r"""
    Calculates the Absolute Salinity of a seawater sample, for given values of
    its density, in situ temperature and sea pressure (in dbar).

    One use for this function is in the laboratory where a measured value of
    the in situ density :math:`\rho` of a seawater sample may have been made at
    the laboratory temperature :math:`t` and at atmospheric pressure :math:`p`.
    The present function will return the Absolute Salinity SA of this seawater
    sample.

    Parameters
    ----------
    rho : array_like
          in situ density [kg m :sup:`-3`]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    This is expressed on the Reference-Composition Salinity Scale of
    Millero et al. (2008).

    After two iterations of a modified Newton-Raphson iteration, the error in SA
    is typically no larger than 2 :math:`^\times` 10 :sup:`-13` [g kg :sup:`-1`]

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> rho = [1021.839, 1022.262, 1024.426, 1027.792, 1029.839, 1032.002]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.SA_from_rho(rho, t, p)
    array([ 34.71022966,  34.89057683,  35.02332421,  34.84952096,
            34.73824809,  34.73188384])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.5.

    .. [2] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008:
    The composition of Standard Seawater and the definition of the
    Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.

    Modifications:
    2010-08-23. Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """
    n0, n1 = 0, 1
    v_lab = np.ones( rho.shape ) / rho
    v_0 = lib._gibbs(n0, n0, n1, 0, t, p)
    v_120 = lib._gibbs(n0, n0, n1, 120, t, p)
    SA = 120 * ( v_lab - v_0 ) / ( v_120 - v_0 ) # initial estimate of SA
    Ior = (SA < 0) | (SA > 120)
    v_SA = ( v_120 - v_0 ) / 120 # initial estimate of v_SA, SA derivative of v

    for iter in range(0,2):
        SA_old = SA
        delta_v = lib._gibbs(n0, n0, n1, SA_old, t, p) - v_lab
        SA = SA_old - delta_v / v_SA # half way through the modified N-R method
        SA_mean = 0.5 * ( SA + SA_old )
        v_SA = lib._gibbs(n1, n0, n1, SA_mean, t, p)
        SA = SA_old - delta_v / v_SA

    SA[Ior] = np.ma.masked

    return SA

if __name__=='__main__':
    import doctest
    doctest.testmod()
