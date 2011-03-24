# -*- coding: utf-8 -*-

"""
Basic thermodynamic properties in terms of (SA, t, p)

Functions:
----------
  rho(SA, t, p)
      in-situ density
  pot_rho(SA, t, p, pr=0)
      potential density
  specvol(SA, t, p)
      specific volume
  specvol_anom(SA, t, p)
      specific volume anomaly
  alpha_wrt_CT(SA, t, p)
      thermal expansion coefficient with respect to Conservative Temperature.
  alpha_wrt_pt(SA, t, p)
      thermal expansion coefficient with respect to potential temperature
  alpha_wrt_t(SA, t, p)
      thermal expansion coefficient with respect to in-situ temperature
  beta_const_CT(SA, t, p)
      saline contraction coefficient at constant Conservative Temperature.
  beta_const_pt(SA, t, p)
      saline contraction coefficient at constant potential temperature
  beta_const_t(SA, t, p)
      saline contraction coefficient at constant in-situ temperature
  entropy(SA, t, p)
      entropy
  internal_energy(SA, t, p)
      internal energy
  enthalpy(SA, t, p)
      enthalpy
  cp(SA, t, p)
      isobaric heat capacity
  isochoric_heat_cap(SA, t, p)
      isochoric heat capacity
  chem_potential_relative(SA, t, p)
      relative chemical potential
  chem_potential_water(SA, t, p)
      chemical potential of water in seawater
  chem_potential_salt(SA, t, p)
      chemical potential of salt in seawater
  Helmholtz_energy(SA, t, p)
      Helmholtz energy
  sound_speed(SA, t, p)
      sound speed
  kappa(SA, t, p)
      isentropic compressibility
  kappa_const_t(SA, t, p)
      isothermal compressibility
  adiabatic_lapse_rate(SA, t, p)
      adiabatic lapse rate
  molality(SA)
      molality of seawater
  ionic_strength(SA)
      ionic strength of seawater
  osmotic_coefficient(SA, t, p)
      osmotic coefficient of seawater
  t_maxdensity(SA, p)
      in_situ temperature of maximum density
  pt_maxdensity(SA, p)
      potential temperature of maximum density
  CT_maxdensity(SA, p)
      Conservative Temperature of maximum density
  temps_maxdensity(SA, p)
      temperatures of maximum density of seawater

This is part of the python Gibbs Sea Water library
http://code.google.com/p/python-seawater/

"""

from __future__ import division

import numpy as np
import seawater.constants as cte
from library import match_args_return
import library as lib
from conversions import *

# ----------------------------

__all__ =  ['rho',
            'pot_rho',
            'specvol',
            'specvol_anom',
            'alpha_wrt_CT',
            'alpha_wrt_pt',
            'alpha_wrt_t',
            'beta_const_CT',
            'beta_const_pt',
            'beta_const_t',
            'entropy',
            'internal_energy',
            'enthalpy',
            'cp',
            'isochoric_heat_cap',
            'chem_potential_relative',
            'chem_potential_water',
            'chem_potential_salt',
            'Helmholtz_energy',
            'sound_speed',
            'kappa',
            'kappa_const_t',
            'adiabatic_lapse_rate',
            'molality',
            'ionic_strength',
            'osmotic_coefficient',
            't_maxdensity',
            'pt_maxdensity',
            'CT_maxdensity',
            'temps_maxdensity']

# -----------------------------

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

# -----------------------

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

    pt = pt_from_t(SA, t, p, pr=pr)

    pot_rho = rho(SA, pt, pr)

    return pot_rho

# ---------------------------------

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

# ------------------------------

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

    t_zero = pt_from_t(SSO, pt_zero, 0, p)

    specvol_anom = ( lib._gibbs(n0, n0, n1, SA, t, p) -
                     lib._gibbs(n0, n0, n1, SSO, t_zero, p) )

    return specvol_anom

# --------------------------------------

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

    # pt0 = pt_from_t(SA, t, p) #NOTE: original call a faster version
    #"pt0_from_t" instead of pt_from_t
    pt0 = pt0_from_t(SA, t, p)

    fac = -cte.cp0 / ( (cte.Kelvin + pt0) * lib._gibbs(n0, n2, n0, SA, t, p ) )

    g011 = lib._gibbs(n0, n1, n1, SA, t, p)
    g001 = lib._gibbs(n0, n0, n1, SA, t, p)
    alpha_wrt_CT = fac * ( g011 / g001 )

    return alpha_wrt_CT

# --------------------------------------

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

    #pt0 = pt_from_t(SA, t, p) #NOTE: original call a faster version
    # (pt0_from_t) instead of pt_from_t
    pt0 = pt0_from_t(SA, t, p)

    fac = lib._gibbs(n0, n2, n0, SA, pt0, 0) / lib._gibbs(n0, n2, n0, SA, t, p)
    g011 = lib._gibbs(n0, n1, n1, SA, t, p)
    g001 = lib._gibbs(n0, n0, n1, SA, t, p)
    alpha_wrt_pt = fac * ( g011 / g001 )

    return alpha_wrt_pt

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

# ------------------------------

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

    # pt0 = pt_from_t(SA, t, p) #NOTE: original call a faster version
    # (pt0_from_t) instead of pt_from_t
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

# -----------------------------------

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

    # pt0 = pt_from_t(SA, t, p) #NOTE: original call a faster version
    # "pt0_from_t" instead of pt_from_t
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

# -------------------------------------------

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

# -----------------------------

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

# ----------------------

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

# -------------------------


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

# -------------------------------------------

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

# -------------------------------------

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

# -------------------------------

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

# --------------------------------

@match_args_return
def Helmholtz_energy(SA, t, p):
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
    >>> gsw.Helmholtz_energy(SA, t, p)
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

# -------------------------------------

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

# ---------------------------------------

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

# ---------------------------------------

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

# ----------------------------------

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
def t_maxdensity(SA, p):
    r"""
    in-situ temperature at maximum density

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
    >>> gsw.t_maxdensity(SA, p)
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
    return np.ma.array(t, mask=mask, copy=False)

def pt_maxdensity(SA, p):
    """potential temperature of maximum density

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    p : array_like
        pressure [dbar]

    Returns
    -------
    pt_maxdensity : array_like, potential temperature [deg C]

    """

    return pt0_from_t(SA, t_maxdensity(SA, p), p)

def CT_maxdensity(SA, p):
    """Conservative Temperature at maximum density

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    p : array_like
        pressure [dbar]

    Returns
    -------
    CT_maxdensity : array_like, Conservative Temperature [deg C]

    """

    return CT_from_pt(SA, pt_maxdensity(SA, p))

def temps_maxdensity(SA, p):
    """temperatures at maximum density

    See individual functions t_maxdensity, pt_maxdensity, CT_maxdensity
    Retained for compatibility with the Matlab toolbox.

    """
    
    t_maxdens  = t_maxdensity(SA, p)
    pt_maxdens = pt0_from_t(SA, t_maxden, p)
    CT_maxdens = CT_from_pt(SA, pt_maxden)
    return t_maxdens, pt_maxdens, CT_maxdens


# -------------------

if __name__=='__main__':
    import doctest
    doctest.testmod()






