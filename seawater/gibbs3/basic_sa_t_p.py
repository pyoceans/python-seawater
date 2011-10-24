# -*- coding: utf-8 -*-

"""
Basic thermodynamic properties in terms of (SA, t, p)

Functions:
----------
  rho_t_exact(SA, t, p)
      in-situ density
  enthalpy_t_exact(SA, t, p)
      enthalpy
  specvol_t_exact(SA, t, p)
      specific volume
  entropy_t_exact(SA, t, p)
      entropy
  cp_t_exact(SA, t, p)
      isobaric heat capacity
  Helmholtz_energy_t_exact(SA, t, p)
      Helmholtz energy
  sound_speed_t_exact(SA, t, p)
      sound speed
  pot_rho__t_exact (SA, t, p, pr=0)
      potential density
  specvol_anom_t_exact(SA, t, p)
      specific volume anomaly
  alpha_wrt_CT_t_exact (SA, t, p)
      thermal expansion coefficient with respect to Conservative Temperature.
  alpha_wrt_pt_t_exact (SA, t, p)
      thermal expansion coefficient with respect to potential temperature
  alpha_wrt_t_exact (SA, t, p)
      thermal expansion coefficient with respect to in-situ temperature
  beta_const_CT_t_exact (SA, t, p)
      saline contraction coefficient at constant Conservative Temperature.
  beta_const_pt_t_exact (SA, t, p)
      saline contraction coefficient at constant potential temperature
  beta_const_t_exact (SA, t, p)
      saline contraction coefficient at constant in-situ temperature
  internal_energy_t_exact (SA, t, p)
      internal energy
  isochoric_heat_cap_t_exact (SA, t, p)
      isochoric heat capacity
  chem_potential_relative_t_exact (SA, t, p)
      relative chemical potential
  chem_potential_water_t_exact (SA, t, p)
      chemical potential of water in seawater
  chem_potential_salt_t_exact (SA, t, p)
      chemical potential of salt in seawater
  kappa_t_exact (SA, t, p)
      isentropic compressibility
  kappa_const_t_exact (SA, t, p)
      isothermal compressibility
  adiabatic_lapse_rate_t_exact (SA, t, p)
      adiabatic lapse rate
  osmotic_coefficient_t_exact (SA, t, p)
      osmotic coefficient of seawater

This is part of the python Gibbs Sea Water library
http://code.google.com/p/python-seawater/
"""

from __future__ import division

import numpy as np

from library import gibbs
from utilities import match_args_return
from constants import Kelvin, db2Pascal, P0, SSO
from conversion import pt_from_CT, pt_from_t

__all__ = [
           'rho_t_exact',
           #'pot_rho_t_exact',
           'specvol_t_exact',
           'specvol_anom_t_exact',
           #'alpha_wrt_CT',
           #'alpha_wrt_pt',
           #'alpha_wrt_t',
           #'beta_const_CT',
           #'beta_const_pt',
           #'beta_const_t',
           'entropy_t_exact',
           #'internal_energy_t_exact',
           'enthalpy_t_exact',
           'cp_t_exact',
           #'isochoric_heat_cap_t_exact',
           #'chem_potential_relative_t_exact',
           #'chem_potential_water_t_exact',
           #'chem_potential_salt_t_exact',
           'Helmholtz_energy_t_exact',
           'sound_speed_t_exact',
           #'kappa_t_exact',
           #'kappa_const_t',
           #'adiabatic_lapse_rate_t_exact',
           #'molality_t_exact',
           #'ionic_strength_t_exact',
           #'osmotic_coefficient_t_exact',
           #'t_maxdensity_t_exact',
           #'pt_maxdensity_t_exact',
           #'CT_maxdensity_t_exact',
           #'temps_maxdensity_t_exact',
          ]


@match_args_return
def Helmholtz_energy_t_exact(SA, t, p):
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
    >>> gsw.Helmholtz_energy_t_exact(SA, t, p)
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
    g000 = gibbs(n0, n0, n0, SA, t, p)
    g001 = gibbs(n0, n0, n1, SA, t, p)
    return (g000 - (db2Pascal * p + P0) * g001)


@match_args_return
def rho_t_exact(SA, t, p):
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
    rho_t_exact : array_like
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
    rho = 1. / gibbs(n0, n0, n1, SA, t, p)

    return rho


@match_args_return
def enthalpy_t_exact(SA, t, p):
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
    g000 = gibbs(n0, n0, n0, SA, t, p)
    g010 = gibbs(n0, n1, n0, SA, t, p)
    return g000 - (t + Kelvin) * g010


@match_args_return
def specvol_t_exact(SA, t, p):
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
    return gibbs(n0, n0, n1, SA, t, p)


@match_args_return
def entropy_t_exact(SA, t, p):
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
    >>> gsw.entropy_t_exact(SA, t, p)
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
    return -gibbs(n0, n1, n0, SA, t, p)


@match_args_return
def cp_t_exact(SA, t, p):
    r"""
    Calculates the isobaric heat capacity of seawater.

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
    cp_t_exact : array_like
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
    >>> gsw.cp_t_exact(SA, t, p)
    array([ 4002.88800396,  4000.98028393,  3995.54646889,  3985.07676902,
            3973.59384348,  3960.18408479])

    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    Modifications:
    2011-03-29. David Jackett, Trevor McDougall and Paul Barker
    """

    n0, n2 = 0, 2

    return -(t + Kelvin) * gibbs(n0, n2, n0, SA, t, p)


@match_args_return
def sound_speed_t_exact(SA, t, p):
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
    >>> gsw.sound_speed_t_exact(SA, t, p)
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
    g020 = gibbs(n0, n2, n0, SA, t, p)
    g011 = gibbs(n0, n1, n1, SA, t, p)
    g001 = gibbs(n0, n0, n1, SA, t, p)
    g002 = gibbs(n0, n0, n2, SA, t, p)
    return g001 * np.sqrt(g020 / (g011 ** 2 - g020 * g002))


@match_args_return
def specvol_anom_t_exact(SA, t, p):
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
    specvol_anom_t_exact : array_like
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
    >>> gsw.specvol_anom_t_exact(SA, t, p)
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

    pt_zero = pt_from_CT(SSO, 0)

    t_zero = pt_from_t(SSO, pt_zero, 0, p)

    return  (gibbs(n0, n0, n1, SA, t, p) - gibbs(n0, n0, n1, SSO, t_zero, p))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
