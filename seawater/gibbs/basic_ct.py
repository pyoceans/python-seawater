# -*- coding: utf-8 -*-

"""
Basic thermodynamic properties in terms of CT and pt

Functions:
  rho_CT(SA, CT, p)
      in-situ density from CT
  alpha_CT(SA, CT, p)
      thermal expansion coefficient
  beta_CT(SA, CT, p)
      saline contraction coefficient
  rho_alpha_beta_CT(SA, CT, p)
      in-situ density, thermal expansion & saline contraction 
  specvol_CT(SA, CT, p)
      specific volume from CT
  specvol_anom_CT
      specific volume anomaly from CT
  sigma0_CT(SA, CT)
      sigma0 in terms of SA & CT with reference pressure of 0 dbar
  sigma1_CT(SA, CT)
      sigma1 in terms of SA & CT with reference pressure of 1000 dbar
  sigma2_CT(SA, CT)
      sigma2 in terms of SA & CT with reference pressure of 2000 dbar
  sigma3_CT(SA, CT)
      sigma3 in terms of SA & CT with reference pressure of 3000 dbar
  sigma4_CT(SA, CT)
      sigma4 in terms of SA & CT with reference pressure of 4000 dbar
  enthalpy_CT(SA, CT, p)
      enthalpy from CT
  enthalpy_diff_CT
      difference of enthalpy from CT between two pressures
  entropy_from_pt(SA, pt)
      entropy from potential temperature
  entropy_from_CT(SA, CT)
      entropy from Conservative Temperature
  pt_from_entropy(SA, entropy)
      potential temperature from entropy
  CT_from_entropy(SA, entropy)
      Conservative Temperature from entropy

This is part of the python Gibbs Sea Water library
http://code.google.com/p/python-seawater/

"""

from __future__ import division

import numpy as np
import seawater.constants as cte
from library import match_args_return
import library as lib
from conversions import *
from basic_sa_t_p import *

# ----------------------------------

__all__ = ['rho_CT',
           'alpha_CT',
           'beta_CT',
           'rho_alpha_beta_CT',
           'specvol_CT',
           #'specvol_anom_CT',
           'sigma0_CT',
           'sigma1_CT',
           'sigma2_CT',
           'sigma3_CT',
           'sigma4_CT',
           'enthalpy_CT',
           #'enthalpy_diff_CT',
           'entropy_from_pt',
           'entropy_from_CT',
           'pt_from_entropy',
           'CT_from_entropy']

# ------------------------------------

@match_args_return
def rho_CT(SA, CT, p):
    """
    density from Absolute Salinity, Conservative Temperature and pressure

    Parameters
    ----------
    SA : array_like
          Absolute Salinity  [g/kg]
    CT : array_like
          Conservative Temperature [deg C]
    p : array_like
          sea pressure  [dbar]

    Returns
    -------
    rho_CT:  array_like
         in-situ density  [kg/m**3]

    """
    
    pt0 = pt_from_CT(SA, CT)
    t = pt_from_t(SA, pt0, 0.0, p)
    return rho(SA, t, p)

# ------------------------

@match_args_return
def alpha_CT(SA, CT, p):
    r"""
    Thermal expansion coefficient of seawater with respect to
    Conservative Temperature as function of Absolute Salinty,
    Conservative Temperature and pressure.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
        in Conservative Temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    alpha_CT : array_like
                   thermal expansion coefficient [K :sup:`-1`]

    """
    pt0 = pt_from_CT(SA, CT)
    t = pt_from_t(SA, pt0, 0, p)
    return alpha_wrt_CT(SA, t, p)

# -------------------------------

@match_args_return
def beta_CT(SA, CT, p):
    r"""
    Saline contraction coefficient of seawater at constant
    Conservative Temperature as function of Absolute Salinty,
    Conservative Temperature and pressure.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
        in Conservative Temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    beta_CT : array_like
                   thermal expansion coefficient [K :sup:`-1`]

    """
    pt0 = pt_from_CT(SA, CT)
    t = pt_from_t(SA, pt0, 0, p)
    return beta_const_CT(SA, t, p)

# -----------------------------------
    
def rho_alpha_beta_CT(SA, CT, p):
    """
    Density, thermal expansion and saline contraction
    as functions of Abs. Salinity, Cons. temp. and pressure.

    See the individual functions rho_CT, alpha_CT, and beta_CT.
    Retained for compatibility with the Matlab GSW toolbox.
    
    """
    
    return (rho_CT(SA, CT, p), 
            alpha_CT(SA, CT, p), 
            beta_CT(SA, CT, p))

# -------------------------

def specvol_CT(SA, CT, p):
    """
    specific volume from Absolute Salinity, Conservative Temperature
    and pressure

    Parameters
    ----------
    SA : array_like
          Absolute Salinity  [g/kg]
    CT : array_like
          Conservative Temperature [deg C]
    p : array_like
          sea pressure  [dbar]

    Returns
    -------
    specvol_CT:  array_like
         in-situ specific volume  [m**3/kg]

    """
    
    pt = pt_from_CT(SA, CT)
    t = pt_from_t(SA, pt, 0.0, p)
    return specvol(SA, t, p)
    #return 1.0 / rho_CT(SA, CT, p)

# ---------------------------------

    # NOT IMPLEMENTED YET
    #def specvol_anom_CT

# --------------------------------

def sigma0_CT(SA, CT):
    """
    potential density anomaly with zero reference pressure as
    function of Absolute Salinity and Conservative Temperature

    Parameters
    ----------
    SA : array_like
          Absolute Salinity  [g/kg]
    CT : array_like
          Conservative Temperature [deg C]

    Returns
    -------
    sigma0_CT:  array_like
         Potential density anomaly  [kg/m**3]

    The anomaly is taken with respect to 1000.0 kg/m**3

    """
    
    return rho_CT(SA, CT, 0) - 1000

# ----------------------

def sigma1_CT(SA, CT):
    return rho_CT(SA, CT, 1000) - 1000
def sigma2_CT(SA, CT):
    return rho_CT(SA, CT, 2000) - 1000
def sigma3_CT(SA, CT):
    return rho_CT(SA, CT, 3000) - 1000
def sigma4_CT(SA, CT):
    return rho_CT(SA, CT, 4000) - 1000

# -------------------------

def enthalpy_CT(SA, CT, p):
    """
    specific enthalpy from Absolute Salinity, Conservative Temperature
    and pressure

    Parameters
    ----------
    SA : array_like
          Absolute Salinity  [g/kg]
    CT : array_like
          Conservative Temperature [deg C]
    p : array_like
          sea pressure  [dbar]

    Returns
    -------
    enthalpy_CT:  array_like
         specific enthalpy  [J/kg]

    """

    pt = pt_from_CT(SA, CT)
    t  = pt_from_t(SA, pt, 0, p)
    return enthalpy(SA, t, p)

# -------------------------

    # NOT IMPLEMENTED YET
    #def enthalpy_diff_CT

# -------------------------

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

# -------------------------------------

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

# --------------------------------
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

