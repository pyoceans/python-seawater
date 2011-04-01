# -*- coding: utf-8 -*-

"""
Neutral and non-linear properties, based on the 25-term expression

Functions
---------
  cabbeling_CT25(SA, CT, p)                  
      cabbeling coefficient
  thermobaric_CT25(SA, CT, p)             
      thermobaric coefficient
  isopycnal_slope_ratio_CT25(SA, CT, p, pr)
      ratio of the slopes of isopycnals on the SA-CT diagram for p & pr
  isopycnal_vs_ntp_CT_ratio_CT25(SA, CT, p, pr)  
      ratio of the gradient of Conservative Temperature in a
      potential density surface to that in the neutral tangent plane
  ntp_pt_vs_CT_ratio_CT25(SA, CT, p)
     ratio of gradients of potential temperature &
     Conservative Temperature in a neutral tangent plane

"""

import numpy as np
from library import match_args_return
from density25 import *
from derivatives import pt_derivative_SA, pt_derivative_CT

# --------

__all__ = ['cabbeling_CT25',
           'thermobaric_CT25',  
           'isopycnal_slope_ratio_CT25',
           'isopycnal_vs_ntp_CT_ratio_CT25',
           'ntp_pt_vs_CT_ratio_CT25'
          ]

# --------

@match_args_return
def cabbeling_CT25(SA, CT, p):
    """Cabbeling coefficient wrt Conservative Temperature by 25 term eq.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g / kg]
    t : array_like
        in situ temperature [deg C (ITS-90)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    cabbeling_CT25 : array_like
        cabbeling coefficient  [K^-2]

    """

    dCT = 1e-3    # increment in Conservative Temperature is 1e-3 deg C.
    CT_l = CT - dCT;
    CT_u = CT + dCT;

    alpha_CT   = alpha_CT25(SA, CT, p)
    alpha_CT_u = alpha_CT25(SA, CT_u, p)
    alpha_CT_l = alpha_CT25(SA, CT_l, p)
    beta_CT    = beta_CT25(SA, CT, p)
    #beta_CT_u  = beta_CT25(SA, CT_u, p)
    #beta_CT_l  = beta_CT25(SA, CT_l, p)
    
    ratio_CT = alpha_CT / beta_CT
    alpha_CT_CT = (alpha_CT_u - alpha_CT_l) / (CT_u-CT_l)

    dSA = 1e-3    # increment in Absolute Salinity is 1e-3 g/kg.
    SA_l = np.clip(SA - dSA, 0.0, np.inf)
    SA_u = SA + dSA

    alpha_CT_u = alpha_CT25(SA_u, CT, p)
    alpha_CT_l = alpha_CT25(SA_l, CT, p)
    beta_CT_u  = beta_CT25(SA_u, CT, p)
    beta_CT_l  = beta_CT25(SA_l, CT, p)

    alpha_CT_SA = (alpha_CT_u - alpha_CT_l) / (SA_u-SA_l)
    beta_CT_SA  = (beta_CT_u - beta_CT_l)   / (SA_u-SA_l)

    return alpha_CT_CT + ratio_CT*(2*alpha_CT_SA - ratio_CT*beta_CT_SA)

# ----------------------------------------------

@match_args_return
def thermobaric_CT25(SA, CT, p):

    """Thermobaric coefficient 25 term equation

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g / kg]
    t : array_like
        in situ temperature [deg C (ITS-90)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    thermobaric_CT25 : array_like
        thermobaric coefficient  [1/(K Pa)]

    """

    db2Pa = 10000.0
    dp = 0.1         # pressure increment is 1e-1 dbar (10 Pa)

    p_l = np.clip(p - dp, 0.0, np.inf)
    p_u = p + dp
    
    alpha_CT   = alpha_CT25(SA, CT, p)
    alpha_CT_u = alpha_CT25(SA, CT, p_u)
    alpha_CT_l = alpha_CT25(SA, CT, p_l)
    beta_CT    = beta_CT25(SA, CT, p)
    beta_CT_u  = beta_CT25(SA, CT, p_u)
    beta_CT_l  = beta_CT25(SA, CT, p_l)

    alpha_CT_p = (alpha_CT_u - alpha_CT_l) / (p_u-p_l)
    beta_CT_p  = (beta_CT_u - beta_CT_l) / (p_u-p_l)

    thermobaric_CT = alpha_CT_p - (alpha_CT / beta_CT) * beta_CT_p
    return thermobaric_CT / db2Pa   # To have units of 1/(K Pa)

# ------------------------------------

def isopycnal_slope_ratio_CT25(SA, CT, p, pr=0):

    alpha_CT    = alpha_CT25(SA, CT, p)
    beta_CT     = beta_CT25(SA, CT, p)
    alpha_CT_pr = alpha_CT25(SA, CT, pr)
    beta_CT_pr  = beta_CT25(SA, CT, pr)

    out = alpha_CT * beta_CT_pr / (alpha_CT_pr * beta_CT)
    return out

# ------------------------------------

def isopycnal_vs_ntp_CT_ratio_CT25(SA, CT, p, pr=0):
    """isopycnal Conservative Temperature gradient ratio

    Ratio of the gradient of Conservative Temperature in a potential
    density surface to that in a neutral tangent plane (i.e. in a locally
    referenced potential density surface) (using 25-term equation)

    parameters:
    -----------
    SA : array-like, Absolute Salinity           [g/kg]
    CT : array-like, Conservative Temperature    [deg C]
    p  : array-like, sea pressure                [dbar]
    pr : scalar, optional,
            reference pressure, default = 0      [dbar]

    """
  
    if not np.isscalar(pr):
        raise ArgumentError, "The reference pressures should be scalar"

    SA, CT, p = np.broadcast_arrays(SA, CT, p)

    p_mid  = 0.5*(p[:-1, ...]  + p[1:, ...])
    SA_mid = 0.5*(SA[:-1, ...] + SA[1:, ...])
    CT_mid = 0.5*(CT[:-1, ...] + CT[1:, ...])


    dSA = SA[:-1,...] - SA[1:,...]
    dCT = CT[:-1,...] - CT[1:,...]

    alpha = alpha_CT25(SA_mid, CT_mid, p_mid)
    beta  = beta_CT25(SA_mid, CT_mid, p_mid)
    alpha_pr = alpha_CT25(SA_mid, CT_mid, pr)
    beta_pr  = beta_CT25(SA_mid, CT_mid, pr)

    anum   = dCT * alpha / beta - dSA
    adenom = dCT * alpha_pr / beta_pr - dSA

    return anum / adenom

# ---------------------------------

def ntp_pt_vs_CT_ratio_CT25(SA, CT, p):
    """ratio of gradients of potential temperature and
    Conservative Temperature in a neutral tangent  plane
    (in a locally-referenced potential density surface)(25-term equation)

    parameters:
    -----------
    SA : array-like, Absolute Salinity           [g/kg]
    CT : array-like, Conservative Temperature    [deg C]
    p  : array-like, sea pressure                [dbar]

    """    

    alpha_CT = alpha_CT25(SA, CT, p)
    beta_CT  = beta_CT25(SA, CT, p)

    pt_SA = pt_derivative_SA(SA, CT)
    pt_CT = pt_derivative_CT(SA, CT)

    return pt_CT + pt_SA * alpha_CT / beta_CT


