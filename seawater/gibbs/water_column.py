# -*- coding: utf-8 -*-

"""
water column properties, based on the 25-term expression for density

Functions
---------
  Nsquared_CT25(SA, CT, p, lat=None)
      buoyancy (Brunt-Vaisala) frequency squared (N^2)
  #Turner_Rsubrho_CT25(SA, CT, p)
  #    Turner angle & Rsubrho
  IPV_vs_fNsquared_ratio_CT25(SA, CT, p, pr)
      ratio of the vertical gradient of potential density to
      the vertical gradient of locally-referenced potential density

  NOTE: Present version works with the tests provided,
        but should be consistent with corner cases

"""

from __future__ import division

import numpy as np
#import seawater.constants as cte
from earth import grav
from density25 import *

# ---------------------

__all__ = ['Nsquared_CT25',
           'Turner_CT25',
           'Rsubrho_CT25',
           #'Turner_Rsubrho_CT25',
           'IPV_vs_fNsquared_ratio_CT25',
          ]

# -----------------------

def Nsquared_CT25(SA, CT, p, lat=None):
    u"""Brunt-Väisälä frequency squared

    parameters
    ----------
      SA  : array-like, Absolute Salinity         [g/kg]
      CT  : array-like, Conservative Temperature [deg C]
      p   : array-like, pressure                  [dbar]
      lat : optional, array-like                  [deg N]

    returns
      N2 : square of Brunt-Väisälä frequency

    TODO: Describe precisely how the dimensions

    """

    db2Pa = 1.0e4

    if lat != None:
        g = grav(lat, p)
    else:
        g = 9.7963 # Standard value from Griffies, 2004

    SA, CT, p, g = np.broadcast_arrays(SA, CT, p, g)

    p_mid = 0.5*(p[1:, ...] + p[:-1, ...])
    drho  = ( rho_CT25(SA[1:, ...], CT[1:, ...], p_mid) -
              rho_CT25(SA[:-1, ...], CT[:-1, ...], p_mid) )
    grav_local = 0.5*(g[1:, ...] + g[:-1, ...])
    dp = p[1:, ...] - p[:-1, ...]

    return grav_local * grav_local * drho / (db2Pa * dp)

# ------------------

def Turner_CT25(SA, CT, p):
    """
    Turner angle and Rsubrho

    DESCRIPTION:
    Calculates the Turner angle and the Rsubrho as a function of pressure 
    down a vertical water column.  These quantities express the relative 
    contributions of the vertical gradients of Conservative Temperature 
    and Absolute Salinity to the vertical stability (the square of the 
    Brunt-Vaisala Frequency squared, N^2).  Tu and Rsubrho are evaluated at 
    the mid pressure between the individual data points in the vertical.  
    This function uses computationally-efficient 25-term expression for 
    density in terms of SA, CT and p (McDougall et al., 2010). 

    parameters
    ----------
    SA : array_like, Absolute Salinity             [g/kg]
    CT : array_like, Conservative Temperature      [deg C]
    p  : array_like, sea pressure                  [dbar]

    returns
    -------
    Tu       :    Turner angle, on the same (M-1)xN grid as p_mid.
                 Turner angle has units of:        [ degrees of rotation ]
    Rsubrho  :    Stability Ratio, on the same (M-1)xN grid as p_mid.
                 Rsubrho is dimensionless.                    [ unitless ]

    """

    SA, CT, p = np.broadcast_arrays(SA, CT, p)

    # Ensure that SA is non-negative.
    SA.clip(0, np.inf)

    p_mid =  0.5 * (p[:-1,...] + p[1:,...])
    SA_mid = 0.5 * (SA[:-1,...] + SA[1:,...])
    CT_mid = 0.5 * (CT[:-1,...] + CT[1:,...])

    dSA = SA[:-1,...] - SA[1:,...]
    dCT = CT[:-1,...] - CT[1:,...]

    alpha = alpha_CT25(SA_mid, CT_mid, p_mid)
    beta  = beta_CT25(SA_mid, CT_mid, p_mid)

    Tu = np.arctan2(alpha*dCT + beta*dSA, alpha*dCT - beta*dSA)
    Tu = Tu * 180 / np.pi

    Rsubrho = np.nan + np.ones_like(dSA)
    I = (dSA != 0)
    Rsubrho[I] = (alpha[I]*dCT[I]) / (beta[I]*dSA[I])

    return Tu

# ------------------------------

def Rsubrho_CT25(SA, CT, p):
    """
    Turner angle and Rsubrho

    DESCRIPTION:
    Calculates the Turner angle and the Rsubrho as a function of pressure 
    down a vertical water column.  These quantities express the relative 
    contributions of the vertical gradients of Conservative Temperature 
    and Absolute Salinity to the vertical stability (the square of the 
    Brunt-Vaisala Frequency squared, N^2).  Tu and Rsubrho are evaluated at 
    the mid pressure between the individual data points in the vertical.  
    This function uses computationally-efficient 25-term expression for 
    density in terms of SA, CT and p (McDougall et al., 2010). 

    parameters
    ----------
    SA : array_like, Absolute Salinity             [g/kg]
    CT : array_like, Conservative Temperature      [deg C]
    p  : array_like, sea pressure                  [dbar]

    returns
    -------
    Tu       :    Turner angle, on the same (M-1)xN grid as p_mid.
                 Turner angle has units of:        [ degrees of rotation ]
    Rsubrho  :    Stability Ratio, on the same (M-1)xN grid as p_mid.
                 Rsubrho is dimensionless.                    [ unitless ]

    """

    SA, CT, p = np.broadcast_arrays(SA, CT, p)

    # Ensure that SA is non-negative.
    SA.clip(0, np.inf)

    p_mid =  0.5 * (p[:-1,...] + p[1:,...])
    SA_mid = 0.5 * (SA[:-1,...] + SA[1:,...])
    CT_mid = 0.5 * (CT[:-1,...] + CT[1:,...])

    dSA = SA[:-1,...] - SA[1:,...]
    dCT = CT[:-1,...] - CT[1:,...]

    alpha = alpha_CT25(SA_mid, CT_mid, p_mid)
    beta  = beta_CT25(SA_mid, CT_mid, p_mid)

    Rsubrho = np.nan + np.ones_like(dSA)
    I = (dSA != 0)
    Rsubrho[I] = (alpha[I]*dCT[I]) / (beta[I]*dSA[I])

    return Rsubrho

# --------------------------------

def IPV_vs_fNsquared_ratio_CT25(SA, CT, p, pr):
    """ratio of gradients of pot. density to local pot. density

   Calculates the ratio of the vertical gradient of potential density to 
   the vertical gradient of locally-referenced potential density.  This 
   ratio is also the ratio of the planetary Isopycnal Potential Vorticity
   (IPV) to f times N^2, hence the name for this variable,
   IPV_vs_fNsquared_ratio_CT25 (see Eqn. (3.20.5) of IOC et al. (2010)). 
   The reference sea pressure of the potential density surface must have a 
   constant value.

   IPV_vs_fNsquared_ratio_CT25 is evaluated at the mid pressure between the 
   individual data points in the vertical.  This function uses the 
   computationally-efficient 25-term expression for density in terms of 
   SA, CT and p (McDougall et al., 2010). 

   parameters
   ----------

    SA : array_like, Absolute Salinity             [g/kg]
    CT : array_like, Conservative Temperature      [deg C]
    p  : array_like, sea pressure                  [dbar]
    pr : array_like, reference pressure            [dbar]
         constant in depth
         
    returns
    -------
    ratio : array-like   [dimensionless]
       ratio of gradients of pot. density to local pot. density  

    """

    SA, CT, p = np.broadcast_arrays(SA, CT, p)

    p_mid =  0.5 * (p[:-1, ...]  + p[1:, ...])
    SA_mid = 0.5 * (SA[:-1, ...] + SA[1:, ...])
    CT_mid = 0.5 * (CT[:-1, ...] + CT[1:, ...])

    dSA = SA[:-1, ...] - SA[1:, ...]
    dCT = CT[:-1, ...] - CT[1:, ...]

    alpha = alpha_CT25(SA_mid, CT_mid, p_mid)
    beta  = beta_CT25(SA_mid, CT_mid, p_mid)
    alpha_pr = alpha_CT25(SA_mid, CT_mid, pr)
    beta_pr  = beta_CT25(SA_mid, CT_mid, pr)
    
    anum   = dCT * alpha_pr - dSA * beta_pr
    adenom = dCT * alpha    - dSA * beta

    I = (adenom != 0)
    ratio = np.nan + np.zeros_like(SA_mid)
    ratio[I] = anum[I] / adenom[I]

    return ratio


