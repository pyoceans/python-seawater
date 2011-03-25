# -*- coding: utf-8 -*-

"""
gibbs25

Functions using the faster and slightly less accurate
25 term equation of state in terms of Absolute Salinity
and Conservative Temperature

Reference:
----------
McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
   Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term 
   expression for the density of seawater in terms of Conservative 
   Temperature, and related properties of seawater.
   To be submitted to Ocean Science Discussions. 

"""

import numpy as np
from library import match_args_return
from density25 import *

def in_funnel(SA, CT, p):
    u"""oceanographic funnel check for the 25-term equation
 
    Parameters
    ----------
    SA : array_like
         Absolute Salinity            [g/kg]
    CT : array_like
         Conservative Temperature     [°C]
    p  : array_like
         sea pressure                 [dbar]
           (ie. absolute pressure - 10.1325 dbar)
          
    Returns
    -------
    in_funnel : boolean ndarray or scalar
        True,  if SA, CT and p are inside the "funnel" 
        False, if SA, CT and p are outside the "funnel",
               or one of the values are NaN or masked

    Note. The term "funnel" describes the range of SA, CT and p over which 
    the error in the fit of the computationally-efficient 25-term 
    expression for density in terms of SA, CT and p was calculated
    (McDougall et al., 2010).

    author: 
    Trevor McDougall and Paul Barker    [ help_gsw@csiro.au ]
    2011-02-27: Bjørn Ådlandsvik, python version

"""

    # Check variables and resize if necessary
    scalar = np.isscalar(SA) and np.isscalar(CT) and np.isscalar(p)
    SA, CT, p = np.broadcast_arrays(SA, CT, p)

    input_nan = np.isnan(SA) | np.isnan(CT) | np.isnan(p)

    infunnel = ((p <= 8000)  &
                (SA >= 0)    &
                (SA <= 42.2) &
                (CT >= (-0.3595467 - 0.0553734*SA)) &
                ((p >= 5500) | (SA >= 0.006028*(p - 500))) &
                ((p >= 5500) | (CT <= (33.0 - 0.003818181818182*p))) &
                ((p <= 5500) | (SA >= 30.14)) &
                ((p <= 5500) | (CT <= 12.0)))

    infunnel = infunnel & np.logical_not(input_nan)

    if scalar:
        infunnel = bool(infunnel)

    return infunnel



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

    #if p==0: p = p+dp
    #inds = find(p>=dp); 
    #p_u = zeros(size(p));
    #p_l = dp*ones(size(p));
    #if ~isempty(inds):
        #p_u(inds) = p(inds)-dp;
        #p_l(inds) = p(inds)+dp;
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



def isopycnal_slope_ratio_CT25(SA, CT, p, pr=0):

    alpha_CT    = alpha_CT25(SA, CT, p)
    beta_CT     = beta_CT25(SA, CT, p)
    alpha_CT_pr = alpha_CT25(SA, CT, pr)
    beta_CT_pr  = beta_CT25(SA, CT, pr)

    out = alpha_CT * beta_CT_pr / (alpha_CT_pr * beta_CT)
    return out


