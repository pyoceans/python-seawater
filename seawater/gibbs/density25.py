# -*- coding: utf-8 -*-

"""
Density and enthalpy, based on the 25-term expression for density

Functions:
----------

  rho_CT25
      in-situ density
  alpha_CT25
      thermal expansion coefficient
  beta_CT25
      saline contraction coefficient
  rho_alpha_beta_CT25
      in-situ density, thermal expansion & saline contraction coefficient
  specvol_CT25
      specific volume
  specvol_anom_CT25
      specific volume anomaly
   enthalpy_CT25                   
      enthalpy
   enthalpy_diff_CT25
      difference of enthalpy between two pressures

This is part of the python Gibbs Sea Water library
http://code.google.com/p/python-seawater/

"""


from __future__ import division

import numpy as np
import seawater.constants as cte
#import library as lib
from library import match_args_return
from basic_ct import *

# -----------------------

__all__ = ['rho_CT25',
           'alpha_CT25',
           'beta_CT25',
           'rho_alpha_beta_CT25',
           'specvol_CT25',
           'specvol_anom_CT25',
           'enthalpy_CT25',              
           'enthalpy_diff_CT25']

# -----------------------------

def _anum(SA, CT, p):
    """Numerator in the 25 term expression for rho"""

    CT2 = CT*CT

    anum =    ( 9.9984380290708214e+002 +
             CT * ( 7.1188090678940910e+000 +
             CT * (-1.9459922513379687e-002 +
             CT * 6.1748404455874641e-004)) +
             SA * ( 2.8925731541277653e+000 +
             CT * 2.1471495493268324e-003 +
             SA * 1.9457531751183059e-003) +
              p * ( 1.1930681818531748e-002 +
            CT2 * 2.6969148011830758e-007 +
             SA * 5.9355685925035653e-006 +
              p * (-2.5943389807429039e-008 +
            CT2 * -7.2734111712822707e-012)) )          

    return anum

def _adenom(SA, CT, p):
    """Denominator of the 25 term expression for rho"""

    CT2 = CT*CT

    adenom =    ( 1.0 +
             CT * ( 7.0547681896071576e-003 +
             CT * (-1.1753695605858647e-005 +
             CT * ( 5.9219809488274903e-007 +
             CT *   3.4887902228012519e-010))) +
             SA * ( 2.0777716085618458e-003 +
             CT * (-2.2210857293722998e-008 +
            CT2 *  -3.6628141067895282e-010) +
    np.sqrt(SA) * ( 3.4688210757917340e-006 +
            CT2 *   8.0190541528070655e-010)) +
              p * ( 6.8314629554123324e-006 +
    (p*CT)*(CT2 * -8.5294794834485446e-017 +
              p * -9.2275325145038070e-018))
                )
    return adenom


@match_args_return
def rho_CT25(SA, CT, p):
    """
    Calculates in situ density of seawater from Absolute Salinity and in situ
    temperature using the computationally-efficient 25-term expression 
    for density in terms of SA, CT and p

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
    rho_CT25 : array_like
          in situ density [kg m^-3]


    """
    return _anum(SA, CT, p) / _adenom(SA, CT, p)



@match_args_return
def alpha_CT25(SA, CT, p):

    CT2 = CT * CT
    sqrtSA = np.sqrt(SA)
    pCT = p * CT


    anum_CT =   (   7.118809067894091e+00 +
             CT * (-3.8919845026759374e-02 +
             CT *   1.8524521336762394e-03) +
             SA *   2.1471495493268324e-03 +
            pCT * ( 5.393829602366152e-07 -
              p *   1.454682234256454e-11)
                )

    adenom_CT =     (   7.0547681896071576e-03 +
                 CT * (-2.35073912117172980e-05 +
                 CT * ( 1.7765942846482467e-06 +
                 CT *   1.3955160891205007e-09)) +
                 SA * (-2.2210857293722998e-08 -
                CT2 *   1.09884423203685860e-09 +
          CT*sqrtSA *   1.6038108305614131e-09) -
           p*p*(CT2 *   2.5588438450345636e-16 +
                  p *   9.227532514503807e-18)
                    )

    return (adenom_CT - anum_CT / rho_CT25(SA,CT,p)) / _adenom(SA,CT,p)
            
@match_args_return
def beta_CT25(SA, CT, p):
    
    CT2 = CT * CT
    sqrtSA = np.sqrt(SA)

    anum_SA =   ( 2.8925731541277653 +
             CT * 2.1471495493268324e-03 +
             SA * 3.8915063502366117e-03 +
              p * 5.935568592503565e-06
                )
    
    adenom_SA = ( 2.077771608561846e-03 +
             CT * (-2.2210857293722998e-08 -
            CT2 *   3.6628141067895287e-10) +
         sqrtSA *  (5.203231613687601e-06 +
            CT2 *   1.2028581229210597e-09)
                )
          
    return (anum_SA/rho_CT25(SA,CT,p) - adenom_SA)/_adenom(SA,CT,p)

def rho_alpha_beta_CT25(SA, CT, p):
    return ( rho_CT25(SA, CT, p),
             alpha_CT25(SA, CT, p),
             beta_CT25(SA, CT, p) )


@match_args_return
def specvol_CT25(SA, CT, p):
    return _adenom(SA, CT, p) / _anum(SA, CT, p)


#### Move to library
@match_args_return
def _specvol_SSO_0_CT25(p):
    
    SSO = 35.16504

    return           ( (1.0 + 
                 SSO * (2.0777716085618458e-003 +
        np.sqrt(SSO) * 3.4688210757917340e-006) +
                   p * 6.8314629554123324e-006) /
                          (9.9984380290708214e+002 +
                 SSO * (2.8925731541277653e+000 +
                 SSO * 1.9457531751183059e-003) +
                   p * ( 1.1930681818531748e-002 +
                 SSO * 5.9355685925035653e-006 +
                   p * -2.5943389807429039e-008))
                     )

@match_args_return
def specvol_anom_CT25(SA, CT, p):
    return _adenom(SA, CT, p) / _anum(SA, CT, p) - _specvol_SSO_0_CT25(p)



@match_args_return
def enthalpy_CT25(SA, CT, p):
    """
   Calculates specific enthalpy of seawater using the computationally-
   efficient 25-term expression for density in terms of SA, CT and p

   Parameters:
   -----------
   SA : array_like
        Absolute Salinity                         [g/kg]
   CT : array_like
        Conservative Temperature                  [deg C]
   p  : arra_like
        sea pressure                              [dbar]

   Returns
   -------
   enthalpy : array_like
              specific enthalpy                   [J/kg]

   """
    
    db2Pa = 1e4 
    cp0 = 3991.86795711963         

    CT2 = CT*CT 
    CT3 = CT*CT2

    a0 =      (  1 +
           CT * (7.0547681896071576e-3 +
           CT * (-1.1753695605858647e-5 + 
           CT * (5.9219809488274903e-7 + 
           CT * 3.4887902228012519e-10))) + 
           SA * ( 2.0777716085618458e-3 + 
           CT * (-2.2210857293722998e-8 + 
          CT2 * -3.6628141067895282e-10) + 
  np.sqrt(SA) * (3.4688210757917340e-6 + 
          CT2 * 8.0190541528070655e-10))
              )   

    a1 = 6.8314629554123324e-6
    a2 = CT3 * -8.5294794834485446e-17
    a3 = CT * -9.2275325145038070e-18

    b0 =     (   9.9984380290708214e2 + 
          CT * (7.1188090678940910e0 + 
          CT * (-1.9459922513379687e-2 + 
          CT * 6.1748404455874641e-4)) + 
          SA * (2.8925731541277653e0 + 
          CT * 2.1471495493268324e-3 + 
          SA * 1.9457531751183059e-3)
             )
    
    b1 =     (  0.5*(1.1930681818531748e-2 + 
         CT2 * 2.6969148011830758e-7 + 
          SA * 5.9355685925035653e-6)
             )
    
    b2 = CT2*-7.2734111712822707e-12 - 2.5943389807429039e-8

    b1sq = b1*b1 
    sqrt_disc = np.sqrt(b1sq - b0*b2)

    N = a0 + (2*a3*b0*b1/b2 - a2*b0)/b2
    
    M = a1 + (4*a3*b1sq/b2 - (a3*b0 + 2*a2*b1))/b2

    A = b1 - sqrt_disc
    B = b1 + sqrt_disc

    part = (N*b2 - M*b1)/(b2*(B - A))


    return ( cp0*CT + 
               db2Pa*(p*(a2 - 2*a3*b1/b2 + 0.5*a3*p)/b2 + 
               (M/(2*b2))*np.log(1 + p*(2*b1 + b2*p)/b0) + 
               part*np.log(1 + (b2*p*(B - A))/(A*(B + b2*p))))
           )

# --------------

@match_args_return
def enthalpy_diff_CT25(SA, CT, p_shallow, p_deep):
    """
    enthalpy_CT25(SA, CT, p_deep) - enthalpy_CT25(SA, CT, p_shallow)

    Hopefully with higher accuracy and faster

    Parameters
    ----------
    SA        : array_like, Absolute Salinity        [g/kg]
    CT        : array_like, Conservative Temperature [deg C]
    p_shallow : array_like, upper sea pressure       [dbar]
    p_deep    : array_like, lower sea temperature    [dbar]

    Returns
    -------
    enthalpy_diff_CT25 : array_like
           difference of specific enthalpy          [ J/kg ]

    """

    SA = SA.clip(0, np.inf)  # Make SA non-negative

    db2Pa = 1e4
    
    CT2 = CT * CT
    CT3 = CT * CT2

    a0 = ( 1 + CT * ( 7.0547681896071576e-3 +
               CT * (-1.1753695605858647e-5 + 
               CT * (5.9219809488274903e-7 + 
               CT * 3.4887902228012519e-10))) + 
               SA * ( 2.0777716085618458e-3 + 
               CT * (-2.2210857293722998e-8 + 
              CT2 * -3.6628141067895282e-10) + 
      np.sqrt(SA) * (3.4688210757917340e-6 + 
              CT2 * 8.0190541528070655e-10)) )
    
    a1 =  6.8314629554123324e-6
    a2 =  CT3 * -8.5294794834485446e-17
    a3 =  CT * -9.2275325145038070e-18

    b0 =   ( 9.9984380290708214e2 + 
        CT * (7.1188090678940910 + 
        CT *   (-1.9459922513379687e-2 + 
        CT *     6.1748404455874641e-4)) + 
        SA * (2.8925731541277653e0 + 
        CT *   2.1471495493268324e-3 + 
        SA *   1.9457531751183059e-3)
           )
              
    b1 = ( 0.5 * (1.1930681818531748e-2 + 
           CT2 * 2.6969148011830758e-7 + 
            SA * 5.9355685925035653e-6) )
    b2 = CT2 * -7.2734111712822707e-12 - 2.5943389807429039e-8
    b1sq = b1*b1
    sqrt_disc = np.sqrt(b1sq - b0*b2)

    N = a0 + (2*a3*b0*b1/b2 - a2*b0)/b2
    M = a1 + (4*a3*b1sq/b2 - (a3*b0 + 2*a2*b1))/b2

    A = b1 - sqrt_disc
    B = b1 + sqrt_disc
    
    delta_p = p_deep - p_shallow;
    p_sum = p_deep + p_shallow;
    part1 = b0 + p_shallow*(2*b1 + b2*p_shallow)
    part2 = (B + b2*p_deep)*(A + b2*p_shallow)
    part3 = (N*b2 - M*b1)/(b2*(B - A))

    return  ( db2Pa*(delta_p*(a2 - 2*a3*b1/b2 + 0.5*a3*p_sum)/b2 + 
              (M/(2*b2))*np.log(1 + delta_p*(2*b1 + b2*p_sum)/part1) +  
              part3*np.log(1 + delta_p*b2*(B - A)/part2)) )

