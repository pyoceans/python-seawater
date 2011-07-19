# -*- coding: utf-8 -*-

"""
Practical Salinity (SP), PSS-78

Functions:
----------

  gsw_SP_from_C
    Practical Salinity from conductivity, C (inc. for SP < 2)
  gsw_C_from_SP
    conductivity, C, from Practical Salinity (inc. for SP < 2)
  gsw_SP_from_R
    Practical Salinity from conductivity ratio, R (inc. for SP < 2)
  gsw_R_from_SP
    conductivity ratio, R, from Practical Salinity (inc. for SP < 2)
  gsw_SP_salinometer
    Practical Salinity from a laboratory salinometer (inc. for SP < 2)


This is part of the python Gibbs Sea Water library
http://code.google.com/p/python-seawater/

"""

from __future__ import division

import numpy as np
from library import Hill_ratio_at_SP2

#import seawater.constants as cte
from utilities import match_args_return
#from conversions import z_from_p
#import seawater.csiro as sw


# ------------------------------

__all__ = ['SP_from_C',
           ]

# -------------------------------

@match_args_return
def SP_from_C(C,t,p):

    """Practical Salinity from conductivity 

 USAGE: 
  SP = gsw_SP_from_C(C,t,p)

 DESCRIPTION:
  Calculates Practical Salinity, SP, from conductivity, C, primarily using
  the PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical 
  Salinity is only valid in the range 2 < SP < 42.  If the PSS-78 
  algorithm produces a Practical Salinity that is less than 2 then the 
  Practical Salinity is recalculated with a modified form of the Hill et 
  al. (1986) formula.  The modification of the Hill et al. (1986)
  expression is to ensure that it is exactly consistent with PSS-78 
  at SP = 2.  Note that the input values of conductivity need to be in 
  units of mS/cm (not S/m). 

 INPUT:
  C  =  conductivity                                             [ mS/cm ]
  t  =  in-situ temperature (ITS-90)                             [ deg C ]
  p  =  sea pressure                                              [ dbar ]
        ( i.e. absolute pressure - 10.1325 dbar )

  OUTPUT:
   SP  =  Practical Salinity on the PSS-78 scale               [ unitless ]

 AUTHOR:  
   Paul Barker, Trevor McDougall and Rich Pawlowicz   [ help_gsw@csiro.au ]

 VERSION NUMBER: 3.0 (1st April, 2010)

  The software is available from http://www.TEOS-10.org

"""

    C, t, p = np.broadcast_arrays(C, t, p)

    a0 =  0.0080
    a1 = -0.1692
    a2 = 25.3851
    a3 = 14.0941
    a4 = -7.0261
    a5 =  2.7081

    b0 =  0.0005
    b1 = -0.0056
    b2 = -0.0066
    b3 = -0.0375
    b4 =  0.0636
    b5 = -0.0144

    c0 =  0.6766097
    c1 =  2.00564e-2
    c2 =  1.104259e-4
    c3 = -6.9698e-7
    c4 =  1.0031e-9
    
    d1 =  3.426e-2
    d2 =  4.464e-4
    d3 =  4.215e-1
    d4 = -3.107e-3

    e1 =  2.070e-5
    e2 = -6.370e-10
    e3 =  3.989e-15

    k  =  0.0162

    t68 = t * 1.00024
    ft68 = (t68 - 15) / (1 + k*(t68 - 15))

    # The dimensionless conductivity ratio, R, is the conductivity input, C,
    # divided by the present estimate of C(SP=35, t_68=15, p=0) which is 
    # 42.9140 mS/cm (=4.29140 S/m), (Culkin and Smith, 1980). 

    R = 0.023302418791070513 * C  # 0.023302418791070513 = 1./42.9140

    # rt_lc corresponds to rt as defined in the UNESCO 44 (1983) routines.  
    rt_lc = c0 + (c1 + (c2 + (c3 + c4*t68)*t68)*t68)*t68
    Rp = ( 1 + (p*(e1 + e2*p + e3*p*p)) / 
          (1 + d1*t68 + d2*t68*t68 + (d3 + d4*t68)*R) )
    Rt = R / (Rp*rt_lc)   

    Rt[Rt < 0] = np.nan
    Rtx = np.sqrt(Rt)

    SP = ( a0 + (a1 + (a2 + (a3 + (a4 + a5*Rtx)*Rtx)*Rtx)*Rtx)*Rtx + 
        ft68*(b0 + (b1 + (b2 + (b3 + (b4 + b5*Rtx)*Rtx)*Rtx)*Rtx)*Rtx) )

    # The following section of the code is designed for SP < 2 based on the
    # Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
    # exactly equal to the PSS-78 algorithm at SP = 2.

    I2, = np.nonzero(np.ravel(SP) < 2)  # find
    if len(I2) > 0:
        Hill_ratio = Hill_ratio_at_SP2(t[I2]) 
        x = 400 * Rt[I2]
        sqrty = 10 * Rtx[I2]
        part1 = 1 + x * (1.5 + x) 
        part2 = 1 + sqrty * (1 + sqrty * (1 + sqrty))
        SP_Hill_raw = SP[I2] - a0/part1 - b0*ft68[I2]/part2
        SP[I2] = Hill_ratio * SP_Hill_raw

    # Ensure that SP is non-negative.
    SP[SP < 0] = 0.0

    return SP

