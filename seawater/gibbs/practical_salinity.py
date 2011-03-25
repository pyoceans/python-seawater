# -*- coding: utf-8 -*-

"""
Practical Salinity (SP), PSS-78

Functions:
----------
  
  SP_from_cndr(R, t, p)
     Practical Salinity from conductivity ratio
  cndr_from_SP(SP, t, p)                    
     conductivity ratio from Practical Salinity

This is part of the python Gibbs Sea Water library
http://code.google.com/p/python-seawater/

"""

from __future__ import division

#import numpy as np
#import seawater.constants as cte
#from library import match_args_return
#from conversions import z_from_p
import seawater.csiro as sw

# ------------------------------

__all__ = ['SP_from_cndr',
           'cndr_from_SP']

# -------------------------------

# Presently imports from old seawater package
# may be reimplemented to make the gibbs package
# independent

SP_from_cndr = sw.salt

cndr_from_SP = sw.cndr
