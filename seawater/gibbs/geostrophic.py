# -*- coding: utf-8 -*-

"""
Geostrophic streamfunctions, based on the 25-term expression for density

Functions
---------

  #geo_strf_dyn_height
  #    dynamic height anomaly
  #geo_strf_dyn_height_pc
  #    dynamic height anomaly for piecewise constant profiles
  #geo_strf_McD_Klocker
  #    McDougall-Klocker geostrophic streamfunction
  #geof_str_McD_Klocker_pc
  #    McDougall-Klocker geostrophic streamfunction for
  #    piecewise constant profiles
  #geo_strf_Montgomery
  #    Montgomery geostrophic streamfunction
  #geo_strf_Cunningham
  #    Cunningham geostrophic streamfunction
  #geostrophic_velocity
  #    geostrophic velocity

"""

from __future__ import division

import numpy as np
#import seawater.constants as cte
from earth import grav
from density25 import *

# ---------------------

__all__ = [#'geo_strf_dyn_height',
           #'geo_strf_dyn_height_pc',
           #'geo_strf_McD_Klocker',
           #'geof_str_McD_Klocker_pc',
           #'geo_strf_Montgomery',
           #'geo_strf_Cunningham',
           #'geostrophic_velocity'
          ]

# -----------------------

