# -*- coding: utf-8 -*-

"""
===============
Python Seawater
===============

Introduction:
-------------

This module contains a translation of the original CSIRO Matlab package (SEAWATER-3.2) for calculating the properties of sea water.
It consists of a self contained library easy to use. The only requirent is NumPy.

The author has no intention to do things in a "pythonic-way", it is just a "work around" from someone that couldn't afford Matlab anymore.

The package uses the formulas from Unesco's joint panel on oceanographic tables and standards, UNESCO 1981 and UNESCO 1983.

The present version is 1.0.4, released 23 August 2010.

User documentation:
-------------------

The functions input must be called as NumPy arrays, returning array of the same shape.

Original seawater functions
---------------------------

+---------------------+-------------+----------------------------------------------------------------------------------------------+
| function            |  units      |  description                                                                                 |
+=====================+=============+==============================================================================================+
| adtg(s, t, p)       | K / dbar    | Calculates adiabatic temperature gradient as per UNESCO 1983 routines.                       |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| alpha(s, t, p)      | 1 / K       | A function to calculate the thermal expansion coefficient.                                   |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| aonb(s, t, p)       | psu / degC  | Calculate alpha/beta.  See alpha and beta.                                                   |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| beta(s, t, p)       | psu**-1     | The saline contraction coefficient as defined by T.J. McDougall.                             |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| bfrq(s, t, p, lat)  | s**-2       | Calculates Brunt-Vaisala Frequency squared (N^2) at the mid depths.                          |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| depth(p, lat)       | m           | Calculates depth in metres from pressure in dbars.                                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| grav(lat, z=0)      | m / s**2    | Calculates acceleration due to gravity as function of latitude.                              |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| cor(lat)            | s**-1       | Calculates the Coriolis factor "f" defined by f = 2*Omega*Sin(lat)                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| cndr(s, t, p)       | no units    | Calculates conductivity ratio from S, T, P                                                   |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| sals(rt, t)         | psu         | Salinity of sea water as a function of Rt and T (PSS-78). UNESCO 1983 polynomial.            |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| salds(rtx, delt)    | no units    | Calculates Salinity differential dS/d(sqrt(Rt)) at constant T. UNESCO 1983 polynomial.       |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| salrt(t)            | no units    | Equation rt(T) = C(35,T,0) / C(35,15(IPTS-68), 0) used in calculating salinity.              |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| salt(cndr, t, p)    | no units    | Calculates Salinity from conductivity ratio. UNESCO 1983 polynomial.                         |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| salrp(r, t, p)      | no units    | Equation Rp(S,T,P) = C(S,T,P)/C(S,T,0) used in calculating salinity.                         |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| fp(s, p)            | degC        | Freezing point of Sea Water using UNESCO 1983 polynomial.                                    |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| svel(s, t, p)       | m / s       | Sound Velocity in sea water using UNESCO 1983 polynomial.                                    |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| pres(depth, lat)    | dbar        | Calculates pressure in dbars from depth in meters.                                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| dist(lon, lat)      | m           | Calculate distance between two positions on globe.                                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| satAr(s, t)         | ml / l      | Solubility (satuaration) of Argon (Ar) in sea water.                                         |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| satN2(s, t)         | ml / l      | Solubility (satuaration) of Nitrogen (N2) in sea water.                                      |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| satO2(s, t)         | ml / l      | Solubility (satuaration) of Oxygen (O2) in sea water.                                        |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| dens0(s, t)         | kg / m**3   | Density of Sea Water at atmospheric pressure.                                                |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| smow(t)             | kg / m**3   | Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980.                            |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| seck(s, t, p=0)     | bars        | Secant Bulk Modulus (K) of Sea Water using Equation of state 1980.                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| dens(s, t, p)       | kg / m**3   | Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.                                  |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| pden(s, t, p, pr=0) | kg / m**3   | Calculates potential density of water mass relative to a reference pressure.                 |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| svan(s, t, p=0)     | m**3 / kg   | Specific Volume Anomaly.                                                                     |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| gpan(s, t, p)       | m**2 s**-2  | Geopotential Anomaly.                                                                        |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| gvel(ga, lon, lat)  | m / s       | Calculates geostrophic velocity given the geopotential anomaly and position of each station. |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| gvel2(ga, dist, lat)| m / s       | Calculates geostrophic velocity given the geopotential anomaly and distance.                 |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| cp(s, t, p)         | J / (kg*K)  | Heat Capacity of Sea Water using UNESCO 1983 polynomial.                                     |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| ptmp(s, t, p, pr=0) | degC        | Calculates potential temperature as per UNESCO 1983 report.                                  |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| temp(s, pt, p, pr)  | degC        | Calculates temperature from potential temperature.                                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| swvel(lenth, depth) | m / s       | Calculates surface wave velocity.                                                            |
+---------------------+-------------+----------------------------------------------------------------------------------------------+

Extra Functions
----------------

+---------------------+-------------+----------------------------------------------------------------------------------------------+
| function            |  units      |  description                                                                                 |
+=====================+=============+==============================================================================================+
| sigma(s, t, p)      | kg / m**3   | Density of Sea Water in sigma.                                                               |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| drhodt(s, t, p)     | kg /(K*m**3)| Temperature derivative of density.                                                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| drhods(s, t, p)     | kg / m**3   | Salinity derivative of density.                                                              |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| test()              |   --        | Execute test routines.                                                                       |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
"""

import numpy as mp
# TODO: only used on the test routine
from os    import uname
from time  import asctime, localtime
from sys   import version

from seawater import *
from extras import *

"""
seawater.py
Translated from matlab CSIRO seawater toolbox Version 3.2

Filipe P. A. Fernandes
e-mail:   ocefpaf@gmail.com
web:      http://ocefpaf.tiddlyspot.com/
date:     14-Jan-2010
modified: 17-Aug-2010
obs:      flag: TODO
"""

# CONSTANTS:

"""
The International Practical Temperature Scale of 1968 (IPTS-68)
#:math:`T68  = 1.00024 * T90`
this linear transformation is accurate within 0.5 m C for conversion between IPTS-68 and ITS-90 over the oceanographic temperature range (Saunders,et al 1991).
"""
T68conv  = 1.00024

"""
0.017453292519943295
"""
DEG2RAD = np.pi/180.

"""
A.E.Gill p.597
..:math:
  \Omega = \frac{2*\\pi}{\\textrm{sidereal day}}

1 sidereal day = 23.9344696 hours
units : radians/sec
"""
OMEGA   = 7.292e-5

"""
Conductivity at S=35 psu , T=15 C [ITPS 68] and P=0 db)
units : mmho cm :sup:`-1` == mS cm :sup:`-1`
Reference: R.C. Millard and K. Yang 1992. "CTD Calibration and Processing Methods used by Woods Hole Oceanographic Institution"  Draft April 14, 1992 (Personal communication)
"""
C3515   = 42.914

"""
acceleration of gravity in m s :sup:`2`
"""
g = 9.8

__authors__    = ['Filipe Fernandes']
__copyright__  = "CSIRO"
__credits__    = ["Filipe Fernandes", "Lindsay Pender","Phil Morgan"]
__license__    = "CSIRO"
__version__    = "1.0.4"
__maintainer__ = "Filipe Fernandes"
__email__      = "ocefpaf@gmail.com"
__status__     = "Production"
__all__        = ['seawater']