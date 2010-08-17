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

The present version is 1.0.4, released 17 August 2010.

User documentation:
-------------------

The functions input must be called as NumPy arrays, returning array of the same shape.

Original seawater functions
---------------------------

+---------------------+-------------+----------------------------------------------------------------------------------------------+
| function            |  units      |  description                                                                                 |
+=====================+=============+==============================================================================================+
| adtg(S, T, P)       | K / dbar    | Calculates adiabatic temperature gradient as per UNESCO 1983 routines.                       |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| alpha(S, PTMP, P)   | 1 / K       | A function to calculate the thermal expansion coefficient.                                   |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| aonb(S, PTMP, P)    | psu / degC  | Calculate alpha/beta.  See alpha and beta.                                                   |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| beta(S, PTMP, P)    | psu**-1     | The saline contraction coefficient as defined by T.J. McDougall.                             |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| bfrq(S, T, P, LAT)  | s**-2       | Calculates Brunt-Vaisala Frequency squared (N^2) at the mid depths.                          |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| depth(P, LAT)       | m           | Calculates depth in metres from pressure in dbars.                                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| grav(LAT, z=0)      | m / s**2    | Calculates acceleration due to gravity as function of latitude.                              |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| cor(LAT)            | s**-1       | Calculates the Coriolis factor "f" defined by f = 2*Omega*Sin(lat)                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| c3515()             | mS / cm     | Returns conductivity at S=35 psu , T=15 C [ITPS 68] and P=0 db).                             |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| cndr(S, T, P)       | no units    | Calculates conductivity ratio from S, T, P                                                   |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| sals(Rt, T)         | psu         | Salinity of sea water as a function of Rt and T (PSS-78). UNESCO 1983 polynomial.            |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| salds(Rtx, delT)    | no units    | Calculates Salinity differential dS/d(sqrt(Rt)) at constant T. UNESCO 1983 polynomial.       |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| salrt(T)            | no units    | Equation rt(T) = C(35,T,0) / C(35,15(IPTS-68), 0) used in calculating salinity.              |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| salt(cndr, T, P)    | no units    | Calculates Salinity from conductivity ratio. UNESCO 1983 polynomial.                         |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| salrp(R, T, P)      | no units    | Equation Rp(S,T,P) = C(S,T,P)/C(S,T,0) used in calculating salinity.                         |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| fp(S, P)            | degC        | Freezing point of Sea Water using UNESCO 1983 polynomial.                                    |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| svel(S, T, P)       | m / s       | Sound Velocity in sea water using UNESCO 1983 polynomial.                                    |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| pres(DEPTH, LAT)    | dbar        | Calculates pressure in dbars from depth in meters.                                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| dist(lon, lat)      | m           | Calculate distance between two positions on globe.                                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| satAr(S, T)         | ml / l      | Solubility (satuaration) of Argon (Ar) in sea water.                                         |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| satN2(S, T)         | ml / l      | Solubility (satuaration) of Nitrogen (N2) in sea water.                                      |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| satO2(S,T)          | ml / l      | Solubility (satuaration) of Oxygen (O2) in sea water.                                        |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| dens0(S,T)          | kg / m**3   | Density of Sea Water at atmospheric pressure.                                                |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| smow(T)             | kg / m**3   | Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980.                            |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| seck(S, T, P=0)     | bars        | Secant Bulk Modulus (K) of Sea Water using Equation of state 1980.                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| dens(S, T, P)       | kg / m**3   | Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.                                  |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| pden(S, T, P, PR=0) | kg / m**3   | Calculates potential density of water mass relative to a reference pressure.                 |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| svan(S, T, P=0)     | m**3 / kg   | Specific Volume Anomaly.                                                                     |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| gpan(S, T, P)       | m**2 s**-2  | Geopotential Anomaly.                                                                        |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| gvel(ga, lon, lat)  | m / s       | Calculates geostrophic velocity given the geopotential anomaly and position of each station. |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| gvel2(ga, dist, lat)| m / s       | Calculates geostrophic velocity given the geopotential anomaly and distance.                 |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| cp(S, T, P)         | J / (kg*K)  | Heat Capacity of Sea Water using UNESCO 1983 polynomial.                                     |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| ptmp(S, T, P, PR=0) | degC        | Calculates potential temperature as per UNESCO 1983 report.                                  |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| temp(S, PTMP, P, PR)| degC        | Calculates temperature from potential temperature.                                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| swvel(lenth, depth) | m / s       | Calculates surface wave velocity.                                                            |
+---------------------+-------------+----------------------------------------------------------------------------------------------+

Extra Functions
----------------

+---------------------+-------------+----------------------------------------------------------------------------------------------+
| function            |  units      |  description                                                                                 |
+=====================+=============+==============================================================================================+
| sigma(S, T, P)      | kg / m**3   | Density of Sea Water in sigma.                                                               |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| drhodt(S, T, P)     | kg /(K*m**3)| Temperature derivative of density.                                                           |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| drhods(S, T, P)     | kg / m**3   | Salinity derivative of density.                                                              |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
| test()              |   --        | Execute test routines.                                                                       |
+---------------------+-------------+----------------------------------------------------------------------------------------------+
"""

from seawater import *
from extras import *

__authors__    = ['Filipe Fernandes']
__copyright__  = "CSIRO"
__credits__    = ["Filipe Fernandes", "Lindsay Pender","Phil Morgan"]
__license__    = "CSIRO"
__version__    = "1.0.4"
__maintainer__ = "Filipe Fernandes"
__email__      = "ocefpaf@gmail.com"
__status__     = "Production"
__all__        = ['seawater']
