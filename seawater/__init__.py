# -*- coding: utf-8 -*-

"""
 Seawater -- Python functions for properties of sea water
    Version 1.0-3.2 30 January 2010

 Functions:
 adtg(S, T, P)          Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
 alpha(S, PTMP, P)      A function to calculate the thermal expansion coefficient.
 aonb(S, PTMP, P)       Calculate alpha/beta.  See alpha and beta.
 beta(S, PTMP, P)       The saline contraction coefficient as defined by T.J. McDougall.
 bfrq(S, T, P, LAT)     Calculates Brunt-Vaisala Frequency squared (N^2) at the mid depths.
 depth(P, LAT)          Calculates depth in metres from pressure in dbars.
 grav(LAT, z=0)         Calculates acceleration due to gravity as function of latitude.
 cor(LAT)               Calculates the Coriolis factor "f" defined by f = 2*Omega*Sin(lat)
 c3515()                Returns conductivity at S=35 psu , T=15 C [ITPS 68] and P=0 db).
 cndr(S, T, P)          Calculates conductivity ratio from S, T, P
 sals(Rt, T)            Salinity of sea water as a function of Rt and T. UNESCO 1983 polynomial.
 salds(Rtx, delT)       Calculates Salinity differential dS/d(sqrt(Rt)) at constant T. UNESCO 1983 polynomial.
 salrt(T)               Equation rt(T) = C(35,T,0) / C(35,15(IPTS-68), 0) used in calculating salinity.
 salt(cndr, T, P)       Calculates Salinity from conductivity ratio. UNESCO 1983 polynomial.
 salrp(R, T, P)         Equation Rp(S,T,P) = C(S,T,P)/C(S,T,0) used in calculating salinity.
 fp(S, P)               Freezing point of Sea Water using UNESCO 1983 polynomial.
 svel(S, T, P)          Sound Velocity in sea water using UNESCO 1983 polynomial.
 pres(DEPTH, LAT)       Calculates pressure in dbars from depth in meters.
 dist(lon, lat)         Calculate distance between two positions on globe.
 satAr(S, T)            Solubility (satuaration) of Argon (Ar) in sea water.
 satN2(S, T)            Solubility (satuaration) of Nitrogen (N2) in sea water.
 satO2(S,T)             Solubility (satuaration) of Oxygen (O2) in sea water.
 dens0(S,T)             Density of Sea Water at atmospheric pressure.
 smow(T)                Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980.
 seck(S, T, P=0)        Secant Bulk Modulus (K) of Sea Water using Equation of state 1980.
 dens(S, T, P)          Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.
 pden(S, T, P, PR=0)    Calculates potential density of water mass relative to a reference pressure.
 svan(S, T, P=0)        Specific Volume Anomaly.
 gpan(S, T, P)          Geopotential Anomaly.
 gvel(ga, lon, lat)     Calculates geostrophic velocity given the geopotential anomaly and position of each station.
 gvel2(ga, dist, lat)   Calculates geostrophic velocity given the geopotential anomaly and distance.
 cp(S, T, P)            Heat Capacity of Sea Water using UNESCO 1983 polynomial.
 ptmp(S, T, P, PR=0)    Calculates potential temperature as per UNESCO 1983 report.
 temp(S, PTMP, P, PR)   Calculates temperature from potential temperature.
 swvel(lenth, depth)    Calculates surface wave velocity.
 test()                 Execute test routines.
"""

# --- Exceptions ---
class OutOfRangeError(Exception): pass

from seawater import *

__all__     = ['seawater']
__authors__ = ['Filipe Fernandes <ocefpaf@gmail.com>']
__version__ = '1.0.3'
