"""
Constants used
==============
"""

from __future__ import division

C3515 = 42.914  # NOTE: reference value, not used.
"""
Conductivity of 42.914 [mmho cm :sup:`-1` == mS cm :sup:`-1`] at
Salinity 35 psu, Temperature 15 :math:`^\\circ` C [ITPS 68] and Pressure 0 db.

References
----------
.. [1] R.C. Millard and K. Yang 1992. "CTD Calibration and Processing Methods
used by Woods Hole Oceanographic Institution" Draft April 14, 1992
(Personal communication)
"""

earth_radius = 6371000.
"""
mean radius of earth  A.E. Gill
"""

OMEGA = 7.292115e-5
"""
:math:`\\Omega = \\frac{2\\pi}{\\textrm{sidereal day}}` = 7.292e-5.radians sec :sup:`-1`

1 sidereal day = 23.9344696 hours

Changed to a more precise value at Groten 2004

References
----------
.. [1] A.E. Gill 1982. p.54  eqn 3.7.15 "Atmosphere-Ocean Dynamics" Academic
Press: New York. ISBN: 0-12-283522-0. page: 597
.. [2] Groten, E., 2004: Fundamental Parameters and Current (2004) Best
Estimates of the Parameters of Common Relevance to Astronomy, Geodesy, and
Geodynamics. Journal of Geodesy, 77, pp. 724-797.
"""

gdef = 9.8
"""
Acceleration of gravity [m s :sup:`2`] used by sw.swvel and bfrq without lat
info.
"""

DEG2NM = 60
NM2KM = 1.8520
"""
1 nm = 1.8520 km
Used by sw.dist() to convert nautical miles to kilometers.

References
----------
.. [1] S. Pond & G.Pickard 2nd Edition 1986 Introductory Dynamical
Oceanography Pergamon Press Sydney. ISBN 0-08-028728-X. page: 303
"""

Kelvin = 273.15
"""
Offset to convert :math:`^\\circ` C to Kelvin. Used in the gas solubility
functions.
"""

db2Pascal = 1e4
"""
Decibar to pascal
"""

#--- only for gsw
gamma = 2.26e-7  # TODO: define
"""
Gamma (A.E. Gill)
"""

M_S = 0.0314038218
"""
Mole-weighted average atomic weight of the elements of sea salt, in units of
kg mol :sup:`-1`
"""

cp0 = 3991.86795711963  # TODO: define
"""
from Eqn. (3.3.3) of IOC et al. (2010).
"""

SSO = 35.16504  # TODO: define
"""
from section 2.4 of IOC et al. (2010)
"""

sfac = 0.0248826675584615  # TODO: define
"""
sfac = 1 / (40 * ( SSO / 35 ) )
"""

R = 8.314472
"""
the molar gas constant = 8.314472 m :sup:`2` kg s:sup:`-21 K :sup:`-1` mol :sup:`-1`
"""

r1 = 0.35
"""
TODO
"""

uPS = 35.16504 / 35
"""
The unit conversion factor for salinities (35.16504/35) g/kg (Millero et
al., 2008). Reference Salinity SR is uPS times Practical Salinity SP.

Ratio, unit conversion factor for salinities [g kg :sup:`-1`]

References
----------
Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: The
composition of Standard Seawater and the definition of the
Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72. See
section 6, Eqn.(6.1).
"""
