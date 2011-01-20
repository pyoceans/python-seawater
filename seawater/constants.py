"""
Constants used
==============
"""

a = 6371000.
""" mean radius of earth  A.E.Gill """

OMEGA = 7.292115e-5
"""
:math:`\\Omega = \\frac{2\\pi}{\\textrm{sidereal day}}` = 7.292e-5.radians sec :sup:`-1`

1 sidereal day = 23.9344696 hours

Changed to a more precise value at Groten 2004

References
----------
.. [1] A.E. Gill 1982. p.54  eqn 3.7.15 "Atmosphere-Ocean Dynamics" Academic Press: New York. ISBN: 0-12-283522-0. page: 597
.. [2] Groten, E., 2004: Fundamental Parameters and Current (2004) Best Estimates of the Parameters of Common Relevance to Astronomy, Geodesy, and Geodynamics. Journal of Geodesy, 77, pp. 724-797.
"""

C3515 = 42.914 #FIXME: not used?
"""
Conductivity of 42.914 [mmho cm :sup:`-1` == mS cm :sup:`-1`] at Salinity 35 psu, Temperature 15 :math:`^\\circ` C [ITPS 68] and Pressure 0 db.

References
----------
.. [1] R.C. Millard and K. Yang 1992. "CTD Calibration and Processing Methods used by Woods Hole Oceanographic Institution" Draft April 14, 1992 (Personal communication)
"""


gdef = 9.8
"""
Acceleration of gravity [m s :sup:`2`] used by sw.swvel and bfrq without lat info.
"""

DEG2NM = 60
NM2KM = 1.8520
"""
1 nm = 1.8520 km
Used by sw.dist() to convert nautical miles to kilometers.

References
----------
.. [1] S. Pond & G.Pickard 2nd Edition 1986 Introductory Dynamical Oceanogrpahy Pergamon Press Sydney. ISBN 0-08-028728-X. page: 303
"""

Kelvin = 273.15
"""
offset to convert :math:`^\\circ` C to Kelvin. Used in the gas solubility functions.
"""

db2Pascal = 1e4
"""
Decibar to pascal
"""

#--- only for gsw
gamma = 2.26e-7
"""
TODO: define Gamma
"""

M_S = 0.0314038218
"""mole-weighted average atomic weight of the elements of sea salt, in units of kg mol :sup:`-1` """

cp0 = 3991.86795711963
""" FIXME: from Eqn. (3.3.3) of IOC et al. (2010). """

SSO = 35.16504
""" from section 2.4 of IOC et al. (2010) """

sfac = 0.0248826675584615
""" TODO: sfac = 1 / (40 * ( 35.16504 / 35 ) ) """

R = 8.314472
""" the molar gas constant = 8.314472 m :sup:`2` kg s:sup:`-21 K :sup:`-1` mol :sup:`-1` """