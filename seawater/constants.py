"""Constants.
=========

"""

from numpy import pi

# dbar to pascal.
db2Pascal = 1e4

# The Celsius zero point.
Kelvin = 273.15

# Acceleration of gravity [m/s**2]
gdef = 9.8

# 1 nm = 1.8520 km
DEG2NM, NM2KM = 60.0, 1.8520

# Sidereal day = 23.9344696 hours.
OMEGA = 7.292e-5  # 7.292115e-5

# Mean radius of earth [m] A.E. Gill.
earth_radius = 6371000.0

# Angle conversions.
deg2rad, rad2deg = pi / 180.0, 180.0 / pi

# Conductivity at S=35 psu , T=15 C [ITPS 68] and P=0 dbar.
c3515 = 42.914
