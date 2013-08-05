# -*- coding: utf-8 -*-
#
# extras.py
#
# purpose:  Non EOS-08 functions.
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.github.io/
# created:  05-Aug-2013
# modified: Mon 05 Aug 2013 10:02:38 AM BRT
#
# obs:
#


from __future__ import division

import numpy as np
from library import T68conv
from constants import OMEGA, DEG2NM, NM2KM, Kelvin, deg2rad, rad2deg, gdef

__all__ = ['dist',
           'f',
           'satAr',
           'satN2',
           'satO2',
           'swvel']


def dist(lon, lat, units='km'):
    """Calculate distance between two positions on globe using the "Plane
    Sailing" method. Also uses simple geometry to calculate the bearing of
    the path between position pairs.

    Parameters
    ----------
    lon : array_like
          decimal degrees (+ve E, -ve W) [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [- 90.. +90]
    units : string, optional
            default kilometers

    Returns
    -------
    dist : array_like
           distance between positions in units
    phaseangle : array_like
                 angle of line between stations with x axis (East).
                 Range of values are -180..+180. (E=0, N=90, S=-90)

    Examples
    --------
    >>> import seawater as sw
    >>> lon = [35, 35]
    >>> lat = [41, 40]
    >>> sw.dist(lon, lat)
    (array([ 111.12]), array([-90.]))
    >>> # Create a distance vector.
    >>> lon = np.arange(30,40,1)
    >>> lat = 35
    >>> np.cumsum(np.append(0, sw.dist(lon, lat, units='km')[0]))
    array([   0.        ,   91.02417516,  182.04835032,  273.07252548,
            364.09670065,  455.12087581,  546.14505097,  637.16922613,
            728.19340129,  819.21757645])

    References
    ----------
    .. [1] The PLANE SAILING method as described in "CELESTIAL NAVIGATION" 1989
    by Dr. P. Gormley. The Australian Antarctic Division.

    Modifications: 92-02-10. Phil Morgan and Steve Rintoul.
                   99-06-25. Lindsay Pender, name change distance to sw_dist.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
    """

    lon, lat = map(np.asanyarray, (lon, lat))

    if lat.size == 1:
        lat = np.repeat(lat, lon.size)
    elif lon.size == 1:
        lon = np.repeat(lon, lat.size)

    npositions = max(lon.shape)

    ind = np.arange(0, npositions - 1, 1)  # Index to first of position pairs.

    dlon = np.diff(lon, axis=0)
    if np.any(np.abs(dlon) > 180):
        flag = abs(dlon) > 180
        dlon[flag] = -np.sign(dlon[flag]) * (360 - np.abs(dlon[flag]))

    latrad = np.abs(lat * deg2rad)
    dep = np.cos((latrad[ind + 1] + latrad[ind]) / 2) * dlon
    dlat = np.diff(lat, axis=0)
    dist = DEG2NM * (dlat ** 2 + dep ** 2) ** 0.5

    if units == 'km':
        dist = dist * NM2KM

    # Calculate angle to x axis.
    phaseangle = np.angle(dep + dlat * 1j) * rad2deg
    return dist, phaseangle


def f(lat):
    """Calculates the Coriolis factor :math:`f` defined by:

    .. math::
        f = 2 \Omega \sin(lat)

    where:

    .. math::
        \Omega = \frac{2 \pi}{\textrm{sidereal day}} = 7.2921150e^{-5}
        \textrm{ radians sec}^{-1}


    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90].

    Returns
    -------
    f : array_like
        Coriolis factor [s :sup:`-1`]

    Examples
    --------
    >>> import seawater as sw
    >>> sw.f(45)
    0.00010312607931384281

    References
    ----------
    .. [1] S. Pond & G.Pickard 2nd Edition 1986 Introductory Dynamical
    Oceanography Pergamon Press Sydney. ISBN 0-08-028728-X

    .. [2] A.E. Gill 1982. p.54  Eqn. 3.7.15 "Atmosphere-Ocean Dynamics"
    Academic Press: New York. ISBN: 0-12-283522-0

    .. [3] Groten, E., 2004: Fundamental Parameters and Current (2004) Best
    Estimates of the Parameters of Common Relevance to Astronomy, Geodesy,
    and Geodynamics. Journal of Geodesy, 77, pp. 724-797.

    Modifications: 93-04-20. Phil Morgan.
    """
    lat = np.asanyarray(lat)
    # Eqn p27.  UNESCO 1983.
    return 2 * OMEGA * np.sin(lat * deg2rad)


def satAr(s, t):
    """Solubility (saturation) of Argon (Ar) in sea water.

    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    satAr : array_like
            solubility of Ar [ml l :sup:`-1`]

    Examples
    --------
    Data from Weiss 1970.
    >>> import seawater as sw
    >>> t = sw.T90conv([[ -1, -1], [ 10, 10], [ 20, 20], [ 40, 40]])
    >>> s = [[ 20, 40], [ 20, 40], [ 20, 40], [ 20, 40]]
    >>> sw.satAr(s, t)
    array([[ 0.4455784 ,  0.38766011],
           [ 0.33970659,  0.29887756],
           [ 0.27660227,  0.24566428],
           [ 0.19861429,  0.17937698]])

    References
    ----------
    .. [1] Weiss, R. F. 1970. The Solubility of Nitrogen, Oxygen and Argon in
    Water and Seawater Deep-Sea Research Vol. 17, p. 721-735.
    doi:10.1016/0011-7471(70)90037-9

    Modifications: 97-11-05. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
    """

    s, t = map(np.asanyarray, (s, t))

    # Convert T to Kelvin.
    t = Kelvin + T68conv(t)

    # Constants for Eqn (4) of Weiss 1970.
    a = [-173.5146, 245.4510, 141.8222, -21.8020]
    b = [-0.034474, 0.014934, -0.0017729]

    # Eqn (4) of Weiss 1970.
    lnC = (a[0] + a[1] * (100 / t) + a[2] * np.log(t / 100) + a[3] *
           (t / 100) + s * (b[0] + b[1] * (t / 100) + b[2] * ((t / 100) ** 2)))

    return np.exp(lnC)


def satN2(s, t):
    """Solubility (saturation) of Nitrogen (N2) in sea water.

    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    satN2 : array_like
            solubility of N2  [ml l :sup:`-1`]

    Examples
    --------
    Data from Weiss 1970.
    >>> import seawater as sw
    >>> t = sw.T90conv([[ -1, -1], [ 10, 10], [ 20, 20], [ 40, 40]])
    >>> s = [[ 20, 40], [ 20, 40], [ 20, 40], [ 20, 40]]
    >>> sw.satN2(s, t)
    array([[ 16.27952432,  14.00784526],
           [ 12.64036196,  11.01277257],
           [ 10.46892822,   9.21126859],
           [  7.78163876,   6.95395099]])


    References
    ----------
    .. [1] Weiss, R. F. 1970. The Solubility of Nitrogen, Oxygen and Argon in
    Water and Seawater Deep-Sea Research Vol. 17, p. 721-735.
    doi:10.1016/0011-7471(70)90037-9

    Modifications: 97-11-05. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
    """

    s, t = map(np.asanyarray, (s, t))

    # Convert T to Kelvin.
    t = Kelvin + T68conv(t)

    # Constants for Eqn (4) of Weiss 1970.
    a = (-172.4965, 248.4262, 143.0738, -21.7120)
    b = (-0.049781, 0.025018, -0.0034861)

    # Eqn (4) of Weiss 1970.
    lnC = (a[0] + a[1] * (100 / t) + a[2] * np.log(t / 100) + a[3] *
           (t / 100) + s * (b[0] + b[1] * (t / 100) + b[2] * ((t / 100) ** 2)))

    return np.exp(lnC)


def satO2(s, t):
    """Solubility (saturation) of Oxygen (O2) in sea water.

    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [:math:`^\circ` C (ITS-68)]

    Returns
    -------
    satO2 : array_like
            solubility of O2  [ml l :sup:`-1` ]

    Examples
    --------
    Data from Weiss 1970
    >>> import seawater as sw
    >>> t = sw.T90conv([[ -1, -1], [ 10, 10], [ 20, 20], [ 40, 40]])
    >>> s = [[ 20, 40], [ 20, 40], [ 20, 40], [ 20, 40]]
    >>> sw.satO2(s, t)
    array([[ 9.162056  ,  7.98404249],
           [ 6.95007741,  6.12101928],
           [ 5.64401453,  5.01531004],
           [ 4.0495115 ,  3.65575811]])

    References
    ----------
    .. [1] Weiss, R. F. 1970. The Solubility of Nitrogen, Oxygen and Argon in
    Water and Seawater Deep-Sea Research Vol. 17, p. 721-735.
    doi:10.1016/0011-7471(70)90037-9

    Modifications: 97-11-05. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
    """

    s, t = map(np.asanyarray, (s, t))

    # Convert T to Kelvin.
    t = Kelvin + T68conv(t)

    # Constants for Eqn (4) of Weiss 1970.
    a = (-173.4292, 249.6339, 143.3483, -21.8492)
    b = (-0.033096, 0.014259, -0.0017000)

    # Eqn (4) of Weiss 1970.
    lnC = (a[0] + a[1] * (100 / t) + a[2] * np.log(t / 100) + a[3] *
           (t / 100) + s * (b[0] + b[1] * (t / 100) + b[2] * ((t / 100) ** 2)))

    return np.exp(lnC)


def swvel(length, depth):
    """Calculates surface wave velocity.

    length : array_like
            wave length
    depth : array_like
            water depth [meters]

    Returns
    -------
    speed : array_like
            surface wave speed [m s :sup:`-1`]

    Examples
    --------
    >>> import seawater as sw
    >>> sw.swvel(10, 100)
    3.9493270848342941

    Modifications: Lindsay Pender 2005
    """
    length, depth = map(np.asanyarray, (length, depth))
    k = 2.0 * np.pi / length
    return np.sqrt(gdef * np.tanh(k * depth) / k)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
