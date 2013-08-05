# -*- coding: utf-8 -*-
#
# geostrophic.py
#
# purpose:  Geostrophic velocity calculation.
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.github.io/
# created:  05-Aug-2013
# modified: Mon 05 Aug 2013 10:35:28 AM BRT
#
# obs:
#

from __future__ import division

import numpy as np

from eos80 import dens
from extras import dist, f
from constants import db2Pascal

__all__ = ['svan',
           'gpan',
           'gvel']


def svan(s, t, p=0):
    """Specific Volume Anomaly calculated as
    svan = 1 / dens(s, t, p) - 1 / dens(35, 0, p).

    Note that it is often quoted in literature as 1e8 * units.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    svan : array_like
           specific volume anomaly  [m :sup:`3` kg :sup:`-1`]

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22.
    >>> import seawater as sw
    >>> s = [[0, 1, 2], [15, 16, 17], [30, 31, 32], [35, 35, 35]]
    >>> t = sw.T90conv([[15]*3]*4)
    >>> p = [[0], [250], [500], [1000]]
    >>> sw.svan(s, t, p)
    array([  2.74953924e-05,   2.28860986e-05,   3.17058231e-05,
             3.14785290e-05,   0.00000000e+00,   0.00000000e+00,
             6.07141523e-06,   9.16336113e-06])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] S. Pond & G.Pickard 2nd Edition 1986 Introductory Dynamical
    Oceanography Pergamon Press Sydney. ISBN 0-08-028728-X

    Modifications: 92-11-05. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
    """
    s, t, p = map(np.asanyarray, (s, t, p))
    return 1 / dens(s, t, p) - 1 / dens(35, 0, p)


def gpan(s, t, p):
    """Geopotential Anomaly calculated as the integral of svan from the
    the sea surface to the bottom. THUS RELATIVE TO SEA SURFACE.

    Adapted method from Pond and Pickard (p76) to calculate gpan relative to
    sea surface whereas P&P calculated relative to the deepest common depth.
    Note that older literature may use units of "dynamic decimeter" for above.


    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    gpan : array_like
           geopotential anomaly
           [m :sup:`3` kg :sup:`-1`
            Pa = m :sup:`2` s :sup:`-2` = J kg :sup:`-1`]

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22.
    >>> import seawater as sw
    >>> s = [[0, 1, 2], [15, 16, 17], [30, 31, 32], [35,35,35]]
    >>> t = [[15]*3]*4
    >>> p = [[0], [250], [500], [1000]]
    >>> sw.gpan(s, t, p)
    array([[   0.        ,    0.        ,    0.        ],
           [  56.35465209,   56.35465209,   56.35465209],
           [  84.67266947,   84.67266947,   84.67266947],
           [ 104.95799186,  104.95799186,  104.95799186]])

    References
    ----------
    .. [1] S. Pond & G.Pickard 2nd Edition 1986 Introductory Dynamical
    Oceanography Pergamon Press Sydney. ISBN 0-08-028728-X

    Modifications: 92-11-05. Phil Morgan.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
    """

    s, t, p = map(np.asanyarray, (s, t, p))
    s, t, p = np.broadcast_arrays(s, t, p)

    if p.ndim > 1:
        m, n = p.shape
    else:
        m, n = p.size, 1

    svn = svan(s, t, p)
    mean_svan = 0.5 * (svn[1:m, ...] + svn[0:-1, ...])

    if n == 1:
        top = svn[0, 0] * p[0, 0] * db2Pascal
    else:
        top = svn[0, :] * p[0, :] * db2Pascal

    delta_ga = (mean_svan * np.diff(p, axis=0)) * db2Pascal
    ga = np.r_[np.atleast_2d(top), delta_ga]
    return np.cumsum(ga, axis=0)


def gvel(ga, lon, lat):
    """Calculates geostrophic velocity given the geopotential anomaly and
    position of each station.

    Parameters
    ----------
    ga : array_like
         geopotential anomaly relative to the sea surface.
    lon : array_like
          longitude of each station (+ve = E, -ve = W) [-180..+180]
    lat : array_like
          latitude  of each station (+ve = N, -ve = S) [ -90.. +90]

    Returns
    -------
    vel : array_like
           geostrophic velocity relative to the sea surface.
           Dimension will be MxN-1 (N: stations)

    Examples
    --------
    >>> import numpy as np
    >>> import seawater as sw
    >>> lon = np.array([-30, -30, -30, -30, -30, -30])
    >>> lat = np.linspace(-22, -21, 6)
    >>> t = np.array([[0,  0,  0,  0,  0,  0],
    ...               [10, 10, 10, 10, 10, 10],
    ...               [20, 20, 20, 20, 20, 20],
    ...               [30, 30, 30, 30, 30, 30],
    ...               [40, 40, 40, 40, 40, 40]])
    >>> s = np.array([[25, 25, 25, 35, 35, 35],
    ...               [25, 25, 25, 35, 35, 35],
    ...               [25, 25, 25, 35, 35, 35],
    ...               [25, 25, 25, 35, 35, 35],
    ...               [25, 25, 25, 35, 35, 35]])
    >>> p = np.array([[0, 5000, 10000, 0, 5000, 10000],
    ...               [0, 5000, 10000, 0, 5000, 10000],
    ...               [0, 5000, 10000, 0, 5000, 10000],
    ...               [0, 5000, 10000, 0, 5000, 10000],
    ...               [0, 5000, 10000, 0, 5000, 10000]])

    >>> ga = sw.gpan(s, t, p)
    >>> sw.gvel(ga, lon, lat)
    array([[-0.        , -0.        ],
           [ 0.11385677,  0.07154215],
           [ 0.22436555,  0.14112761],
           [ 0.33366412,  0.20996272]])

    References
    ----------
    .. [1] S. Pond & G.Pickard 2nd Edition 1986 Introductory Dynamical
    Oceanography Pergamon Press Sydney. ISBN 0-08-028728-X

    Modifications: 92-03-26. Phil Morgan.
    """

    ga, lon, lat = map(np.asanyarray, (ga, lon, lat))
    distm = dist(lon, lat, units='km')[0] * 1e3
    lf = f((lat[0:-1] + lat[1:]) / 2) * distm
    return -np.diff(ga, axis=1) / lf
