# -*- coding: utf-8 -*-

"""
Planet Earth properties

Functions:
----------
  f
      Coriolis parameter (f) 
  grav                            
      Gravitational acceleration
  distance                        
      spherical earth distance between points
      in lon, lat coordinates at a given pressure

This is part of the python Gibbs Sea Water library
http://code.google.com/p/python-seawater/

"""

from __future__ import division

import numpy as np
import seawater.constants as cte
from library import match_args_return
from conversions import z_from_p

# -----------

__all__ = ['f', 
           'grav', 
           'distance']

# -----------

def f(lat):
    """Coriolis parameter

    parameter
    ---------
    lat : array_like, latitude      [degrees north]

    returns
    -------
    f : array_like, Coriolis paramter  [1/s]

    """

    lat = np.asanyarray(lat)
    rad = np.pi / 180.0
    omega = 7.292115e-5
    return 2 * omega * np.sin(lat*rad)


@match_args_return
def grav(lat, p=0):
    r"""
    Calculates acceleration due to gravity as a function of latitude and as a
    function of pressure in the ocean.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    p : number or array_like. Default p = 0
        pressure [dbar]

    Returns
    -------
    g : array_like
        gravity [m s :sup:`2`]

    See Also
    --------
    TODO

    Notes
    -----
    In the ocean z is negative.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> lat = [-90, -60, -30, 0]
    >>> p = 0
    >>> gsw.grav(lat, p)
    array([ 9.83218621,  9.81917886,  9.79324926,  9.780327  ])
    >>> gsw.grav(45)
    9.8061998770458008

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    .. [2] Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.

    .. [3] Saunders, P.M., and N.P. Fofonoff (1976) Conversion of pressure to
    depth in the ocean. Deep-Sea Res.,pp. 109 - 111.

    Modifications:
    2010-07-23. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    X = np.sin( np.deg2rad(lat) )
    sin2 = X**2
    gs = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2)
    z = z_from_p(p, lat)
    grav = gs * (1 - cte.gamma * z) # z is the height corresponding to p
    return grav

@match_args_return
def distance(lon, lat, p=0):
    r"""
    Calculates the distance in met res between successive points in the vectors
    lon and lat, computed using the Haversine formula on a spherical earth of
    radius 6,371 km, being the radius of a sphere having the same volume as
    Earth. For a spherical Earth of radius 6,371,000 m, one nautical mile is
    1,853.2488 m, thus one degree of latitude is 111,194.93 m.

    Haversine formula:
        R = earth's radius (mean radius = 6,371 km)

    .. math::
        a = \sin^2(\delta \text{lat}/2) + \cos(\text{lat}_1) \cos(\text{lat}_2) \sin^2(\delta \text{lon}/2)

        c = 2 \times \text{atan2}(\sqrt{a}, \sqrt{(1-a)})

        d = R \times c

    Parameters
    ----------
    lon : array_like
          decimal degrees east [0..+360] or [-180 ... +180]
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    p : number or array_like. Default p = 0
        pressure [dbar]

    Returns
    -------
    dist: array_like
          distance between points on a spherical Earth at pressure (p) [m]

    See Also
    --------
    TODO

    Notes
    -----
    z is height and is negative in the oceanographic.

    Distances are probably good to better than 1\% of the "true" distance on the
    ellipsoidal earth.

    The check value below differ from the original online docs at
    "http://www.teos-10.org/pubs/gsw/html/gsw_distance.html" but agree with the
    result.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> lon = [159, 220]
    >>> lat = [-35, 35]
    >>> gsw.distance(lon, lat)
    array([[ 10030974.652916]])
    >>> p = [200, 1000]
    >>> gsw.distance(lon, lat, p)
    array([[ 10030661.63878009]])
    >>> p = [[200], [1000]]
    >>> gsw.distance(lon, lat, p)
    array([[ 10030661.63878009],
           [ 10029412.58776001]])

    References
    ----------
    .. [1] http://www.eos.ubc.ca/~rich/map.html

    Modifications:
    2000-11-06. Rich Pawlowicz
    2010-07-28. Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """
    #FIXME? The argument handling seems much too complicated.
    # Maybe we can come up with some simple specifications of
    # what argument combinations are permitted, and handle everything
    # with broadcasting. - EF

    #FIXME: Eric what do you think? This assume p(stations, depth)
    lon, lat, = np.atleast_2d(lon), np.atleast_2d(lat)

    if (lon.size == 1) & (lat.size == 1):
        raise ValueError('more than one point is needed to compute distance')
    elif lon.ndim != lat.ndim:
        raise ValueError('lon, lat must have the same dimension')

    lon, lat, p = np.broadcast_arrays(lon, lat, p)

    dlon = np.deg2rad( np.diff(lon) )
    dlat = np.deg2rad( np.diff(lat) )

    a = ( ( np.sin(dlat/2.) )**2 + np.cos( np.deg2rad( lat[:,:-1] ) ) *
    np.cos( np.deg2rad( lat[:,1:] ) ) * ( np.sin(dlon/2.) )**2 )

    angles = 2. * np.arctan2( np.ma.sqrt(a), np.ma.sqrt(1-a) )

    p_mid = 0.5 * (   p[:,0:-1] +   p[:,0:-1] )
    lat_mid = 0.5 * ( lat[:,:-1] + lat[:,1:] )

    z = z_from_p(p_mid, lat_mid)

    distance = (cte.a + z) * angles

    return distance



