#TODO: Go over PDFs to improve documentation

import numpy as np
from seawater import constants as cte
from seawater import library as lib

def  z_from_p(p, lat):
    """
    Calculates height from sea pressure using the computationally-efficient 25-term expression for density in terms of SA, CT and p.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90].
    p : array_like
        pressure [db].

    Returns
    -------
    z : array_like
        height [m]

    See Also
    --------
    TODO


    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> p = [0., 15., 100., 550., 1500., 2000., 3000., 5000., 10000.]
    >>> lat = 32.
    >>> gsw.z_from_p(p, lat)
    array([   -0.        ,   -14.89499448,   -99.27948265,  -545.44412444,
           -1484.209721  , -1976.61994868, -2958.05761312, -4907.87772419,
           -9712.16369644])
    >>> lat = [0., 15., 20., 35., 42., 63., 77., 85., 90.]
    >>> gsw.z_from_p(p, lat)
    array([   -0.        ,   -14.9118282 ,   -99.36544813,  -545.30528098,
           -1482.90095076, -1971.26442738, -2947.61650675, -4889.44474273,
           -9675.31755921])

    Notes
    -----
    At sea level z = 0, and since z (HEIGHT) is defined to be positive upwards, it follows that while z is positive in the atmosphere, it is NEGATIVE in the ocean.

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. Available from http://www.TEOS-10.org
    ,, [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater. FIXME: To be submitted to Ocean Science Discussions.
    .. [3] Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.

    Modifications:
    2010-08-26. Trevor McDougall, Claire Roberts-Thomson & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    p, lat = np.asarray(p), np.asarray(lat)

    X     = np.sin( np.deg2rad(lat) )
    sin2  = X**2
    B     = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2 )
    A     = -0.5 * cte.gamma * B
    C     = lib._enthalpy_SSO_0_CT25(p)
    z     = -2 * C / ( B + np.sqrt( B**2 - 4 * A * C ) )

    return z

def grav(lat, p=0):
    """
    Calculates acceleration due to gravity as a function of latitude and as a function of pressure in the ocean.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90].
    p : number or array_like. Default p = 0
        pressure [db]

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
    >>> lat = [0., 15., 20., 35., 42., 63., 77., 85., 90.]
    >>> gsw.grav(lat)
    array([ 9.780327  ,  9.78378673,  9.78636994,  9.79733807,  9.80349012,
            9.82146051,  9.82955107,  9.83179057,  9.83218621])
    >>> p = [0., 15., 100., 550., 1500., 2000., 3000., 5000., 10000.]
    >>> gsw.grav(lat, p)
    array([ 9.780327  ,  9.7838197 ,  9.78658971,  9.79854548,  9.80677561,
            9.82583603,  9.83609914,  9.84265484,  9.85368548])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. Available from http://www.TEOS-10.org. Appendix D of this TEOS-10 Manual.
    .. [2] Moritz (2000) Goedetic reference system 1980. J. Geodesy,74,128-133.
    .. [3] Saunders, P.M., and N.P. Fofonoff (1976) Conversion of pressure to depth in the ocean. Deep-Sea Res.,pp. 109 - 111.

    Modifications:
    2010-07-23. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    lat, p, = np.asarray(lat), np.asarray(p)

    X = np.sin( np.deg2rad(lat) )
    sin2 = X**2
    gs = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2)
    z = z_from_p(p, lat)
    grav = gs * (1 - cte.gamma * z) # z is the height corresponding to p
    return grav

def molality(SA):
    """
    Calculates the molality of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]

    Returns
    -------
    molality : array_like
        seawater molality [mol kg :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
    >>> gsw.molality(SA)
    array([[ 1.78214644,  0.98484303,  0.32164907,  0.64986241],
           [ 0.32164907,         nan,  0.48492271,  0.25680047]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. Available from http://www.TEOS-10.org. Appendix D of this TEOS-10 Manual.

    Modifications:
    2010-09-28. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA = np.asarray(SA)

    Isalty = (SA >= 0).nonzero()
    molality = np.ones( SA.shape )*np.nan
    # molality of seawater in mol kg :sup:`-1`
    molality[Isalty] = SA[Isalty] / (cte.M_S * ( 1000 - SA[Isalty] ) )

    return molality

def ionic_strength(SA):
    """
    Calculates the ionic strength of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]

    Returns
    -------
    ionic_strength : array_like
        ionic strength of seawater [mol kg :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
    >>> gsw.ionic_strength(SA)
    array([[ 1.10964439,  0.61320749,  0.20027315,  0.40463351],
           [ 0.20027315,         nan,  0.30193465,  0.1598955 ]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
    see Table L.1 of this TEOS-10 Manual.

    .. [2] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: The composition of Standard Seawater and the definition of the Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72. see Eqns. (5.9) and (5.12) of this paper.

    Modifications:
    2010-09-28. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA = np.asarray(SA)

    Z_2 = 1.2452898 # the valence factor of sea salt


    molal = molality(SA) # molality of seawater in mol kg

    ionic_strength = 0.5*Z_2*molal

    return ionic_strength