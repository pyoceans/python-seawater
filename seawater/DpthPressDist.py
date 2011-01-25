import numpy as np
from seawater import constants as cte
from seawater import library as lib

def  z_from_p(p, lat):
    """
    Calculates height from sea pressure using the computationally-efficient 25-term expression for density in terms of SA, CT and p.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    p : array_like
        pressure [db]

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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

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

def  p_from_z(z, lat): # Also in DpthPressDist.py
    """
    Calculates sea pressure from height using computationally-efficient 25-term expression for density, in terms of SA, CT and p.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    z : array_like
        height [m]

    Returns
    -------
    p : array_like
        pressure [db]

    See Also
    --------
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> z = [-0., -14.89499448, -99.27948265, -545.44412444, -1484.209721, -1976.61994868, -2958.05761312, -4907.87772419, -9712.16369644]
    >>> lat = 32.
    >>> gsw.p_from_z(z, lat)
    array([     0.,     15.,    100.,    550.,   1500.,   2000.,   3000.,
             5000.,  10000.])

    Notes
    -----
    Height (z) is NEGATIVE in the ocean.  Depth is -z. Depth is not used in the gibbs library.

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    ,, [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater. FIXME: To be submitted to Ocean Science Discussions.

    .. [3] Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.

    .. [4] Saunders, P. M., 1981: Practical conversion of pressure to depth. Journal of Physical Oceanography, 11, 573-574.

    Modifications:
    2010-08-26. Trevor McDougall, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    z, lat = np.asarray(z), np.asarray(lat)

    X     = np.sin( np.deg2rad(lat) )
    sin2  = X**2
    gs    = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2 )
    # get the first estimate of p from Saunders (1981)
    c1 =  5.25e-3 * sin2 + 5.92e-3
    p  = -2 * z / ( (1-c1) + np.sqrt( (1-c1) * (1-c1) + 8.84e-6 * z ) )
    # end of the first estimate from Saunders (1981)
    df_dp = cte.db2Pascal * lib._specvol_SSO_0_CT25(p) # initial value of the derivative of f
    f     = lib._enthalpy_SSO_0_CT25(p) + gs * ( z - 0.5 * cte.gamma * ( z**2 ) )
    p_old = p
    p     = p_old - f / df_dp
    pm    = 0.5 * (p + p_old)
    df_dp = cte.db2Pascal * lib._specvol_SSO_0_CT25(pm)
    p     = p_old - f / df_dp

    return p

def grav(lat, p=0):
    """
    Calculates acceleration due to gravity as a function of latitude and as a function of pressure in the ocean.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90]
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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

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

def distance(lon, lat, p=None):
    """
    Calculates the distance in metres between successive points in the vectors lon and lat, computed using the Haversine formula on a spherical earth of radius 6,371 km, being the radius of a sphere having the same volume as Earth. For a sperical Earth of radius 6,371,000 m, one nautical mile is 1,853.2488 m, thus one degree of latitude is 111,194.93 m.

    Parameters
    ----------
    lon : array_like
          decimal degrees east [0..+360] or [-180 ... +180]
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    p : number or array_like. Default p = 0
        pressure [db]

    Returns
    -------
    dist: array_like
          distance between points on a spherical Earth at pressure (p) [m]

    See Also
    --------
    TODO

    Notes
    -----
    Distances are probably good to better than 1% of the "true" distance on the ellipsoidal earth.

    Examples
    --------
    >>> import seawater.DpthPressDist as gsw
    lon = [142, 183, 20]
    lat = [11, 9.5, 59]
    p = [[0., 0., 0.,],[500., 500., 500.,]]
    gsw.distance(lon, lat, p)
    TODO

    References
    ----------
    .. [1] http://www.eos.ubc.ca/~rich/map.html

    Modifications:
    2000-11-06. Rich Pawlowicz
    2010-07-28. Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    #NOTE: dtype=np.float64 is need because the test data is uint8
    lon, lat = np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64)

    if (lon.size == 1) & (lat.size == 1):
        raise NameError('more than one point is needed to compute distance')
    elif lon.ndim != lat.ndim:
        raise NameError('lon, lat must have the same dimension')

    if (lon.size == 1) & (lat.size != 1): # fill lon with lat size
        lon = lib.check_dim(lon, lat)
    elif (lat.size == 1) & (lon.size != 1): # fill lat with lon size
        lat = lib.check_dim(lat, lon)

    if (lon.ndim == 1) & (lat.ndim == 1): #NOTE: Ugly way to matlabsize it
        lon = lon[np.newaxis,:]
        lat = lat[np.newaxis,:]

    # check for lon/lat and p dimensions
    if p == None:
        p = np.zeros( lon.shape )
    else:
        p = np.asarray(p)

    if p.ndim > lat.ndim:
        lon = lib.check_dim(lon, p)
        lat = lib.check_dim(lat, p)
    elif p.ndim == 1:
        p = lib.check_dim(p, lon)

    dlon = np.deg2rad( np.diff(lon, axis=1) )
    dlat = np.deg2rad( np.diff(lat, axis=1) )

    a = ( np.sin(dlat/2.) )**2 + np.cos( np.deg2rad( lat[:,:-1] ) ) * np.cos( np.deg2rad( lat[:,1:] ) ) * ( np.sin(dlon/2.) )**2

    angles = 2. * np.arctan2( np.sqrt(a), np.sqrt(1-a) )

    p_mid = 0.5 * (   p[:,0:-1] +   p[:,0:-1] )
    lat_mid = 0.5 * ( lat[:,:-1] + lat[:,1:] )

    z = z_from_p(p_mid, lat_mid) #NOTE: z is height and is negative in th ocean

    distance = (cte.a + z) * angles

    return distance

if __name__=='__main__':
    try:
        import cPickle as pickle
    except:
        import pickle
    import numpy as np
    import scipy.io as sio

    """ load test data """
    class Dict2Struc(object):
        """ all the variables from a dict in a "matlab-like-structure" """
        def __init__(self, adict):
            self.__dict__.update(adict)

    data = pickle.load( open('gsw_cv.pkl','rb') )
    gsw_cv = Dict2Struc(data) # then type dat.<tab> to navigate through your variables

    def test_print(method, comp_value=None):
        """
        Run a function test mimicking the original logic. This is done to allow for a direct comparison of the result from the Matlab to the python package.
        """

        if comp_value is None:
            comp_value = method

        # test for floating differences with: computed - check_value >= defined_precision
        try:
            exec( "unequal = (gsw_cv." +comp_value+ " - " +method+ " ) >= gsw_cv." +comp_value+ "_ca")
        except:
            print "unequal = (gsw_cv." +comp_value+ " - " +method+ " ) >= gsw_cv." +method+ "_ca"
            exec( "unequal = (gsw_cv." +comp_value+ " - " +method+ " ) >= gsw_cv." +method+ "_ca")

        width = 23
        if unequal.any():
            print "%s: Failed" % method.rjust(width)
        else:
            # test if check value is identical to computed value
            if eval( "( gsw_cv." +comp_value+ "[~np.isnan(gsw_cv."+comp_value+")] == " +method+ "[~np.isnan("+method+")] ).all()" ):
                print "%s: Passed" % method.rjust(width)
            else:
                # test for differences in case their aren't equal. This is an attempt to place all tests together (i.e. term25 and small float differences that will appear)
                exec("nmax = ( gsw_cv."+comp_value+" - "+method+" )[~np.isnan(gsw_cv."+comp_value+")].max()")
                exec("nmin = ( gsw_cv."+comp_value+" - "+method+" )[~np.isnan(gsw_cv."+comp_value+")].min()")
                print "%s: Passed, but small diff ranging from: %s to %s" % ( method.rjust(width), nmax, nmin)

    import seawater.DpthPressDist as gsw

    """ z_from_p """
    z_from_p = gsw.z_from_p(gsw_cv.p_chck_cast, gsw_cv.lat_chck_cast)
    test_print("z_from_p")

    """ p_from_z """ #NOTE: show diff not present in the original
    p_from_z = gsw.p_from_z( z_from_p, gsw_cv.lat_chck_cast )
    test_print("p_from_z")

    """ grav """
    grav = gsw.grav(gsw_cv.lat_chck_cast, gsw_cv.p_chck_cast )
    test_print("grav")

    """ distance """
    distance = gsw.distance(gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast, gsw_cv.p_chck_cast)
    test_print("distance")

#class DpthPressDist: #TODO: an attempt to create a class for Depth pressure conversions
    #"""
    #Class to convert pressure to depth and vice-versa:

    #This class will (maybe) replace z_from_p, p_from_z, and grav functions
    #TODO: Copy documentation
    #"""
    #def __init__(self, lat, p=None, z=None):
        #self.lat, self.p, self.z = np.asarray(lat), np.asarray(p), np.asarray(z)

        #X     = np.sin( np.deg2rad(self.lat) )
        #sin2   = X**2
        #gs = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2)

        #if p:
            #A      = -0.5 * cte.gamma * gs
            #C      = lib._enthalpy_SSO_0_CT25(self.p)
            #self.z = -2 * C / ( gs + np.sqrt( gs**2 - 4 * A * C ) )
        #elif z:
            ## get the first estimate of p from Saunders (1981)
            #c1 =  5.25e-3 * sin2 + 5.92e-3
            #p  = -2 * self.z / ( ( 1-c1 ) + np.sqrt( ( 1-c1 ) * ( 1-c1 ) + 8.84e-6 * self.z ) )
            ## end of the first estimate from Saunders (1981)
            ## initial value of the derivative of f
            #df_dp = cte.db2Pascal * lib._specvol_SSO_0_CT25(p)
            #f     = lib._enthalpy_SSO_0_CT25(p) + gs * ( self.z - 0.5 * cte.gamma * ( self.z**2 ) )
            #p_old = p
            #p     = p_old - f / df_dp
            #pm    = 0.5 * ( p + p_old )
            #df_dp = cte.db2Pascal * lib._specvol_SSO_0_CT25(pm)
            #self.p     = p_old - f / df_dp
        #else:
            #raise NameError('need latitude (mandatory) and pressure or depth')

        #self.grav = gs * (1 - cte.gamma * self.z)


#if __name__=='__main__':
    #lat = 32.
    #p = [0., 15., 100., 550., 1500., 2000., 3000., 5000., 10000.]
    #dp = DpthPressDist(lat=lat, p=p)
    #z = [-0., -14.89499448, -99.27948265, -545.44412444, -1484.209721, -1976.61994868, -2958.05761312, -4907.87772419, -9712.16369644]
    #pd = DpthPressDist(lat=lat, z=z)