#TODO: study the possibility of transforming this into a DataSet like class with "properties" defining the different salinities.

import numpy as np
from seawater import constants as cte
from seawater import library as lib

def check_dim(prop1, prop2):
    """
    Broadcast prop1 to the shape prop2. Prop1 can be scalar, row equal or column equal to prop2.
    TODO: Needs lots of improvements and cleanups...
    """
    if prop1.ndim == 1:
        prop1 = prop1.flatten()

    if (prop1.ndim == 1) & (prop1.size == 1):
        prop1 = prop1 * np.ones( prop2.shape )
    elif (prop1.ndim == 1) & (prop2.ndim != 1):
        if prop1.size == prop2.shape[1]:
            prop1 = prop1 * np.ones(prop2.shape)
            #prop1 = np.repeat(prop1[np.newaxis,:], prop2.shape[1], axis=1).reshape(prop2.shape)
        elif prop1.size == prop2.shape[0]:
            prop1 = prop1[:,np.newaxis] * np.ones(prop2.shape)
            #prop1 = np.repeat(prop1[np.newaxis,:], prop2.shape[0], axis=0).reshape(prop2.shape)
        else:
            raise NameError('Blahrg')

    if prop1.ndim == 0:
        prop1 = prop1 * np.ones(prop2.shape)

    return prop1

def SA_from_SP(SP, p, lon, lat):
    """
    Calculates Absolute Salinity from Practical Salinity.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    p : array_like
        pressure [db]
    lon : array_like
          decimal degrees east [0..+360]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    Since SP is non-negative by definition, this function changes any negative input values of SP to be zero.

    Examples
    --------
    >>> import seawater.Salinity as gsw
    >>> SP = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
    >>> p = [0., 500., 1500., 2000.]
    >>> lon, lat = -69., 42.
    >>> gsw.SA_from_SP(SP, p, lon, lat)[0]
    array([[  5.32510274e+01,   3.01448066e+01,   1.00503768e+01,
              2.00980483e+01],
           [  1.00482640e+01,   3.34377888e-03,   1.50739539e+01,
              8.04146315e+00]])


    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SP, p, lon, lat = np.asarray(SP), np.asarray(p), np.asarray(lon), np.asarray(lat)

    p = check_dim(p, SP)
    lat = check_dim(lat, SP)
    lon = check_dim(lon, SP)

    SP[SP < 0] = 0
    lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(SP) # pythonic

    SA = np.NaN * np.zeros( SP.shape )
    dSA = np.NaN * np.zeros( SP.shape )
    SA_baltic = np.NaN * np.zeros( SP.shape )
    in_ocean = np.NaN * np.zeros( SP.shape ) #FIXME: change to boolean

    dSA[inds], in_ocean[inds] = lib._delta_SA( p[inds], lon[inds], lat[inds] )
    SA[inds] = ( cte.SSO / 35 ) * SP[inds] + dSA[inds]

    SA_baltic = lib._SA_from_SP_Baltic( SP, lon, lat )

    indsbaltic = ~np.isnan(SA_baltic)

    SA[indsbaltic] = SA_baltic[indsbaltic]

    return SA, in_ocean

def SA_from_rho(rho, t, p):
    """
    Calculates the Absolute Salinity of a seawater sample, for given values of its density, in-situ temperature and sea pressure (in dbar).

    Parameters
    ----------
    rho : array_like
          in-situ density [kg m :sup:`-3`]
    t : array_like
        in-situ temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]

    Returns
    -------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    This is expressed on the Reference-Composition Salinity Scale of Millero et al. (2008).
    # After two iterations of a modified Newton-Raphson iteration, the error in SA is typically no larger than 2 :math:`^\times` 10 :sup:`-13` [g kg :sup:`-1`]

    Examples
    --------
    >>> import seawater.Salinity as gsw
    >>> rho = [[1041.77425464, 1024.2413978, 1011.923534, 1018.28328036],[1006.74841976, 1002.37206267, 1014.78353156, 1010.8696052]]
    >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    >>>
    >>> p = [0., 500., 1500., 2000.]
    >>> gsw.SA_from_rho(rho, t, p)
    FIXME: [[53., 30, 10., 20.],[10., -5., 15., 8.]]

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.5.

    .. [2] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: The composition of Standard Seawater and the definition of the Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.

    Modifications:
    2010-08-23. Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    rho, t, p = np.asarray(rho), np.asarray(t), np.asarray(p)

    n0, n1 = 0, 1
    v_lab = np.ones( rho.shape ) / rho
    v_0 = lib._gibbs(n0, n0, n1, 0, t, p)
    v_120 = lib._gibbs(n0, n0, n1, 120, t, p)
    SA = 120 * ( v_lab - v_0 ) / ( v_120 - v_0 ) # initial estimate of SA
    Ior = (SA < 0) | (SA > 120)
    v_SA = ( v_120 - v_0 ) / 120 # initial estimate of v_SA, the SA derivative of v

    for iter in range(0,2):
        SA_old = SA
        delta_v = lib._gibbs(n0, n0, n1, SA_old, t, p) - v_lab
        SA = SA_old - delta_v / v_SA # this is half way through the modified N-R method
        SA_mean = 0.5 * ( SA + SA_old )
        v_SA = lib._gibbs(n1, n0, n1, SA_mean, t, p)
        SA = SA_old - delta_v / v_SA

    SA[Ior] = np.NaN

    return SA

if __name__=='__main__':
    try:
        import cPickle as pickle
    except:
        import pickle
    import numpy as np
    import seawater.Salinity as gsw
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
        exec( "unequal = (gsw_cv." +comp_value+ " - " +method+ " ) >= gsw_cv." +comp_value+ "_ca")

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

    """ derived values (for comparison/test) """
    rho_chck_cast = sio.loadmat("derived_prop.mat", squeeze_me=True)['rho']

    """ SA_from_SP """ #TODO: test the second output (in funnel)
    SA_from_SP  = gsw.SA_from_SP(gsw_cv.SP_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
    test_print("SA_from_SP")

    """ SA_from_rho """
    SA_from_rho = gsw.SA_from_rho(rho_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
    test_print("SA_from_rho")

