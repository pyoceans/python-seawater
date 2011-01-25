#TODO: study the possibility of transforming this into a DataSet like class with "properties" defining the different salinities.

import numpy as np
from seawater import constants as cte
from seawater import library as lib


""" cndr_from_SP == sw.cndr """
from  seawater.csiro import cndr as cndr_from_SP
""" SP_from_cndr == sw.salt """
from  seawater.csiro import salt as SP_from_cndr

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
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

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

    p = lib.check_dim(p, SP)
    lat = lib.check_dim(lat, SP)
    lon = lib.check_dim(lon, SP)

    SP[SP < 0] = 0
    lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(SP) # pythonic

    SA = np.NaN * np.zeros( SP.shape )
    dSA = np.NaN * np.zeros( SP.shape )
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

    After two iterations of a modified Newton-Raphson iteration, the error in SA is typically no larger than 2 :math:`^\times` 10 :sup:`-13` [g kg :sup:`-1`]

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

def SA_from_Sstar(Sstar, p, lon, lat):
    """
    Calculates Absolute Salinity from Preformed Salinity.

    Parameters
    ----------
    Sstar : array_like
            Preformed Salinity [g kg :sup:`-1`]
    p : array_like
        pressure [db]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    TODO

    Notes
    -----
    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometres inland from the coast.

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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    Sstar, p, lon, lat = np.asarray(Sstar), np.asarray(p), np.asarray(lon), np.asarray(lat)

    p = lib.check_dim(p, Sstar)
    lon = lib.check_dim(lon, Sstar)
    lat = lib.check_dim(lat, Sstar)

    lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(Sstar)
    SA = np.NaN * np.zeros( Sstar.shape )
    dSA = np.NaN * np.zeros( Sstar.shape )
    in_ocean = np.NaN * np.zeros( Sstar.shape ) #FIXME: change to boolean

    dSA[inds], in_ocean[inds] = lib._delta_SA( p[inds], lon[inds], lat[inds] )

    SA[inds] = Sstar[inds] + ( 1 + cte.r1 ) * dSA[inds]

    #NOTE: In the Baltic Sea, SA = Sstar, and note that gsw_delta_SA returns zero for dSA in the Baltic.

    return SA, in_ocean

def SP_from_SA(SA, p, lon, lat):
    """
    Calculates Practical Salinity from Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    p : array_like
        pressure [db]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    TODO

    Notes
    -----
    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometres inland from the coast.

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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, p, lon, lat = np.asarray(SA), np.asarray(p), np.asarray(lon), np.asarray(lat)

    p = lib.check_dim(p, SA)
    lat = lib.check_dim(lat, SA)
    lon = lib.check_dim(lon, SA)

    lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(SA)
    SP = np.NaN * np.zeros( SA.shape )
    dSA = np.NaN * np.zeros( SA.shape )
    in_ocean = np.NaN * np.zeros( SA.shape ) #FIXME: change to boolean

    dSA[inds], in_ocean[inds] = lib._delta_SA( p[inds], lon[inds], lat[inds] )

    SP[inds] = (35./cte.SSO) * ( SA[inds] - dSA[inds] )

    SP_baltic = lib._SP_from_SA_Baltic( SA, lon, lat )

    #SP_baltic[inds] = lib._SP_from_SA_Baltic(SA[inds], lon[inds], lat[inds] )

    indsbaltic = ~np.isnan(SP_baltic)

    SP[indsbaltic] = SP_baltic[indsbaltic]

    return SP, in_ocean

def Sstar_from_SA(SA, p, lon, lat):
    """
    Converts Preformed Salinity from Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    p : array_like
        pressure [db]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    Sstar : array_like
            Preformed Salinity [g kg :sup:`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    TODO

    Notes
    -----
    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometres inland from the coast.

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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, p, lon, lat = np.asarray(SA), np.asarray(p), np.asarray(lon), np.asarray(lat)

    p = lib.check_dim(p, SA)
    lon = lib.check_dim(lon, SA)
    lat = lib.check_dim(lat, SA)

    lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(SA)
    Sstar = np.NaN * np.zeros( SA.shape )
    dSA = np.NaN * np.zeros( SA.shape )
    in_ocean = np.NaN * np.zeros( SA.shape ) #FIXME: change to boolean

    dSA[inds], in_ocean[inds] = lib._delta_SA( p[inds], lon[inds], lat[inds] )

    Sstar[inds] =  SA[inds] - ( 1 + cte.r1 ) * dSA[inds]
    #NOTE: In the Baltic Sea, SA = Sstar, and note that gsw_delta_SA returns zero for dSA in the Baltic.

    return Sstar, in_ocean

def SP_from_Sstar(Sstar, p, lon, lat):
    """
    Calculates Practical Salinity from Preformed Salinity.

    Parameters
    ----------
    Sstar : array_like
            Preformed Salinity [g kg :sup:`-1`]
    p : array_like
        pressure [db]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    TODO

    Notes
    -----
    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometres inland from the coast.

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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    Sstar, p, lon, lat = np.asarray(Sstar), np.asarray(p), np.asarray(lon), np.asarray(lat)

    p = lib.check_dim(p, Sstar)
    lon = lib.check_dim(lon, Sstar)
    lat = lib.check_dim(lat, Sstar)

    lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(Sstar)
    SP = np.NaN * np.zeros( Sstar.shape )
    dSA = np.NaN * np.zeros( Sstar.shape )
    in_ocean = np.NaN * np.zeros( Sstar.shape ) #FIXME: change to boolean

    dSA[inds], in_ocean[inds] = lib._delta_SA( p[inds], lon[inds], lat[inds] )
    SP[inds] = (35/cte.SSO) * ( Sstar[inds] + cte.r1 * dSA[inds] )

    # In the Baltic Sea, SA = Sstar.
    SP_baltic = lib._SP_from_SA_Baltic( Sstar, lon, lat )

    indsbaltic = ~np.isnan(SP_baltic)

    SP[indsbaltic] = SP_baltic[indsbaltic]

    return SP, in_ocean

def Sstar_from_SP(SP, p, lon, lat):
    """
    Calculates Preformed Salinity from Absolute Salinity.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    p : array_like
        pressure [db]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    Sstar : array_like
            Preformed Salinity [g kg :sup:`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    TODO

    Notes
    -----
    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometres inland from the coast.

    Since SP is non-negative by definition, this function changes any negative input values of SP to be zero.

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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    ,, [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SP, p, lon, lat = np.asarray(SP), np.asarray(p), np.asarray(lon), np.asarray(lat)

    p = lib.check_dim(p, SP)
    lon = lib.check_dim(lon, SP)
    lat = lib.check_dim(lat, SP)

    SP[SP < 0] = 0
    lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(SP)
    Sstar = np.NaN * np.zeros( SP.shape )
    dSA = np.NaN * np.zeros( SP.shape )
    in_ocean = np.NaN * np.zeros( SP.shape ) #FIXME: change to boolean

    dSA[inds], in_ocean[inds] = lib._delta_SA( p[inds], lon[inds], lat[inds] )
    Sstar[inds] = (cte.SSO/35.) * SP[inds] - cte.r1 * dSA[inds]

    # In the Baltic Sea, SA = Sstar.
    Sstar_baltic = lib._SA_from_SP_Baltic( SP, lon, lat )

    indsbaltic = ~np.isnan(Sstar_baltic)

    Sstar[indsbaltic] = Sstar_baltic[indsbaltic]

    return Sstar, in_ocean

def SA_Sstar_from_SP(SP, p, lon, lat):
    """
    Calculates Absolute Salinity and Preformed Salinity from Practical Salinity.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    p : array_like
        pressure [db]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    Sstar : array_like
            Preformed Salinity [g kg :sup:`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    TODO

    Notes
    -----
    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometres inland from the coast.

    Since SP is non-negative by definition, this function changes any negative input values of SP to be zero.

    Examples
    --------
    >>> import seawater.Salinity as gsw
    >>> SP = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
    >>> p = [0., 500., 1500., 2000.]
    >>> lon, lat = -69., 42.
    >>> gsw.SA_from_SP(SP, p, lon, lat)[0]
    TODO


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

    p = lib.check_dim(p, SP)
    lat = lib.check_dim(lat, SP)
    lon = lib.check_dim(lon, SP)

    SP[SP < 0] = 0
    lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(SP)

    SA = np.NaN * np.zeros( SP.shape )
    Sstar = np.NaN * np.zeros( SP.shape )
    dSA = np.NaN * np.zeros( SP.shape )
    in_ocean = np.NaN * np.zeros( SP.shape ) #FIXME: change to boolean

    dSA[inds], in_ocean[inds] = lib._delta_SA( p[inds], lon[inds], lat[inds] )

    SA[inds] = ( cte.SSO / 35 ) * SP[inds] + dSA[inds]
    Sstar[inds] = ( cte.SSO / 35 ) * SP[inds] - cte.r1 * dSA[inds]

    SA_baltic = lib._SA_from_SP_Baltic( SP, lon, lat )

    indsbaltic = ~np.isnan(SA_baltic)

    SA[indsbaltic] = SA_baltic[indsbaltic]
    Sstar[indsbaltic] = SA_baltic[indsbaltic] # In the Baltic Sea, Sstar = SA.

    return SA, Sstar, in_ocean

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
        try:
            exec( "unequal = (gsw_cv." +comp_value+ " - " +method+ " ) >= gsw_cv." +comp_value+ "_ca")
        except:
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

    """ SA_from_SP """ #TODO: test the second output (in_ocean)
    SA_from_SP  = gsw.SA_from_SP(gsw_cv.SP_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
    test_print("SA_from_SP")

    """ SP_from_SA """ #TODO: test the second output (in_ocean)
    SA_chck_cast = sio.loadmat("derived_prop.mat", squeeze_me=True)['SA_chck_cast']
    SP_from_SA = gsw.SP_from_SA(SA_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
    test_print("SP_from_SA") #NOTE: the original make the comparison with SP_chck_cast, why?

    """ SA_from_rho """
    rho_chck_cast = sio.loadmat("derived_prop.mat", squeeze_me=True)['rho']
    SA_from_rho = gsw.SA_from_rho(rho_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
    test_print("SA_from_rho")

    """ SA_from_Sstar """
    Sstar_from_SA = gsw.Sstar_from_SA(SA_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
    test_print("Sstar_from_SA")

    """ Sstar_from_SA """
    #FIXME: why SA_chck_cast instead of Sstar?
    SA_from_Sstar = gsw.SA_from_Sstar(SA_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
    test_print("SA_from_Sstar")

    """ SP_from_Sstar """
    #FIXME: original run SP_from_SA instead of SP_from_Sstar
    SP_from_Sstar = gsw.SP_from_SA(SA_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
    test_print("SP_from_Sstar")

    """ Sstar_from_SP """
    Sstar_from_SP = gsw.Sstar_from_SP(gsw_cv.SP_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
    test_print("Sstar_from_SP")

    """ cndr_from_SP """ #NOTE: show diff not present in the original
    cndr = gsw.cndr_from_SP(gsw_cv.SP_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
    test_print("cndr")

    """ SP_from_cndr """
    cndr = sio.loadmat("derived_prop.mat", squeeze_me=True)['cndr']
    SP_from_cndr = gsw.SP_from_cndr(cndr, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
    test_print("SP_from_cndr")

    """ SA_Sstar_from_SP """
    SA_SA_Sstar_from_SP, Sstar_SA_Sstar_from_SP = gsw.SA_Sstar_from_SP(gsw_cv.SP_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0:2]
    test_print("SA_SA_Sstar_from_SP")
    test_print("Sstar_SA_Sstar_from_SP")