#TODO: Go over PDFs to improve documentation
#TODO: Examples: simple with the data range (copy-and-paste numbers) and complex (real data)
#TODO: Check original authors and dates
#TODO: csiro vs gibbs (table?)
#TODO: check_dim for p in all "p" functions
#TODO: def abc(): in alphabetical order
#FIXME: some function return values even with nan in the input, check this behaviour (also present in the original).

import numpy as np
from seawater import constants as cte
from seawater import library as lib

def check_dim(prop1, prop2):
    """
    Broadcast prop1 to the shape prop2. Prop1 can be scalar, row equal or column equal to prop2.
    TODO: Needs lots of improvement and cleanups...
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

    prop1.dtype = 'float64' # FIXME: somehow lon is been loaded as uint8
    return prop1

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

def entropy_from_t(SA, t, t_type='pt'):
    """
    Calculates specific entropy of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
        temperature [:math:`^\\circ` C]
    t_type : str, optional
             'pt' for potential temperature [:math:`^\\circ` C (ITS-90)]
             'CT' for Conservative Temperature [:math:`^\\circ` C (TEOS-10)] FIXME: the orginal has (ITS-90) Copy and paste issue?

    Returns
    -------
    entropy : array_like
              specific entropy [J kg :sup:`-1` K :sup:`-1`]

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
    >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    >>> gsw.entropy_from_t(SA, t)
    array([[  6.36913727e+01,   2.16184242e+02,   3.23536391e+02,
              4.54589274e+02],
           [  2.24455426e+02,  -1.47644587e-01,   3.62859557e+02,
              4.07493891e+02]])
    >>> gsw.entropy_from_t(SA, t, 'CT')
    array([[  6.71618495e+01,   2.14564404e+02,   3.12442303e+02,
              4.45351241e+02],
           [  2.16404703e+02,  -3.71020898e-01,   3.52932205e+02,
              3.92869342e+02]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See appendix A.10.

    Modifications:
    2010-10-13. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, t = np.asarray(SA), np.asarray(t)

    SA[SA < 0] = 0

    n0 = 0
    n1 = 1
    pr0 = np.zeros( SA.shape )

    if t_type == 'pt':
        entropy = -lib._gibbs(n0, n1, n0, SA, t, pr0)
    elif t_type == 'CT':
        pt0 = pt_from_CT(SA, t)
        entropy = -lib._gibbs(n0, n1, n0, SA, pt0, pr0)

    return entropy

def pt_from_CT(SA, CT):
    """
    Calculates potential temperature (with a reference sea pressure of zero dbar) from Conservative Temperature.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\\circ` C (TEOS-10)] FIXME: the orginal has (ITS-90) Copy and paste issue?

    Returns
    -------
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar [:math:`^\\circ` C (ITS-90)]

    See Also
    --------
    TODO

    Notes
    -----
    This function uses 1.5 iterations through a modified Newton-Raphson (N-R) iterative solution proceedure, starting from a rational-function-based initial condition for both pt and dCT_dpt.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
    >>> CT = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    >>> gsw.pt_from_CT(SA, CT)
    array([[  5.24810333e+00,   1.48839096e+01,   2.12076507e+01,
              3.13091326e+01],
           [  1.44387776e+01,  -1.44601365e-02,   2.42789834e+01,
              2.69372570e+01]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See sections 3.1 and 3.3.

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater.  To be submitted to Ocean Science Discussions.

    Modifications:
    2010-08-26. McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, CT = np.asarray(SA), np.asarray(CT)

    SA[SA < 0] = 0

    s1 = SA * 35. / 35.16504

    a0 = -1.446013646344788e-2
    a1 = -3.305308995852924e-3
    a2 =  1.062415929128982e-4
    a3 =  9.477566673794488e-1
    a4 =  2.166591947736613e-3
    a5 =  3.828842955039902e-3

    b0 =  1.000000000000000e+0
    b1 =  6.506097115635800e-4
    b2 =  3.830289486850898e-3
    b3 =  1.247811760368034e-6

    a5CT = a5 * CT
    b3CT = b3 * CT
    CT_factor = ( a3 + a4 * s1 + a5CT )
    pt_num = a0 + s1 * ( a1 + a2 * s1 ) + CT * CT_factor
    pt_den = b0 + b1 * s1 + CT * ( b2 + b3CT )
    pt = pt_num / pt_den

    dCT_dpt = pt_den / (CT_factor + a5CT - (b2 + b3CT + b3CT) * pt)

    # start the 1.5 iterations through the modified Newton-Rapshon iterative method
    CT_diff = CT_from_pt(SA, pt) - CT
    pt_old = pt
    pt = pt_old - CT_diff /dCT_dpt # 1/2-way through the 1st modified N-R loop
    ptm = 0.5 * (pt + pt_old)

    #This routine calls gibbs_pt0_pt0(SA,pt0) to get the second derivative of the Gibbs function with respect to temperature at zero sea pressure.

    dCT_dpt = -(ptm + 273.15) * lib._gibbs_pt0_pt0(SA, ptm) / cte.cp0
    pt = pt_old - CT_diff / dCT_dpt # end of 1st full modified N-R iteration
    CT_diff = CT_from_pt(SA, pt) - CT
    pt_old = pt
    pt = pt_old - CT_diff / dCT_dpt # 1.5 iterations of the modified N-R method

    return pt

def CT_from_pt(SA, pt):
    """
    Calculates Conservative Temperature of seawater from potential temperature (whose reference sea pressure is zero dbar).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar [:math:`^\\circ` C (ITS-90)]

    Returns
    -------
    CT : array_like
         Conservative Temperature [:math:`^\\circ` C (TEOS-10)] FIXME: the orginal has (ITS-90) Copy and paste issue?

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
    >>> pt = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    >>> gsw.CT_from_pt(SA, pt)
    array([[  4.75807226e+00,   1.51169032e+01,   2.28191712e+01,
              3.27053824e+01],
           [  1.55805693e+01,   1.52844796e-02,   2.57405705e+01,
              2.91013409e+01]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See sections 3.1 and 3.3.

    Modifications:
    2010-08-05. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, pt = np.asarray(SA), np.asarray(pt)

    SA[SA < 0] = 0

    sfac = 0.0248826675584615 # sfac = 1 / (40 * ( 35.16504 / 35 ) )

    x2 = sfac * SA
    x = np.sqrt(x2)
    y = pt * 0.025 # normalize for F03 and F08

    pot_enthalpy =  61.01362420681071 + y * ( 168776.46138048015 + \
    y * ( -2735.2785605119625 + y * ( 2574.2164453821433 + \
    y * ( -1536.6644434977543 + y * ( 545.7340497931629 + \
    ( -50.91091728474331 - 18.30489878927802 * y ) * y ) ) ) ) ) + \
    x2 * ( 268.5520265845071 + y * ( -12019.028203559312 + \
    y * ( 3734.858026725145 + y * ( -2046.7671145057618 + \
    y * ( 465.28655623826234 + ( -0.6370820302376359 - \
    10.650848542359153 * y ) * y ) ) ) ) + \
    x * ( 937.2099110620707 + y * ( 588.1802812170108 + \
    y * ( 248.39476522971285 + ( -3.871557904936333 - \
    2.6268019854268356 * y ) * y ) ) + \
    x * ( -1687.914374187449 + x * ( 246.9598888781377 + \
    x * ( 123.59576582457964 - 48.5891069025409 * x ) ) + \
    y * ( 936.3206544460336 + \
    y * ( -942.7827304544439 + y * ( 369.4389437509002 + \
    ( -33.83664947895248 - 9.987880382780322 * y ) * y ) ) ) ) ) )

    #The above polynomial for pot_enthalpy is the full expression for potential enthalpy in terms of SA and pt, obtained from the Gibbs function as below.
    #The above polynomial has simply collected like powers of x and y so that it is computationally faster than calling the Gibbs function twice as is done in the commented code below. When this code below is run, the results are identical to calculating pot_enthalpy as above, to machine precision.
    #n0 = 0
    #n1 = 1
    #pr0 = np.zeros( SA.shape )
    #pot_enthalpy = lib._gibbs(n0, n0, n0, SA, pt, pr0) - (273.15 + pt) * lib._gibbs(n0, n1, n0, SA, pt, pr0)

    #----------------This is the end of the alternative code------------------
    #%timeit gsw.CT_from_pt(SA, pt)
    #1000 loops, best of 3: 1.34 ms per loop <- calling gibbs
    #1000 loops, best of 3: 254 us per loop <- standard

    CT = pot_enthalpy / cte.cp0

    return CT

def t_from_entropy(SA, entropy, t_type='pt'):
    """
    Calculates potential temperature with reference pressure pr = 0 dbar or Conservative temperature from entropy.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    entropy : array_like
              specific entropy [J kg :sup:`-1` K :sup:`-1`]
    t_type : str, optional
             'pt' for potential temperature [:math:`^\\circ` C (ITS-90)]
             'CT' for Conservative Temperature [:math:`^\\circ` C (TEOS-10)] FIXME: the orginal has (ITS-90) Copy and paste issue?

    Returns
    -------
    t : array_like
        potential temperature [:math:`^\\circ` C (ITS-90)] (Default) or Conservative Temperature [:math:`^\\circ` C (TEOS-10)]

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
    >>> entropy = [[63.6913727, 216.184242, 323.536391, 454.589274], [224.455426, -0.147644587, 362.859557, 407.493891]]
    >>> gsw.t_from_entropy(SA, entropy)
    array([[  5.00000000e+00,   1.50000000e+01,   2.20000000e+01,
              3.20000000e+01],
           [  1.50000000e+01,  -1.46624844e-12,   2.50000000e+01,
              2.80000000e+01]])
    >>> entropy = [[67.1618495, 214.564404, 312.442303, 445.351241], [216.404703, -0.371020898, 352.932205, 392.869342]]
    >>> gsw.t_from_entropy(SA, entropy, t_type='CT')
    array([[  5.00000000e+00,   1.50000000e+01,   2.20000000e+01,
              3.20000000e+01],
           [  1.50000000e+01,  -2.64654543e-11,   2.50000000e+01,
              2.80000000e+01]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See appendix  A.10.

    Modifications:
    2010-08-13. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, entropy = np.asarray(SA), np.asarray(entropy)

    SA[SA < 0] = 0

    n0 = 0
    n1 = 1
    pr0 = np.zeros( SA.shape )

    part1 = 1 - SA / cte.SSO
    part2 = 1 - 0.05 * part1
    ent_SA = (cte.cp0 / 273.15) * part1 * ( 1 - 1.01 * part1)
    c = (entropy - ent_SA) * part2 / cte.cp0
    pt = 273.15 * (np.exp(c) - 1)
    dentropy_dt = cte.cp0 / ( (273.15 + pt) * part2) # this is the intial value of dentropy_dt

    for Number_of_iterations in range(0,2,1):
        pt_old = pt
        dentropy = entropy_from_t(SA, pt_old, t_type='pt') - entropy
        pt = pt_old - dentropy / dentropy_dt # this is half way through the modified method
        ptm = 0.5 * (pt + pt_old)
        dentropy_dt = -lib._gibbs_pt0_pt0(SA, ptm)
        t = pt_old - dentropy / dentropy_dt

    if t_type == 'CT':
        t = CT_from_pt(SA, t)

    return t

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
    >>> import seawater.gibbs as gsw
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

    SA = np.nan * np.zeros( SP.shape )
    dSA = np.nan * np.zeros( SP.shape )
    SA_baltic = np.nan * np.zeros( SP.shape )
    in_ocean = np.nan * np.zeros( SP.shape ) #FIXME: change to boolean

    dSA[inds], in_ocean[inds] = lib._delta_SA( p[inds], lon[inds], lat[inds] )
    SA[inds] = ( 35.16504 / 35 ) * SP[inds] + dSA[inds]

    SA_baltic = lib._SA_from_SP_Baltic( SP, lon, lat )

    indsbaltic = ~np.isnan(SA_baltic)

    SA[indsbaltic] = SA_baltic[indsbaltic]

    return SA, in_ocean

def sigma_CT(SA, CT, p=0):
    """
    Calculates potential density anomaly with reference pressure (default is 0 dbar). Returns potential density minus 1000 kg m :sup:`-3`.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
        Conservative Temperature [:math:`^\\circ` C (TEOS-10)]
    p : array_like
        pressure [db], default = 0 db

    Returns
    -------
    sigma_CT : array_like
               potential density anomaly [kg m :sup:`-3`]

    See Also
    --------
    TODO

    Notes
    -----
    Original has 5 versions for this (gsw_sigma0_CT, gsw_sigma1_CT, gsw_sigma2_CT, gsw_sigma3_CT and gsw_sigma4_CT). Here just changed the pressure to the desireded reference.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
    >>> CT = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    >>> gsw.sigma_CT(SA, CT)
    array([[ 41.73409047,  22.04181414,   5.48105772,  10.02188228],
           [  6.84287855,  -0.15791025,   8.44540164,   2.49335766]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. A.30.1.

    Modifications:
    2010-08-26. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, CT, p = np.asarray(SA), np.asarray(CT), np.asarray(p)

    pr0 = np.zeros( SA.shape )
    pt0 = pt_from_CT(SA, CT)

    pref = p + np.zeros( SA.shape )
    tref = pt_from_t(SA, pt0, pr0, pref)

    n0 = 0
    n1 = 1

    rho = np.ones( SA.shape ) / lib._gibbs(n0, n0, n1, SA, tref, pref)

    sigma_CT = rho - 1000

    return sigma_CT

def enthalpy(SA, t, p, t_type='CT', term25=False):
    """
    Calculates the specific enthalpy of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
        temperature [:math:`^\\circ` C]
    p : array_like
        pressure [db]
    t_type : str, optional
            'CT' for Conservative Temperature [:math:`^\\circ` C (TEOS-10)], default
    term25 : bool
             using the computationally-efficient 25-term expression for density in terms of SA, CT and p, default is False

    Returns
    -------
    enthalpy : array_like
               specific enthalpy [ J kg :sup:`-1`]
    in_funnel : bool
                False, if SA, CT and p are outside the "funnel"
                True, if SA, CT and p are inside the "funnel"
                "funnel" is the range of SA, CT and p where the fit error for density was calculated (McDougall et al., 2010).

    See Also
    --------
    TODO

    Notes
    -----
    TODO:
    gsw_enthalpy_CT: Conservative Temperature
    gsw_enthalpy_CT25: Calculates specific enthalpy of seawater using the computationally- efficient 25-term expression for density in terms of SA, CT and p (McDougall et al., 2010)

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
    >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    >>> p = [0., 500., 1500., 2000.]
    >>> gsw.enthalpy(SA, t, p)
    array([[  19959.3397856 ,   64764.8981313 ,  102690.44744487,
             147459.53882186],
           [  59878.01935679,    4994.46567716,  114623.36652738,
             131637.09809679]])
    >>> gsw.enthalpy(SA, t, p, t_type='CT', term25=True)[0]
    array([[  19959.3397856 ,   64764.89794958,  102690.46525807,
             147459.50132288],
           [  59878.01935679,    4994.43975554,  114623.36824859,
             131637.09114752]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See apendix A.11.

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater.  To be submitted to Ocean Science Discussions.

    Modifications:
    2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, t, p = np.asarray(SA), np.asarray(t), np.asarray(p)

    if (t_type == 'CT') & (term25 == True):
        SA[SA<0] = 0
        CT = t
        in_funnel = lib._infunnel(SA, CT, p)

        CT2 = CT**2
        CT3 = CT**3

        a0 = 1 + CT * ( 7.0547681896071576e-3 + \
                 CT * (-1.1753695605858647e-5 + \
                 CT * ( 5.9219809488274903e-7 + \
                 CT * 3.4887902228012519e-10 ) ) ) + \
                 SA * ( 2.0777716085618458e-3 + \
                 CT * ( -2.2210857293722998e-8 + \
                 CT2 * -3.6628141067895282e-10 ) + \
                 np.sqrt(SA) * ( 3.4688210757917340e-6 + \
                 CT2 * 8.0190541528070655e-10 ) )
        a1 = 6.8314629554123324e-6
        a2 = CT3 * -8.5294794834485446e-17
        a3 = CT * -9.2275325145038070e-18

        b0 = 9.9984380290708214e2 + \
             CT * ( 7.1188090678940910e0 + \
             CT * ( -1.9459922513379687e-2 + \
             CT * 6.1748404455874641e-4 ) ) + \
             SA * ( 2.8925731541277653e0 + \
             CT * 2.1471495493268324e-3 + \
             SA * 1.9457531751183059e-3 )
        b1 = 0.5 * ( 1.1930681818531748e-2 + \
             CT2 * 2.6969148011830758e-7 + \
             SA * 5.9355685925035653e-6 )
        b2 = CT2 * -7.2734111712822707e-12 - 2.5943389807429039e-8

        sqrt_disc = np.sqrt( b1**2 - b0 * b2)
        N = a0 + ( 2 * a3 * b0 * b1 / b2 - a2 * b0 ) / b2
        M = a1 + ( 4 * a3 * b1**2 / b2 - ( a3 * b0 + 2 * a2 * b1 ) ) / b2
        A = b1 - sqrt_disc
        B = b1 + sqrt_disc
        part = ( N * b2 - M * b1 ) / ( b2 * (B - A) )

        enthalpy = cte.cp0 * CT + \
                   cte.db2Pascal * ( p * ( a2 - 2 * a3 * b1 / b2 + 0.5 * a3 * p ) / b2 + \
                   ( M / ( 2 * b2 ) ) * np.log( 1 + p * ( 2 * b1 + b2 * p ) / b0 ) + \
                   part * np.log( 1 + ( b2 * p * (B - A) ) / ( A * ( B + b2 * p ) ) ) )

        return enthalpy, in_funnel

    elif (t_type == 'CT') & (term25 == False):
        pt = pt_from_CT(SA, t)
        t = pt_from_t(SA, pt, 0, p)
    else:
        raise NameError('Wrong combination. Read help for mor info')

    n0 = 0
    n1 = 1

    enthalpy = lib._gibbs(n0, n0, n0, SA, t, p) - ( t + cte.Kelvin ) * lib._gibbs(n0, n1, n0, SA, t, p)
    return enthalpy

def t_from_CT(SA, CT, p):
    """
    Calculates in-situ temperature from Conservative Temperature of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
        Conservative Temperature [:math:`^\\circ` C (TEOS-10)]
    p : array_like
        pressure [db]

    Returns
    -------
    t : array_like
        in-situ temperature [:math:`^\\circ` C (ITS-90)]

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
    >>> CT = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    >>> p = 900
    >>> gsw.t_from_CT(SA, CT, p)
    array([[  5.35055483,  15.01761692,  21.36145947,  31.53232787],
           [ 14.54639292,  -0.04443084,  24.45692118,  27.12600316]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See sections 3.1 and 3.3.

    Modifications:
    2010-08-26. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, CT, p = np.asarray(SA), np.asarray(CT), np.asarray(p)

    pr0 = np.zeros( SA.shape )
    pt0 = pt_from_CT(SA, CT)
    t = pt_from_t(SA, pt0, pr0, p)

    return t