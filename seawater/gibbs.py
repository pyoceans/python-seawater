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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-09-28. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA = np.asarray(SA)

    Isalty = np.where(SA >= 0)
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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    See Table L.1 of this TEOS-10 Manual.

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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    See appendix A.10 of this TEOS-10 Manual.

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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    See sections 3.1 and 3.3 of this TEOS-10 Manual.

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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    See sections 3.1 and 3.3 of this TEOS-10 Manual.

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

def pt0_from_t(SA, t, p):
    """
    Calculates potential temperature with reference pressure, pr = 0 dbar. The present routine is computationally faster than the more general function "pt_from_t(SA, t, p, pr)" which can be used for any reference pressure value.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
         in-situ temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]

    Returns
    -------
    pt0 : array_like
          potential temperature relative to 0 db [:math:`^\\circ` C (FIXME: ?ITS-90)]

    See Also
    --------
    entropy_part, entropy_part_zero, gibbs_pt0_pt0

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
    >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    >>> p = 900
    >>> gsw.pt0_from_t(SA, t, p)
    array([[  4.89971486e+00,   1.48664023e+01,   2.18420392e+01,
              3.17741959e+01],
           [  1.48891940e+01,   2.95267636e-02,   2.48187231e+01,
              2.78058513e+01]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    See section 3.1 of this TEOS-10 Manual.

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater.  To be submitted to Ocean Science Discussions.

    Modifications:
    2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, t, p = np.asarray(SA), np.asarray(t), np.asarray(p)

    SA[SA < 0] = 0

    s1 = SA*(35. / cte.SSO)

    pt0 = t + p * ( 8.65483913395442e-6  - \
    s1 * 1.41636299744881e-6 - \
    p * 7.38286467135737e-9 + \
    t * ( -8.38241357039698e-6 + \
    s1 * 2.83933368585534e-8 + \
    t * 1.77803965218656e-8 + \
    p * 1.71155619208233e-10 ) )

    dentropy_dt = cte.cp0 / ( (273.15 + pt0) * ( 1 - 0.05 * ( 1 - SA / cte.SSO ) ) )

    true_entropy_part = lib._entropy_part(SA, t, p)

    for Number_of_iterations in range(0,2,1):
        pt0_old = pt0
        dentropy = lib._entropy_part_zerop(SA, pt0_old) - true_entropy_part
        pt0 = pt0_old - dentropy / dentropy_dt # this is half way through the modified method
        pt0m = 0.5 * (pt0 + pt0_old);
        dentropy_dt = -lib._gibbs_pt0_pt0(SA, pt0m)
        pt0 = pt0_old - dentropy / dentropy_dt

    # maximum error of 6.3x10^-9 degrees C for one iteration.
    # maximum error is 1.8x10^-14 degrees C for two iterations (two iterations is the default, "for Number_of_iterations = 1:2")
    # These errors are over the full "oceanographic funnel" of McDougall et al. (2010), which reaches down to p = 8000 dbar.

    return pt0

def CT_from_t(SA, t, p):
    """
    Calculates Conservative Temperature of seawater from in-situ temperature.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
         in-situ temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]

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
    >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    >>> p = 900
    >>> gsw.CT_from_t(SA, t, p)
    array([[  4.66028901,  14.98237022,  22.6558658 ,  32.47483113],
           [ 15.46594688,   0.04649395,  25.55437701,  28.90014276]])


    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    See section 3.3 of this TEOS-10 Manual.

    Modifications:
    2010-08-26. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, t, p = np.asarray(SA), np.asarray(t), np.asarray(p)

    pt0 = pt0_from_t(SA, t, p)
    CT = CT_from_pt(SA, pt0)

    return CT

def pt_from_t(SA, t, p, pr=0):
    """
    Calculates potential temperature with the general reference pressure, pr, from in-situ temperature.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
         in-situ temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]
    pr : int, float, optional
         reference pressure, default = 0

    Returns
    -------
    pt : array_like
         potential temperature [:math:`^\\circ` C (ITS-90)]

    See Also
    --------
    entropy_part

    Notes
    -----
    This function calls "entropy_part" which evaluates entropy except for the parts which are a function of Absolute Salinity alone. A faster routine exists pt0_from_t(SA,t,p) if pr is indeed zero dbar.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
    >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    >>> p = 900
    >>> gsw.pt_from_t(SA, t, p)
    array([[  4.89971486e+00,   1.48664023e+01,   2.18420392e+01,
              3.17741959e+01],
           [  1.48891940e+01,   2.95267636e-02,   2.48187231e+01,
              2.78058513e+01]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    See section 3.1 of this TEOS-10 Manual.

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010: A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater. To be submitted to Ocean Science Discussions.

    Modifications:
    2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, t, p, pr = np.asarray(SA), np.asarray(t), np.asarray(p), np.asarray(pr)

    SA[SA < 0] = 0

    n0 = 0
    n2 = 2

    s1 = SA * 35. / cte.SSO

    pt = t + (p-pr) * ( 8.65483913395442e-6  - \
    s1 * 1.41636299744881e-6 - \
    (p+pr) * 7.38286467135737e-9 + \
    t * ( -8.38241357039698e-6 + \
    s1 * 2.83933368585534e-8 + \
    t * 1.77803965218656e-8 + \
    (p+pr) * 1.71155619208233e-10 ) )

    dentropy_dt = cte.cp0 / ( (273.15 + pt) * ( 1 - 0.05 * ( 1 - SA / cte.SSO ) ) )

    true_entropy_part = lib._entropy_part(SA, t, p)

    for Number_of_iterations in range(0,2,1):
        pt_old = pt
        dentropy = lib._entropy_part(SA, pt_old, pr) - true_entropy_part
        pt = pt_old - dentropy / dentropy_dt # this is half way through the modified method
        ptm = 0.5 * (pt + pt_old)
        dentropy_dt = -lib._gibbs(n0, n2, n0, SA, ptm, pr)
        pt = pt_old - dentropy / dentropy_dt


    # maximum error of 6.3x10^-9 degrees C for one iteration.
    # maximum error is 1.8x10^-14 degrees C for two iterations (two iterations is the default, "for Number_of_iterations = 1:2).
    # These errors are over the full "oceanographic funnel" of McDougall et al. (2010), which reaches down to p = 8000 dbar.

    return pt

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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    See appendix  A.10 of this TEOS-10 Manual.

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

def entropy(SA, t, p):
    """
    Calculates potential temperature with reference pressure pr = 0 dbar or Conservative temperature from entropy.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
         in-situ temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]

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
    >>> p = [[0., 500., 1500., 2000.], [0., 500., 1500., 2000.]]
    >>> gsw.entropy(SA, t, p)
    array([[  6.36913727e+01,   2.15161921e+02,   3.19806445e+02,
              4.47838663e+02],
           [  2.24455426e+02,   1.43185159e-01,   3.58666432e+02,
              4.01487857e+02]])

    TODO

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, t, p = np.asarray(SA), np.asarray(t), np.asarray(p)

    n0 = 0
    n1 = 1
    entropy = -lib._gibbs(n0, n1, n0, SA, t, p)

    return entropy

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
    >>> p = [[0., 500., 1500., 2000.], [0., 500., 1500., 2000.]]
    >>> gsw.SA_from_SP


    TODO

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    See section 2.5 and appendices A.4 and A.5 of this TEOS-10 Manual.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SP, p, lon, lat = np.asarray(SP), np.asarray(p), np.asarray(lon), np.asarray(lat)

    #TODO: maybe make this check a function?
    p = check_dim(p, SP)
    lat = check_dim(p, SP)
    lon = check_dim(p, SP)


    SP[SP < 0] = 0

    inds = np.where( np.isfinite(SP) )

    SA = np.nan * np.zeros( SP.shape )
    dSA = SA
    in_ocean = SA #FIXME

    dSA[inds], in_ocean[inds] = lib._delta_SA( p[inds], lon[inds], lat[inds] )
    SA[inds] = ( 35.16504 / 35 ) * SP[inds] + dSA[inds]

    SA_baltic[inds] = lib._SA_from_SP_Baltic( SP[inds], lon[inds], lat[inds] )

    indsbaltic = np.where( ~np.isnan( SA_baltic[inds] ) )

    SA[inds[indsbaltic]] = SA_baltic[inds[indsbaltic]]

    return SA, in_ocean

def check_dim(prop1, prop2):
    """
    Broadcast prop1 to the shape prop2. Prop1 can be scalar, row equal or column equal to prop2.
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
            print "add a proper error msg"

    prop1.dtype = 'float64' # FIXME: somehow lon is been loaded as uint8
    return prop1