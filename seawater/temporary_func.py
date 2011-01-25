#FIXME: temporary_func.py is a place for functions without a home yet...

import numpy as np
from seawater import constants as cte
from seawater import library as lib

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

# TEST:
""" enthalpy_CT """
#enthalpy_CT =  gsw.enthalpy(SA_chck_cast,CT_chck_cast, gsw_cv.p_chck_cast, t_type='CT')
""" enthalpy_CT25 """
#enthalpy_CT25 =  gsw.enthalpy(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast, t_type='CT', erm25=True)[0]
""" sigma_CT """
#sigma0_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast)
#sigma1_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast, 1000)
#sigma2_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast, 2000)
#sigma3_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast, 3000)
#sigma4_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast, 4000)