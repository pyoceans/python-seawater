#FIXME: temporary_func.py is a place for functions without a home yet...
#TODO: create a temperature and salinity class that are aware of its type
#TODO: compare SIGMAS
#TODO: create term25 class
#TODO: add pt0_from_t as an option for potential_t in SaTePr

#def sigma_CT(SA, CT, p=0):
    #"""
    #Calculates potential density anomaly with reference pressure (default is 0 dbar). Returns potential density minus 1000 kg m :sup:`-3`.

    #Parameters
    #----------
    #SA : array_like
         #Absolute salinity [g kg :sup:`-1`]
    #CT : array_like
        #Conservative Temperature [:math:`^\\circ` C (TEOS-10)]
    #p : array_like
        #pressure [db], default = 0 db

    #Returns
    #-------
    #sigma_CT : array_like
               #potential density anomaly [kg m :sup:`-3`]

    #See Also
    #--------
    #TODO

    #Notes
    #-----
    #Original has 5 versions for this (gsw_sigma0_CT, gsw_sigma1_CT, gsw_sigma2_CT, gsw_sigma3_CT and gsw_sigma4_CT). Here just changed the pressure to the desireded reference.

    #Examples
    #--------
    #>>> import seawater.gibbs as gsw
    #>>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
    #>>> CT = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    #>>> gsw.sigma_CT(SA, CT)
    #array([[ 41.73409047,  22.04181414,   5.48105772,  10.02188228],
           #[  6.84287855,  -0.15791025,   8.44540164,   2.49335766]])

    #References
    #----------
    #.. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. A.30.1.

    #Modifications:
    #2010-08-26. Trevor McDougall & Paul Barker
    #2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    #"""

    ## Convert input to numpy arrays
    #SA, CT, p = np.asarray(SA), np.asarray(CT), np.asarray(p)

    #pr0 = np.zeros( SA.shape )
    #pt0 = pt_from_CT(SA, CT)

    #pref = p + np.zeros( SA.shape )
    #tref = pt_from_t(SA, pt0, pr0, pref)

    #n0 = 0
    #n1 = 1

    #rho = np.ones( SA.shape ) / lib._gibbs(n0, n0, n1, SA, tref, pref)

    #sigma_CT = rho - 1000

    #return sigma_CT

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


    X = np.sin( np.deg2rad(lat) )
    sin2 = X**2
    gs = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2)
    z = z_from_p(p, lat)
    grav = gs * (1 - cte.gamma * z) # z is the height corresponding to p
    return grav

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

if __name__=='__main__':
    """ entropy_from_t 'CT'"""
    entropy_from_CT =  gsw.entropy_from_t(SA_chck_cast, CT_chck_cast, t_type='CT')
    Ientropy_from_CT = (gsw_cv.entropy_from_CT - entropy_from_CT) >= gsw_cv.entropy_from_CT_ca

    if Ientropy_from_CT.any():
        print "entropy_from_CT:   Failed"
    else:
        print "entropy_from_CT:   Passed"

    """ entropy_from_t 'pt' """ #FIXME: pass, but small float are detected, investigate further
    entropy_from_pt =  gsw.entropy_from_t(SA_chck_cast, pt)
    Ientropy_from_pt = (gsw_cv.entropy_from_pt - entropy_from_pt) >= gsw_cv.entropy_from_pt_ca

    if Ientropy_from_pt.any():
        print "entropy_from_pt:   Failed"
    else:
        print "entropy_from_pt:   Passed"

    """ pt_from_CT """
    pt_from_CT = gsw.pt_from_CT(SA_chck_cast, CT_chck_cast)
    Ipt_from_CT = (gsw_cv.pt - pt_from_CT) >= gsw_cv.pt_ca

    if Ipt_from_CT.any():
        print "pt_from_CT:   Failed"
    else:
        print "pt_from_CT:   Passed"

    """ CT_from_pt """
    CT_from_pt = gsw.CT_from_pt(SA_chck_cast, pt) #FIXME: pass, but small float are detected, investigate further
    ICT_from_pt = (gsw_cv.CT_from_pt - CT_from_pt) >= gsw_cv.CT_from_pt_ca

    if ICT_from_pt.any():
        print "CT_from_pt:   Failed"
    else:
        print "CT_from_pt:   Passed"

    """ CT_from_entropy """ #FIXME: pass, but small float are detected, investigate further
    CT_from_entropy =  gsw.t_from_entropy(SA_chck_cast, entropy, 'CT')
    ICT_from_entropy = (gsw_cv.CT_from_entropy - CT_from_entropy) >= gsw_cv.CT_from_entropy_ca

    if ICT_from_entropy.any():
        print "CT_from_entropy:   Failed"
    else:
        print "CT_from_entropy:   Passed"

    """ t_from_CT """ #FIXME: pass, but small float are detected, investigate further
    t_from_CT =  gsw.t_from_CT(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
    It_from_CT = (gsw_cv.t_chck_cast - t_from_CT) >= gsw_cv.t_from_CT_ca

    if It_from_CT.any():
        print "t_from_CT:   Failed"
    else:
        print "t_from_CT:   Passed"

    """ enthalpy_CT """
    enthalpy_CT =  gsw.enthalpy(SA_chck_cast,CT_chck_cast, gsw_cv.p_chck_cast, t_type='CT')
    Ienthalpy_CT = (gsw_cv.enthalpy_CT - enthalpy_CT) >= gsw_cv.enthalpy_CT_ca

    if Ienthalpy_CT.any():
        print "enthalpy_CT:   Failed"
    else:
        print "enthalpy_CT:   Passed"

    """ enthalpy_CT25 """ #FIXME: pass, but small float are detected, investigate further
    enthalpy_CT25 =  gsw.enthalpy(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast, t_type='CT', term25=True)[0]
    Ienthalpy_CT25 = (gsw_cv.enthalpy_CT25 - enthalpy_CT25) >= gsw_cv.enthalpy_CT25_ca

    if Ienthalpy_CT25.any():
        print "enthalpy_CT25:   Failed"
    else:
        print "enthalpy_CT25:   Passed"

#""" sigma_CT """
#sigma0_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast)
#Isigma_CT = (gsw_cv.sigma0_CT - sigma0_CT) >= gsw_cv.sigma0_CT_ca

#if Isigma_CT.any():
    #print "sigma_CT at 0 db:   Failed"
#else:
    #print "sigma_CT at 0 db:   Passed"

#sigma1_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast, 1000)
#Isigma_CT = (gsw_cv.sigma1_CT - sigma1_CT) >= gsw_cv.sigma1_CT_ca

#if Isigma_CT.any():
    #print "sigma_CT at 1000 db:   Failed"
#else:
    #print "sigma_CT at 1000 db:   Passed"

#sigma2_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast, 2000)
#Isigma_CT = (gsw_cv.sigma2_CT - sigma2_CT) >= gsw_cv.sigma2_CT_ca

#if Isigma_CT.any():
    #print "sigma_CT at 2000 db:   Failed"
#else:
    #print "sigma_CT at 2000 db:   Passed"

#sigma3_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast, 3000)
#Isigma_CT = (gsw_cv.sigma3_CT - sigma3_CT) >= gsw_cv.sigma3_CT_ca

#if Isigma_CT.any():
    #print "sigma_CT at 3000 db:   Failed"
#else:
    #print "sigma_CT at 3000 db:   Passed"

#sigma4_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast, 4000)
#Isigma_CT = (gsw_cv.sigma4_CT - sigma4_CT) >= gsw_cv.sigma4_CT_ca

#if Isigma_CT.any():
    #print "sigma_CT at 4000 db:   Failed"
#else:
    #print "sigma_CT at 4000 db:   Passed"

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

    #prop1.dtype = 'float64' # FIXME: somehow lon is been loaded as uint8
    return prop1