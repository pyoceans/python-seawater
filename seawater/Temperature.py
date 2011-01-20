#TODO: study the possibility of transforming this into a DataSet like class with "properties" defining the different temperatures.

import numpy as np
from seawater import constants as cte
from seawater import library as lib

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

    x2 = cte.sfac * SA
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

def pt_from_CT(SA, CT): #NOTE: used at specvol_anom
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

def potential_t(SA, t, p, pr=0): # FIXME: repeated inside SaTePr due to different call (see spec_anom)
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
    TODO

    Notes
    -----
    This function calls "entropy_part" which evaluates entropy except for the parts which are a function of Absolute Salinity alone. A faster routine exists pt0_from_t(SA,t,p) if pr is indeed zero dbar.

    Examples
    --------
    >>> from seawater.gibbs import SaTePr
    >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
    >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    >>> p = 900
    >>> STP = SaTePr(SA, t, p)
    >>> STP.potential_t()
    TODO
    >>> STP.potential_t(pr = 900)
    TODO

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 3.1.

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010: A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater. To be submitted to Ocean Science Discussions.

    Modifications:
    2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, t, p = np.asarray(SA), np.asarray(t), np.asarray(p)
    pr = np.asarray(pr)

    SA[SA < 0] = 0
    n0, n2 = 0, 2

    s1 = SA * 35. / cte.SSO

    pt = t + ( p - pr ) * ( 8.65483913395442e-6  - \
    s1 * 1.41636299744881e-6 - \
    ( p + pr ) * 7.38286467135737e-9 + \
    t * ( -8.38241357039698e-6 + \
    s1 * 2.83933368585534e-8 + \
    t * 1.77803965218656e-8 + \
    ( p + pr ) * 1.71155619208233e-10 ) )

    dentropy_dt = cte.cp0 / ( (cte.Kelvin + pt) * ( 1 - 0.05 * ( 1 - SA / cte.SSO ) ) )

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

#def pt0_from_t(self): #NOTE: Not necessary at the moment, since potential_t does this
    #"""
    #Calculates potential temperature with reference pressure, pr = 0 dbar. The present routine is computationally faster than the more general function "pt_from_t(SA, t, p, pr)" which can be used for any reference pressure value.

    #Returns
    #-------
    #pt0 : array_like
            #potential temperature relative to 0 db [:math:`^\\circ` C (ITS-90)]

    #See Also
    #--------
    #TODO

    #Notes
    #-----
    #TODO

    #Examples
    #--------
    #>>> from seawater.SaTePr import SaTePr
    #>>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
    #>>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
    #>>> p = 900
    #>>> STP = SaTePr(SA, t, p)
    #>>> STP.pt0_from_t()
    #array([[  4.89971486e+00,   1.48664023e+01,   2.18420392e+01,
                #3.17741959e+01],
            #[  1.48891940e+01,   2.95267636e-02,   2.48187231e+01,
                #2.78058513e+01]])

    #References
    #----------
    #.. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 3.1.

    #.. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater.  To be submitted to Ocean Science Discussions.

    #Modifications:
    #2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
    #2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    #"""

    #SA = self.masked_SA.filled() # ensure that SA is non-negative

    #s1 = SA * (35. / cte.SSO)

    #pt0 = self.t + self.p * ( 8.65483913395442e-6  - \
            #s1 * 1.41636299744881e-6 - \
            #self.p * 7.38286467135737e-9 + \
            #self.t * ( -8.38241357039698e-6 + \
            #s1 * 2.83933368585534e-8 + \
            #self.t * 1.77803965218656e-8 + \
            #self.p * 1.71155619208233e-10 ) )

    #dentropy_dt = cte.cp0 / ( (cte.Kelvin + pt0) * ( 1 - 0.05 * ( 1 - SA / cte.SSO ) ) )

    #true_entropy_part = lib._entropy_part(SA, self.t, self.p)

    #for Number_of_iterations in range(0,2,1):
        #pt0_old = pt0
        #dentropy = lib._entropy_part_zerop(SA, pt0_old) - true_entropy_part
        #pt0 = pt0_old - dentropy / dentropy_dt # this is half way through the modified method
        #pt0m = 0.5 * (pt0 + pt0_old);
        #dentropy_dt = -lib._gibbs_pt0_pt0(SA, pt0m)
        #pt0 = pt0_old - dentropy / dentropy_dt

    ## maximum error of 6.3x10^-9 degrees C for one iteration.
    ## maximum error is 1.8x10^-14 degrees C for two iterations (two iterations is the default, "for Number_of_iterations = 1:2")
    ## These errors are over the full "oceanographic funnel" of McDougall et al. (2010), which reaches down to p = 8000 dbar.

    #return pt0