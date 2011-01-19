#TODO: Go over PDFs to improve documentation
#TODO: Examples: simple with the data range (copy-and-paste numbers) and complex (real data)
#TODO: Check original authors and dates
#TODO: csiro vs gibbs (table?)
#TODO: check_dim for p in all "p" functions
#FIXME: some function return values even with nan in the input, check this behaviour (also present in the original).

import numpy as np
from seawater import constants as cte
from seawater import library as lib

import numpy.ma as ma

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

def  z_from_p(p, lat): # Also in DepthPressure.py
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

def  p_from_z(z, lat): # Also in DepthPressure.py
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

def grav(lat, p=0): # Also in DepthPressure.py
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

class SaTePr: #TODO: find a better name!
    """
    Class that agregatte all Sa, t, p functions

    Parameters
    ----------
    SA : array_like
        Absolute salinity [g kg :sup:`-1`]
    t : array_like
        in-situ temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]

    """
    def __init__(self, SA=None, t=None, p=None):
        # Convert input to numpy arrays
        self.SA, self.t, self.p = np.asarray(SA), np.asarray(t), np.asarray(p)

        self.masked_SA = ma.masked_less(self.SA, 0)
        self.masked_SA.fill_value = 0

        # order for the gibbs function
        self.n0, self.n1, self.n2 = 0, 1, 2

    def entropy(self):
        """
        Calculates potential temperature with reference pressure pr = 0 dbar or Conservative temperature from entropy.

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
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = [0., 500., 1500., 2000.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.entropy()
        array([[  6.36913727e+01,   2.15161921e+02,   3.19806445e+02,
                  4.47838663e+02],
               [  2.24455426e+02,   1.43185159e-01,   3.58666432e+02,
                  4.01487857e+02]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

        Modifications:
        2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        entropy = -lib._gibbs(self.n0, self.n1, self.n0, self.SA, self.t, self.p)

        return entropy

    def rho(self):
        """
        Calculates in-situ density of seawater from Absolute Salinity and in-situ temperature.

        Returns
        -------
        rho : array_like
            in-situ density [kg m :sup:`-3`]

        See Also
        --------
        TODO

        Notes
        -----
        TODO

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = [0., 500., 1500., 2000.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.rho()
        array([[ 1041.77425464,  1024.2413978 ,  1011.923534  ,  1018.28328036],
               [ 1006.74841976,  1002.37206267,  1014.78353156,  1010.8696052 ]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.8.

        Modifications:
        2010-07-23. David Jackett, Trevor McDougall & Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        rho = 1. / lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p)

        return rho

    def cp(self):
        """
        Calculates the isobaric heat capacity of seawater.

        Returns
        -------
        cp : array_like
             heat capacity of seawater [ J kg :sup:`-1` K:sup:`-1`]

        See Also
        --------
        TODO

        Notes
        -----
        TODO

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = 900
        >>> STP = SaTePr(SA, t, p)
        >>> STP.cp()
        array([[ 3869.46487578,  3996.62909658,  4102.39689639,  4056.09090058],
               [ 4102.00085198,  4176.72470928,  4077.47206662,  4114.01189933]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

        Modifications:
        2010-07-23. David Jackett, Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        cp = -( self.t + cte.Kelvin ) * lib._gibbs(self.n0, self.n2, self.n0, self.SA, self.t, self.p)

        return cp

    def helmholtz_energy(self):
        """
        Calculates the Helmholtz energy of seawater.

        Returns
        -------
        Helmholtz_energy : array_like
                           Helmholtz energy [ J kg :sup:`-1`]

        See Also
        --------
        TODO

        Notes
        -----
        TODO

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = [0., 500., 1500., 2000.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.helmholtz_energy()
        array([[  1.18057894e+03,  -2.04243623e+03,  -4.45224072e+03,
                 -8.18003196e+03],
               [ -2.58190138e+03,   6.54845497e+00,  -5.48590282e+03,
                 -6.56341929e+03]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.13.

        Modifications:
        2010-08-26. Trevor McDougall
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        helmholtz_energy = lib._gibbs(self.n0, self.n0, self.n0, self.SA, self.t, self.p) - \
                        ( cte.db2Pascal * self.p + 101325 ) * lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p)

        return helmholtz_energy

    def internal_energy(self):
        """
        Calculates the Helmholtz energy of seawater.

        Returns
        -------
        internal_energy(u) : array_like #FIXME: function of "u" !? What is "u"?
                             specific internal energy [ J kg :sup:`-1`]

        See Also
        --------
        TODO

        Notes
        -----
        TODO

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = [0., 500., 1500., 2000.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.internal_energy()
        array([[  1.88963342e+04,   5.99564714e+04,   8.99386314e+04,
                  1.28477936e+05],
               [  6.20949295e+04,   4.56594812e+01,   1.01450494e+05,
                  1.14344649e+05]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.13.

        Modifications:
        2010-08-26. Trevor McDougall
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        internal_energy = lib._gibbs(self.n0, self.n0, self.n0, self.SA, self.t, self.p) - \
                        (cte.Kelvin + self.t) * lib._gibbs(self.n0, self.n1, self.n0, self.SA, self.t, self.p) - \
                        (cte.db2Pascal * self.p + 101325) * lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p)

        return internal_energy

    def sound_speed(self):
        """
        Calculates the speed of sound in seawater.

        Returns
        -------
        sound_speed : array_like
                      speed of sound in seawater [ m s :sup:`-1`]

        See Also
        --------
        TODO

        Notes
        -----
        TODO

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = [0., 500., 1500., 2000.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.sound_speed()
        array([[ 1493.5609568 ,  1508.86141015,  1524.04873089,  1567.35919386],
               [ 1477.63190763,  1410.40969768,  1537.60287636,  1546.09128039]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.17.1)

        Modifications:
        2010-07-23. Trevor McDougall
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        g_tt = lib._gibbs(self.n0, self.n2, self.n0, self.SA, self.t, self.p)
        g_tp = lib._gibbs(self.n0, self.n1, self.n1, self.SA, self.t, self.p)

        sound_speed = lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p) * \
        np.sqrt( g_tt / ( g_tp * g_tp - g_tt * lib._gibbs(self.n0, self.n0, self.n2, self.SA, self.t, self.p ) ) )

        return sound_speed

    def kappa(self):
        """
        Calculates the isentropic compressibility of seawater.

        Returns
        -------
        kappa : array_like
                Isentropic compressibility [ Pa :sup:`-1`]

        See Also
        --------
        TODO

        Notes
        -----
        The output is Pascal and not dbar.

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = [0., 500., 1500., 2000.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.kappa()
        array([[  4.30309045e-10,   4.28843638e-10,   4.25455945e-10,
                  3.99755378e-10],
               [  4.54932038e-10,   5.01511014e-10,   4.16810090e-10,
                  4.13842034e-10]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqns. (2.16.1) and the row for kappa in Table P.1 of appendix P

        Modifications:
        2010-07-23. David Jackett, Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        g_tt = lib._gibbs(self.n0, self.n2, self.n0, self.SA, self.t, self.p)
        g_tp = lib._gibbs(self.n0, self.n1, self.n1, self.SA, self.t, self.p)

        kappa = ( g_tp * g_tp - g_tt * lib._gibbs(self.n0, self.n0, self.n2, self.SA, self.t, self.p) ) / ( lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p ) * g_tt)

        return kappa

    def adiabatic_lapse_rate(self):
        """
        Calculates the adiabatic lapse rate of sea water.

        Returns
        -------
        adiabatic_lapse_rate : array_like
                               Adiabatic lapse rate [ K Pa :sup:`-1`]

        See Also
        --------
        TODO

        Notes
        -----
        The output is in unit of degress Celsius per Pa, (or equivilently K/Pa) not in units of K/dbar

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = [0., 500., 1500., 2000.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.adiabatic_lapse_rate()
        array([[  1.05756574e-08,   1.49457941e-08,   1.85280735e-08,
                  2.58480453e-08],
               [  1.18016760e-08,  -3.17131249e-09,   2.09612644e-08,
                  2.26342914e-08]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.22.1).

        Modifications:
        2010-07-23. Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        adiabatic_lapse_rate = - lib._gibbs(self.n0, self.n1, self.n1, self.SA, self.t, self.p) / ( lib._gibbs(self.n0, self.n2, self.n0, self.SA, self.t, self.p ) )

        return adiabatic_lapse_rate

    def chem_potential_relative(self):
        """
        Calculates the adiabatic lapse rate of sea water.

        Returns
        -------
        chem_potential_relative : array_like
                                  relative chemical potential [ J kg :sup:`-1`]

        See Also
        --------
        TODO

        Notes
        -----
        TODO

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = [0., 500., 1500., 2000.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.chem_potential_relative()
        array([[ 96.16894372,  55.28852987, -27.68689914,  25.15417648],
               [-18.85739517,          nan,   2.85509781, -44.75926356]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

        Modifications:
        2010-08-26. Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        chem_potential_relative = lib._gibbs(self.n1, self.n0, self.n0, self.SA, self.t, self.p)

        return chem_potential_relative

    def specvol(self):
        """
        Calculates the specific volume of seawater.

        Returns
        -------
        specvol : array_like
                  specific volume [ m :sup:`-3` kg :sup:`-1`] #TODO: original has a typo [ kg m :sup:`-3`]

        See Also
        --------
        TODO

        Notes
        -----
        TODO

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = [0., 500., 1500., 2000.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.specvol()
        array([[ 0.0009599 ,  0.00097633,  0.00098822,  0.00098204],
               [ 0.0009933 ,  0.00099763,  0.00098543,  0.00098925]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.7.

        Modifications:
        2010-08-26. David Jackett & Paul Barker.
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        specvol = lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p)

        return specvol

    def conservative_t(self):
        """
        Calculates Conservative Temperature of seawater from in-situ temperature.

        Returns
        -------
        CT : array_like
             Conservative Temperature [:math:`^\\circ` C (TEOS-10)]

        See Also
        --------
        TODO

        Notes
        -----
        TODO

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = 900
        >>> STP = SaTePr(SA, t, p)
        >>> STP.conservative_t()
        array([[  4.66028901,  14.98237022,  22.6558658 ,  32.47483113],
               [ 15.46594688,   0.04649395,  25.55437701,  28.90014276]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 3.3.

        Modifications:
        2010-08-26. David Jackett, Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        #pt0 = self.pt0_from_t() #TODO: pt0_from_t
        pt0 = self.potential_t()

        #CT = gsw.CT_from_pt(self.SA, pt0) NOTE: incorporated below

        SA = self.masked_SA.filled() # ensure that SA is non-negative

        sfac = 0.0248826675584615 # sfac = 1 / (40 * ( 35.16504 / 35 ) )

        x2 = sfac * SA
        x = np.sqrt(x2)
        y = pt0 * 0.025 # normalize for F03 and F08

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

        #The above polynomial for pot_enthalpy is the full expression for potential enthalpy in terms of SA and pt, obtained from the Gibbs function as below. The above polynomial has simply collected like powers of x and y so that it is computationally faster than calling the Gibbs function twice as is done in the commented code below. When this code below is run, the results are identical to calculating pot_enthalpy as above, to machine precision.

        #pr0 = np.zeros( self.SA.shape )
        #pot_enthalpy = lib._gibbs(self.n0, self.n0, self.n0, self.SA, pt0, pr0) - (cte.Kelvin + pt0) * lib._gibbs(self.n0, self.n1, self.n0, self.SA, pt0, pr0)

        #----------------This is the end of the alternative code------------------
        #%timeit results
        #1000 loops, best of 3: 1.34 ms per loop <- calling gibbs
        #1000 loops, best of 3: 254 us per loop <- standard

        CT = pot_enthalpy / cte.cp0

        return CT

    def molality(self):
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
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> STP = SaTePr(SA)
        >>> STP.molality()
        array([[ 1.78214644,  0.98484303,  0.32164907,  0.64986241],
               [ 0.32164907,         nan,  0.48492271,  0.25680047]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

        Modifications:
        2010-09-28. Trevor McDougall & Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        # only >= than zero
        self.masked_SA.fill_value = np.nan
        SA = self.masked_SA

        # molality of seawater in mol kg :sup:`-1`
        molality = SA / (cte.M_S * ( 1000 - SA ) )

        return molality.filled()

    def ionic_strength(self):
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
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> STP = SaTePr(SA)
        >>> STP.ionic_strength()
        array([[ 1.10964439,  0.61320749,  0.20027315,  0.40463351],
               [ 0.20027315,         nan,  0.30193465,  0.1598955 ]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Table L.1.

        .. [2] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: The composition of Standard Seawater and the definition of the Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72. see Eqns. 5.9 and 5.12.

        Modifications:
        2010-09-28. Trevor McDougall & Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        Z_2 = 1.2452898 # the valence factor of sea salt

        molal = self.molality() # molality of seawater in mol kg

        ionic_strength = 0.5 * Z_2 * molal

        return ionic_strength

    def potential_t(self, pr=0):
        """
        Calculates potential temperature with the general reference pressure, pr, from in-situ temperature.

        Parameters
        ----------
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
        array([[  4.89971486e+00,   1.48664023e+01,   2.18420392e+01,
                  3.17741959e+01],
               [  1.48891940e+01,   2.95267636e-02,   2.48187231e+01,
                  2.78058513e+01]])
        >>> STP.potential_t(pr = 900)
        array([[  5.,  15.,  22.,  32.],
               [ 15.,   0.,  25.,  28.]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 3.1.

        .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010: A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater. To be submitted to Ocean Science Discussions.

        Modifications:
        2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        # Convert input to numpy arrays
        pr = np.asarray(pr)
        SA = self.masked_SA.filled() # ensure that SA is non-negative

        s1 = SA * 35. / cte.SSO

        pt = self.t + ( self.p - pr ) * ( 8.65483913395442e-6  - \
        s1 * 1.41636299744881e-6 - \
        ( self.p + pr ) * 7.38286467135737e-9 + \
        self.t * ( -8.38241357039698e-6 + \
        s1 * 2.83933368585534e-8 + \
        self.t * 1.77803965218656e-8 + \
        ( self.p + pr ) * 1.71155619208233e-10 ) )

        dentropy_dt = cte.cp0 / ( (cte.Kelvin + pt) * ( 1 - 0.05 * ( 1 - SA / cte.SSO ) ) )

        true_entropy_part = lib._entropy_part(SA, self.t, self.p)

        for Number_of_iterations in range(0,2,1):
            pt_old = pt
            dentropy = lib._entropy_part(SA, pt_old, pr) - true_entropy_part
            pt = pt_old - dentropy / dentropy_dt # this is half way through the modified method
            ptm = 0.5 * (pt + pt_old)
            dentropy_dt = -lib._gibbs(self.n0, self.n2, self.n0, SA, ptm, pr)
            pt = pt_old - dentropy / dentropy_dt

        # maximum error of 6.3x10^-9 degrees C for one iteration.
        # maximum error is 1.8x10^-14 degrees C for two iterations (two iterations is the default, "for Number_of_iterations = 1:2).
        # These errors are over the full "oceanographic funnel" of McDougall et al. (2010), which reaches down to p = 8000 dbar.

        return pt

    def enthalpy(self):
        """
        Calculates the specific enthalpy of seawater.

        Returns
        -------
        enthalpy : array_like
                   specific enthalpy [ J kg :sup:`-1`]

        See Also
        --------
        TODO

        Notes
        -----
        TODO

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = [0., 500., 1500., 2000.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.enthalpy()
        array([[  18993.59620275,   64937.05999321,  104862.01693673,
                 148218.3415969 ],
               [  62195.57534579,    5134.91245416,  116331.82020187,
                 134229.82985461]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See apendix A.11.

        Modifications:
        2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        enthalpy = lib._gibbs(self.n0, self.n0, self.n0, self.SA, self.t, self.p) - ( self.t + cte.Kelvin ) * lib._gibbs(self.n0, self.n1, self.n0, self.SA, self.t, self.p)

        return enthalpy

    def alpha(self):
        """
        Calculates the thermal expansion coefficient of seawater with respect to in-situ temperature.

        Returns
        -------
        alpha : array_like
                thermal expansion coefficient [ K :sup:`-1`]

        See Also
        --------
        TODO

        Notes
        -----
        original name: gsw_alpha_wrt_t (A.K.A with respect to in-situ temperature)

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = [0., 500., 1500., 2000.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.alpha()
        array([[  1.54174741e-04,   2.12859667e-04,   2.59617457e-04,
                  3.47907236e-04],
               [  1.70265060e-04,  -4.88225022e-05,   2.89880704e-04,
                  3.10594834e-04]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.18.1)

        .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.

        Modifications:
        2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        alpha = lib._gibbs(self.n0, self.n1, self.n1, self.SA, self.t, self.p) / lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p)

        return alpha

if __name__=='__main__':
    try:
        import cPickle as pickle
    except:
        import pickle
    import numpy as np

    """ load test data """
    class Dict2Struc(object):
        """ all the variables from a dict in a "matlab-like-structure" """
        def __init__(self, adict):
            self.__dict__.update(adict)


    data = pickle.load( open('gsw_cv.pkl','rb') )
    gsw_cv = Dict2Struc(data) # then type dat.<tab> to navigate through your variables

    def test_print(STP, method, comp_value):
        """
        TODO
        """
        exec( "inequal = (gsw_cv." +comp_value+ " - STP." +method+ "() ) >= gsw_cv." +comp_value+ "_ca")

        width = 23
        if inequal.any():
            print "%s: Failed" % method.rjust(width)
        else:
            if eval( "( gsw_cv." +comp_value+ "[~np.isnan(gsw_cv."+comp_value+")] == STP." +method+ "()[~np.isnan(STP."+method+"())] ).all()" ):
                print "%s: Passed, equal" % method.rjust(width)
            else:
                exec("nmax = ( gsw_cv."+comp_value+" - STP."+method+"() )[~np.isnan(gsw_cv."+comp_value+")].max()")
                exec("nmin = ( gsw_cv."+comp_value+" - STP."+method+"() )[~np.isnan(gsw_cv."+comp_value+")].min()")
                print "%s: Passed, but small diff ranging from: %s to %s" % ( method.rjust(width), nmax, nmin)

    STP = SaTePr(gsw_cv.SA_from_SP, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)

    test_print(STP, "entropy", "entropy")
    test_print(STP, "rho", "rho")
    test_print(STP, "cp", "cp")
    test_print(STP, "helmholtz_energy", "Helmholtz_energy")
    test_print(STP, "internal_energy", "internal_energy")
    test_print(STP, "sound_speed", "sound_speed")
    test_print(STP, "kappa", "kappa")
    test_print(STP, "adiabatic_lapse_rate", "adiabatic_lapse_rate")
    test_print(STP, "chem_potential_relative", "chem_potential")
    test_print(STP, "specvol", "specvol")
    test_print(STP, "molality", "molality")
    test_print(STP, "ionic_strength", "ionic_strength")
    test_print(STP, "potential_t", "pt_from_t")
    test_print(STP, "potential_t", "pt0") #TODO: pt0_from_t
    test_print(STP, "conservative_t", "CT_from_t")
    test_print(STP, "enthalpy", "enthalpy")
    test_print(STP, "alpha", "alpha_wrt_t")