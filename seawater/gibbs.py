#TODO: Go over PDFs to improve documentation
#TODO: Examples: simple with the data range (copy-and-paste numbers) and complex (real data)
#A short demonstation of the GSW Oceanographic Toolbox now follows. The following vertical profile, from the North Pacific, is of Practical Salinity, SP, and in situ temperature, t, as a function of pressure, p,
#SP = [ 34.3454  34.5427  34.6289  34.6663  34.6839  34.6915  34.6914 ]
#t  = [ 27.9620   4.4726   2.1178   1.6031   1.4601   1.4753   1.5998 ]
#p  = [       0     1010     2025     3045     4069     5098     6131 ]
#TODO: Check original authors and dates
#TODO: csiro vs gibbs (table?)
#TODO: check_dim for p in all "p" functions
#FIXME: some function return values even with NaN in the input, check this behaviour (also present in the original).

import numpy as np
from seawater import constants as cte
from seawater import library as lib
from seawater import Temperature as temp

import numpy.ma as ma

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

        # broadcast p to SA dimension TODO
        #self.p = check_dim(self.p, self.SA)

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
             heat capacity of seawater [J kg :sup:`-1` K:sup:`-1`]

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
                           Helmholtz energy [J kg :sup:`-1`]

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
                             specific internal energy [J kg :sup:`-1`]

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
                      speed of sound in seawater [m s :sup:`-1`]

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

    def adiabatic_lapse_rate(self):
        """
        Calculates the adiabatic lapse rate of sea water.

        Returns
        -------
        adiabatic_lapse_rate : array_like
                               Adiabatic lapse rate [K Pa :sup:`-1`]

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
                                  relative chemical potential [J kg :sup:`-1`]

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
                  specific volume [m :sup:`3` kg :sup:`-1`]
                  #TODO: the original has this reversed

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

        pt0 = self.potential_t() # NOTE: pt0_from_t

        CT = temp.CT_from_pt(self.SA, pt0)

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
        self.masked_SA.fill_value = np.NaN #TODO: check if this affect other masks
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
                   specific enthalpy [J kg :sup:`-1`]

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

    def chem_potential_water(self):
        """
        Calculates the chemical potential of water in seawater.

        Returns
        -------
        chem_potential_water : array_like
                              chemical potential of water in seawater [J kg :sup:`-1`]

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
        >>> STP.chem_potential_water()
        array([[ -3819.11311835,   1279.4964416 ,  10748.01377942,  11057.29014031],
               [ -2292.68162515,   5095.80142795,   9352.59715652,  13679.8357573 ]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

        Modifications:
        2010-09-28. Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        SA = self.masked_SA.filled() # ensure that SA is non-negative

        x2 = cte.sfac * SA
        x = np.sqrt(x2)
        y = self.t * 0.025
        z = self.p * 1e-4 # Note that the input pressure (p) is sea pressure in units of dbar

        g03_g = 101.342743139674 + z * ( 100015.695367145 + \
            z * ( -2544.5765420363 + z * ( 284.517778446287 + \
            z * ( -33.3146754253611 + ( 4.20263108803084 - 0.546428511471039 * z ) * z ) ) ) ) + \
            y * ( 5.90578347909402 + z * ( -270.983805184062 + \
            z * ( 776.153611613101 + z * ( -196.51255088122 + ( 28.9796526294175 - 2.13290083518327 * z ) * z ) ) ) + \
            y * ( -12357.785933039 + z * ( 1455.0364540468 + \
            z * ( -756.558385769359 + z * ( 273.479662323528 + z * ( -55.5604063817218 + 4.34420671917197 * z ) ) ) ) + \
            y * ( 736.741204151612 + z * ( -672.50778314507 + \
            z * ( 499.360390819152 + z * ( -239.545330654412 + ( 48.8012518593872 - 1.66307106208905 * z ) * z ) ) ) + \
            y * ( -148.185936433658 + z * ( 397.968445406972 + \
            z * ( -301.815380621876 + ( 152.196371733841 - 26.3748377232802 * z ) * z ) ) + \
            y * ( 58.0259125842571 + z * ( -194.618310617595 + \
            z * ( 120.520654902025 + z * ( -55.2723052340152 + 6.48190668077221 * z ) ) ) + \
            y * ( -18.9843846514172 + y * ( 3.05081646487967 - 9.63108119393062 * z ) + \
            z * ( 63.5113936641785 + z * ( -22.2897317140459 + 8.17060541818112 * z ) ) ) ) ) ) ) )

        g08_g = x2 * ( 1416.27648484197 + \
            x * ( -2432.14662381794 + x * ( 2025.80115603697 + \
            y * ( 543.835333000098 + y * ( -68.5572509204491 + \
            y * ( 49.3667694856254 + y * ( -17.1397577419788 + 2.49697009569508 * y ) ) ) - 22.6683558512829 * z ) + \
            x * ( -1091.66841042967 - 196.028306689776 * y + \
            x * ( 374.60123787784 - 48.5891069025409 * x + 36.7571622995805 * y ) + 36.0284195611086 * z ) + \
            z * ( -54.7919133532887 + ( -4.08193978912261 - 30.1755111971161 * z ) * z ) ) + \
            z * ( 199.459603073901 + z * ( -52.2940909281335 + ( 68.0444942726459 - 3.41251932441282 * z ) * z ) ) + \
            y * ( -493.407510141682 + z * ( -175.292041186547 + ( 83.1923927801819 - 29.483064349429 * z ) * z ) + \
            y * ( -43.0664675978042 + z * ( 383.058066002476 + z * ( -54.1917262517112 + 25.6398487389914 * z ) ) + \
            y * ( -10.0227370861875 - 460.319931801257 * z + y * ( 0.875600661808945 + 234.565187611355 * z ) ) ) ) ) + \
            y * ( 168.072408311545 ) )

        g_SA_part = 8645.36753595126 + \
            x * ( -7296.43987145382 + x * ( 8103.20462414788 + \
            y * ( 2175.341332000392 + y * ( -274.2290036817964 + \
            y * ( 197.4670779425016 + y * ( -68.5590309679152 + 9.98788038278032 * y ) ) ) - 90.6734234051316 * z ) + \
            x * ( -5458.34205214835 - 980.14153344888 * y + \
            x * ( 2247.60742726704 - 340.1237483177863 * x + 220.542973797483 * y ) + 180.142097805543 * z ) + \
            z * ( -219.1676534131548 + ( -16.32775915649044 - 120.7020447884644 * z ) * z ) ) + \
            z * ( 598.378809221703 + z * ( -156.8822727844005 + ( 204.1334828179377 - 10.23755797323846 * z ) * z ) ) + \
            y * ( -1480.222530425046 + z * ( -525.876123559641 + ( 249.57717834054571 - 88.449193048287 * z ) * z ) + \
            y * ( -129.1994027934126 + z * ( 1149.174198007428 + z * ( -162.5751787551336 + 76.9195462169742 * z ) ) + \
            y * ( -30.0682112585625 - 1380.9597954037708 * z + y * ( 2.626801985426835 + 703.695562834065 * z ) ) ) ) ) + \
            y * ( 1187.3715515697959)

        chem_potential_water =  g03_g + g08_g  - 0.5 * cte.sfac * SA * g_SA_part

        return chem_potential_water

    def chem_potential_salt(self):
        """
        Calculates the chemical potential of salt in seawater.

        Returns
        -------
        chem_potential_salt : array_like
                              chemical potential of salt in seawater [J kg :sup:`-1`]

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
        >>> STP.chem_potential_salt()
        array([[ -3722.94417463,   1334.78497147,  10720.32688028,  11082.44431678],
               [ -2311.53902032,             nan,   9355.45225433,  13635.07649374]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

        Modifications:
        2010-09-28. Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        chem_potential_salt = self.chem_potential_relative() + \
                              self.chem_potential_water()

        return chem_potential_salt

    def isochoric_heat_cap(self):
        """
        Calculates the isochoric heat capacity of seawater.

        Returns
        -------
        isochoric_heat_cap : array_like
                             isochoric heat capacity [J kg :sup:`-1` K :sup:`-1`]

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
        >>> STP.isochoric_heat_cap()
        array([[ 3877.6517887 ,  3977.23066819,  4041.25381871,  3944.74636445],
               [ 4111.17043765,  4193.91078826,  4004.77189377,  4019.75411647]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.21.

        Modifications:
        2010-08-26. Trevor McDougall
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        g_tt = lib._gibbs(self.n0, self.n2, self.n0, self.SA, self.t, self.p)
        g_tp = lib._gibbs(self.n0, self.n1, self.n1, self.SA, self.t, self.p)
        g_pp = lib._gibbs(self.n0, self.n0, self.n2, self.SA, self.t, self.p)

        isochoric_heat_cap = -(cte.Kelvin + self.t) * (g_tt - g_tp * g_tp / g_pp)

        return isochoric_heat_cap

    def kappa(self):
        """
        Calculates the isentropic compressibility of seawater.

        Returns
        -------
        kappa : array_like
                Isentropic compressibility [Pa :sup:`-1`]

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

    def kappa_const_t(self):
        """
        Calculates isothermal compressibility of seawater at constant in-situ temperature.

        Returns
        -------
        kappa : array_like
                Isothermal compressibility [Pa :sup:`-1`]

        See Also
        --------
        TODO

        Notes
        -----
        This is the compressibility of seawater at constant in-situ temperature.

        Examples
        --------
        >>> from seawater.gibbs import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = [0., 500., 1500., 2000.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.kappa_const_t()
        array([[  4.31939544e-10,   4.32024995e-10,   4.30266156e-10,
                  4.08748100e-10],
               [  4.56941451e-10,   5.01665845e-10,   4.22886356e-10,
                  4.20872128e-10]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.15.1)

        Modifications:
        2010-07-23. David Jackett, Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        kappa = - lib._gibbs(self.n0, self.n0, self.n2, self.SA, self.t, self.p) / lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p)

        return kappa

    def osmotic_coefficient(self):
        """
        Calculates the osmotic coefficient of seawater.

        Returns
        -------
        osmotic_coefficient : array_like
                              osmotic coefficient of seawater [unitless]

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
        >>> STP.osmotic_coefficient()
        array([[ 0.90488718,  0.89901313,  0.90280557,  0.89943715],
               [ 0.90152697,         nan,  0.89931561,  0.90564689]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

        Modifications:
        2010-09-28. Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        M_S = 0.0314038218 # mole-weighted average atomic weight of the elements of sea salt, in units of kg mol :sup:`-1`

        # only >= than zero
        self.masked_SA.fill_value = np.NaN #TODO: check if this affect other masks
        SA = self.masked_SA.filled()

        molality = SA / ( M_S * ( 1000 - SA ) ) # molality of seawater in mol kg :sup:`-1`
        part = molality * cte.R * ( cte.Kelvin + self.t )

        #FIXME: the SAzero is needed to fixlib._gibbs: ValueError: shape mismatch: objects cannot be broadcast to a single shape
        SAzero = np.zeros( SA.shape )

        osmotic_coefficient = ( lib._gibbs(self.n0, self.n0, self.n0, SAzero, self.t, self.p) - \
        self.chem_potential_water() ) / part

        return osmotic_coefficient

    def pot_rho(self, pr=0):
        """
        Calculates potential density of seawater.

        Parameters
        ----------
        pr : int, float, optional
            reference pressure, default = 0

        Returns
        -------
        pot_rho : array_like
                  potential density  [kg m :sup:`-3`]

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
        >>> STP.pot_rho()
        array([[ 1041.77425464,  1022.03286026,  1005.35590628,  1009.95952733],
               [ 1006.74841976,   999.84434287,  1008.33162777,  1002.31311402]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 3.4.

        Modifications:
        2010-07-23. David Jackett, Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        pt = self.potential_t(pr=pr)
        #TODO: the line below recompute rho with pr=0 and pt
        pot_rho = 1. / lib._gibbs(self.n0, self.n0, self.n1, self.SA, pt, pr)

        return pot_rho

    def specvol_anom(self):
        """
        Calculates specific volume anomaly from Absolute Salinity, in-situ temperature and pressure, using the full TEOS-10 Gibbs function.

        The reference value of Absolute Salinity is SSO and the reference value of Conservative Temperature is equal to 0 :math:`^\\circ` C.

        Returns
        -------
        specvol_anom : array_like
                       specific volume anomaly  [m :sup:`3` kg :sup:`-1`]
                       #TODO: the original has this reversed

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
        >>> STP.pot_rho()
        array([[ 1041.77425464,  1022.03286026,  1005.35590628,  1009.95952733],
               [ 1006.74841976,   999.84434287,  1008.33162777,  1002.31311402]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (3.7.3)

        Modifications:
        2010-08-26. Trevor McDougall & Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        SSO = cte.SSO * np.ones( self.SA.shape )
        CT0 = np.zeros( SSO.shape )
        pr0 = np.zeros( SSO.shape )
        pt_zero = temp.pt_from_CT(SSO, CT0)
        #TODO: If I figure out a way to recalculate from self.potential_t() this call to Temperature.py will be unecessary
        t_zero = temp.potential_t(SSO, pt_zero, pr0, self.p)

        specvol_anom = lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p) - \
        lib._gibbs(self.n0, self.n0, self.n1, SSO, t_zero, self.p)

        return specvol_anom

    def alpha_wrt_t(self):
        """
        Calculates the thermal expansion coefficient of seawater with respect to in-situ temperature.

        Returns
        -------
        alpha_wrt_t : array_like
                      thermal expansion coefficient [K :sup:`-1`]

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
        >>> STP.alpha_wrt_t()
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

        alpha_wrt_t = lib._gibbs(self.n0, self.n1, self.n1, self.SA, self.t, self.p) / lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p)

        return alpha_wrt_t

    def alpha_wrt_CT(self):
        """
        Calculates the thermal expansion coefficient of seawater with respect to Conservative Temperature.

        Returns
        -------
        alpha_wrt_CT : array_like
                       thermal expansion coefficient [K :sup:`-1`]

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
        >>> STP.alpha_wrt_CT()
        array([[  1.58116830e-04,   2.12123759e-04,   2.53807153e-04,
                  3.44888273e-04],
               [  1.64597103e-04,  -4.64529201e-05,   2.85086325e-04,
                  3.03728728e-04]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.18.3).

        Modifications:
        2010-07-23. David Jackett, Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        pt0 = self.potential_t() # NOTE: pt0_from_t
        factor = -cte.cp0 / ( (cte.Kelvin + pt0) * lib._gibbs(self.n0, self.n2, self.n0, self.SA, self.t, self.p ) )
        alpha_wrt_CT = factor  * ( lib._gibbs(self.n0, self.n1, self.n1, self.SA, self.t, self.p) / lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p ) )

        return alpha_wrt_CT

    def alpha_wrt_pt(self):
        """
        Calculates the thermal expansion coefficient of seawater with respect to potential temperature, with a reference pressure of zero.

        Returns
        -------
        alpha_wrt_pt : array_like
                       thermal expansion coefficient [K :sup:`-1`]

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
        >>> STP.alpha_wrt_pt()
        array([[  1.54174741e-04,   2.13608621e-04,   2.62397019e-04,
                  3.52131126e-04],
               [  1.70265060e-04,  -4.91000706e-05,   2.92817943e-04,
                  3.14759131e-04]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.18.2).

        Modifications:
        2010-07-23. David Jackett, Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        pr0 = np.zeros( self.p.shape )
        pt0 = self.potential_t() # NOTE: pt0_from_t
        factor = lib._gibbs(self.n0, self.n2, self.n0, self.SA, pt0, pr0) / lib._gibbs(self.n0, self.n2, self.n0, self.SA, self.t, self.p)

        alpha_wrt_pt = factor * ( lib._gibbs(self.n0, self.n1, self.n1, self.SA, self.t, self.p) / lib._gibbs( self.n0, self.n0, self.n1, self.SA, self.t, self.p ) )

        return alpha_wrt_pt

    def beta_const_t(self):
        """
        Calculates the saline (i.e. haline) contraction coefficient of seawater at constant in-situ temperature.

        Returns
        -------
        beta_const_t : array_like
                       saline contraction coefficient [kg g :sup:`-1`]

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
        >>> STP.beta_const_t()
        array([[ 0.00076014,  0.00074453,  0.0007323 ,  0.0007157 ],
               [ 0.00075704,  0.00081627,  0.00072689,  0.00072291]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.19.1)

        Modifications:
        2010-07-23. David Jackett, Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        beta_const_t = -lib._gibbs(self.n1, self.n0, self.n1, self.SA, self.t, self.p) / lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p)

        return beta_const_t

    def beta_const_CT(self):
        """
        Calculates the saline (i.e. haline) contraction coefficient of seawater at constant Conservative Temperature.

        Returns
        -------
        beta_const_CT : array_like
                        saline contraction coefficient [kg g :sup:`-1`]

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
        >>> STP.beta_const_CT()
        array([[ 0.0007578 ,  0.00073925,  0.0007239 ,  0.00069986],
               [ 0.0007534 ,         nan,  0.00071632,  0.00071045]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.19.3)

        Modifications:
        2010-07-23. David Jackett, Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        pr0 = np.zeros( self.p.shape )

        pt0 = self.potential_t() # NOTE: pt0_from_t

        gp = lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p)

        factora = lib._gibbs(self.n1, self.n1, self.n0, self.SA, self.t, self.p) - \
        lib._gibbs(self.n1, self.n0, self.n0, self.SA, pt0, pr0) / (cte.Kelvin + pt0)

        factor = factora / ( gp * lib._gibbs(self.n0, self.n2, self.n0, self.SA, self.t, self.p) )

        beta_const_CT = lib._gibbs(self.n0, self.n1, self.n1, self.SA, self.t, self.p) * factor - \
        lib._gibbs(self.n1, self.n0, self.n1, self.SA, self.t, self.p) / gp

        return beta_const_CT

    def beta_const_pt(self):
        """
        Calculates the saline (i.e. haline) contraction coefficient of seawater at constant potential temperature with a reference pressure of 0 dbar.

        Returns
        -------
        beta_const_pt : array_like
                        saline contraction coefficient [kg g :sup:`-1`]

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
        >>> STP.beta_const_pt()
        array([[ 0.00076014,  0.0007444 ,  0.0007319 ,  0.00071523],
               [ 0.00075704,         nan,  0.00072649,  0.0007224 ]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.19.2)

        Modifications:
        2010-07-23. David Jackett, Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        pr0 = np.zeros( self.p.shape )

        pt0 = self.potential_t() # NOTE: pt0_from_t

        gp = lib._gibbs(self.n0, self.n0, self.n1, self.SA, self.t, self.p)

        factora = lib._gibbs(self.n1, self.n1, self.n0, self.SA, self.t, self.p) - \
        lib._gibbs(self.n1, self.n1, self.n0, self.SA, pt0, pr0)

        factor = factora / ( gp * lib._gibbs(self.n0, self.n2, self.n0, self.SA, self.t, self.p) )

        beta_const_pt = lib._gibbs(self.n0, self.n1, self.n1, self.SA, self.t, self.p) * factor - \
        lib._gibbs(self.n1, self.n0, self.n1, self.SA, self.t, self.p) / gp

        return beta_const_pt

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

    def test_print(STP, method, comp_value=None):
        """
        Run a function test mimicking the original logic. This is done to allow for a direct comparison of the result from the Matlab to the python package.
        """

        if comp_value is None:
            comp_value = method

        # test for floating differences with: computed - check_value >= defined_precision
        exec( "unequal = (gsw_cv." +comp_value+ " - STP." +method+ "() ) >= gsw_cv." +comp_value+ "_ca")

        width = 23
        if unequal.any():
            print "%s: Failed" % method.rjust(width)
        else:
            # test if check value is identical to computed value
            if eval( "( gsw_cv." +comp_value+ "[~np.isnan(gsw_cv."+comp_value+")] == STP." +method+ "()[~np.isnan(STP."+method+"())] ).all()" ):
                print "%s: Passed" % method.rjust(width)
            else:
                # test for differences in case their aren't equal. This is an attempt to place all tests together (i.e. term25 and small float differences that will appear)
                exec("nmax = ( gsw_cv."+comp_value+" - STP."+method+"() )[~np.isnan(gsw_cv."+comp_value+")].max()")
                exec("nmin = ( gsw_cv."+comp_value+" - STP."+method+"() )[~np.isnan(gsw_cv."+comp_value+")].min()")
                print "%s: Passed, but small diff ranging from: %s to %s" % ( method.rjust(width), nmax, nmin)

    STP = SaTePr(gsw_cv.SA_from_SP, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)

    test_print(STP, "entropy")
    test_print(STP, "rho")
    test_print(STP, "cp")
    test_print(STP, "helmholtz_energy", "Helmholtz_energy")
    test_print(STP, "internal_energy")
    test_print(STP, "sound_speed")
    test_print(STP, "adiabatic_lapse_rate")
    test_print(STP, "chem_potential_relative", "chem_potential")
    test_print(STP, "specvol")
    test_print(STP, "molality")
    test_print(STP, "ionic_strength")
    test_print(STP, "potential_t", "pt_from_t")
    test_print(STP, "potential_t", "pt0") #NOTE: pt0_from_t
    test_print(STP, "conservative_t", "CT_from_t")
    test_print(STP, "enthalpy")
    test_print(STP, "alpha_wrt_t")
    test_print(STP, "alpha_wrt_CT")
    test_print(STP, "alpha_wrt_pt")
    test_print(STP, "beta_const_t")
    test_print(STP, "beta_const_CT")
    test_print(STP, "beta_const_pt")
    test_print(STP, "chem_potential_water")
    test_print(STP, "chem_potential_salt")
    test_print(STP, "isochoric_heat_cap")
    test_print(STP, "kappa")
    test_print(STP, "kappa_const_t")
    test_print(STP, "osmotic_coefficient")
    test_print(STP, "pot_rho")
    test_print(STP, "specvol_anom")