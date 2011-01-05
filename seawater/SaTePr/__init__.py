import numpy as np
import numpy.ma as ma
from seawater import constants as cte
from seawater import library as lib
import seawater.gibbs as gsw

class SaTePr:
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

    def pt0_from_t(self):
        """
        Calculates potential temperature with reference pressure, pr = 0 dbar. The present routine is computationally faster than the more general function "pt_from_t(SA, t, p, pr)" which can be used for any reference pressure value.

        Returns
        -------
        pt0 : array_like
              potential temperature relative to 0 db [:math:`^\\circ` C (ITS-90)]

        See Also
        --------
        TODO

        Notes
        -----
        TODO

        Examples
        --------
        >>> from seawater.SaTePr import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = 900
        >>> STP = SaTePr(SA, t, p)
        >>> STP.pt0_from_t()
        array([[  4.89971486e+00,   1.48664023e+01,   2.18420392e+01,
                  3.17741959e+01],
               [  1.48891940e+01,   2.95267636e-02,   2.48187231e+01,
                  2.78058513e+01]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 3.1.

        .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater.  To be submitted to Ocean Science Discussions.

        Modifications:
        2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        SA = self.masked_SA.filled() # ensure that SA is non-negative

        s1 = SA * (35. / cte.SSO)

        pt0 = self.t + self.p * ( 8.65483913395442e-6  - \
              s1 * 1.41636299744881e-6 - \
              self.p * 7.38286467135737e-9 + \
              self.t * ( -8.38241357039698e-6 + \
              s1 * 2.83933368585534e-8 + \
              self.t * 1.77803965218656e-8 + \
              self.p * 1.71155619208233e-10 ) )

        dentropy_dt = cte.cp0 / ( (273.15 + pt0) * ( 1 - 0.05 * ( 1 - SA / cte.SSO ) ) )

        true_entropy_part = lib._entropy_part(SA, self.t, self.p)

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
        >>> from seawater.SaTePr import SaTePr
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

        n0, n1 = 0, 1

        entropy = -lib._gibbs(n0, n1, n0, self.SA, self.t, self.p)

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
        >>> from seawater.SaTePr import SaTePr
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

        n0, n1 = 0, 1

        rho = 1. / lib._gibbs(n0, n0, n1, self.SA, self.t, self.p)

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
        >>> from seawater.SaTePr import SaTePr
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

        n0, n2 = 0, 2

        cp = -( self.t + cte.Kelvin ) * lib._gibbs(n0, n2, n0, self.SA, self.t, self.p)

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
        >>> from seawater.SaTePr import SaTePr
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

        n0, n1 = 0, 1

        helmholtz_energy = lib._gibbs(n0, n0, n0, self.SA, self.t, self.p) - \
                        ( cte.db2Pascal * self.p + 101325 ) * lib._gibbs(n0, n0, n1, self.SA, self.t, self.p)

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
        >>> from seawater.SaTePr import SaTePr
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

        n0, n1 = 0, 1

        internal_energy = lib._gibbs(n0, n0, n0, self.SA, self.t, self.p) - \
                        (cte.Kelvin + self.t) * lib._gibbs(n0, n1, n0, self.SA, self.t, self.p) - \
                        (cte.db2Pascal * self.p + 101325) * lib._gibbs(n0, n0, n1, self.SA, self.t, self.p)

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
        >>> from seawater.SaTePr import SaTePr
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

        n0, n1, n2  = 0, 1, 2

        g_tt = lib._gibbs(n0, n2, n0, self.SA, self.t, self.p)
        g_tp = lib._gibbs(n0, n1, n1, self.SA, self.t, self.p)

        sound_speed = lib._gibbs(n0, n0, n1, self.SA, self.t, self.p) * \
        np.sqrt( g_tt / ( g_tp * g_tp - g_tt * lib._gibbs(n0, n0, n2, self.SA, self.t, self.p ) ) )

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
        >>> from seawater.SaTePr import SaTePr
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

        n0, n1, n2 = 0, 1, 2

        g_tt = lib._gibbs(n0, n2, n0, self.SA, self.t, self.p)
        g_tp = lib._gibbs(n0, n1, n1, self.SA, self.t, self.p)

        kappa = ( g_tp * g_tp - g_tt * lib._gibbs(n0, n0, n2, self.SA, self.t, self.p) ) / ( lib._gibbs(n0, n0, n1, self.SA, self.t, self.p ) * g_tt)

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
        >>> from seawater.SaTePr import SaTePr
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

        n0, n1, n2 = 0, 1, 2

        adiabatic_lapse_rate = - lib._gibbs(n0, n1, n1, self.SA, self.t, self.p) / ( lib._gibbs(n0, n2, n0, self.SA, self.t, self.p ) )

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
        >>> from seawater.SaTePr import SaTePr
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

        n0, n1 = 0, 1

        chem_potential_relative = lib._gibbs(n1, n0, n0, self.SA, self.t, self.p)

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
        >>> from seawater.SaTePr import SaTePr
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

        n0, n1= 0, 1

        specvol = lib._gibbs(n0, n0, n1, self.SA, self.t, self.p)

        return specvol

    def CT_from_t(self):
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
        >>> from seawater.SaTePr import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = 900
        >>> STP = SaTePr(SA, t, p)
        >>> STP.CT_from_t()
        array([[  4.66028901,  14.98237022,  22.6558658 ,  32.47483113],
               [ 15.46594688,   0.04649395,  25.55437701,  28.90014276]])

        References
        ----------
        .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 3.3.

        Modifications:
        2010-08-26. David Jackett, Trevor McDougall and Paul Barker
        2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
        """

        pt0 = self.pt0_from_t()

        CT = gsw.CT_from_pt(self.SA, pt0)

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
        >>> from seawater.SaTePr import SaTePr
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
        >>> from seawater.SaTePr import SaTePr
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

    def pt_from_t(self, pr=0):
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
        >>> from seawater.SaTePr import SaTePr
        >>> SA = [[53., 30, 10., 20.],[10., -5., 15., 8.]]
        >>> t = [[5., 15., 22., 32.],[15., 0., 25., 28.]]
        >>> p = 900
        >>> STP = SaTePr(SA, t, p)
        >>> STP.pt_from_t()
        array([[  4.89971486e+00,   1.48664023e+01,   2.18420392e+01,
                  3.17741959e+01],
               [  1.48891940e+01,   2.95267636e-02,   2.48187231e+01,
                  2.78058513e+01]])
        >>> STP.pt_from_t(pr = 900)
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

        n0, n2 = 0, 2

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
            dentropy_dt = -lib._gibbs(n0, n2, n0, SA, ptm, pr)
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
        >>> from seawater.SaTePr import SaTePr
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

        n0, n1 = 0, 1

        enthalpy = lib._gibbs(n0, n0, n0, self.SA, self.t, self.p) - ( self.t + cte.Kelvin ) * lib._gibbs(n0, n1, n0, self.SA, self.t, self.p)

        return enthalpy

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
        """
        exec( "inequal = (gsw_cv." +comp_value+ " - STP." +method+ "() ) >= gsw_cv." +comp_value+ "_ca")

        width = 23
        if inequal.any():
            print "%s: Failed" % method.rjust(width)
        else:
            if eval( "( gsw_cv." +comp_value+ "[~np.isnan(gsw_cv."+comp_value+")] == STP." +method+ "()[~np.isnan(STP."+method+"())] ).all()" ):
                print "%s: Passed, equal" % method.rjust(width)
            else:
                exec("nmax = STP." +method+ "()[~np.isnan(STP."+method+"())].max()")
                print "%s: Passed, but small diff %s" % ( method.rjust(width), nmax)

    STP = SaTePr(gsw_cv.SA_from_SP, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)

    test_print(STP, "pt0_from_t", "pt0")
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
    test_print(STP, "CT_from_t", "CT_from_t")
    test_print(STP, "molality", "molality")
    test_print(STP, "ionic_strength", "ionic_strength")
    test_print(STP, "pt_from_t", "pt_from_t")
    test_print(STP, "enthalpy", "enthalpy")