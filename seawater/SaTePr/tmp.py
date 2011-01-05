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

        n0 = 0
        n1 = 1

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

        n0 = 0
        n1 = 1
        n2 = 2

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

        n0 = 0
        n1 = 1
        n2 = 2

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

        n0 = 0
        n1 = 1
        n2 = 2

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

        n0 = 0
        n1 = 1

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

        n0 = 0
        n1 = 1

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

        pt0 = pt0_from_t()
        CT = CT_from_pt(self.SA, pt0)

        return CT