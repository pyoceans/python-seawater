import numpy as np
from seawater import constants as cte
from seawater import library as lib

class SaTePr():
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
    def __init__(self, SA, t, p): #TODO: pr=0
        # Convert input to numpy arrays
        self.SA, self.t, self.p = np.asarray(SA), np.asarray(t), np.asarray(p)

        if self.SA.shape:
            self.SA[self.SA < 0] = 0
        elif self.SA < 0:
            self.SA = 0

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

        s1 = self.SA*(35. / cte.SSO)

        pt0 = self.t + self.p * ( 8.65483913395442e-6  - \
              s1 * 1.41636299744881e-6 - \
              self.p * 7.38286467135737e-9 + \
              self.t * ( -8.38241357039698e-6 + \
              s1 * 2.83933368585534e-8 + \
              self.t * 1.77803965218656e-8 + \
              self.p * 1.71155619208233e-10 ) )

        dentropy_dt = cte.cp0 / ( (273.15 + pt0) * ( 1 - 0.05 * ( 1 - self.SA / cte.SSO ) ) )

        true_entropy_part = lib._entropy_part(self.SA, self.t, self.p)

        for Number_of_iterations in range(0,2,1):
            pt0_old = pt0
            dentropy = lib._entropy_part_zerop(self.SA, pt0_old) - true_entropy_part
            pt0 = pt0_old - dentropy / dentropy_dt # this is half way through the modified method
            pt0m = 0.5 * (pt0 + pt0_old);
            dentropy_dt = -lib._gibbs_pt0_pt0(self.SA, pt0m)
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

        n0 = 0
        n1 = 1
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

        n0 = 0
        n1 = 1
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

        n0 = 0
        n2 = 2
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

        n0 = 0
        n1 = 1

        helmholtz_energy = lib._gibbs(n0, n0, n0, self.SA, self.t, self.p) - \
                        ( cte.db2Pascal * self.p + 101325 ) * lib._gibbs(n0, n0, n1, self.SA, self.t, self.p)

        return helmholtz_energy

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

        width = 17
        if inequal.any():
            print "%s: Failed" % method.rjust(width)
        else:
            if eval( "( gsw_cv." +comp_value+ "[~np.isnan(gsw_cv."+comp_value+")] == STP." +method+ "()[~np.isnan(STP."+method+"())] ).all()" ):
                print "%s: Passed, equal" % method.rjust(width)
            else:
                exec("nmax = STP." +method+ "()[~np.isnan(STP."+method+"())].max()")
                print "%s: Passed, but small diff %s" % ( method.rjust(width), nmax)

    STP = SaTePr(gsw_cv.SA_from_SP, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)

    #methods_list = ["pt0_from_t", "entropy", "rho", "cp" , "helmholtz_energy"]
    test_print(STP, "pt0_from_t", "pt0")
    test_print(STP, "entropy", "entropy") #FIXME: pass, but small float are detected, investigate further
    test_print(STP, "rho", "rho")
    test_print(STP, "cp", "cp")
    test_print(STP, "helmholtz_energy", "Helmholtz_energy")