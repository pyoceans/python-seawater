import numpy as np
from seawater import constants as cte
from seawater import library as lib

class DepthPressure:
    """
    """
    def __init__(self, lat, p=None, z=None):
        self.lat, self.p, self.z = np.asarray(lat), np.asarray(p), np.asarray(z)

        if p:
            self.z = z_from_p()
        elif z:
            self.p = p_from_z()
        else:
            raise NameError('need latitude plus pressure or depth')

    def  z_from_p(self):
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

        X     = np.sin( np.deg2rad(self.lat) )
        sin2  = X**2
        B     = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2 )
        A     = -0.5 * cte.gamma * B
        C     = lib._enthalpy_SSO_0_CT25(self.p)
        self.z     = -2 * C / ( B + np.sqrt( B**2 - 4 * A * C ) )

        return self.z

if __name__=='__main__':
    lat = 32.
    p = [0., 15., 100., 550., 1500., 2000., 3000., 5000., 10000.]
    dp = DepthPressure(lat=lat, p=p)