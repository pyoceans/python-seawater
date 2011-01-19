import numpy as np
from seawater import constants as cte
from seawater import library as lib

class DepthPressure:
    """
    Class to convert pressure to depth and vice-versa:

    This class will (maybe) replace z_from_p, p_from_z, and grav functions
    """
    def __init__(self, lat, p=None, z=None):
        self.lat, self.p, self.z = np.asarray(lat), np.asarray(p), np.asarray(z)

        X     = np.sin( np.deg2rad(self.lat) )
        sin2   = X**2
        gs = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2)

        if p:
            A      = -0.5 * cte.gamma * gs
            C      = lib._enthalpy_SSO_0_CT25(self.p)
            self.z = -2 * C / ( gs + np.sqrt( gs**2 - 4 * A * C ) )
        elif z:
            # get the first estimate of p from Saunders (1981)
            c1 =  5.25e-3 * sin2 + 5.92e-3
            p  = -2 * self.z / ( ( 1-c1 ) + np.sqrt( ( 1-c1 ) * ( 1-c1 ) + 8.84e-6 * self.z ) )
            # end of the first estimate from Saunders (1981)
            # initial value of the derivative of f
            df_dp = cte.db2Pascal * lib._specvol_SSO_0_CT25(p)
            f     = lib._enthalpy_SSO_0_CT25(p) + gs * ( self.z - 0.5 * cte.gamma * ( self.z**2 ) )
            p_old = p
            p     = p_old - f / df_dp
            pm    = 0.5 * ( p + p_old )
            df_dp = cte.db2Pascal * lib._specvol_SSO_0_CT25(pm)
            self.p     = p_old - f / df_dp
        else:
            raise NameError('need latitude (mandatory) and pressure or depth')

        self.grav = gs * (1 - cte.gamma * self.z)


if __name__=='__main__':
    lat = 32.
    p = [0., 15., 100., 550., 1500., 2000., 3000., 5000., 10000.]
    dp = DepthPressure(lat=lat, p=p)
    z = [-0., -14.89499448, -99.27948265, -545.44412444, -1484.209721, -1976.61994868, -2958.05761312, -4907.87772419, -9712.16369644]
    pd = DepthPressure(lat=lat, z=z)