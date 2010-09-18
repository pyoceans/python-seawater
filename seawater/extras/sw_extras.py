# -*- coding: utf-8 -*-

"""
Extra seawater functions
===========================
"""

import numpy as np
import seawater.csiro as sw

def sigma_t(s, t, p):
    """
    :math:`\\sigma_{t}` is the remainder of subtracting 1000 kg m :sup:`-3` from the density of a sea water sample at atmospheric pressure.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"

    Returns
    -------
    sgmt : array_like
           density  [kg m :sup:`3`]

    See Also
    --------
    dens, sigmatheta

    Notes
    -----
    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22

    >>> import numpy as np
    >>> import seawater.csiro as sw
    >>> from seawater.csiro import T68conv
    >>> import seawater.extras.sw_extras as swe
    >>> s = np.array([0, 0, 0, 0, 35, 35, 35, 35])
    >>> t = np.array([0, 0, 30, 30, 0, 0, 30, 30]) / T68conv
    >>> p = np.array([0, 10000, 0, 10000, 0, 10000, 0, 10000])
    >>> swe.sigma_t(s, t, p)
    array([ -0.157406  ,  45.33710972,  -4.34886626,  36.03148891,
            28.10633141,  70.95838408,  21.72863949,  60.55058771])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for computation of fundamental properties of seawater. UNESCO Tech. Pap. in Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39. http://www.scor-int.org/Publications.htm

    .. [2] Millero, F.J., Chen, C.T., Bradshaw, A., and Schleicher, K. A new high pressure equation of state for seawater. Deap-Sea Research., 1980, Vol27A, pp255-264. doi:10.1016/0198-0149(80)90016-3

    Modifications: Filipe Fernandes, 2010
                   10-01-26. Filipe Fernandes, first version.
    """

    sgmt = sw.dens(s, t, p) - 1000.0
    return sgmt

def sigmatheta(s, t, p, pr=0):
    """
    :math:`\\sigma_{\\theta}` is a measure of the density of ocean water where the quantity :math:`\\sigma_{t}` is calculated using the potential temperature (:math:`\\theta`) rather than the in situ temperature and potential density of water mass relative to the specified reference pressure.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"
    pr : number
         reference pressure [db], default = 0

    Returns
    -------
    sgmte : array_like
           density  [kg m :sup:`3`]

    See Also
    --------
    dens, sigma_t

    Notes
    -----
    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22

    >>> import numpy as np
    >>> import seawater.extras.sw_extras as swe
    >>> from seawater.csiro import T68conv
    >>> s = np.array([0, 0, 0, 0, 35, 35, 35, 35])
    >>> t = np.array([0, 0, 30, 30, 0, 0, 30, 30]) / T68conv
    >>> p = np.array([0, 10000, 0, 10000, 0, 10000, 0, 10000])
    >>> swe.sigmatheta(s, t, p)
    array([ -0.157406  ,  -0.20476006,  -4.34886626,  -3.63884068,
            28.10633141,  28.15738545,  21.72863949,  22.59634627])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for computation of fundamental properties of seawater. UNESCO Tech. Pap. in Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39. http://www.scor-int.org/Publications.htm

    .. [2] Millero, F.J., Chen, C.T., Bradshaw, A., and Schleicher, K. A new high pressure equation of state for seawater. Deap-Sea Research., 1980, Vol27A, pp255-264. doi:10.1016/0198-0149(80)90016-3

    Modifications: Filipe Fernandes, 2010
                   10-01-26. Filipe Fernandes, first version.
    """

    sgmte = sw.pden(s, t, p, pr) - 1000.0
    return sgmte

def N(bvfr2):
    """
    Buoyancy frequency is the frequency with which a parcel or particle of fluid displaced a small vertical distance from its equilibrium position in a stable environment will oscillate. It will oscillate in simple harmonic motion with an angular frequency defined by

    .. math:: N = \\left(\\frac{-g}{\\sigma_{\\theta}} \\frac{d\\sigma_{\\theta}}{dz}\\right)^{2}

    Parameters
    ----------
    n2 : array_like
         Brünt-Väisälä Frequency squared [s :sup:`-2`]

    Returns
    -------
    n : array_like
        Brünt-Väisälä Frequency not-squared [s :sup:`-1`]

    Examples
    --------
    >>> import numpy as np
    >>> import seawater.csiro as sw
    >>> import seawater.extras.sw_extras as swe
    >>> s = np.array([[0, 0, 0], [15, 15, 15], [30, 30, 30],[35,35,35]])
    >>> t = np.repeat(15, s.size).reshape(s.shape)
    >>> p = np.array([0, 250, 500, 1000])
    >>> lat = np.array([30,32,35])
    >>> swe.N(sw.bfrq(s, t, p, lat)[0])
    array([[ 0.02124956,  0.02125302,  0.02125843],
           [ 0.02110919,  0.02111263,  0.02111801],
           [ 0.00860812,  0.00860952,  0.00861171]])

    References
    ----------
    .. [1] A.E. Gill 1982. p.54  eqn 3.7.15 "Atmosphere-Ocean Dynamics" Academic Press: New York. ISBN: 0-12-283522-0

    .. [2] Jackett, David R., Trevor J. Mcdougall, 1995: Minimal Adjustment of Hydrographic Profiles to Achieve Static Stability. J. Atmos. Oceanic Technol., 12, 381-389. doi: 10.1175/1520-0426(1995)012<0381:MAOHPT>2.0.CO;2

    Modifications: Filipe Fernandes, 2010
                   10-01-26. Filipe Fernandes, first version.
    """

    bvfr  = np.sqrt( np.abs( bvfr2 ) ) * np.sign( bvfr2 )
    return bvfr

def shear(p, u, v=0):
    """
    Calculates the vertical shear for u, v velocity section.

    .. math::
        \\textrm{shear} = \\frac{\\partial (u^2 + v^2)^{0.5}}{\partial z}

    Parameters
    ----------
    p : array_like
        pressure [db]. The shape can be "broadcasted"
    u(p) : array_like
           Eastward velocity [m s :sup:`-1`]
    v(p) : array_like
           Northward velocity [m s :sup:`-1`]

    Returns
    -------
    shr : array_like
          frequency [s :sup:`-1`]
    p_ave : array_like
            mid pressure between p grid (M-1xN)  [db]

    See Also
    --------
    TODO

    Notes
    -----
    TODO check where depth increases to find dimension

    Examples
    --------
    >>> import numpy as np
    >>> import seawater.csiro as sw
    >>> import seawater.extras.sw_extras as swe
    >>> p = np.array([0, 250, 500, 1000])
    >>> vel = np.array([[0.5, 0.5, 0.5], [0.15, 0.15, 0.15], [0.03, 0.03, .03],[0.,0.,0.]])
    >>> swe.shear(p, vel)[0]
    array([[ -1.40000000e-03,  -1.40000000e-03,  -1.40000000e-03],
           [ -4.80000000e-04,  -4.80000000e-04,  -4.80000000e-04],
           [ -6.00000000e-05,  -6.00000000e-05,  -6.00000000e-05]])

    AUTHOR:  Filipe Fernandes, 2010

    MODIFICATIONS:
    10-01-28. Filipe Fernandes, first version.
    """
    # if pressure is a vector make it a array of the same size as t/s
    if p.ndim == 1:
        p = np.repeat(p[np.newaxis,:], u.shape[1], axis=1).reshape(u.shape)

    m,n      = p.shape
    iup      = np.arange(0,m-1)
    ilo      = np.arange(1,m)
    p_ave    = ( p[iup,:] + p[ilo,:] )/2.
    vel      = np.sqrt( u**2 + v**2 )
    diff_vel = np.diff( vel, axis=0 )
    diff_z   = np.diff(   p, axis=0 ) # TODO to Z ?
    shr      = diff_vel / diff_z

    return shr, p_ave

def richnumb(n, s):
    """
    Calculates  the ratio of buoyancy to inertial forces which measures the stability of a fluid layer.
    this functions computes the gradient Richardson number in the form of:

    .. math::
        Ri = \\frac{N^2}{S^2}

    Representing a dimensionless number that expresses the ratio of the energy extracted by buoyancy forces to the energy gained from the shear of the large-scale velocity field.

    Parameters
    ----------
    n : array_like
        Brünt-Väisälä [s :sup:`-1`]
    shr : array_like
          shear [s :sup:`-1`]

    Returns
    -------
    ri : array_like
         non-dimensional

    Examples
    --------
    TODO: check the example and add real values
    >>> import numpy as np
    >>> import seawater.extras.sw_extras as swe
    >>> import seawater.csiro as sw
    >>> s   = np.array([[0, 0, 0], [15, 15, 15], [30, 30, 30],[35,35,35]])
    >>> t   = np.repeat(15, s.size).reshape(s.shape)
    >>> p   = np.array([0, 250, 500, 1000])
    >>> lat = np.array([30,32,35])
    >>> n   = swe.N(sw.bfrq(s, t, p, lat)[0])
    >>> vel = np.array([[0.5, 0.5, 0.5], [0.15, 0.15, 0.15], [0.03, 0.03, .03],[0.,0.,0.]])
    >>> s   = swe.shear(p, vel)[0]
    >>> swe.richnumb(n, s)
    array([[   230.37941215,    230.45444299,    230.57181258],
           [  1934.01949759,   1934.64933431,   1935.63457818],
           [ 20583.24410868,  20589.94661835,  20600.43125069]])

    Modifications: Filipe Fernandes, 2010
                   10-01-26. Filipe Fernandes, first version.
    """

    n2 = n**2 * np.sign(n)
    s2 = s**2
    ri = n2 / s2
    return ri

def inertial_period(lat):
    """
    Calculate the inertial period as:

    .. math::
        Ti = \\frac{2\\pi}{f} = \\frac{T_{sd}}{2\\sin\\phi}

    Parameters
    ----------
    lat : array_like
          latitude in decimal degress north [-90..+90]

    Returns
    -------
    Ti : array_like
         period in hours

    Examples
    --------
    >>> import seawater.extras.sw_extras as swe
    >>> lat = 30
    >>> swe.inertial_period(lat)
    23.934849862785651

    Modifications: Filipe Fernandes, 2010
                   10-01-26. Filipe Fernandes, first version.
    """

    Ti = 2*np.pi / sw.cor(lat)/3600

    return Ti


def get_wavenum(T, h, L, thetao, Ho):

    """
    Solves the wave dispersion relationship

    .. math::
        \Omega = ...

    Parameters
    ----------
    T : array_like
        Wave period in seconds
    TODO: h -> meters
    TODO: L -> meters
    TODO: thetao
    TODO: Ho

    Returns
    -------
    omega : array_like
            Wave frequency
    TODO: hoLo, hoL, Lo, L, k, T, Co, C, Cg, G, Ks, Kr, theta, H

    Examples
    --------
    >>> import seawater.extras.sw_extras as swe
    >>> TODO: use [hoLo,hoL,Lo,L,k,omega,T,Co,C,Cg,G,Ks,Kr,theta, H]= get_wavenum(T,h,L,thetao, Ho)

    Modifications: Filipe Fernandes, 2010
                   10-01-26. Filipe Fernandes, first version.
    """

    if isnan(L)
        omega = 2*pi/T
        Lo    = (g.*T.^2)./2./pi;
        # returns wavenumber of the gravity wave dispersion relation using newtons method
        k = omega / sqrt(g) # the initial guess will be the shallow water wavenumber
        f = g *k*tanh(k*h) - omega**2

        while max(abs(f))>1e-10
            dfdk = g*k*h*( sech(k*h) )**2 + g*tanh(k*h)
            k    = k - f/dfdk
            f    = g*k*tanh(k*h) - omega**2

        L = 2*pi/k
    else
        Lo    = L/tanh(2*pi*h/L)
        k     = 2*pi/L
        T     = sqrt(2*pi*Lo/g)
        omega = 2*pi/T

    hoL   = h/L
    hoLo  = h/Lo
    C     = omega/k
    Co    = Lo/T
    G     = 2*k*h/sinh(2*k*h)
    n     = (1+G)/2
    Cg    = n*C
    Ks    = sqrt( 1/(1+G) / tanh(k*h) )

    if isnan(thetao)
        theta = NaN
        Kr    = NaN
    else
        theta = asin(C/Co * sin(thetao*pi/180) ) * 180/pi
        Kr = sqrt( cos(thetao*pi/180) / cos(theta*pi/180) )

    if isnan(Ho)
        H = NaN
    else
        H = Ho*Ks*Kr

    return hoLo, hoL, Lo, L, k, omega, T, Co, C, Cg, G, Ks, Kr, theta, H