# -*- coding: utf-8 -*-
import numpy as np

"""
Sea presure is absolute pressure - 10.1325 dbar (or minus atmospheric pressure)
"""

"""
The Gibbs function approach allows the calculation of internal energy, entropy, enthalpy, potential enthalpy and the chemical potentials of seawater as well as the freezing temperature, and the latent heats of freezing and of evaporation. These quantities were not available from EOS-80 but are essential for the accurate accounting of heat in the ocean and for the consistent and accurate treatment of air-sea and ice-sea heat fluxes.
"""

"""
For the first time the influence of the spatially varying composition of seawater can be systematically taken into account through the use of Absolute Salinity. In the open ocean, this has a non-trivial effect on the horizontal density gradient, and thereby on the ocean velocities and heat transports calculated via the thermal wind relationship.
"""

"""
The new salinity variable, Absolute Salinity, is measured in SI units (e.g. g/kg).
"""

"""
The thermodynamic quantities available from TEOS-10 are totally consistent with each other, while this was not the case with EOS-80
"""

"""
The new equation incorporates a more accurate representation of salinity known as Absolute Salinity. The main justification for preferring Absolute Salinity over the most recent salinity definition, Practical Salinity, is that the thermodynamic properties of seawater are directly influenced by the total mass of dissolved constituents (Absolute Salinity). However, the mass of dissolved constituents are regionally variable and are not always accurately represented when using conductivity measurements of seawater, the key parameter used in the calculation of Practical Salinity.
"""

"""
In addition, it would be useful to have a conservative tracer for salinity. The Absolute Salinity S_A is not conservative, because it slowly increases in the ocean due to biogeochemical processes, which remineralize carbon and nutrients. The Reference Salinity S_R is also not conservative, since it is measured with conductivity which will also slowly increase due to these increases in concentrations of nutrient and carbon system ions.

However, we can "remove" these effects and construct a conservative salinity tracer S_*..

Typically this PREFORMED SALINITY S_* will equal S_R (and S_A) in Standard Seawater (SSW), but as ions are added it will become smaller than either S_R and S_A (as illustrated by the example in Table 1).
"""

"""
CT == \BigTheta = frac{h_o}{C^o_p}.

The scale factor C^o_p is a constant carefully chosen so that potential temperature \smalltheta and Conservative Temperature \BigTheta will be numerically similar for typical seawaters at SP = 35, or near t = 0degC. However, the difference between the two can exceed 1◦C when salinities are low and temperatures high (for details, see IOC et al., 2010).
"""

def _specvol_SSO_0_CT25(p):
    """
    Calcualtes specifc volume at the Standard Ocean Salinty (SSO) and Conservative Temperature of zero degrees C (CT=0), as a function of pressure (p [db]). It uses a streamlined version of the 25-term CT version of specific volume ( _rho_alpha_beta_CT25(SA,CT,p) ).

    Parameters
    ----------
    p : array_like
        pressure [db]

    Returns
    -------
    specvol_SSO_0_CT25 : array_like
                         Specific volume at (SSO, CT=0, p), 25-term equation.

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    TODO

    References
    ----------
    TODO

    Modifications
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    p = np.asarray(p)

    SSO = 35.16504 * np.ones( p.shape )
    specvol_SSO_0_CT25 = (1.00000000e+00 + SSO * ( 2.0777716085618458e-003 +np.sqrt(SSO) * 3.4688210757917340e-006) + p * 6.8314629554123324e-006) / (9.9984380290708214e+002 + SSO * ( 2.8925731541277653e+000 + SSO * 1.9457531751183059e-003) + p * ( 1.1930681818531748e-002 + SSO * 5.9355685925035653e-006 + p * -2.5943389807429039e-008) )

    return specvol_SSO_0_CT25

def  _enthalpy_SSO_0_CT25(p):
    """
     Calcualtes enthalpy at the Standard Ocean Salinty (SSO) and at a Conservative Temperature of zero degrees C (CT=0), as a function of pressure (p [db]). It Uses a streamlined version of the 25-term CT version of the Gibbs function ( _enthalpy_CT25(SA,CT,p) )

    Parameters
    ----------
    p : array_like
        pressure [db]

    Returns
    -------
    enthalpy_CT25 : array_like
                    enthalpy_CT25 at (SSO, CT = 0, p), 25-term equation.

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    TODO

    References
    ----------
    TODO

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    p = np.asarray(p)

    SSO = 35.16504 * np.ones( p.shape )

    a0 = 1 + SSO * (2.0777716085618458e-3 + np.sqrt(SSO) * 3.4688210757917340e-6)
    a1 = 6.8314629554123324e-6
    b0 = 9.9984380290708214e2 + SSO * (2.8925731541277653e0 + SSO * 1.9457531751183059e-3)
    b1 = 0.5 * (1.1930681818531748e-2 + SSO * 5.9355685925035653e-6)
    b2 = -2.5943389807429039e-8
    A = b1 - np.sqrt(b1**2 - b0 * b2)
    B = b1 + np.sqrt(b1**2 - b0 * b2)

    part = ( a0 * b2 - a1 * b1) / (b2 * (B - A) )

    enthalpy_SSO_0_CT25 = db2Pascal * ( ( a1 / (2*b2) ) * np.log( 1 + p * ( 2 * b1 + b2 * p ) / b0 ) + part * np.log( 1 + ( b2 * p * (B - A) ) / (A * (B + b2 * p ) ) ) )

    return enthalpy_SSO_0_CT25

def  _gibbs_pt0_pt0(SA, pt0):
    """
    Calculates the second derivative of the specific Gibbs function with respect to temperature at zero sea pressure.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :math:`-1`]
    pt0 : array_like
          potential temperature relative to 0 db [:math:`^\\circ` C (ITS-90)]

    Returns
    -------
    gibbs_pt0_pt0 : array_like
                    TODO: write the eq for the second derivative of the specific Gibbs function.

    See Also
    --------
    pt_from_CT, pt0_from_t

    Notes
    -----
    TODO

    Examples
    --------
    TODO

    References
    ----------
    TODO

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, pt0 = np.asarray(SA), np.asarray(pt0)

    # Ensure that SP is non-negative. FIXME: Shouldn't it be SA?
    SA[SA < 0] = 0
    sfac = 0.0248826675584615 # sfac = 1 / ( 40 * ( 35.16504 / 35 ) )

    x2 = sfac * SA
    x = np.sqrt(x2)
    y = pt0 * 0.025 # FIXME: 0.025d0 -> decimal.Decimal('0.025') ?

    g03 = -24715.571866078 + \
    y * ( 4420.4472249096725 + \
    y * ( -1778.231237203896 + \
    y * ( 1160.5182516851419 + \
    y * ( -569.531539542516 + y * 128.13429152494615) ) ) )

    g08 = x2 * ( 1760.062705994408 + x * ( -86.1329351956084 + \
    x * ( -137.1145018408982 + y * ( 296.20061691375236 + \
    y * ( -205.67709290374563 + 49.9394019139016 * y ) ) ) + \
    y * ( -60.136422517125 + y * 10.50720794170734 ) ) + \
    y * ( -1351.605895580406 + y * ( 1097.1125373015109 +  \
    y * ( -433.20648175062206 + 63.905091254154904 * y ) ) ) )

    gibbs_pt0_pt0 = ( g03 + g08 ) * 0.000625

    return gibbs_pt0_pt0

def  _entropy_part_zerop(SA, pt0):
    """
    Calculates entropy at a sea surface (p = 0 db), except that it does not evaluate any terms that are functions of Absolute Salinity alone. By not calculating these terms, which are a function only of Absolute Salinity, several unnecessary computations are avoided (including saving the computation of a natural logarithm). These terms are a necessary part of entropy, but are not needed when calculating potential temperature from in-situ temperature.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :math:`-1`]
    pt0 : array_like
          potential temperature relative to 0 db [:math:`^\\circ` C (ITS-90)]

    Returns
    -------
    entropy_part_zerop : array_like
                    TODO: [J/(K.mol)]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    TODO

    References
    ----------
    TODO

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, pt0 = np.asarray(SA), np.asarray(pt0)

    # Ensure that SA is non-negative
    SA[SA < 0] = 0
    sfac = 0.0248826675584615 # sfac = 1 / ( 40 * ( 35.16504 / 35 ) )

    x2 = sfac * SA
    x = np.sqrt(x2)
    y = pt0 * 0.025 # FIXME: 0.025d0 -> decimal.Decimal('0.025') ?

    g03 = y * ( -24715.571866078 + y * ( 2210.2236124548363 + \
    y * ( -592.743745734632 + y * ( 290.12956292128547 + \
    y * ( -113.90630790850321 + y * 21.35571525415769) ) ) ) )

    g08 = x2 * ( x * ( x * ( y * ( -137.1145018408982 + y * ( 148.10030845687618 + \
    y * ( -68.5590309679152 + 12.4848504784754 * y ) ) ) ) + \
    y * ( -86.1329351956084 + y * ( -30.0682112585625 + y * 3.50240264723578 ) ) ) + \
    y * ( 1760.062705994408 + y * ( -675.802947790203 + \
    y * ( 365.7041791005036 + y * ( -108.30162043765552 + 12.78101825083098 * y ) ) ) ) )

    entropy_part_zerop = -( g03 + g08 ) * 0.025 # FIXME: 0.025d0 -> decimal.Decimal('0.025') ?
    return entropy_part_zerop

def  _SP_from_SA_Baltic(SA, lon, lat):
    """
    Calculates Practical Salinity (SP) for the Baltic Sea, from a value computed analytically from Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :math:`-1`]
    lon : array_like
          decimal degrees east [0..+360]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SP_baltic : array_like
                salinity [psu (PSS-78)], unitless

    See Also
    --------
    TODO: explain why these two need Baltic salinity computations
    SP_from_SA, SP_from_Sstar

    Notes
    -----
    This program will only produce Practical Salinty values for the Baltic Sea. Calculates entropy at a sea surface (p = 0 db), except that it does not evaluate any terms that are functions of Absolute Salinity alone. By not calculating these terms, which are a function only of Absolute Salinity, several unnecessary computations are avoided (including saving the computation of a natural logarithm). These terms are a necessary part of entropy, but are not needed when calculating potential temperature from in-situ temperature.

    Examples
    --------
    TODO

    References
    ----------
    .. [1] Feistel, R., S. Weinreben, H. Wolf, S. Seitz, P. Spitzer, B. Adel, G. Nausch, B. Schneider and D. G. Wright, 2010c: Density and Absolute Salinity of the Baltic Sea 2006-2009.  Ocean Science, 6, 3-24.
    http://www.ocean-sci.net/6/3/2010/os-6-3-2010.pdf
    .. [2] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    Available from http://www.TEOS-10.org
    .. [3] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, lon, lat = np.asarray(SA), np.asarray(lon), np.asarray(lat)

    xb1, xb2, xb3 = 12.6, 7., 26.
    xb1a, xb3a = 45., 26.
    yb1, yb2, yb3 = 50., 59., 69.

    inds = (xb2 < lon) & (lon < xb1a) & (yb1 < lat) & (lat < yb3)

    SP_baltic = np.ones( SA.shape )*np.nan

    if np.any(inds):
        xx_left = np.interp( lat[inds], [yb1,yb2,yb3], [xb1,xb2,xb3])
        xx_right = np.interp( lat[inds], [yb1,yb3], [xb1a,xb3a] )
        inds1 = (xx_left <= lon[inds]) & (lon[inds] <= xx_right)

        if np.any(inds1):
            SP_baltic[inds[inds1]] = ( 35 / ( 35.16504 - 0.087 ) ) * ( SA[inds[inds1]] - 0.087)


        SP_baltic = np.reshape( SP_baltic, lon.shape )

    return SP_baltic

def  _SA_from_SP_Baltic(SP, lon, lat):
    """
    Calculates Absolute Salinity in the Baltic Sea, from Practical Salinity.
    Since SP is non-negative by definition, this function changes any negative input values of SP to be zero.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    lon : array_like
          decimal degrees east [0..+360]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SA_baltic : array_like
                Absolute salinity [g kg :math:`-1`]

    See Also
    --------
    TODO: explain why these three need Baltic salinity computations
    SA_from_SP, Sstar_from_SP, SA_Sstar_from_SP

    Notes
    -----
    This programme will only produce Absolute Salinity values for the Baltic Sea.

    Examples
    --------
    TODO

    References
    ----------
    .. [1] Feistel, R., S. Weinreben, H. Wolf, S. Seitz, P. Spitzer, B. Adel, G. Nausch, B. Schneider and D. G. Wright, 2010c: Density and Absolute Salinity of the Baltic Sea 2006-2009.  Ocean Science, 6, 3-24.
    http://www.ocean-sci.net/6/3/2010/os-6-3-2010.pdf
    .. [2] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    Available from http://www.TEOS-10.org
    .. [3] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SP, lon, lat = np.asarray(SP), np.asarray(lon), np.asarray(lat)

    xb1, xb2, xb3 = 12.6, 7., 26.
    xb1a, xb3a = 45., 26.
    yb1, yb2, yb3 = 50., 59., 69.

    inds = (xb2 < lon) & (lon < xb1a) & (yb1 < lat) & (lat < yb3)

    SA_baltic = np.ones( SP.shape )*np.nan

    if np.any(inds):
        xx_left = np.interp( lat[inds], [yb1,yb2,yb3], [xb1,xb2,xb3])
        xx_right = np.interp( lat[inds], [yb1,yb3], [xb1a,xb3a] )
        inds1 = (xx_left <= lon[inds]) & (lon[inds] <= xx_right)

        SA_baltic[inds[inds1]] = ( ( 35.16504 - 0.087 ) / 35 ) * SP[inds[inds1]] + 0.087
        SA_baltic = np.reshape( SA_baltic, lon.shape )

    return SA_baltic

def  _infunnel(SA, CT, p):
    """
    Calculates Absolute Salinity in the Baltic Sea, from Practical Salinity.
    Since SP is non-negative by definition, this function changes any negative input values of SP to be zero.

    Parameters
    ----------
    SA(p) : array_like
         Absolute salinity [g kg :math:`-1`]
    CT(p) : array_like
         Conservative Temperature [:math:`^\\circ` C (TEOS-10)]
    p : array_like
        pressure [db]

    Returns
    -------
    in_funnel : bool
                False, if SA, CT and p are outside the "funnel"
                True, if SA, CT and p are inside the "funnel"

    See Also
    --------
    TODO

    Notes
    -----
    The term "funnel" describes the range of SA, CT and p over which the error in the fit of the computationally-efficient 25-term expression for density in terms of SA, CT and p was calculated (McDougall et al., 2010).

    Examples
    --------
    TODO

    References
    ----------
    .. [1] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, CT, p = np.asarray(SA), np.asarray(CT), np.asarray(p)

    in_funnel = np.ones( SA.shape )
    Inan = ( np.isnan(SA) | np.isnan(CT) | np.isnan(p) )

    Ifunnel = (p > 8000) | (SA < 0) | (SA > 42.2) | \
        ( CT < ( -0.3595467 - 0.0553734 * SA ) ) | \
        ( (p < 5500) & ( SA < 0.006028 * ( p - 500 ) ) ) | \
        ( (p < 5500) & ( CT > ( 33.0 - 0.003818181818182 * p ) ) ) | \
        ( (p > 5500) & ( SA < 30.14 ) ) | \
        ( (p > 5500) & ( CT > 12.0 ) )

    Ifunnel = (Ifunnel == False) # reverse True <-> False
    Ifunnel[Inan] = np.NaN; # TEST

    return Ifunnel