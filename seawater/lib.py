# -*- coding: utf-8 -*-
import numpy as np
try:
    import cPickle as pickle
except:
    import pickle

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
    Calcualtes specifc volume at the Standard Ocean Salinty (SSO) and Conservative Temperature of zero degrees C (CT=0), as a function of pressure (p [db]).

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
    It uses a streamlined version of the 25-term CT version of specific volume ( _rho_alpha_beta_CT25(SA,CT,p) ).

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
     Calcualtes enthalpy at the Standard Ocean Salinty (SSO) and at a Conservative Temperature of zero degrees C (CT=0), as a function of pressure (p [db]).

    Parameters
    ----------
    p : array_like
        pressure [db]

    Returns
    -------
    enthalpy_CT25 : array_like
                    enthalpy_CT25 at (SSO, CT = 0, p), 25-term equation.
                    FIXME: [J]

    See Also
    --------
    TODO

    Notes
    -----
    It Uses a streamlined version of the 25-term CT version of the Gibbs function ( _enthalpy_CT25(SA,CT,p) )

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
                    TODO: write the eq for the second derivative of the specific Gibbs function. FIXME: [units]

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
    Calculates entropy at a sea surface (p = 0 db), except that it does not evaluate any terms that are functions of Absolute Salinity alone.

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
    By not calculating these terms, which are a function only of Absolute Salinity, several unnecessary computations are avoided (including saving the computation of a natural logarithm). These terms are a necessary part of entropy, but are not needed when calculating potential temperature from in-situ temperature.

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
    Inan = np.nonzero( np.isnan(SA) | np.isnan(CT) | np.isnan(p) )

    Ifunnel = (p > 8000) | (SA < 0) | (SA > 42.2) | \
        ( CT < ( -0.3595467 - 0.0553734 * SA ) ) | \
        ( (p < 5500) & ( SA < 0.006028 * ( p - 500 ) ) ) | \
        ( (p < 5500) & ( CT > ( 33.0 - 0.003818181818182 * p ) ) ) | \
        ( (p > 5500) & ( SA < 30.14 ) ) | \
        ( (p > 5500) & ( CT > 12.0 ) )

    Ifunnel = (Ifunnel == False) # reverse True <-> False
    Ifunnel[Inan] = False; # TODO: Nans will become False, change to mask array

    return Ifunnel

def  _entropy_part(SA, t, p):
    """
    Calculates entropy, except that it does not evaluate any terms that are functions of Absolute Salinity alone.

    Parameters
    ----------
    SA(p) : array_like
         Absolute salinity [g kg :math:`-1`]
    t(p) : array_like
         FIXME?: in situ [:math:`^\\circ` C (TEOS-10)]
    p : array_like
        pressure [db]

    Returns
    -------
    entropy_part : array_like
                   entropy minus the terms that due to SA alone TODO: [J/(K.mol)]

    See Also
    --------
    _entropy_part_zerop, TODO

    Notes
    -----
    By not calculating these terms, which are a function only of Absolute Salinity, several unnecessary computations are avoided (including saving the computation of a natural logarithm). These terms are a necessary part of entropy, but are not needed when calculating potential temperature from in-situ temperature.

    Examples
    --------
    TODO

    References
    ----------
    TODO

    Modifications:
    2010-07-23. Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, t, p = np.asarray(SA), np.asarray(t), np.asarray(p)

    SA[SA < 0] = 0

    sfac = 0.0248826675584615 # sfac = 1 / ( 40 * ( 35.16504 / 35 ) )
    x2 = sfac * SA
    x = np.sqrt(x2)
    y = t * 0.025 #FIXME 0.025d0 use Decimal?
    z = p * 1e-4 #TODO: matlab original 1d-4

    g03 = z * ( -270.983805184062 + \
    z * ( 776.153611613101 + z * ( -196.51255088122 + ( 28.9796526294175 - 2.13290083518327 * z ) * z ) ) ) + \
    y * ( -24715.571866078 + z * ( 2910.0729080936 + \
    z * ( -1513.116771538718 + z * ( 546.959324647056 + z * ( -111.1208127634436 + 8.68841343834394 * z ) ) ) ) + \
    y * ( 2210.2236124548363 + z * ( -2017.52334943521 + \
    z * ( 1498.081172457456 + z * ( -718.6359919632359 + ( 146.4037555781616 - 4.9892131862671505 * z ) * z ) ) ) + \
    y * ( -592.743745734632 + z * ( 1591.873781627888 + \
    z * ( -1207.261522487504 + ( 608.785486935364 - 105.4993508931208 * z ) * z ) ) + \
    y * ( 290.12956292128547 + z * ( -973.091553087975 + \
    z * ( 602.603274510125 + z * ( -276.361526170076 + 32.40953340386105 * z ) ) ) + \
    y * ( -113.90630790850321 + y * ( 21.35571525415769 - 67.41756835751434 * z ) + \
    z * ( 381.06836198507096 + z * ( -133.7383902842754 + 49.023632509086724 * z ) ) ) ) ) ) )

    g08 = x2 * ( z * ( 729.116529735046 + \
    z * ( -343.956902961561 + z * ( 124.687671116248 + z * ( -31.656964386073 + 7.04658803315449 * z ) ) ) ) + \
    x * ( x * ( y * ( -137.1145018408982 + y * ( 148.10030845687618 + y * ( -68.5590309679152 + 12.4848504784754 * y ) ) ) - \
    22.6683558512829 * z ) + z * ( -175.292041186547 + ( 83.1923927801819 - 29.483064349429 * z ) * z ) + \
    y * ( -86.1329351956084 + z * ( 766.116132004952 + z * ( -108.3834525034224 + 51.2796974779828 * z ) ) + \
    y * ( -30.0682112585625 - 1380.9597954037708 * z + y * ( 3.50240264723578 + 938.26075044542 * z ) ) ) ) + \
    y * ( 1760.062705994408 + y * ( -675.802947790203 + \
    y * ( 365.7041791005036 + y * ( -108.30162043765552 + 12.78101825083098 * y ) + \
    z * ( -1190.914967948748 + ( 298.904564555024 - 145.9491676006352 * z ) * z ) ) + \
    z * ( 2082.7344423998043 + z * ( -614.668925894709 + ( 340.685093521782 - 33.3848202979239 * z ) * z ) ) ) + \
    z * ( -1721.528607567954 + z * ( 674.819060538734 + \
    z * ( -356.629112415276 + ( 88.4080716616 - 15.84003094423364 * z ) * z ) ) ) ) )

    entropy_part = -( g03 + g08 )  * 0.025 #FIXME: 0.025d0

    return entropy_part

def  _dsa_add_barrier(dsa, lon, lat, longs_ref, lats_ref, dlongs_ref, dlats_ref):
    """
    Adds a barrier through Central America (Panama) and then averages over the appropriate side of the barrier.

    Parameters
    ----------
    dsa : array_like
          Absolute Salinity anomaly of the 4 adjacent neighbours  [g kg :math:`-1`]
    lon : array_like
          decimal degrees east [0..+360]
    lat : array_like
          decimal degrees [-90..+90]
    longs_ref : array_like
          longitudes of regular grid in decimal degrees east [0..+360]
    lats_ref : array_like
          latitudes of regular grid in decimal degrees north [-90..+90]
    dlongs_ref : array_like
          longitudes difference of regular grid in decimal degrees east [0..+360]
    dlats_ref : array_like
          latitudes difference of regular grid in decimal degrees north [-90..+90]

    Returns
    -------
    delta_SA : array_like
          Absolute Salinity anomaly of the 4 adjacent neighbours  [g kg :math:`-1`]

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
    2010-07-23. David Jackett, Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """
    #FIXME: I do not understand the logix for this...

    # Convert input to numpy arrays
    dsa = np.asarray(dsa)
    lon, lat = np.asarray(lon), np.asarray(lat)
    longs_ref, lats_ref = np.asarray(longs_ref), np.asarray(lats_ref)
    dlongs_ref, dlats_ref = np.asarray(dlongs_ref), np.asarray(dlats_ref)

    longs_pan = np.array([260.0000, 272.5900, 276.5000, 278.6500, 280.7300, 295.2170])
    lats_pan = np.array([19.5500, 13.9700, 9.6000, 8.1000, 9.3300, 0])

    lats_lines0 = interp1 ( longs_pan,  lats_pan, lon)
    lats_lines0 = np.interp(lon, longs_pan, lats_pan)

    lats_lines1 = np.interp( longs_ref, lats_pan, longs_pan)
    lats_lines2 = np.interp( (longs_ref+dlongs_ref), lats_pan, longs_pan)

    for k0 in range(0, len(lon.shape) )
        if lats_lines0[k0] <= lat[k0]:
            above_line0 = True
        else:
            above_line0 = False

        if lats_lines1[k0] <= lats_ref[k0]:
            above_line[0] = True
        else:
            above_line[0] = False

        if lats_lines1[k0] <= (lats_ref[k0] + dlats_ref):
            above_line[3] = True
        else:
            above_line[4] = False

        if lats_lines2[k0] <= lats_ref[k0]:
            above_line[1] = True
        else:
            above_line[1] = False

        if lats_lines2[k0] <= (lats_ref[k0] + dlats_ref):
            above_line[2] = True
        else:
            above_line[2] = False

        inds = np.nonzero( above_line ~= above_line0 ) # indices of different sides of CA line
        dsa[inds,k0]) = np.nan

    dsa_mean = dsa.mean()
    inds_nan = np.nonzero( np.isnan( dsa_mean ) )
    no_nan = len(inds_nan)

    for kk in range(0,no_nan):
        col = inds_nan[kk]
        inds_kk = np.nonzero( np.isnan( dsa[:,col] ) )
        Inn = np.nonzero( ~np.isnan( dsa[:,col] ) )
        if Inn.size == 0:
            dsa[inds_kk,col] = dsa[Inn,col].mean()


    delta_SA = dsa
    return delta_SA

def  _dsa_add_mean(dsa):
    """
    Replaces NaN's with nanmean of the 4 adjacent neighbours

    Parameters
    ----------
    dsa : array_like
          Absolute Salinity anomaly of the 4 adjacent neighbours  [g kg :math:`-1`]

    Returns
    -------
    delta_SA : array_like
          Absolute Salinity anomaly of the 4 adjacent neighbours  [g kg :math:`-1`]

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
    2010-07-23. David Jackett, Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """
    #FIXME: I do not understand the logix for this...

    # Convert input to numpy arrays
    dsa = np.asarray(dsa)

    #FIXME: there must be a better way
    dsa_mean = dsa.mean() #FIXME: should be nanmean here in the original...
    inds_nan = np.nonzero( np.isnan(dsa_mean) )
    no_nan = len(inds_nan)

    for kk in range(0, no_nan):
        col = inds_nan[kk]
        inds_kk = np.nonzero( np.isnan( dsa[:,col] ) )
        Inn = np.nonzero(~np.isnan( dsa[:,col] ) )
        if Inn.size == 0:
            dsa[inds_kk, col] = dsa[Inn,col].mean()

    delta_SA = dsa

    return delta_SA

def  _delta_SA(p, lon, lat):
    """
    Calculates the Absolute Salinity anomaly, SA - SR, in the open ocean by spatially interpolating the global reference data set of delta_SA to the location of the seawater sample.

    Parameters
    ----------
    p : array_like
        pressure [db]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    delta_SA : array_like
               Absolute Salinity anomaly [g kg :math:`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    TODO

    Notes
    -----
    The Absolute Salinity Anomaly in the Baltic Sea is evaluated separately, since it is a function of Practical Salinity, not of space.   The present function returns a delta_SA of zero for data in the Baltic Sea.  The correct way of calculating Absolute Salinity in the Baltic Sea is by calling _SA_from_SP.

    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometres inland from the coast.

    Examples
    --------
    TODO

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    Available from http://www.TEOS-10.org
    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean.  Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    p = np.asarray(p), np.asarray(lon), np.asarray(lat)

    data = pickle.load( open('gsw_data_v2_0.pkl', 'rb') )
    delta_SA_ref = data['delta_SA_ref']
    lats_ref = data['lats_ref']
    longs_ref = data['longs_ref']
    p_ref = data['p_ref']
    ndepth_ref = data['ndepth_ref']

    nx = len(longs_ref)
    ny = len(lats_ref)
    nz = len(p_ref)
    n0 = len(p)

    dlongs_ref = longs_ref[1] - longs_ref[0]
    dlats_ref = lats_ref[1] - lats_ref[0]

    indsx0 = np.floor( 1 + (nx-1) * (lon - longs_ref[0] ) / (longs_ref[nx]- longs_ref[0] ) )
    #indsx0 = indsx0(:) #FIXME: need to flat it?
    inds = np.nonzero(indsx0 == nx)
    indsx0[inds] = nx - 1

    indsy0 = np.floor( 1 + (ny-1) * ( lat - lats_ref[0] ) / ( lats_ref[ny]- lats_ref[0] ) )
    #indsy0 = indsy0(:) #FIXME: need to flat it?
    inds = np.nonzero(indsy0 == ny)
    indsy0[inds] = ny - 1

    indsz0 = np.sum( np.ones(nz) * p >= p_ref * np.ones(n0) #FIXME: will break matlab matrix mult
    #indsz0 = indsz0(:) #adjust in the vertical FIXME

    indsn1 = sub2ind([ny,nx],indsy0,indsx0);        % casts containing data
    indsn2 = sub2ind([ny,nx],indsy0,indsx0+1);
    indsn3 = sub2ind([ny,nx],indsy0+1,indsx0+1);
    indsn4 = sub2ind([ny,nx],indsy0+1,indsx0);

    nmax = max([ndepth_ref(indsn1)';ndepth_ref(indsn2)';ndepth_ref(indsn3)';ndepth_ref(indsn4)']);

    inds1 = find(indsz0(:)' > nmax);                % casts deeper than GK maximum

    p(inds1) = p_ref(nmax(inds1));                  % have reset p here so have to reset indsz0

    indsz0 = sum(ones(nz,1)*p(:)' >= p_ref(:)*ones(1,n0));
    indsz0 = indsz0(:);
    inds = find(indsz0 == nz);
    indsz0(inds) = nz - 1;

    inds0 = sub2ind([nz,ny,nx],indsz0,indsy0,indsx0);

    data_indices = [indsx0,indsy0,indsz0,inds0];
    data_inds = data_indices(:,3);

    r1 = (long(:) - longs_ref(indsx0))./(longs_ref(indsx0+1) - longs_ref(indsx0));
    s1 = (lat(:) - lats_ref(indsy0))./(lats_ref(indsy0+1) - lats_ref(indsy0));
    t1 = (p(:) - p_ref(indsz0))./(p_ref(indsz0+1) - p_ref(indsz0));

    nksum = 0;
    no_levels_missing = 0;

    sa_upper = nan(size(data_inds));
    sa_lower = nan(size(data_inds));
    delta_SA = nan(size(data_inds));
    in_ocean = ones(size(delta_SA));

    for k = 1:nz-1

        inds_k = find(indsz0 == k);
        nk = length(inds_k);

        if nk>0
            nksum = nksum+nk;
            indsx = indsx0(inds_k);
            indsy = indsy0(inds_k);
            indsz = k*ones(size(indsx));
            inds_di = find(data_inds == k);             % level k interpolation
            dsa = nan(4,n0);
            inds1 = sub2ind([nz,ny,nx], indsz, indsy, indsx);
            dsa(1,inds_k) = delta_SA_ref(inds1);
            inds2 = sub2ind([nz,ny,nx], indsz, indsy, indsx+1);
            dsa(2,inds_k) = delta_SA_ref(inds2);                % inds0 + ny*nz
            inds3 = sub2ind([nz,ny,nx], indsz, indsy+1, indsx+1);
            dsa(3,inds_k) = delta_SA_ref(inds3);           % inds0 + ny*nz + nz
            inds4 = sub2ind([nz ny,nx], indsz, indsy+1, indsx);
            dsa(4,inds_k) = delta_SA_ref(inds4);                   % inds0 + nz

            inds = find(260<=long(:) & long(:)<=295.217 & ...
                0<=lat(:) & lat(:)<=19.55 & indsz0(:)==k);
            if ~isempty(inds)
                dsa(:,inds) = gsw_dsa_add_barrier(dsa(:,inds),long(inds), ...
                    lat(inds),longs_ref(indsx0(inds)),lats_ref(indsy0(inds)),dlongs_ref,dlats_ref);
            end

            inds = find(isnan(sum(dsa))' & indsz0==k);
            if ~isempty(inds)
                dsa(:,inds) = gsw_dsa_add_mean(dsa(:,inds));
            end

            sa_upper(inds_di) = (1-s1(inds_di)).*(dsa(1,inds_k)' + ...
                r1(inds_di).*(dsa(2,inds_k)'-dsa(1,inds_k)')) + ...
                s1(inds_di).*(dsa(4,inds_k)' + ...
                r1(inds_di).*(dsa(3,inds_k)'-dsa(4,inds_k)'));  % level k+1 interpolation

            dsa = nan(4,n0);
            inds1 = sub2ind([nz,ny,nx], indsz+1, indsy, indsx);
            dsa(1,inds_k) = delta_SA_ref(inds1);
            inds2 = sub2ind([nz,ny,nx], indsz+1, indsy, indsx+1);
            dsa(2,inds_k) = delta_SA_ref(inds2);                % inds1 + ny*nz
            inds3 = sub2ind([nz,ny,nx], indsz+1, indsy+1, indsx+1);
            dsa(3,inds_k) = delta_SA_ref(inds3);           % inds1 + ny*nz + nz
            inds4 = sub2ind([nz ny,nx], indsz+1, indsy+1, indsx);
            dsa(4,inds_k) = delta_SA_ref(inds4);                   % inds1 + nz

            inds = find(260<=long(:) & long(:)<=295.217 & ...
                0<=lat(:) & lat(:)<=19.55 & indsz0(:)==k);
            if ~isempty(inds)
                dsa(:,inds) = gsw_dsa_add_barrier(dsa(:,inds),long(inds), ...
                    lat(inds),longs_ref(indsx0(inds)),lats_ref(indsy0(inds)),dlongs_ref,dlats_ref);
            end

            inds = find(isnan(sum(dsa))' & indsz0==k);

            if ~isempty(inds)
                dsa(:,inds) = gsw_dsa_add_mean(dsa(:,inds));
            end

            sa_lower(inds_di) = (1-s1(inds_di)).*(dsa(1,inds_k)' + ...
                r1(inds_di).*(dsa(2,inds_k)'-dsa(1,inds_k)')) + ...
                s1(inds_di).*(dsa(4,inds_k)' + ...
                r1(inds_di).*(dsa(3,inds_k)'-dsa(4,inds_k)'));

            inds_different = find(isfinite(sa_upper(inds_di)) & isnan(sa_lower(inds_di)));

            if ~isempty(inds_different)
                sa_lower(inds_di(inds_different)) = sa_upper(inds_di(inds_different));
            end

            delta_SA(inds_di) = sa_upper(inds_di) + t1(inds_di).*(sa_lower(inds_di) - sa_upper(inds_di));

        else
            no_levels_missing = no_levels_missing + 1;
        end
    end

    inds = find(~isfinite(delta_SA));
    delta_SA(inds) = 0;

    in_ocean(inds) = 0;

    returns delta_SA, in_ocean

#def  _interp_McD_Klocker(spycnl, A='gn'):
    #"""
    #Interpolates the reference cast with respect to the interpolating variable "spycnl". This reference cast is located at (188 E, 4 N) from the reference data set which underlies the Jackett & McDougall (1997) Neutral Density computer code.

    #Parameters
    #----------
    #spycnl : array_like
             #TODO
    #A : array_like, default 'gn'
        #TODO: isopycnal ?
        #'s2'

    #Returns
    #-------
    #TODO: SA_iref_cast, CT_iref_cast, p_iref_cast

    #See Also
    #--------
    #TODO

    #Notes
    #-----
    #This function finds the values of SA, CT and p on this reference cast which correspond to the value of isopycnal which is passed to this function from the function "_geo_str_McD_Klocker". The isopycnal could be either gamma_n or sigma_2. if A is set to any of the following 's2','S2','sigma2','sigma_2' the interpolation will take place in sigma 2 space, any other input will result in the programme working in gamma_n space.

    #Examples
    #--------
    #TODO

    #References
    #----------
    #Jackett, D. R. and T. J. McDougall, 1997: A neutral density variable for the world's oceans. Journal of Physical Oceanography, 27, 237-263.

    #Modifications:
    #2010-07-23. Trevor McDougall and Paul Barker
    #2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    #"""

    ## Convert input to numpy arrays
    #spycnl= np.asarray(spycnl)

    #data = pickle.load( open('gsw_data_v2_0.pkl', 'rb') )
    #SA_ref_cast = data['SA_ref_cast']
    #CT_ref_cast = data['CT_ref_cast']
    #p_ref_cast = data['p_ref_cast']

    #if A == 's2':
        #spycnl_ref_cast = data['sigma_2_ref_cast']
    #elif A == 'gn':
        #spycnl_ref_cast = data['gamma_n_ref_cast']
    #else:
        #print "unknown method" #FIXME: add a proper python error

    #min_spycnl_ref_cast, Imin_spycnl_ref_cast = spycnl_ref_cast.min(), spycnl_ref_cast.argmin()

    #Ishallow = np.nonzero( spycnl <= min_spycnl_ref_cast ) # Set equal to the shallowest bottle.

    #SA_iref_cast[Ishallow] = SA_ref_cast[Imin_spycnl_ref_cast]
    #CT_iref_cast[Ishallow] = CT_ref_cast[Imin_spycnl_ref_cast]
    #p_iref_cast[Ishallow] = p_ref_cast[Imin_spycnl_ref_cast]

    #max_spycnl_ref_cast, Imax_spycnl_ref_cast = spycnl_ref_cast.max(), spycnl_ref_cast.argmax()

    #Ideep = np.nonzero( spycnl >= max_spycnl_ref_cast ) # Set equal to the deepest bottle.

    #SA_iref_cast[Ideep] = SA_ref_cast[Imax_spycnl_ref_cast]
    #CT_iref_cast[Ideep] = CT_ref_cast[Imax_spycnl_ref_cast]
    #p_iref_cast[Ideep] = p_ref_cast[Imax_spycnl_ref_cast]

    #I = np.nonzero(spycnl >= 21.805 & spycnl <= 28.3614)

    #xi = spycnl[I]
    #x = spycnl_ref_cast

    #siz = xi.shape # FIXME: unfinished
    #if xi.shape > 1:
        #[xxi, k] = sort(xi)
        #[dum, j] = sort([x;xxi])
        #r(j) = 1:length(j)
        #r = r(length(x)+1:end) - (1:length(xxi))
        #r(k) = r
        #r(xi==x(end)) = length(x)-1
        #ind = find((r>0) & (r<length(x)))
        #ind = ind(:)
        #SA_ref_casti = NaN(length(xxi),size(SA_ref_cast,2),superiorfloat(x,SA_ref_cast,xi))
        #CT_ref_casti = NaN(length(xxi),size(CT_ref_cast,2),superiorfloat(x,CT_ref_cast,xi))
        #p_ref_casti = NaN(length(xxi),size(p_ref_cast,2),superiorfloat(x,p_ref_cast,xi))
        #rind = r(ind)
        #xrind = x(rind)
        #u = (xi(ind)-xrind)./(x(rind+1)-xrind)
        #SArind = SA_ref_cast(rind,:)
        #CTrind = CT_ref_cast(rind,:)
        #prind = p_ref_cast(rind,:)
        #SA_ref_casti(ind,:) = SArind + bsxfun(@times,SA_ref_cast(rind+1,:)-SArind,u)
        #CT_ref_casti(ind,:) = CTrind + bsxfun(@times,CT_ref_cast(rind+1,:)-CTrind,u)
        #p_ref_casti(ind,:) = prind + bsxfun(@times,p_ref_cast(rind+1,:)-prind,u)
    #else: # Special scalar xi case
        #r = find(x <= xi,1,'last')
        #r(xi==x(end)) = length(x)-1
        #if isempty(r) || r<=0 || r>=length(x):
            #SA_ref_casti = NaN(1,size(SA,2),superiorfloat(x,SA,xi))
            #CT_ref_casti = NaN(1,size(CT,2),superiorfloat(x,CT,xi))
            #p_ref_casti = NaN(1,size(p,2),superiorfloat(x,p,xi))
        #else:
            #u = (xi-x(r))./(x(r+1)-x(r))
            #SAr_ref_cast = SA_ref_cast(r,:)
            #CTr_ref_cast = CT_ref_cast(r,:)
            #pr_ref_cast = p_ref_cast(r,:)
            #SA_ref_casti = SAr_ref_cast + bsxfun(@times,SA_ref_cast(r+1,:)-SAr_ref_cast,u)
            #CT_ref_casti = CTr_ref_cast + bsxfun(@times,CT_ref_cast(r+1,:)-CTr_ref_cast,u)
            #p_ref_casti = pr_ref_cast + bsxfun(@times,p_ref_cast(r+1,:)-pr_ref_cast,u)


    #if min(size(SA_ref_casti)) == 1 && numel(xi) > 1:
        #SA_ref_casti = reshape(SA_ref_casti,siz)
        #CT_ref_casti = reshape(CT_ref_casti,siz)
        #p_ref_casti = reshape(p_ref_casti,siz)


    #SA_iref_cast[I] = SA_ref_casti
    #CT_iref_cast[I] = CT_ref_casti
    #p_iref_cast[I] = p_ref_casti

    #return SA_iref_cast, CT_iref_cast, p_iref_cast

def sub2ind(shape, I, J, row_major=True):
    if row_major:
        ind = (I % shape[0]) * shape[1] + (J % shape[1])
    else:
        ind = (J % shape[1]) * shape[0] + (I % shape[0])
    return ind