# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
from seawater import constants as cte

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

    enthalpy_SSO_0_CT25 = cte.db2Pascal * ( ( a1 / (2*b2) ) * np.log( 1 + p * ( 2 * b1 + b2 * p ) / b0 ) + part * np.log( 1 + ( b2 * p * (B - A) ) / (A * (B + b2 * p ) ) ) )

    return enthalpy_SSO_0_CT25

def  _gibbs_pt0_pt0(SA, pt0):
    """
    Calculates the second derivative of the specific Gibbs function with respect to temperature at zero sea pressure.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup::`-1`]
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
    if SA.shape:
        SA[SA < 0] = 0
    elif SA < 0:
        SA = 0

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
         Absolute salinity [g kg :sup::`-1`]
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
    if SA.shape:
        SA[SA < 0] = 0
    elif SA < 0:
        SA = 0

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
         Absolute salinity [g kg :sup::`-1`]
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
                Absolute salinity [g kg :sup::`-1`]

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

    inds_baltic = (xb2 < lon) & (lon < xb1a) & (yb1 < lat) & (lat < yb3)
    SA_baltic = np.ones( SP.shape )*np.nan

    #if inds_itrp.size != 0:
    if list(inds_baltic): #FIXME: find a equivalent for numpy arrays
        xx_left = np.interp( lat[inds_baltic], [yb1,yb2,yb3], [xb1,xb2,xb3])
        xx_right = np.interp( lat[inds_baltic], [yb1,yb3], [xb1a,xb3a] )
        inds_baltic1 = (xx_left <= lon[inds_baltic]) & (lon[inds_baltic] <= xx_right)
        #SA_baltic.flatten('F')[inds[inds1]] = ( ( 35.16504 - 0.087 ) / 35 ) * SP.flatten('F')[inds[inds1]] + 0.087
        SA_baltic[inds_baltic[inds_baltic1]] = ( ( 35.16504 - 0.087 ) / 35 ) * SP[inds_baltic[inds_baltic1]] + 0.087
        #SA_baltic = np.reshape( SA_baltic, lon.shape )

    return SA_baltic

def  _infunnel(SA, CT, p):
    """
    Calculates Absolute Salinity in the Baltic Sea, from Practical Salinity.
    Since SP is non-negative by definition, this function changes any negative input values of SP to be zero.

    Parameters
    ----------
    SA(p) : array_like
         Absolute salinity [g kg :sup::`-1`]
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
    Ifunnel[Inan] = False; # TODO: Nans will become False, change to mask array

    return Ifunnel

def  _entropy_part(SA, t, p):
    """
    Calculates entropy, except that it does not evaluate any terms that are functions of Absolute Salinity alone.

    Parameters
    ----------
    SA(p) : array_like
         Absolute salinity [g kg :sup::`-1`]
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

    if SA.shape:
        SA[SA < 0] = 0
    elif SA < 0:
        SA = 0

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

def _gibbs(ns, nt, npr, SA, t, p):
    """
    Calculates specific Gibbs energy and its derivatives up to order 2 for seawater.

    Parameters
    ----------
    ns : int
         order of SA derivative [0, 1 or 2 ]
    nt : int
          order of t derivative [0, 1 or 2 ]
    npr : int
         order of p derivative [0, 1 or 2 ]
    SA : array_like
         Absolute salinity [g kg :sup::`-1`]
    t : array_like
         in situ temperature [:math:`^\\circ` C (ITS-90)]
    p : array_like
        pressure [db]

    Returns
    -------
    gibbs : array_like
           Specific Gibbs energy or its derivatives.
           The Gibbs energy (when ns = nt = npr = 0) has units of: [ J kg :sup::`-1` ]
           The Absolute Salinity derivatives are output in units of: [ (J kg :sup::`-1`) (g kg :sup::`-1`) :sup::`-ns` ]
           The temperature derivatives are output in units of: [ (J kg :sup::`-1`) K :sup::`-nt` ]
           The pressure derivatives are output in units of: [ (J kg :sup::`-1`) Pa :sup::`-npr` ]
           The mixed derivatives are output in units of: [ (J kg :sup::`-1`) (g kg :sup::`-1`) :sup::`-ns` K :sup::`-nt` Pa :sup::`-npr` ]
           Note. The derivatives are taken with respect to pressure in Pa, not withstanding that the pressure input into this routine is in dbar.


    See Also
    --------
    TODO

    Notes
    -----
    The Gibbs function for seawater is that of TEOS-10 (IOC et al., 2010), being the sum of IAPWS-08 for the saline part and IAPWS-09 for the pure water part. These IAPWS releases are the officially blessed IAPWS descriptions of Feistel (2008) and the pure water part of Feistel (2003). Absolute Salinity, SA, in all of the GSW routines is expressed on the Reference-Composition Salinity Scale of 2008 (RCSS-08) of Millero et al. (2008).

    Examples
    --------
    TODO

    References
    ----------
    .. [1] Feistel, R., 2003: A new extended Gibbs thermodynamic potential of seawater, Progr. Oceanogr., 58, 43-114.

    .. [2] Feistel, R., 2008: A Gibbs function for seawater thermodynamics for -6 to 80 :math:`^\\circ` C and salinity up to 120 g kg :sup::`-1`, Deep-Sea Res. I, 55, 1639-1671.

    .. [3] IAPWS, 2008: Release on the IAPWS Formulation 2008 for the Thermodynamic Properties of Seawater. The International Association for the Properties of Water and Steam. Berlin, Germany, September 2008, available from http://www.iapws.org.  This Release is referred to as IAPWS-08.

    .. [4] IAPWS, 2009: Supplementary Release on a Computationally Efficient Thermodynamic Formulation for Liquid Water for Oceanographic Use. The International Association for the Properties of Water and Steam. Doorwerth, The Netherlands, September 2009, available from http://www.iapws.org.  This Release is referred to as IAPWS-09.

    .. [5] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org See section 2.6 and appendices A.6,  G and H of this TEOS-10 Manual.

    .. [6] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: The composition of Standard Seawater and the definition of the Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.

    Modifications:
    2010-09-24. David Jackett, Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    SA, t, p = np.asarray(SA), np.asarray(t), np.asarray(p)

    # ensure that SA is non-negative.
    SA[SA < 0] = 0

    sfac = 0.0248826675584615 # sfac = 1 / ( 40 * ( 35.16504 / 35 ) )

    x2 = sfac * SA
    x = np.sqrt(x2)
    y = t * 0.025 # FIXME: 0.025d0
    z = p * 1e-4 # The input pressure (p) is sea pressure in units of dbar.
    #FIXME: check if 1d-4 is OK

    if (ns==0) & (nt==0) & (npr==0):
        g03 = 101.342743139674 + z * ( 100015.695367145 + \
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

        g08 = x2 * ( 1416.27648484197 + z * ( -3310.49154044839 + \
        z * ( 384.794152978599 + z * ( -96.5324320107458 + ( 15.8408172766824 - 2.62480156590992 * z ) * z ) ) ) + \
        x * ( -2432.14662381794 + x * ( 2025.80115603697 + \
        y * ( 543.835333000098 + y * ( -68.5572509204491 + \
        y * ( 49.3667694856254 + y * ( -17.1397577419788 + 2.49697009569508 * y ) ) ) - 22.6683558512829 * z ) + \
        x * ( -1091.66841042967 - 196.028306689776 * y + \
        x * ( 374.60123787784 - 48.5891069025409 * x + 36.7571622995805 * y ) + 36.0284195611086 * z ) + \
        z * ( -54.7919133532887 + ( -4.08193978912261 - 30.1755111971161 * z ) * z ) ) + \
        z * ( 199.459603073901 + z * ( -52.2940909281335 + ( 68.0444942726459 - 3.41251932441282 * z ) * z ) ) + \
        y * ( -493.407510141682 + z * ( -175.292041186547 + ( 83.1923927801819 - 29.483064349429 * z ) * z ) + \
        y * ( -43.0664675978042 + z * ( 383.058066002476 + z * ( -54.1917262517112 + 25.6398487389914 * z ) ) + \
        y * ( -10.0227370861875 - 460.319931801257 * z + y * ( 0.875600661808945 + 234.565187611355 * z ) ) ) ) )  + \
        y * ( 168.072408311545 + z * ( 729.116529735046 + \
        z * ( -343.956902961561 + z * ( 124.687671116248 + z * ( -31.656964386073 + 7.04658803315449 * z ) ) ) ) + \
        y * ( 880.031352997204 + y * ( -225.267649263401 + \
        y * ( 91.4260447751259 + y * ( -21.6603240875311 + 2.13016970847183 * y ) + \
        z * ( -297.728741987187 + ( 74.726141138756 - 36.4872919001588 * z ) * z ) ) + \
        z * ( 694.244814133268 + z * ( -204.889641964903 + ( 113.561697840594 - 11.1282734326413 * z ) * z ) ) ) + \
        z * ( -860.764303783977 + z * ( 337.409530269367 + \
        z * ( -178.314556207638 + ( 44.2040358308 - 7.92001547211682 * z ) * z ) ) ) ) ) )

        g08[x>0] = g08[x>0] + x2[x>0] * ( 5812.81456626732 + 851.226734946706 * y[x>0] ) * np.log( x[x>0] )

        gibbs = g03 + g08

    elif (ns==1) & (nt==0) & (npr==0):
        g08 = 8645.36753595126 + z * ( -6620.98308089678 + \
        z * ( 769.588305957198 + z * ( -193.0648640214916 + ( 31.6816345533648 - 5.24960313181984 * z ) * z ) ) ) + \
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
        y * ( 1187.3715515697959 + z * ( 1458.233059470092 + \
        z * ( -687.913805923122 + z * ( 249.375342232496 + z * ( -63.313928772146 + 14.09317606630898 * z ) ) ) ) + \
        y * ( 1760.062705994408 + y * ( -450.535298526802 + \
        y * ( 182.8520895502518 + y * ( -43.3206481750622 + 4.26033941694366 * y ) + \
        z * ( -595.457483974374 + ( 149.452282277512 - 72.9745838003176 * z ) * z ) ) + \
        z * ( 1388.489628266536 + z * ( -409.779283929806 + ( 227.123395681188 - 22.2565468652826 * z ) * z ) ) ) + \
        z * ( -1721.528607567954 + z * ( 674.819060538734 + \
        z * ( -356.629112415276 + ( 88.4080716616 - 15.84003094423364 * z ) * z ) ) ) ) )

        g08[x>0] = g08[x>0] + ( 11625.62913253464 + 1702.453469893412 * y[x>0] ) * np.log( x[x>0] )
        g08[x==0] = np.nan

        gibbs = 0.5 * sfac * g08

    elif (ns==0) & (nt==1) & (npr==0):
        g03 = 5.90578347909402 + z * ( -270.983805184062 + \
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

        g08 = x2 * ( 168.072408311545 + z * ( 729.116529735046 + \
        z * ( -343.956902961561 + z * ( 124.687671116248 + z * ( -31.656964386073 + 7.04658803315449 * z ) ) ) ) + \
        x * ( -493.407510141682 + x * ( 543.835333000098 + x * ( -196.028306689776 + 36.7571622995805 * x ) + \
        y * ( -137.1145018408982 + y * ( 148.10030845687618 + y * ( -68.5590309679152 + 12.4848504784754 * y ) ) ) - \
        22.6683558512829 * z ) + z * ( -175.292041186547 + ( 83.1923927801819 - 29.483064349429 * z ) * z ) + \
        y * ( -86.1329351956084 + z * ( 766.116132004952 + z * ( -108.3834525034224 + 51.2796974779828 * z ) ) + \
        y * ( -30.0682112585625 - 1380.9597954037708 * z + y * ( 3.50240264723578 + 938.26075044542 * z ) ) ) ) + \
        y * ( 1760.062705994408 + y * ( -675.802947790203 + \
        y * ( 365.7041791005036 + y * ( -108.30162043765552 + 12.78101825083098 * y ) + \
        z * ( -1190.914967948748 + ( 298.904564555024 - 145.9491676006352 * z ) * z ) ) + \
        z * ( 2082.7344423998043 + z * ( -614.668925894709 + ( 340.685093521782 - 33.3848202979239 * z ) * z ) ) ) + \
        z * ( -1721.528607567954 + z * ( 674.819060538734 + \
        z * ( -356.629112415276 + ( 88.4080716616 - 15.84003094423364 * z ) * z ) ) ) ) )

        g08[x>0] = g08[x>0] + 851.226734946706 * x2[x>0] * np.log( x[x>0] )

        gibbs = (g03 + g08) * 0.025 # FIXME: 0.025d0

    elif (ns==0) & (nt==0) & (npr==1):
        g03 = 100015.695367145 + z * ( -5089.1530840726 + \
        z * ( 853.5533353388611 + z * ( -133.2587017014444 + ( 21.0131554401542 - 3.278571068826234 * z ) * z ) ) ) + \
        y * ( -270.983805184062 + z * ( 1552.307223226202 + \
        z * ( -589.53765264366 + ( 115.91861051767 - 10.664504175916349 * z ) * z ) ) + \
        y * ( 1455.0364540468 + z * ( -1513.116771538718 + \
        z * ( 820.438986970584 + z * ( -222.2416255268872 + 21.72103359585985 * z ) ) ) + \
        y * ( -672.50778314507 + z * ( 998.720781638304 + \
        z * ( -718.6359919632359 + ( 195.2050074375488 - 8.31535531044525 * z ) * z ) ) + \
        y * ( 397.968445406972 + z * ( -603.630761243752 + ( 456.589115201523 - 105.4993508931208 * z ) * z ) + \
        y * ( -194.618310617595 + y * ( 63.5113936641785 - 9.63108119393062 * y + \
        z * ( -44.5794634280918 + 24.511816254543362 * z ) ) + \
        z * ( 241.04130980405 + z * ( -165.8169157020456 + 25.92762672308884 * z ) ) ) ) ) ) )

        g08 = x2 * ( -3310.49154044839 + z * ( 769.588305957198 + \
        z * ( -289.5972960322374 + ( 63.3632691067296 - 13.1240078295496 * z ) * z ) ) + \
        x * ( 199.459603073901 + x * ( -54.7919133532887 + 36.0284195611086 * x - 22.6683558512829 * y + \
        ( -8.16387957824522 - 90.52653359134831 * z ) * z ) + \
        z * ( -104.588181856267 + ( 204.1334828179377 - 13.65007729765128 * z ) * z ) + \
        y * ( -175.292041186547 + ( 166.3847855603638 - 88.449193048287 * z ) * z + \
        y * ( 383.058066002476 + y * ( -460.319931801257 + 234.565187611355 * y ) + \
        z * ( -108.3834525034224 + 76.9195462169742 * z ) ) ) ) + \
        y * ( 729.116529735046 + z * ( -687.913805923122 + \
        z * ( 374.063013348744 + z * ( -126.627857544292 + 35.23294016577245 * z ) ) )  + \
        y * ( -860.764303783977 + y * ( 694.244814133268 + \
        y * ( -297.728741987187 + ( 149.452282277512 - 109.46187570047641 * z ) * z ) + \
        z * ( -409.779283929806 + ( 340.685093521782 - 44.5130937305652 * z ) * z ) ) + \
        z * ( 674.819060538734 + z * ( -534.943668622914 + ( 176.8161433232 - 39.600077360584095 * z ) * z ) ) ) ) )

        gibbs = (g03 + g08) * 1e-8 #This pressure derivative of the gibbs function is in units of (J kg :sup::`-1`) (Pa :sup::`-1`) = m :sup::`3` kg :sup::`-1`

    elif (ns==1) & (nt==1) & (npr==0):
        g08 = 1187.3715515697959 + z * ( 1458.233059470092 + \
        z * ( -687.913805923122 + z * ( 249.375342232496 + z * ( -63.313928772146 + 14.09317606630898 * z ) ) ) ) + \
        x * ( -1480.222530425046 + x * ( 2175.341332000392 + x * ( -980.14153344888 + 220.542973797483 * x ) + \
        y * ( -548.4580073635929 + y * ( 592.4012338275047 + y * ( -274.2361238716608 + 49.9394019139016 * y ) ) ) - \
        90.6734234051316 * z ) + z * ( -525.876123559641 + ( 249.57717834054571 - 88.449193048287 * z ) * z ) + \
        y * ( -258.3988055868252 + z * ( 2298.348396014856 + z * ( -325.1503575102672 + 153.8390924339484 * z ) ) + \
        y * ( -90.2046337756875 - 4142.8793862113125 * z + y * ( 10.50720794170734 + 2814.78225133626 * z ) ) ) ) + \
        y * ( 3520.125411988816 + y * ( -1351.605895580406 + \
        y * ( 731.4083582010072 + y * ( -216.60324087531103 + 25.56203650166196 * y ) + \
        z * ( -2381.829935897496 + ( 597.809129110048 - 291.8983352012704 * z ) * z ) ) + \
        z * ( 4165.4688847996085 + z * ( -1229.337851789418 + ( 681.370187043564 - 66.7696405958478 * z ) * z ) ) ) + \
        z * ( -3443.057215135908 + z * ( 1349.638121077468 + \
        z * ( -713.258224830552 + ( 176.8161433232 - 31.68006188846728 * z ) * z ) ) ) )

        g08[x>0] = g08[x>0] + 1702.453469893412 * np.log( x[x>0] )
        g08[SA==0] = np.nan
        gibbs = 0.5 * sfac * 0.025 * g08 # FIXME: 0.025d0

    elif (ns==1) & (nt==0) & (npr==1):
        g08 = -6620.98308089678 + z * ( 1539.176611914396 + \
        z * ( -579.1945920644748 + ( 126.7265382134592 - 26.2480156590992 * z ) * z ) ) + \
        x * ( 598.378809221703 + x * ( -219.1676534131548 + 180.142097805543 * x - 90.6734234051316 * y + \
        (-32.65551831298088 - 362.10613436539325 * z ) * z ) + \
        z * ( -313.764545568801 + ( 612.4004484538132 - 40.95023189295384 * z ) * z ) + \
        y * ( -525.876123559641 + ( 499.15435668109143 - 265.347579144861 * z ) * z + \
        y * ( 1149.174198007428 + y * ( -1380.9597954037708 + 703.695562834065 * y ) + \
        z * ( -325.1503575102672 + 230.7586386509226 * z ) ) ) ) + \
        y * ( 1458.233059470092 + z * ( -1375.827611846244 + \
        z * ( 748.126026697488 + z * ( -253.255715088584 + 70.4658803315449 * z ) ) )  + \
        y * ( -1721.528607567954 + y * ( 1388.489628266536 + \
        y * ( -595.457483974374 + ( 298.904564555024 - 218.92375140095282 * z ) * z ) + \
        z * ( -819.558567859612 + ( 681.370187043564 - 89.0261874611304 * z ) * z ) ) + \
        z * ( 1349.638121077468 + z * ( -1069.887337245828 + ( 353.6322866464 - 79.20015472116819 * z ) * z ) ) ) )

        gibbs = g08 * sfac * 0.5e-8 # FIXME: 0.5d
        # This derivative of the Gibbs function is in units of (m :sup::`3` kg :sup::`-1`) / (g kg :sup::`-1`) = m :sup::`3` g :sup::`-1`
        # that is, it is the derivative of specific volume with respect to Absolute Salinity measured in g kg :sup::`-1`.

    elif (ns==0) & (nt==1) & (npr==1):
        g03 = -270.983805184062 + z * ( 1552.307223226202 + z * ( -589.53765264366 + \
        ( 115.91861051767 - 10.664504175916349 * z ) * z ) ) + \
        y * ( 2910.0729080936 + z * ( -3026.233543077436 + \
        z * ( 1640.877973941168 + z * ( -444.4832510537744 + 43.4420671917197 * z ) ) ) + \
        y * ( -2017.52334943521 + z * ( 2996.162344914912 + \
        z * ( -2155.907975889708 + ( 585.6150223126464 - 24.946065931335752 * z ) * z ) ) + \
        y * ( 1591.873781627888 + z * ( -2414.523044975008 + ( 1826.356460806092 - 421.9974035724832 * z ) * z ) + \
        y * ( -973.091553087975 + z * ( 1205.20654902025 + z * ( -829.084578510228 + 129.6381336154442 * z ) ) + \
        y * ( 381.06836198507096 - 67.41756835751434 * y + z * ( -267.4767805685508 + 147.07089752726017 * z ) ) ) ) ) )

        g08 = x2 * ( 729.116529735046 + z * ( -687.913805923122 + \
        z * ( 374.063013348744 + z * ( -126.627857544292 + 35.23294016577245 * z ) ) ) + \
        x * ( -175.292041186547 - 22.6683558512829 * x + ( 166.3847855603638 - 88.449193048287 * z ) * z + \
        y * ( 766.116132004952 + y * ( -1380.9597954037708 + 938.26075044542 * y ) + \
        z * ( -216.7669050068448 + 153.8390924339484 * z ) ) ) + \
        y * ( -1721.528607567954 + y * ( 2082.7344423998043 + \
        y * ( -1190.914967948748 + ( 597.809129110048 - 437.84750280190565 * z ) * z ) + \
        z * ( -1229.337851789418 + ( 1022.055280565346 - 133.5392811916956 * z ) * z ) ) + \
        z * ( 1349.638121077468 + z * ( -1069.887337245828 + ( 353.6322866464 - 79.20015472116819 * z ) * z ) ) ) )

        gibbs = (g03 + g08) * 2.5e-10
        # This derivative of the Gibbs function is in units of (m :sup::`3` (K kg) ) that is, the pressure of the derivative in Pa.

    elif (ns==2) & (nt==0) & (npr==0):
        g08 = 2.0 * ( 8103.20462414788 + \
        y * ( 2175.341332000392 + y * ( -274.2290036817964 + \
        y * ( 197.4670779425016 + y * ( -68.5590309679152 + 9.98788038278032 * y ) ) ) - 90.6734234051316 * z ) + \
        1.5 * x * ( -5458.34205214835 - 980.14153344888 * y + \
        ( 4.0 / 3.0 ) * x * ( 2247.60742726704 - 340.1237483177863 * 1.25 * x + 220.542973797483 * y ) + \
        180.142097805543 * z ) + \
        z * ( -219.1676534131548 + ( -16.32775915649044 - 120.7020447884644 * z ) * z ) )

        g08[x>0] = g08[x>0] + ( -7296.43987145382 + z[x>0] * ( 598.378809221703 + \
        z[x>0] * ( -156.8822727844005 + ( 204.1334828179377 - 10.23755797323846 * z[x>0] ) * z[x>0] ) ) + \
        y[x>0] * ( -1480.222530425046 + z[x>0] * ( -525.876123559641 + \
        ( 249.57717834054571 - 88.449193048287 * z[x>0] ) * z[x>0] ) + \
        y[x>0] * ( -129.1994027934126 + z[x>0] * ( 1149.174198007428 + \
        z[x>0] * ( -162.5751787551336 + 76.9195462169742 * z[x>0] ) ) + \
        y[x>0] * ( -30.0682112585625 - 1380.9597954037708 * z[x>0] + \
        y[x>0] * ( 2.626801985426835 + 703.695562834065 * z[x>0] ) ) ) ) ) / x[x>0] + \
        ( 11625.62913253464 + 1702.453469893412 * y[x>0] ) / x2[x>0]

        g08[x==0] = np.nan

        gibbs = 0.25 * sfac**2 * g08

    elif (ns==0) & (nt==2) & (npr==0):
        g03 = -24715.571866078 + z * ( 2910.0729080936 + z * \
        ( -1513.116771538718 + z * ( 546.959324647056 + z * ( -111.1208127634436 + 8.68841343834394 * z ) ) ) ) + \
        y * ( 4420.4472249096725 + z * ( -4035.04669887042 + \
        z * ( 2996.162344914912 + z * ( -1437.2719839264719 + ( 292.8075111563232 - 9.978426372534301 * z ) * z ) ) ) + \
        y * ( -1778.231237203896 + z * ( 4775.621344883664 + \
        z * ( -3621.784567462512 + ( 1826.356460806092 - 316.49805267936244 * z ) * z ) ) + \
        y * ( 1160.5182516851419 + z * ( -3892.3662123519 + \
        z * ( 2410.4130980405 + z * ( -1105.446104680304 + 129.6381336154442 * z ) ) ) + \
        y * ( -569.531539542516 + y * ( 128.13429152494615 - 404.50541014508605 * z ) + \
        z * ( 1905.341809925355 + z * ( -668.691951421377 + 245.11816254543362 * z ) ) ) ) ) )

        g08 = x2 * ( 1760.062705994408 + x * ( -86.1329351956084 + \
        x * ( -137.1145018408982 + y * ( 296.20061691375236 + y * ( -205.67709290374563 + 49.9394019139016 * y ) ) )  + \
        z * ( 766.116132004952 + z * ( -108.3834525034224 + 51.2796974779828 * z ) ) + \
        y * ( -60.136422517125 - 2761.9195908075417 * z + y * ( 10.50720794170734 + 2814.78225133626 * z ) ) ) + \
        y * ( -1351.605895580406 + y * ( 1097.1125373015109 + y * ( -433.20648175062206 + 63.905091254154904 * y ) + \
        z * ( -3572.7449038462437 + ( 896.713693665072 - 437.84750280190565 * z ) * z ) ) + \
        z * ( 4165.4688847996085 + z * ( -1229.337851789418 + ( 681.370187043564 - 66.7696405958478 * z ) * z ) ) ) + \
        z * ( -1721.528607567954 + z * ( 674.819060538734 + \
        z * ( -356.629112415276 + ( 88.4080716616 - 15.84003094423364 * z ) * z ) ) ) )

        gibbs = (g03 + g08) * 0.000625

    elif (ns==0) & (nt==0) & (npr==2):
        g03 = -5089.1530840726 + z * ( 1707.1066706777221 + \
        z * ( -399.7761051043332 + ( 84.0526217606168 - 16.39285534413117 * z ) * z ) ) + \
        y * ( 1552.307223226202 + z * ( -1179.07530528732 + ( 347.75583155301 - 42.658016703665396 * z ) * z ) + \
        y * ( -1513.116771538718 + z * ( 1640.877973941168 + z * ( -666.7248765806615 + 86.8841343834394 * z ) ) + \
        y * ( 998.720781638304 + z * ( -1437.2719839264719 + ( 585.6150223126464 - 33.261421241781 * z ) * z ) + \
        y * ( -603.630761243752 + ( 913.178230403046 - 316.49805267936244 * z ) * z + \
        y * ( 241.04130980405 + y * ( -44.5794634280918 + 49.023632509086724 * z ) + \
        z * ( -331.6338314040912 + 77.78288016926652 * z ) ) ) ) ) )

        g08 = x2 * ( 769.588305957198 + z * ( -579.1945920644748 + ( 190.08980732018878 - 52.4960313181984 * z ) * z ) + \
        x * ( -104.588181856267 + x * ( -8.16387957824522 - 181.05306718269662 * z ) + \
        ( 408.2669656358754 - 40.95023189295384 * z ) * z + \
        y * ( 166.3847855603638 - 176.898386096574 * z + y * ( -108.3834525034224 + 153.8390924339484 * z ) ) ) + \
        y * ( -687.913805923122 + z * ( 748.126026697488 + z * ( -379.883572632876 + 140.9317606630898 * z ) ) + \
        y * ( 674.819060538734 + z * ( -1069.887337245828 + ( 530.4484299696 - 158.40030944233638 * z ) * z ) + \
        y * ( -409.779283929806 + y * ( 149.452282277512 - 218.92375140095282 * z ) + \
        ( 681.370187043564 - 133.5392811916956 * z ) * z ) ) ) )

        gibbs = (g03 + g08) * 1e-16
        # This is the second derivative of the Gibbs function with respect to pressure, measured in Pa.  This derivative has units of (J kg :sup::`-1`) (Pa :sup::`-2`).
    else:
        raise NameError('Wrong Combination of order/variables')

    return gibbs

def  _dsa_add_barrier(dsa, lon, lat, longs_ref, lats_ref, dlongs_ref, dlats_ref):
    """
    Adds a barrier through Central America (Panama) and then averages over the appropriate side of the barrier.

    Parameters
    ----------
    dsa : array_like
          Absolute Salinity anomaly of the 4 adjacent neighbours  [g kg :sup::`-1`]
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
          Absolute Salinity anomaly of the 4 adjacent neighbours  [g kg :sup::`-1`]

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

    for k0 in range(0, len(lon.shape) ):
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

        inds = ( above_line != above_line0 ) # indices of different sides of CA line
        dsa[inds,k0] = np.nan

    dsa_mean = dsa.mean()
    inds_nan = np.where( np.isnan( dsa_mean ) )[0]
    no_nan = len(inds_nan)

    for kk in range(0,no_nan):
        col = inds_nan[kk]
        inds_kk = np.where( np.isnan( dsa[:,col] ) )[0]
        Inn = np.where( ~np.isnan( dsa[:,col] ) )[0]
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
          Absolute Salinity anomaly of the 4 adjacent neighbours  [g kg :sup::`-1`]

    Returns
    -------
    delta_SA : array_like
          Absolute Salinity anomaly of the 4 adjacent neighbours  [g kg :sup::`-1`]

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
    #FIXME: I do not understand the logic for this function

    # Convert input to numpy arrays
    dsa = np.asarray(dsa)

    #FIXME: there must be a better way
    dsa_mean = dsa.mean(axis = 0) #FIXME: should be nanmean here in the original...
    inds_nan = np.where( np.isnan(dsa_mean) )[0]
    no_nan = len(inds_nan)

    for kk in range(0, no_nan):
        col = inds_nan[kk]
        inds_kk = np.where( np.isnan( dsa[:,col] ) )[0]
        Inn = np.where(~np.isnan( dsa[:,col] ) )[0]
        if Inn.size != 0:
            dsa[inds_kk, col] = dsa[Inn,col].mean()

    delta_SA = dsa

    return delta_SA

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

    #Ishallow = np.where( spycnl <= min_spycnl_ref_cast ) # Set equal to the shallowest bottle.

    #SA_iref_cast[Ishallow] = SA_ref_cast[Imin_spycnl_ref_cast]
    #CT_iref_cast[Ishallow] = CT_ref_cast[Imin_spycnl_ref_cast]
    #p_iref_cast[Ishallow] = p_ref_cast[Imin_spycnl_ref_cast]

    #max_spycnl_ref_cast, Imax_spycnl_ref_cast = spycnl_ref_cast.max(), spycnl_ref_cast.argmax()

    #Ideep = np.where( spycnl >= max_spycnl_ref_cast ) # Set equal to the deepest bottle.

    #SA_iref_cast[Ideep] = SA_ref_cast[Imax_spycnl_ref_cast]
    #CT_iref_cast[Ideep] = CT_ref_cast[Imax_spycnl_ref_cast]
    #p_iref_cast[Ideep] = p_ref_cast[Imax_spycnl_ref_cast]

    #I = np.where(spycnl >= 21.805 & spycnl <= 28.3614)

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
               Absolute Salinity anomaly [g kg :sup::`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    TODO

    Notes
    -----
    The Absolute Salinity Anomaly in the Baltic Sea is evaluated separately, since it is a function of Practical Salinity, not of space. The present function returns a delta_SA of zero for data in the Baltic Sea. The correct way of calculating Absolute Salinity in the Baltic Sea is by calling _SA_from_SP.

    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometres inland from the coast.

    Examples
    --------
    TODO

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean.  Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    # Convert input to numpy arrays
    p, lon, lat = np.asarray(p), np.asarray(lon), np.asarray(lat)

    data = pickle.load( open('gsw_data_v2_0.pkl','rb') )

    delta_SA_ref = data['delta_SA_ref']
    lats_ref = data['lats_ref']
    longs_ref = data['longs_ref']
    p_ref = data['p_ref']
    ndepth_ref = data['ndepth_ref']

    dlongs_ref = longs_ref[1] - longs_ref[0]
    dlats_ref = lats_ref[1] - lats_ref[0]

    indsx0 = np.floor( float( longs_ref.size - 1 ) * (lon - longs_ref[0] ) / (longs_ref[-1]- longs_ref[0] ) )
    indsx0 = np.int64(indsx0)
    indsx0[indsx0 == longs_ref.size] = longs_ref.size - 1

    indsy0 = np.floor( float(lats_ref.size - 1 ) * ( lat - lats_ref[0] ) / (lats_ref[-1]- lats_ref[0] ) )
    indsy0 = np.int64(indsy0)
    indsy0[indsy0 == lats_ref.size] = lats_ref.size - 1

    #FIXME: Ugly matlab matrix dot multiplication, there must be a better way...
    P_REF = np.dot( np.ones(p_ref.size)[:,np.newaxis], p[np.newaxis,:] )
    P = np.dot( p_ref[:,np.newaxis], np.ones(p.size)[np.newaxis,:] )
    indsz0 = np.sum( (P_REF >= P), axis=0 )

    nmax = np.c_[ ndepth_ref[indsy0, indsx0], \
                  ndepth_ref[indsy0, indsx0+1], \
                  ndepth_ref[indsy0+1, indsx0+1], \
                  ndepth_ref[indsy0+1, indsx0] ].max(axis=1)


    #FIXME: check if this work when inds1 !=0
    inds1 = np.where(indsz0 > nmax)[0] # casts deeper than GK maximum
    if inds1.size != 0:
        p[inds1] = p_ref[nmax[inds1]] # have reset p here so have to reset indsz0


    #FIXME: Ugly matlab matrix dot multiplication, there must be a better way...
    P_REF = np.dot( np.ones(p_ref.size)[:,np.newaxis], p[np.newaxis,:] )
    P = np.dot( p_ref[:,np.newaxis], np.ones(p.size)[np.newaxis,:] )
    indsz0 = np.sum( (P_REF >= P), axis=0 ) - 1

    inds = (indsz0 == p_ref.size-1)
    indsz0[inds] = p_ref.size - 2

    inds0 = indsz0 + indsy0 * delta_SA_ref.shape[0] + indsx0 * delta_SA_ref.shape[0] * delta_SA_ref.shape[1]

    data_indices = np.c_[indsx0, indsy0, indsz0, inds0]
    data_inds = data_indices[:,2]

    r1 = ( lon - longs_ref[indsx0] ) / np.float64( longs_ref[indsx0+1] - longs_ref[indsx0] )
    s1 = ( lat - lats_ref[indsy0] ) / np.float64( lats_ref[indsy0+1] - lats_ref[indsy0] )
    t1 = ( p - p_ref[indsz0] ) / np.float64( p_ref[indsz0+1] - p_ref[indsz0] )

    nksum = 0
    no_levels_missing = 0

    sa_upper = np.nan * ( np.ones(data_inds.shape) )
    sa_lower = np.nan * ( np.ones(data_inds.shape) )
    delta_SA = np.nan * ( np.ones(data_inds.shape) )
    in_ocean = np.ones( delta_SA.shape )

    for k in range(0, p_ref.size-1):
        inds_k = (indsz0 == k)
        nk = len(inds_k)

        if nk > 0:
            nksum = nksum + nk
            indsx = indsx0[inds_k]
            indsy = indsy0[inds_k]
            indsz = k * np.ones( indsx.shape, dtype='int64' )
            inds_di = (data_inds == k) # level k interpolation
            dsa = np.nan * np.ones( (4, p.size) )

            dsa[0, inds_k] = delta_SA_ref[indsz, indsy, indsx]
            dsa[1, inds_k] = delta_SA_ref[indsz, indsy, indsx+1]   # inds0 + ny*nz
            dsa[2, inds_k] = delta_SA_ref[indsz, indsy+1, indsx+1] # inds0 + ny*nz + nz
            dsa[3, inds_k] = delta_SA_ref[indsz, indsy+1, indsx] #  inds0 + nz

            inds = np.where( (260. <= lon) & (lon <= 295.217) & (0. <= lat) & (lat <= 19.55) & (indsz0 == k) )[0]

            if inds.size !=0: #FIXME: test case when this is True
                dsa[:,inds] = _dsa_add_barrier( dsa[:,inds], lon[inds], \
                lat[inds], longs_ref[indsx0[0][inds[0]]], lats_ref[indsy0[0][inds[0]]], dlongs_ref, dlats_ref)

            inds = np.where( ( np.isnan( np.sum(dsa, axis=0) ) ) & (indsz0==k))[0]

            if inds.size !=0:  #FIXME: Not working !!!
                dsa[:,inds] = _dsa_add_mean(dsa[:,inds])

            sa_upper[inds_di] = ( 1 - s1[inds_di] ) * ( dsa[0, inds_k] + \
            r1[inds_di] * ( dsa[1, inds_k] - dsa[0, inds_k] ) ) + \
            s1[inds_di] * ( dsa[3, inds_k] + \
            r1[inds_di] * ( dsa[2, inds_k] - dsa[3,inds_k] ) ) # level k+1 interpolation

            dsa = np.nan * np.ones( (4, p.size) )
            dsa[0, inds_k] = delta_SA_ref[indsz+1, indsy, indsx]
            dsa[1, inds_k] = delta_SA_ref[indsz+1, indsy, indsx+1] # inds1 + ny*nz
            dsa[2, inds_k] = delta_SA_ref[indsz+1, indsy+1, indsx+1] # inds1 + ny*nz + nz
            dsa[3, inds_k] = delta_SA_ref[indsz+1, indsy+1, indsx] # inds1 + nz

            inds = np.where( (260. <= lon) & (lon <= 295.217) & (0 <= lat) & (lat <= 19.55) & (indsz0 == k) )[0]

            """ TODO: describe add_barrier """
            if inds.size != 0:  #FIXME: test case when this is True
                dsa[:,inds] = _dsa_add_barrier( dsa[:,inds], lon[inds], \
                lat[inds], longs_ref[ndsx0[0][inds[0]]], lats_ref[indsy0[0][inds[0]]], dlongs_ref, dlats_ref)

            inds = ( np.isnan( np.sum(dsa, axis=0) ) ) & (indsz0==k)

            """ TODO: describe add_mean """
            dsa[:,inds] = _dsa_add_mean(dsa[:,inds])

            sa_lower[inds_di] = ( 1 - s1[inds_di] ) * ( dsa[0, inds_k] + \
            r1[inds_di] * ( dsa[1, inds_k] - dsa[0,inds_k] ) ) + \
            s1[inds_di] * ( dsa[3, inds_k] + \
            r1[inds_di] * ( dsa[2, inds_k] - dsa[3, inds_k] ) )

            inds_different = np.isfinite(sa_upper[inds_di]) & np.isnan(sa_lower[inds_di])
            sa_lower[inds_di[inds_different]] = sa_upper[inds_di[inds_different]]

            delta_SA[inds_di] = sa_upper[inds_di] + t1[inds_di] * ( sa_lower[inds_di] - sa_upper[inds_di] )

        else:
            no_levels_missing = no_levels_missing + 1

    inds = ~np.isfinite(delta_SA)
    delta_SA[inds] = 0
    in_ocean[inds] = False # TODO: change to boolean

    return delta_SA, in_ocean