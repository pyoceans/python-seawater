# -*- coding: utf-8 -*-

"""
Derivatives of enthalpy, entropy, CT and pt

Functions:
----------
  
  CT_derivative_SA(SA, pt)
  CT_derivative_pt(SA, pt)
  CT_first_derivatives(SA, pt)  
     first derivatives of Conservative Temperature
  CT_derivative_SA_SA(SA, pt)
  CT_derivative_SA_pt(SA, pt)
  CT_derivative_pt_pt(SA, pt)
  CT_second_derivatives(SA, pt)           
     second derivatives of Conservative Temperature
  enthalpy_derivative_SA(SA, CT, p)
  enthalpy_derivative_CT(SA, CT, p)
  enthalpy_derivative_p(SA, CT, p)
  enthalpy_first_derivatives(SA, CT, p)
     first derivatives of enthalpy
  enthalpy_derivative_SA_SA(SA, CT, p)
  enthalpy_derivative_SA_CT(SA, CT, p)
  enthalpy_derivative_CT_CT(SA, CT, p)
  enthalpy_second_derivatives(SA, CT, p)
     second derivatives of enthalpy
  entropy_derivative_SA(SA, CT)
  entropy_derivative_CT(SA, CT)
  entropy_first_derivatives(SA, CT)      
     first derivatives of entropy
  entropy_derivative_CT_CT(SA, CT)
  entropy_derivative_SA_CT(SA, CT)
  entropy_derivative_SA_SA(SA, CT)
  entropy_second_derivatives(SA, CT)
     second derivatives of entropy
  pt_derivative_SA(SA, CT)
  pt_derivative_CT(SA, CT)
  pt_first_derivatives(SA, CT)            
     first derivatives of potential temperature
  pt_derivative_SA_SA(SA, CT)
  pt_derivative_SA_CT(SA, CT)
  pt_derivative_CT_CT(SA, CT)
  pt_second_derivatives(SA, CT)           
     second derivatives of potential temperature

"""
from __future__ import division


import numpy as np
import seawater.constants as cte
import library as lib
from library import match_args_return
from conversions import *

# -----------

__all__ =  ['CT_derivative_SA',
            'CT_derivative_pt',
            'CT_first_derivatives',
            'CT_derivative_SA_SA',
            'CT_derivative_SA_pt',
            'CT_derivative_pt_pt',
            'CT_second_derivatives',
            'enthalpy_derivative_SA',
            'enthalpy_derivative_CT',
            'enthalpy_derivative_p',
            'enthalpy_first_derivatives',
            'enthalpy_derivative_SA_SA',
            'enthalpy_derivative_SA_CT',
            'enthalpy_derivative_CT_CT',
            'enthalpy_second_derivatives',
            'entropy_derivative_SA',
            'entropy_derivative_CT',
            'entropy_first_derivatives',
            'entropy_derivative_CT_CT',
            'entropy_derivative_SA_CT',
            'entropy_derivative_SA_SA',
            'entropy_second_derivatives',
            'pt_derivative_SA',
            'pt_derivative_CT',
            'pt_first_derivatives',
            'pt_derivative_SA_SA',
            'pt_derivative_SA_CT',
            'pt_derivative_CT_CT',
            'pt_second_derivatives',
           ]

# ----------------

@match_args_return
def CT_derivative_SA(SA, pt):
    r"""
    Calculates the derivatives of Conservative Temperature (CT_SA) with respect
    to Absolute Salinity at constant potential temperature (with pr = 0 dbar).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    CT_SA : array_like
            The derivative of CT with respect to SA at constant potential
            temperature reference sea pressure of 0 dbar.
            [K (g kg :sup:`-1`) :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> pt = [28.7832, 28.4209, 22.7850, 10.2305, 6.8292, 4.3245]
    >>> gsw.CT_derivative_SA(SA, pt)
    array([-0.04198109, -0.04155814, -0.03473921, -0.0187111 , -0.01407594,
           -0.01057172])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.3) and (A.12.9a,b).

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
    Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    Modifications:
    2010-08-05. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1
    abs_pt = cte.Kelvin + pt

    g100 = lib._gibbs(n1, n0, n0, SA, pt, 0)
    g110 = lib._gibbs(n1, n1, n0, SA, pt, 0)
    CT_SA = ( g100 - abs_pt * g110 ) / cte.cp0

    return CT_SA

@match_args_return
def CT_derivative_pt(SA, pt):
    r"""
    Calculates the derivatives of Conservative Temperature (CT_pt) with respect
    to potential temperature (the regular potential temperature which is
    referenced to 0 dbar) at constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    CT_pt : array_like
            The derivative of CT with respect to pt at constant SA.
            [ unitless ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> pt = [28.7832, 28.4209, 22.7850, 10.2305, 6.8292, 4.3245]
    >>> gsw.CT_derivative_pt(SA, pt)
    array([1.00281494,  1.00255482,  1.00164514,  1.00000377,  0.99971636,
           0.99947433])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.3) and (A.12.9a,b).

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
    Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    Modifications:
    2010-08-05. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n2 = 0, 2
    abs_pt = cte.Kelvin + pt

    g020 = lib._gibbs(n0, n2, n0, SA, pt, 0)
    CT_pt = - (abs_pt * g020 ) / cte.cp0

    return CT_pt

def CT_first_derivatives(SA, pt):
    r"""
    Calculates the following two derivatives of Conservative Temperature
    (1) CT_SA, the derivative with respect to Absolute Salinity at constant
        potential temperature (with pr = 0 dbar), and
    (2) CT_pt, the derivative with respect to potential temperature (the
        regular potential temperature which is referenced to 0 dbar) at
        constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar
         [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    CT_SA : array_like
            The derivative of CT with respect to SA at constant potential
            temperature reference sea pressure of 0 dbar.
            [K (g kg :sup:`-1`) :sup:`-1`]

    CT_pt : array_like
            The derivative of CT with respect to pt at constant SA.
            [ unitless ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> pt = [28.7832, 28.4209, 22.7850, 10.2305, 6.8292, 4.3245]
    >>> gsw.CT_first_derivatives(SA, pt)
    array([[-0.04198109, -0.04155814, -0.03473921, -0.0187111 , -0.01407594,
            -0.01057172],
           [ 1.00281494,  1.00255482,  1.00164514,  1.00000377,  0.99971636,
             0.99947433]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.3) and (A.12.9a,b).

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
    Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
    expression for the density of seawater in terms of Conservative Temperature,
    and related properties of seawater.

    Modifications:
    2010-08-05. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    return CT_derivative_SA(SA, pt), CT_derivative_pt(SA, pt)

# ----------------------------------------------

@match_args_return
def CT_derivative_SA_SA(SA, pt):

    dSA = 1e-3
    SA_l = SA - dSA
    SA_l = SA_l.clip(0.0, np.inf)
    SA_u = SA + dSA

    CT_SA_l = CT_derivative_SA(SA_l, pt)
    CT_SA_u = CT_derivative_SA(SA_u, pt)

    return (CT_SA_u - CT_SA_l) / (SA_u - SA_l)

# ------------------------

@match_args_return
def CT_derivative_SA_pt(SA, pt):

    dpt = 1e-2
    pt_l = pt - dpt
    pt_u = pt + dpt

    CT_SA_l = CT_derivative_SA(SA, pt_l)
    CT_SA_u = CT_derivative_SA(SA, pt_u)

    return (CT_SA_u - CT_SA_l) / (pt_u - pt_l)

# -------------------------------

@match_args_return
def CT_derivative_pt_pt(SA, pt):

    dpt = 1e-2
    pt_l = pt - dpt
    pt_u = pt + dpt

    CT_pt_l = CT_derivative_pt(SA, pt_l)
    CT_pt_u = CT_derivative_pt(SA, pt_u)

    return (CT_pt_u - CT_pt_l) / (pt_u - pt_l)

# ----------------------------------

def CT_second_derivatives(SA, pt):

    return ( CT_derivative_SA_SA(SA, pt),
             CT_derivative_SA_pt(SA, pt),
             CT_derivative_pt_pt(SA, pt) )

# -----------------------------------------------

@match_args_return
def enthalpy_derivative_SA(SA, CT, p):

    pt0 = pt_from_CT(SA, CT)
    t   = pt_from_t(SA, pt0, 0, p)
    temp_ratio = (273.15 + t) / (273.15 + pt0)
    return ( lib._gibbs(1, 0, 0, SA, t, p) - 
             temp_ratio * lib._gibbs(1, 0, 0, SA, pt0, 0) )


@match_args_return
def enthalpy_derivative_CT(SA, CT, p):

    cp0 = 3991.86795711963
    pt0 = pt_from_CT(SA, CT)
    t   = pt_from_t(SA, pt0, 0, p)
    temp_ratio = (273.15 + t) / (273.15 + pt0)
    return cp0 * temp_ratio

@match_args_return
def enthalpy_derivative_p(SA, CT, p):
    pt0 = pt_from_CT(SA, CT)
    t   = pt_from_t(SA, pt0, 0, p)
    return lib._gibbs(0, 0, 1, SA, t, p)

def enthalpy_first_derivatives(SA, CT, p):

    return ( enthalpy_derivative_SA(SA, CT, p),
             enthalpy_derivative_CT(SA, CT, p),
             enthalpy_derivative_p(SA, CT, p) )

# ----------------------------------------------

@match_args_return
def enthalpy_derivative_SA_SA(SA, CT, p):

    cp0 = 3991.86795711963

    pt0 = pt_from_CT(SA, CT)
    abs_pt0 = 273.15 + pt0
    t = pt_from_t(SA, pt0, 0, p)
    temp_ratio = (273.15 + t) / abs_pt0

    rec_gTT_pt0 = 1.0 / lib._gibbs(0, 2, 0, SA, pt0, 0)
    rec_gTT_t   = 1.0 / lib._gibbs(0, 2, 0, SA, t, p)
    gST_pt0  = lib._gibbs(1, 1, 0, SA, pt0, 0)
    gST_t    = lib._gibbs(1, 1, 0, SA, t, p)
    gS_pt0   = lib._gibbs(1, 0, 0, SA, pt0, 0)

    part = (temp_ratio * gST_pt0 * rec_gTT_pt0 - gST_t * rec_gTT_t) / (abs_pt0)
    factor = gS_pt0 / cp0

    h_CT_CT = cp0 * cp0 * ( (temp_ratio * rec_gTT_pt0 - rec_gTT_t)
                            / (abs_pt0 * abs_pt0) )

    return ( lib._gibbs(2, 0, 0, SA, t, p) -
             temp_ratio * lib._gibbs(2, 0, 0, SA, pt0, 0) +
             temp_ratio * gST_pt0 * gST_pt0 * rec_gTT_pt0 -
             gST_t * gST_t * rec_gTT_t - 
             2.0 * gS_pt0 * part +
             factor * factor * h_CT_CT )

# --------

@match_args_return
def enthalpy_derivative_SA_CT(SA, CT, p):

    cp0 = 3991.86795711963

    pt0 = pt_from_CT(SA, CT)
    abs_pt0 = 273.15 + pt0
    t = pt_from_t(SA, pt0, 0, p)
    temp_ratio = (273.15 + t) / abs_pt0

    rec_gTT_pt0 = 1.0 / lib._gibbs(0, 2, 0, SA, pt0, 0)
    rec_gTT_t   = 1.0 / lib._gibbs(0, 2, 0, SA, t, p)
    gST_pt0  = lib._gibbs(1, 1, 0, SA, pt0, 0)
    gST_t    = lib._gibbs(1, 1, 0, SA, t, p)
    gS_pt0   = lib._gibbs(1, 0, 0, SA, pt0, 0)

    part = (temp_ratio * gST_pt0 * rec_gTT_pt0 - gST_t * rec_gTT_t) / (abs_pt0)
    factor = gS_pt0 / cp0

    h_CT_CT = cp0 * cp0 * ( (temp_ratio * rec_gTT_pt0 - rec_gTT_t)
                             / (abs_pt0 * abs_pt0) )

    return cp0 * part - factor * h_CT_CT

# -------------

@match_args_return
def enthalpy_derivative_CT_CT(SA, CT, p):

    cp0 = 3991.86795711963

    pt0 = pt_from_CT(SA, CT)
    abs_pt0 = 273.15 + pt0
    t = pt_from_t(SA, pt0, 0, p)
    temp_ratio = (273.15 + t) / abs_pt0

    rec_gTT_pt0 = 1.0 / lib._gibbs(0, 2, 0, SA, pt0, 0)
    rec_gTT_t   = 1.0 / lib._gibbs(0, 2, 0, SA, t, p)
    #gST_pt0  = lib._gibbs(1, 1, 0, SA, pt0, 0)
    #gST_t    = lib._gibbs(1, 1, 0, SA, t, p)
    #gS_pt0   = lib._gibbs(1, 0, 0, SA, pt0, 0)

    h_CT_CT = cp0 * cp0 * ( (temp_ratio * rec_gTT_pt0 - rec_gTT_t)
             / (abs_pt0 * abs_pt0) )

    return h_CT_CT

# ---------------
    
def enthalpy_second_derivatives(self):
    return ( enthalpy_derivative_SA_SA(SA, CT, p),
             enthalpy_derivative_SA_CT(SA, CT, p),
             enthalpy_derivative_CT_CT(SA, CT, p) )

# -------------------------------------------------


@match_args_return
def entropy_derivative_SA(SA, CT):
    r"""
    Calculates the derivatives of specific entropy (eta_SA) with respect to
    Absolute Salinity at constant Conservative Temperature.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_SA : array_like
             The derivative of specific entropy with respect to SA at constant
             CT [J g :sup:`-1` K :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_derivative_SA(SA, CT)
    array([ -0.2632868 ,  -0.26397728,  -0.2553675 ,  -0.23806659,
            -0.23443826,  -0.23282068])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.8) and (P.14a,c).

    Modifications:
    2010-08-21. Trevor McDougall.
    2011-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1

    pt = pt_from_CT(SA, CT)

    eta_SA = -(lib._gibbs(n1, n0, n0, SA, pt, 0) ) / (cte.Kelvin + pt)

    return eta_SA


@match_args_return
def entropy_derivative_CT(SA, CT):
    r"""
    Calculates the derivatives of specific entropy (eta_CT) with respect to
    Conservative Temperature at constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_CT : array_like
             The derivative of specific entropy with respect to CT at constant
             SA [ J (kg K :sup:`-2`) :sup:`-1` ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_derivative_CT(SA, CT)
    array([13.22103121,  13.23691119,  13.48900463,  14.08659902,
           14.25772958,  14.38642995])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.8) and (P.14a,c).

    Modifications:
    2010-08-21. Trevor McDougall.
    2011-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    pt = pt_from_CT(SA, CT)

    eta_CT = cte.cp0 / (cte.Kelvin + pt)

    return eta_CT


def entropy_first_derivatives(SA, CT):
    r"""
    Calculates the following two partial derivatives of specific entropy (eta)
    (1) eta_SA, the derivative with respect to Absolute Salinity at constant
        Conservative Temperature, and
    (2) eta_CT, the derivative with respect to Conservative Temperature at
        constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_SA : array_like
             The derivative of specific entropy with respect to SA at constant
             CT [J g :sup:`-1` K :sup:`-1`]
    eta_CT : array_like
             The derivative of specific entropy with respect to CT at constant
             SA [ J (kg K :sup:`-2`) :sup:`-1` ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_first_derivatives(SA, CT)
    array([[ -0.2632868 ,  -0.26397728,  -0.2553675 ,  -0.23806659,
             -0.23443826,  -0.23282068],
           [ 13.22103121,  13.23691119,  13.48900463,  14.08659902,
             14.25772958,  14.38642995]])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.12.8) and (P.14a,c).

    Modifications:
    2010-08-21. Trevor McDougall.
    2011-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    return ( entropy_derivative_SA(SA, CT),
             entropy_derivative_CT(SA, CT) )


@match_args_return
def entropy_derivative_CT_CT(SA, CT):
    r"""
    Calculates the the second-order partial derivatives of specific entropy
    with respect to Conservative Temperature at constant Absolute
    Salinity (eta_CT_CT).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_CT_CT : array_like
                The second derivative of specific entropy with respect to CT at
                constant SA [J (kg K :sup:`3`) :sup:`-1` ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_derivative_CT_CT(SA, CT)
    array([-0.04366502, -0.04378134, -0.04550611, -0.04970894, -0.05093869,
           -0.05187502])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (P.14b) and (P.15a,b).

    Modifications:
    2010-08-23. Trevor McDougall and Paul Barker.
    2011-02-28. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n2 = 0, 2

    pt = pt_from_CT(SA, CT)
    abs_pt = 273.15 + pt

    CT_pt = -(abs_pt * lib._gibbs(n0, n2, n0, SA, pt, 0) ) / cte.cp0

    eta_CT_CT = -cte.cp0 / ( CT_pt * abs_pt**2 )

    return eta_CT_CT


@match_args_return
def entropy_derivative_SA_CT(SA, CT):
    r"""
    Calculates the the second-order partial derivatives of specific entropy
    with to Absolute Salinity and Conservative Temperature (eta_SA_CT).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_SA_CT : array_like
                The second derivative of specific entropy with respect to
                SA and CT [J (kg (g kg :sup:`-1` ) K :sup:`2`) :sup:`-1` ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_derivative_SA_CT(SA, CT)
    array([-0.0018331 , -0.00181947, -0.00158084, -0.00093011, -0.00071701,
           -0.00054841])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (P.14b) and (P.15a,b).

    Modifications:
    2010-08-23. Trevor McDougall and Paul Barker.
    2011-02-28. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1 = 0, 1

    pt = pt_from_CT(SA, CT)
    abs_pt = 273.15 + pt

    CT_SA = (lib._gibbs(n1, n0, n0, SA, pt, 0) -
                    ( abs_pt * lib._gibbs(n1, n1, n0, SA, pt, 0 ) ) ) / cte.cp0

    eta_SA_CT = - CT_SA * entropy_derivative_CT_CT(SA, CT)

    return eta_SA_CT


@match_args_return
def entropy_derivative_SA_SA(SA, CT):
    r"""
    Calculates the the second-order partial derivatives of specific entropy
    with respect to Absolute Salinity at constant Conservative
    Temperature (eta_SA_SA).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_SA_SA : array_like
                The second derivative of specific entropy with respect to SA at
                constant CT [J (kg K (g kg :sup:`-1` ) :sup:`2`) :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_derivative_SA_SA(SA, CT)
    TODO

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (P.14b) and (P.15a,b).

    Modifications:
    2010-08-23. Trevor McDougall and Paul Barker.
    2011-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    n0, n1, n2 = 0, 1, 2

    pt = pt_from_CT(SA, CT)
    abs_pt = 273.15 + pt

    CT_SA = (lib._gibbs(n1, n0, n0, SA, pt, 0) -
                    ( abs_pt * lib._gibbs(n1, n1, n0, SA, pt, 0 ) ) ) / cte.cp0

    eta_SA_CT = entropy_derivative_SA_CT(SA, CT)

    eta_SA_SA = -lib._gibbs(n2, n0, n0, SA, pt, 0) / abs_pt - CT_SA * eta_SA_CT

    return eta_SA_SA


def entropy_second_derivatives(SA, CT):
    r"""
    Calculates the following three second-order partial derivatives of specific
    entropy (eta)
    (1) eta_SA_SA, the second derivative with respect to Absolute Salinity at
        constant Conservative Temperature, and
    (2) eta_SA_CT, the derivative with respect to Absolute Salinity and
        Conservative Temperature.
    (3) eta_CT_CT, the second derivative with respect to Conservative
        Temperature at constant Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    eta_SA_SA : array_like
                The second derivative of specific entropy with respect to SA at
                constant CT [J (kg K (g kg :sup:`-1` ) :sup:`2`) :sup:`-1`]
    eta_SA_CT : array_like
                The second derivative of specific entropy with respect to
                SA and CT [J (kg (g kg :sup:`-1` ) K :sup:`2`) :sup:`-1` ]
    eta_CT_CT : array_like
                The second derivative of specific entropy with respect to CT at
                constant SA [J (kg K :sup:`3`) :sup:`-1` ]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_second_derivatives(SA, CT)
    TODO

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (P.14b) and (P.15a,b).

    Modifications:
    2010-08-23. Trevor McDougall and Paul Barker.
    2011-02-25. Filipe Fernandes, Python translation from gsw toolbox.
    """

    return ( entropy_derivative_SA_SA(SA, CT),
             entropy_derivative_SA_CT(SA, CT),
             entropy_derivative_CT_CT(SA, CT),
           )



# -------------------------------------

@match_args_return
def pt_derivative_SA(SA, CT):
    cp0 = 3991.86795711963
    pt = pt_from_CT(SA, CT)
    abs_pt = 273.15 + pt

    CT_SA = (lib._gibbs(1, 0, 0, SA, pt, 0) -
             abs_pt * lib._gibbs(1, 1, 0, SA, pt, 0)) / cp0
    CT_pt = - (abs_pt * lib._gibbs(0, 2, 0,SA, pt, 0)) / cp0

    return - CT_SA / CT_pt

@match_args_return
def pt_derivative_CT(SA, CT):
    cp0 = 3991.86795711963
    pt = pt_from_CT(SA, CT)
    abs_pt = 273.15 + pt

    CT_pt = - (abs_pt * lib._gibbs(0, 2, 0,SA, pt, 0)) / cp0

    return 1.0 / CT_pt

def pt_first_derivatives(SA, CT):
    return ( pt_derivative_SA(SA, CT),
             pt_derivative_CT(SA, CT) )

# -----------------------------------

@match_args_return
def pt_derivative_SA_SA(SA, CT):

    dSA = 1e-3
    SA_l = SA - dSA
    SA_l = SA_l.clip(0.0, np.inf)
    SA_u = SA + dSA
    pt_SA_l = pt_derivative_SA(SA_l, CT)
    pt_SA_u = pt_derivative_SA(SA_u, CT)
 
    return (pt_SA_u - pt_SA_l) / (SA_u - SA_l)


@match_args_return
def pt_derivative_SA_CT(SA, CT):

    dCT  = 1e-2
    CT_l = CT - dCT
    CT_u = CT + dCT
    pt_SA_l = pt_derivative_SA(SA, CT_l)
    pt_SA_u = pt_derivative_SA(SA, CT_u)
    return (pt_SA_u - pt_SA_l) / (CT_u - CT_l) 

    # Alternative
    #dSA = 1e-3
    #SA_l = SA - dSA
    #SA_l = SA_l.clip(0.0, np.inf)
    #SA_u = SA + dSA
    #pt_CT_l = pt_derivative_CT(SA_l, CT)
    #pt_CT_u = pt_derivative_CT(SA_u, CT)
    #return (pt_CT_u - pt_CT_l) / (SA_u - SA_l)

@match_args_return
def pt_derivative_CT_CT(SA, CT):

    dCT  = 1e-2
    CT_l = CT - dCT
    CT_u = CT + dCT
    pt_CT_l = pt_derivative_CT(SA, CT_l)
    pt_CT_u = pt_derivative_CT(SA, CT_u)
    
    return (pt_CT_u - pt_CT_l) / (CT_u - CT_l)

def pt_second_derivatives(SA, CT):
    return ( pt_derivative_SA_SA(SA, CT),
             pt_derivative_SA_CT(SA, CT),
             pt_derivative_CT_CT(SA, CT) )


# -----------------------------------

if __name__=='__main__':
    import doctest
    doctest.testmod()
