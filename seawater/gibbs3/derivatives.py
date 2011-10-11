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


from library import gibbs
from  seawater.constants import cp0, Kelvin
from utilities import match_args_return
from conversions import pt_from_CT, pt_from_t

__all__ = [
           #'CT_derivative_SA',
           #'CT_derivative_pt',
           #'CT_first_derivatives',
           #'CT_derivative_SA_SA',
           #'CT_derivative_SA_pt',
           #'CT_derivative_pt_pt',
           #'CT_second_derivatives',
           'enthalpy_derivative_SA',
           'enthalpy_derivative_CT',
           'enthalpy_derivative_p',
           'enthalpy_first_derivatives',
           #'enthalpy_derivative_SA_SA',
           #'enthalpy_derivative_SA_CT',
           #'enthalpy_derivative_CT_CT',
           #'enthalpy_second_derivatives',
           #'entropy_derivative_SA',
           #'entropy_derivative_CT',
           #'entropy_first_derivatives',
           'entropy_derivative_CT_CT',
           'entropy_derivative_SA_CT',
           'entropy_derivative_SA_SA',
           'entropy_second_derivatives',
           #'pt_derivative_SA',
           #'pt_derivative_CT',
           #'pt_first_derivatives',
           #'pt_derivative_SA_SA',
           #'pt_derivative_SA_CT',
           #'pt_derivative_CT_CT',
           #'pt_second_derivatives',
          ]


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
    abs_pt = Kelvin + pt

    CT_pt = -(abs_pt * gibbs(n0, n2, n0, SA, pt, 0)) / cp0

    eta_CT_CT = - cp0 / (CT_pt * abs_pt ** 2)

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
    abs_pt = Kelvin + pt

    CT_SA = (gibbs(n1, n0, n0, SA, pt, 0) -
                    (abs_pt * gibbs(n1, n1, n0, SA, pt, 0))) / cp0

    eta_SA_CT = -CT_SA * entropy_derivative_CT_CT(SA, CT)

    return eta_SA_CT


@match_args_return
def entropy_derivative_SA_SA(SA, CT):
    r"""
    Calculates the the second-order partial derivatives of specific entropy
    with respect to Absolute Salinity at constant Conservative
    Temperature (eta_SA_SA.)

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
    abs_pt = Kelvin + pt

    CT_SA = (gibbs(n1, n0, n0, SA, pt, 0) -
                    (abs_pt * gibbs(n1, n1, n0, SA, pt, 0))) / cp0

    eta_SA_CT = entropy_derivative_SA_CT(SA, CT)

    eta_SA_SA = -gibbs(n2, n0, n0, SA, pt, 0) / abs_pt - CT_SA * eta_SA_CT

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

    return (entropy_derivative_SA_SA(SA, CT),
            entropy_derivative_SA_CT(SA, CT),
            entropy_derivative_CT_CT(SA, CT),)


@match_args_return
def enthalpy_derivative_SA(SA, CT, p):
    """Part of enthalpy_first_derivatives."""
    # FIXME: The gsw 3.0 has the gibbs derivatives "copy-and-pasted" here
    # instead of the calls to the library! Why?
    n0, n1 = 0, 1
    pt0 = pt_from_CT(SA, CT)
    t = pt_from_t(SA, pt0, 0, p)
    temp_ratio = (Kelvin + t) / (Kelvin + pt0)
    return (gibbs(n1, n0, n0, SA, t, p) -
             temp_ratio * gibbs(n1, n0, n0, SA, pt0, 0))


@match_args_return
def enthalpy_derivative_CT(SA, CT, p):
    """Part of enthalpy_first_derivatives."""
    # FIXME: The gsw 3.0 has the gibbs derivatives "copy-and-pasted" here
    # instead of the calls to the library! Why?
    pt0 = pt_from_CT(SA, CT)
    t = pt_from_t(SA, pt0, 0, p)
    temp_ratio = (Kelvin + t) / (Kelvin + pt0)
    h_CT = cp0 * temp_ratio

    return h_CT


@match_args_return
def enthalpy_derivative_p(SA, CT, p):
    """Part of enthalpy_first_derivatives."""
    # FIXME: The gsw 3.0 has the gibbs derivatives "copy-and-pasted" here
    # instead of the calls to the library! Why?
    n0, n1 = 0, 1
    pt0 = pt_from_CT(SA, CT)
    t = pt_from_t(SA, pt0, 0, p)
    return gibbs(n0, n0, n1, SA, t, p)


def enthalpy_first_derivatives(SA, CT, p):
    r"""
    Calculates the following three derivatives of specific enthalpy (h)
    (1) h_SA, the derivative with respect to Absolute Salinity at
        constant CT and p, and
    (2) h_CT, derivative with respect to CT at constant SA and p.
    (3) h_P, derivative with respect to pressure (in Pa) at constant SA and CT.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    h_SA : array_like
           The first derivative of specific enthalpy with respect to Absolute
           Salinity at constant CT and p. [J/(kg (g/kg))]  i.e. [J/g]
    h_CT : array_like
           The first derivative of specific enthalpy with respect to CT at
           constant SA and p. [J/(kg K)]
    h_P : array_like
          The first partial derivative of specific enthalpy with respect to
          pressure (in Pa) at fixed SA and CT.  Note that h_P is specific
          volume (1/rho.)

    See Also
    --------
    TODO

    Notes
    -----
    TODO


    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See Eqns. (A.11.18), (A.11.15) and (A.11.12.)

    Modifications:
    2010-09-24. Trevor McDougall.
    """

    return (enthalpy_derivative_SA(SA, CT, p),
            enthalpy_derivative_CT(SA, CT, p),
            enthalpy_derivative_p(SA, CT, p),)
