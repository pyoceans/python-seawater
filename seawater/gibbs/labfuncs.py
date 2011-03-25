# -*- coding: utf-8 -*-

"""
Functions for laboratory use

Functions:
----------

  SA_from_rho
     Absolute Salinity from density measurements
  #sigma0_pt                       
  #   sigma_0 expressed in terms of SA & pt0 with ref. press. at 0 dbar

"""

from __future__ import division

import numpy as np
#import seawater.constants as cte
import library as lib
from library import match_args_return
#from conversions import *

# ---------------------------

__all__ = ['SA_from_rho',
           #'sigma0_pt',
          ]

# ----------------------------

@match_args_return
def SA_from_rho(rho, t, p):
    r"""
    Calculates the Absolute Salinity of a seawater sample, for given values of
    its density, in situ temperature and sea pressure (in dbar).

    One use for this function is in the laboratory where a measured value of
    the in situ density :math:`\rho` of a seawater sample may have been made at
    the laboratory temperature :math:`t` and at atmospheric pressure :math:`p`.
    The present function will return the Absolute Salinity SA of this seawater
    sample.

    Parameters
    ----------
    rho : array_like
          in situ density [kg m :sup:`-3`]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    This is expressed on the Reference-Composition Salinity Scale of
    Millero et al. (2008).

    After two iterations of a modified Newton-Raphson iteration, the error in SA
    is typically no larger than 2 :math:`^\times` 10 :sup:`-13` [g kg :sup:`-1`]

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> rho = [1021.839, 1022.262, 1024.426, 1027.792, 1029.839, 1032.002]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.SA_from_rho(rho, t, p)
    array([ 34.71022966,  34.89057683,  35.02332421,  34.84952096,
            34.73824809,  34.73188384])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.5.

    .. [2] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008:
    The composition of Standard Seawater and the definition of the
    Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.

    Modifications:
    2010-08-23. Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """
    n0, n1 = 0, 1
    v_lab = np.ones( rho.shape ) / rho
    v_0 = lib._gibbs(n0, n0, n1, 0, t, p)
    v_120 = lib._gibbs(n0, n0, n1, 120, t, p)
    SA = 120 * ( v_lab - v_0 ) / ( v_120 - v_0 ) # initial estimate of SA
    Ior = (SA < 0) | (SA > 120)
    v_SA = ( v_120 - v_0 ) / 120 # initial estimate of v_SA, SA derivative of v

    for iter in range(0,2):
        SA_old = SA
        delta_v = lib._gibbs(n0, n0, n1, SA_old, t, p) - v_lab
        SA = SA_old - delta_v / v_SA # half way through the modified N-R method
        SA_mean = 0.5 * ( SA + SA_old )
        v_SA = lib._gibbs(n1, n0, n1, SA_mean, t, p)
        SA = SA_old - delta_v / v_SA

    SA[Ior] = np.ma.masked

    return SA

# ----------------------------


if __name__=='__main__':
    import doctest
    doctest.testmod()
