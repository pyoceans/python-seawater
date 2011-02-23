# -*- coding: utf-8 -*-
"""
SEAWATER - originally developed by Phil Morgan, CSIRO

Translated from matlab CSIRO seawater toolbox Version 3.3

Filipe P. A. Fernandes
e-mail:   ocefpaf@gmail.com
web:      http://ocefpaf.tiddlyspot.com/
date:     14-Jan-2010
modified: Mon 21 Feb 2011 01:49:02 PM EST
obs:      fixme (flags what needs attention)
          some keywords and default values are hardcoded!!!
          create proper docstrings

SOFTWARE LICENCE AGREEMENT

1.0 Grant of Licence

1.1  The CSIRO Division of  Oceanography (herein referred to as
     "CSIRO") hereby grants you (hereinafter  referred to  as
     the  "Licencee"),  subject  to  the  Licencee agreeing  to
     comply with the terms and  conditions of this Agreement, a
     non-transferable,   non-exclusive  licence   to  use   the
     computer  programs described in this document (hereinafter
     referred to  as the   "Software")  for  the  purpose   of
     the  Licencee's computing activity.

1.2  CSIRO hereby grants the Licencee the right to  make copies
     of  the  Software  for   the  purpose  of  the  Licencee's
     computing activity only.

1.3  The benefit of  the rights granted to the Licencee  by the
     Licence and this Agreement  generally shall be personal to
     the Licencee and the  Licencee shall not mortgage, charge,
     assign,  rent,  lease,  sell or  otherwise  dispose  of or
     transfer the same or any part to any third party.

1.4  Unless  otherwise agreed  in  writing or  provided  for in
     this  Agreement, CSIRO  shall  be under  no  obligation or
     responsibility  to provide the Licencee with any training,
     maintenance  services,  enhancements  or  updates  of  the
     Software or any services whatsoever.

2.0 Acknowledgment by the Licencee

2.1  The Licencee acknowledges and agrees that it shall not:

     (i)   sell,  let for  hire or  by way  of  trade, offer  or
           exhibit  or  expose  for sale  or  hire  or otherwise
       distribute the Software for  the purposes of trade or
           any other purpose;

     (ii)  authorise  or assist any third person  to do any
           of the acts set out in (i)  above;

     (iii) modify  the  Software source  code  without  advising
           CSIRO.

2.2 The Licencee agrees that:

     (a)  CSIRO  is  the  owner  of  all  copyright  and  other
          Intellectual  Property   Rights  subsisting  in   the
          Software;

     (b)  this   document  must  be   properly   cited  in  any
          publication  reporting  results  derived   from  this
          document or obtained from application and use of this
          software. Any  of   the  Licencee's  documentation
          describing  results  generated  by  the  Licencee's
          use  of  the Software will contain an acknowledgement
          of CSIRO's  ownership of the Software;

     (c)  CSIRO reserves all rights  in the Software other than
          the  rights   granted  to   the   Licencee  by   this
          Agreement;

     (d)  each  item  of  the Software  will  display  a  banner
          summarising   the  terms   of   this   Agreement  and
          acknowledging  the source  of  the Software,  and the
          contents of  a banner  will not be  modified and  its
          display  will  not  be  inactivated  by  the Licencee
          without the approval of CSIRO.

3.0 Indemnity

3.1  To the full  extent permitted by  law, CSIRO  excludes any
     and  all liability  in  respect  of any  loss  or  damage,
     whether  personal  (includes  death  or  illness)  or  of
     property  and  whether direct,  consequential  or  special
     (including  consequential  financial  loss or  damage) of
     the Licencee, its  officers, agents and  employees or  any
     third party  howsoever caused,  which may  be suffered  or
     incurred  or which  may  arise directly  or  indirectly in
     respect of  or  arising  out  of  the  Licencee's  use  or
     inability to use  the Software or the failure or  omission
     on the  part of CSIRO  to comply  with the conditions  and
     warranties under  this  Licence  Agreement.    Insofar  as
     liability for  loss or damages  under or  pursuant to such
     legislation cannot  be  excluded,  CSIRO's  liability  for
     loss or  damages shall  be limited  to the  amount of  One
     Dollar ($1.00).

3.2  CSIRO  make  no  warranties,  expressed  or  implied,  and
     excludes all  other warranties  representations, terms  or
     conditions, whether  express or implied,  oral or written,
     statutory  or  otherwise,  relating  in  any  way  to  the
     Software, or  to  this  Agreement, including  any  implied
     warranty of merchantability  or of fitness for  particular
     purpose.   To the full extent permitted  by the law of the
     Commonwealth of  Australia  or the  laws of  any State  or
     Territory  of  Australia,  any  conditions  or  warranties
     imposed by such  legislation are hereby  excluded.   In so
     far as  liability under  or pursuant  to such  legislation
     may not  be excluded,  CSIRO's liability  to the  Licencee
     pursuant to this Agreement shall be limited as  set out in
     clause 3.1 hereof.

3.3  The  Licencee acknowledges  and agrees  that  the Software
     was developed  for CSIRO  research purposes  and may  have
     inherent  defects, errors or deficiencies, and  that it is
     the  responsibility  of  the   Licencee  to  make  its  own
     assessment  of the  suitability  of the  Software  for the
     purpose  of  the   Licencee's  computing  activity.    The
     Licencee will use  the Software, and  advice, opinions  or
     information supplied by CSIRO,  its officers, employees or
     agents  concerning  the  Software  at the  Licencee's  own
     risk.

3.4  The  Licencee hereby  releases and  indemnifies and  shall
     continue to  release and  indemnify  CSIRO, its  officers,
     employees  and  agents  from  and  against  all   actions,
     claims, proceedings  or demands  (including those  brought
     by third parties) which may be bought against  it or them,
     whether  on their  own or  jointly  with the  Licencee and
     whether at common law, in equity or pursuant to statute or
     otherwise, in respect of  any loss, death, injury, illness
     or  damage  (whether  personal  or  property,  and whether
     direct   or    consequential,   including    consequential
     financial  loss)   and  any  infringement  of   copyright,
     patents,   trade  marks,  designs  or  other  Intellectual
     Property Rights, howsoever  arising out of the  Licencee's
     exercise of its  rights under this  Agreement and from and
     against  all  damages,  costs  and  expenses  incurred  in
     defending  or  settling  any  such  claim,  proceeding  or
     demand.

3.5  The  Licencee's  obligation to  indemnify  CSIRO  and  its
     officers,  employees  and  agents set  out  in  clause 3.4
     hereof   is   a   continuing   obligation   separate  from
     and independent of the Licencee's other obligations  under
     this  Agreement,  and  shall  survive  all  expiration  or
     termination of this Agreement.

4.0 Termination

4.1  The Licence shall terminate  immediately upon the Licencee
     breaching any term or  condition of this Agreement whether
     or  not CSIRO is aware of  the occurrence of the breach at
     the time that it happens.

4.2  CSIRO  may terminate the Licence on  reasonable grounds by
     notice in  writing  to the  Licencee, and  such notice  of
     termination shall  be effective  immediately upon  receipt
     by  the  Licencee.

DISCLAIMER:
  This software is provided "as is" without warranty of any kind.
  See the file copy for conditions of use and licence.


  Computational routines for the properties of sea water

SEAWATER - developed by Phil Morgan and
           Lindsay Pender (Lindsay.Pender@csiro.au) CSIRO


DESCRIPTION:
   SEAWATER is a toolkit of PYTHON routines for calculating the
   properties of sea water. They are a self contained library and
   are extremely easy to use and will run on all computers that
   support NUMPY.
"""

import numpy as np
from seawater import constants as cte

# only used on the test routine
import sys
from platform import uname
from time  import asctime, localtime

__version__ = '3.3' # matlab version

def T68conv(T90):
    r"""
    Convert ITS-90 temperature to IPTS-68

    :math:`T68  = T90 * 1.00024`

    Parameters
    ----------
    t : array_like
           temperature [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    t : array_like
           temperature [:math:`^\circ` C (IPTS-68)]

    See Also
    --------
    TODO

    Notes
    -----
    The International Practical Temperature Scale of 1968 (IPTS-68) need to be
    correct to the ITS-90. This linear transformation is accurate within
    0.5 :math:`^\circ` C for conversion between IPTS-68 and ITS-90 over the
    oceanographic temperature range.

    Examples
    --------
    >>> import seawater.csiro as sw
    >>> sw.T68conv(19.995201151723585)
    20.0

    References
    ----------
    .. [1] Saunders, P. M., 1991: The International Temperature Scale of 1990,
    ITS-90. WOCE Newsletter, No. 10, WOCE International Project Office,
    Southampton, United Kingdom, 10.

    Modifications: Filipe Fernandes, 2010
                   #10-11-24. Filipe Fernandes, first version.
    """

    T90 = np.asanyarray(T90)

    T68 = T90 * 1.00024
    return T68

def T90conv(t, t_type='T68'):
    r"""
    Convert IPTS-68 or IPTS-48 to temperature to ITS-90.

    T48 apply to all data collected prior to 31/12/1967.
    T68 apply to all data collected between 01/10/1968 and 31/12/1989.


    ..math:
        T90 = T68 / 1.00024

        T90 = T48 - (4.4e-6) * T48 * (100-T48) ) / 1.00024

    Parameters
    ----------
    t : array_like
           temperature [:math:`^\circ` C (IPTS-68) or (IPTS-48)]
    t_type : string, optional
            'T68' (default) or 'T48'

    Returns
    -------
    T90 : array_like
           temperature [:math:`^\circ` C (ITS-90)]

    See Also
    --------
    TODO

    Notes
    -----
    The International Practical Temperature Scale of 1968 (IPTS-68) need to be
    correct to the ITS-90. This linear transformation is accurate within
    0.5 :math:`^\circ` C for conversion between IPTS-68 and ITS-90 over the
    oceanographic temperature range.

    Examples
    --------
    >>> import seawater.csiro as sw
    >>> sw.T90conv(20.004799999999999)
    20.0
    >>> sw.T90conv(20., t_type='T48')
    19.988162840918179

    References
    ----------
    .. [1] Saunders, P. M., 1991: The International Temperature Scale of 1990,
    ITS-90. WOCE Newsletter, No. 10, WOCE International Project Office,
    Southampton, United Kingdom, 10.

    .. [2] International Temperature Scales of 1948, 1968 and 1990, an ICES
    note, available from http://www.ices.dk/ocean/procedures/its.htm

    Modifications: Filipe Fernandes, 2010
                   2010-11-24. Filipe Fernandes, first version.
                   2010-12-24. Filipe Fernandes, added T48.
    """

    t = np.asanyarray(t)

    if t_type == 'T68':
        #T90 = t * 0.999760057586179 #NOTE: gsw way is less precise
        T90 = t / 1.00024
    elif t_type == 'T48':
        T90 = (t - 4.4e-6 * t * (100 - t) ) / 1.00024
    else:
        raise NameError('Wrong t_type')

    return T90

# Original seawater functions
def adtg(s, t, p):
    r"""
    Calculates adiabatic temperature gradient as per UNESCO 1983 routines.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db]

    Returns
    -------
    adtg : array_like
           adiabatic temperature gradient [:math:`^\circ` C db :sup:`-1`]

    See Also
    --------
    ptmp

    Notes
    -----


    Examples
    --------
    Data from UNESCO 1983 p45

    >>> import seawater.csiro as sw
    >>> t = sw.T90conv([[ 0,  0,  0,  0,  0,  0], [10, 10, 10, 10, 10, 10], [20, 20, 20, 20, 20, 20], [30, 30, 30, 30, 30, 30], [40, 40, 40, 40, 40, 40]])
    >>> s = [[25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35]]
    >>> p = [0, 5000, 10000, 0, 5000, 10000]
    >>> sw.adtg(s, t, p)
    array([[  1.68710000e-05,   1.04700000e-04,   1.69426000e-04,
              3.58030000e-05,   1.17956500e-04,   1.77007000e-04],
           [  1.00194580e-04,   1.60959050e-04,   2.06874170e-04,
              1.14887280e-04,   1.71364200e-04,   2.12991770e-04],
           [  1.73819840e-04,   2.13534000e-04,   2.44483760e-04,
              1.84273240e-04,   2.21087800e-04,   2.49137960e-04],
           [  2.41720460e-04,   2.64764100e-04,   2.82959590e-04,
              2.47934560e-04,   2.69466550e-04,   2.86150390e-04],
           [  3.07870120e-04,   3.16988600e-04,   3.23006480e-04,
              3.09844920e-04,   3.18839700e-04,   3.24733880e-04]])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Bryden, H. 1973. New Polynomials for thermal expansion, adiabatic
    temperature gradient and potential temperature of sea water. Deep-Sea Res.
    Vol20,401-408. doi:10.1016/0011-7471(73)90063-6

    Modifications: 93-04-22. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-16. Filipe Fernandes, Reformulated docstring.
    """

    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)

    T68 = T68conv(t)

    a0 =  3.5803E-5
    a1 =  8.5258E-6
    a2 = -6.836E-8
    a3 =  6.6228E-10

    b0 =  1.8932E-6
    b1 = -4.2393E-8

    c0 =  1.8741E-8
    c1 = -6.7795E-10
    c2 =  8.733E-12
    c3 = -5.4481E-14

    d0 = -1.1351E-10
    d1 =  2.7759E-12

    e0 = -4.6206E-13
    e1 =  1.8676E-14
    e2 = -2.1687E-16

    adtg = ( a0 + ( a1 + ( a2 + a3 * T68 ) * T68) * T68
            + ( b0 + b1 * T68 ) * ( s-35 )
            + ( ( c0 + ( c1 + ( c2 + c3 * T68 ) * T68 ) * T68 )
            + ( d0 + d1 * T68 ) * ( s-35 ) ) * p
            + (  e0 + (e1 + e2 * T68) * T68 )*p*p )

    return adtg

def alpha(s, t, p, pt=False):
    r"""
    Calculate the thermal expansion coefficient.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature or potential temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].
    pt : bool
         True if temperature is potential, default is False

    Returns
    -------
    alpha : array_like
            thermal expansion coeff :math:`\alpha` [:math:`^\circ` C :sup:`-1`]

    See Also
    --------
    beta, aonb

    Notes
    -----


    Examples
    --------
    Data from McDougall 1987

    >>> import seawater.csiro as sw
    >>> s, t, p = 40, 10, 4000
    >>> sw.alpha(s, t, p, pt=True)
    0.00025061316481624323

    References
    ----------
    .. [1] McDougall, Trevor J., 1987: Neutral Surfaces. J. Phys.
    Oceanogr., 17, 1950-1964. doi: 10.1175/1520-0485(1987)017<1950:NS>2.0.CO;2

    Modifications: 93-xx-xx. N.L. Bindoff.
                   93-04-22. Phil Morgan, Help display modified to suit library.
                   93-04-23. Phil Morgan, Input argument checking.
                   94-10-15. Phil Morgan, Pass S,T,P and keyword for 'ptmp'.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-16. Filipe Fernandes, Reformulated docstring.
    """


    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)
    pt = np.asanyarray(pt)

    alpha = aonb(s, t, p, pt) * beta(s, t, p, pt)
    return alpha

def aonb(s, t, p, pt=False):
    r"""
    Calculate :math:`\alpha/\beta`.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature or potential temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].
    pt : bool
         True if temperature is potential, default is False

    Returns
    -------
    aonb : array_like
           :math:`\alpha/\beta` [psu :math:`^\circ` C :sup:`-1`]

    See Also
    --------
    alpha, beta

    Notes
    -----


    TODO: Test pt=False

    Examples
    --------
    Data from McDouogall 1987

    >>> import seawater.csiro as sw
    >>> s, t, p = 40, 10, 4000
    >>> sw.aonb(s, t, p, pt=True)
    0.347650567047807

    References
    ----------
    .. [1] McDougall, Trevor J., 1987: Neutral Surfaces. J. Phys.
    Oceanogr., 17, 1950-1964. doi: 10.1175/1520-0485(1987)017<1950:NS>2.0.CO;2

    Modifications: 93-xx-xx. N.L. Bindoff.
                   93-04-22. Phil Morgan, Help display modified to suit library.
                   93-04-23. Phil Morgan, Input argument checking.
                   94-10-15. Phil Morgan, Pass S,T,P and keyword for 'ptmp'.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-16. Filipe Fernandes, Reformulated docstring.
    """

    # Ensure we use ptmp in calculations
    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)
    pt = np.asanyarray(pt)

    if not pt:
        t = ptmp(s, t, p, 0) # now we have ptmp

    p = np.float32(p)
    t = T68conv(t)

    c1  = np.array([-0.255019e-7, 0.298357e-5, -0.203814e-3,
                    0.170907e-1, 0.665157e-1])
    c2  = np.array([-0.846960e-4, 0.378110e-2])
    c2a = np.array([-0.251520e-11, -0.164759e-6, 0.0])
    c3  = -0.678662e-5
    c4  = np.array([0.791325e-8, -0.933746e-6, 0.380374e-4])
    c5  =  0.512857e-12
    c6  = -0.302285e-13

    # Now calculate the thermal expansion saline contraction ratio aonb
    sm35  = s - 35.0
    aonb  = ( np.polyval(c1, t) + sm35 * ( np.polyval(c2, t)
            + np.polyval(c2a, p) )
            + sm35**2 * c3 + p * np.polyval(c4, t)
            + c5 * (p**2) * (t**2) + c6 * p**3 )

    return aonb

def beta(s, t, p, pt=False):
    r"""
    Calculate the saline contraction coefficient :math:`\beta` as defined by
    T.J. McDougall.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature or potential temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].
    pt : bool
         True if temperature is potential, default is False

    Returns
    -------
    beta : array_like
           saline Contraction Coefficient [psu :sup:`-1`]

    Examples
    --------
    Data from McDouogall 1987

    >>> import seawater.csiro as sw
    >>> s, t, p = 40, 10, 4000
    >>> sw.beta(s, t, p, pt=True)
    0.00072087661741618932

    Notes
    -----


    TODO: Test pt=False for alpha, beta and aonb.

    References
    ----------
    .. [1] McDougall, Trevor J., 1987: Neutral Surfaces. J. Phys.
    Oceanogr., 17, 1950-1964. doi: 10.1175/1520-0485(1987)017<1950:NS>2.0.CO;2

    Modifications: 93-xx-xx. N.L. Bindoff.
                   93-04-22. Phil Morgan, Help display modified to suit library.
                   93-04-23. Phil Morgan, Input argument checking.
                   94-10-15. Phil Morgan, Pass S,T,P and keyword for 'ptmp'.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-16. Filipe Fernandes, Reformulated docstring.
    """

    s, t, p, = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)
    pt = np.asanyarray(pt)

    # Ensure we use ptmp in calculations
    if not pt:
        t = ptmp(s, t, p, 0) # now we have ptmp

    p = np.float32(p)
    t = T68conv(t)

    c1 = np.array([-0.415613e-9, 0.555579e-7, -0.301985e-5, 0.785567e-3])
    c2 = np.array([0.788212e-8, -0.356603e-6])
    c3 = np.array([-0.602281e-15, 0.408195e-10, 0.0])
    c4 = 0.515032e-8
    c5 = np.array([-0.213127e-11, 0.192867e-9, -0.121555e-7])
    c6 = np.array([-0.175379e-14, 0.176621e-12])
    c7 = 0.121551e-17

    # Now calculate the thermal expansion saline contraction ratio adb
    sm35 = s - 35
    beta = ( np.polyval(c1, t) + sm35 * (np.polyval(c2, t) +
            np.polyval(c3, p) ) + c4 * (sm35**2) +
            p * np.polyval(c5, t) + (p**2) * np.polyval(c6, t)
            + c7 * (p**3) )

    return beta

def bfrq(s, t, p, lat=None):
    r"""
    Calculates Brünt-Väisälä Frequency squared (N :sup:`2`) at the mid depths
    from the equation:

    .. math::
        N^{2} = \frac{-g}{\sigma_{\theta}} \frac{d\sigma_{\theta}}{dz}

    Also calculates Potential Vorticity from:

    .. math::
        q=f \frac{N^2}{g}

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature or potential temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].

    lat : number or array_like, optional
          latitude in decimal degrees north [-90..+90].
          Will grav instead of the default g = 9.8 m :sup:`2` s :sup:`-1`) and
          d(z) instead of d(p)

    Returns
    -------
    n2 : array_like
           Brünt-Väisälä Frequency squared (M-1xN)  [rad s :sup:`-2`]
    q : array_like
           planetary potential vorticity (M-1xN)  [ m s :sup:`-1`]
    p_ave : array_like
            mid pressure between P grid (M-1xN) [db]

    See Also
    --------
    pden, dens

    Notes
    -----
    The value of gravity is a global constant if lat is not provided.



    Examples
    --------
    >>> import numpy as np
    >>> import seawater.csiro as sw
    >>> s = [[0, 0, 0], [15, 15, 15], [30, 30, 30],[35,35,35]]
    >>> t = [[15]*3]*4
    >>> p = [[0], [250], [500], [1000]]
    >>> lat = [30,32,35]
    >>> sw.bfrq(s, t, p, lat)[0]
    array([[  4.51543648e-04,   4.51690708e-04,   4.51920753e-04],
           [  4.45598092e-04,   4.45743207e-04,   4.45970207e-04],
           [  7.40996788e-05,   7.41238078e-05,   7.41615525e-05]])

    References
    ----------
    .. [1] A.E. Gill 1982. p.54  eqn 3.7.15 "Atmosphere-Ocean Dynamics"
    Academic Press: New York. ISBN: 0-12-283522-0

    .. [2] Jackett, David R., Trevor J. Mcdougall, 1995: Minimal Adjustment of
    Hydrographic Profiles to Achieve Static Stability. J. Atmos. Oceanic
    Technol., 12, 381-389. doi: 10.1175/1520-0426(1995)012<0381:MAOHPT>2.0.CO;2

    Modifications: 93-06-24. Phil Morgan.
                   Greg Johnson (gjohnson@pmel.noaa.gov) pot. vort. calculation.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   06-04-19. Lindsay Pender, Corrected sign of PV.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-17. Filipe Fernandes, Reformulated docstring.
    """

    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)
    s, t, p = np.broadcast_arrays(s, t, p)

    if (s.ndim != 2)  and (t.ndim != 2):
        raise ValueError('Arguments must be 2D arrays: n_depths, n_profiles')

    if lat is None:
        z = p
        f = np.nan
        g = cte.gdef*np.ones(p.shape)
    else:
        lat = np.asanyarray(lat)
        z = depth(p, lat)
        g = grav(lat, -z) # note that grav expects height as argument
        f = cor(lat)

    m   = p.shape[0]
    iup = np.arange(0, m-1)
    ilo = np.arange(1, m)

    p_ave    = ( p[iup,:] + p[ilo,:] )/2.
    pden_up  = pden( s[iup,:], t[iup,:], p[iup,:], p_ave )
    pden_lo  = pden( s[ilo,:], t[ilo,:], p[ilo,:], p_ave )
    mid_pden = ( pden_up + pden_lo )/2
    dif_pden = pden_up - pden_lo
    mid_g    = ( g[iup,:] + g[ilo,:] )/2
    dif_z    = np.diff(z, axis=0)
    n2       = -mid_g * dif_pden / ( dif_z * mid_pden )
    q        = -f * dif_pden / ( dif_z * mid_pden )

    return n2, q, p_ave

def depth(p, lat):
    r"""
    Calculates depth in meters from pressure in dbars.

    Parameters
    ----------
    p : array_like
        pressure [db].
    lat : number or array_like
          latitude in decimal degrees north [-90..+90].

    Returns
    -------
    z : array_like
        depth [meters]

    Examples
    --------
    UNESCO 1983 data p30

    >>> import seawater.csiro as sw
    >>> lat = [0, 30, 45, 90]
    >>> p   = [[  500,   500,   500,  500], [ 5000,  5000,  5000, 5000], [10000, 10000, 10000, 10000]]
    >>> sw.depth(p, lat)
    array([[  496.65299239,   495.99772917,   495.3427354 ,   494.03357499],
           [ 4915.04099112,  4908.55954332,  4902.08075214,  4889.13132561],
           [ 9725.47087508,  9712.6530721 ,  9699.84050403,  9674.23144056]])

    Notes
    -----
    Original matlab seawater name is dpth and not depth.

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech.
    Pap. in Mar. Sci., No. 44, 53 pp.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    Modifications: 92-04-06. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    p, lat = np.asanyarray(p), np.asanyarray(lat)

    # Eqn 25, p26.  UNESCO 1983.
    c1 =  9.72659
    c2 = -2.2512E-5
    c3 =  2.279E-10
    c4 = -1.82E-15

    gam_dash = 2.184e-6

    lat = abs(lat)
    X   = np.sin( np.deg2rad(lat) )
    X   = X * X

    bot_line = ( 9.780318 * ( 1.0 + ( 5.2788E-3 + 2.36E-5 * X ) * X ) +
               gam_dash * 0.5 * p )
    top_line = ( ( ( c4 * p + c3 ) * p + c2 ) * p + c1 ) * p
    depthm   = top_line / bot_line
    return depthm

def grav(lat, z=0):
    r"""
    Calculates acceleration due to gravity as function of latitude.

    Parameters
    ----------
    lat : array_like
         latitude in decimal degrees north [-90..+90].

    z : number or array_like. Default z = 0
        height in meters (+ve above sea surface, -ve below).

    Returns
    -------
    g : array_like
        gravity [m s :sup:`2`]

    See Also
    --------
    bfrq

    Notes
    -----
    Original matlab name is g and not grav.

    Examples
    --------
    >>> import seawater.csiro as sw
    >>> sw.grav(45, z=0)
    9.8061898752053995

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap.
    in Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] A.E. Gill 1982. p.54  eqn 3.7.15 "Atmosphere-Ocean Dynamics"
    Academic Press: New York. ISBN: 0-12-283522-0

    Modifications: 93-04-20. Phil Morgan.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    lat, z = np.asanyarray(lat), np.asanyarray(z)

    # Eqn p27.  UNESCO 1983.
    lat     = np.abs(lat)
    X       = np.sin( np.deg2rad(lat) )
    sin2    = X * X
    grav    = 9.780318 * ( 1.0 + ( 5.2788E-3 + 2.36E-5 * sin2 ) * sin2 )
    grav    = grav / ( ( 1 + z / cte.a )**2 )    # from A.E.Gill p.597
    return  grav

def cor(lat):
    r"""
    Calculates the Coriolis factor :math:`f` defined by:

    .. math::
        f = 2 \Omega \sin(lat)

    where:

    .. math::
        \Omega = \frac{2 \pi}{\textrm{sidereal day}} = 7.2921150e^{-5} \textrm{ radians sec}^{-1}


    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90].

    Returns
    -------
    f : array_like
        Coriolis factor [s :sup:`-1`]

    See Also
    --------
    inertial_period

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.csiro as sw
    >>> sw.cor(45)
    0.00010312607931384281

    References
    ----------
    .. [1] S. Pond & G.Pickard 2nd Edition 1986 Introductory Dynamical
    Oceanography Pergamon Press Sydney. ISBN 0-08-028728-X

    .. [2] A.E. Gill 1982. p.54  eqn 3.7.15 "Atmosphere-Ocean Dynamics"
    Academic Press: New York. ISBN: 0-12-283522-0

    .. [3] Groten, E., 2004: Fundamental Parameters and Current (2004) Best
    Estimates of the Parameters of Common Relevance to Astronomy, Geodesy,
    and Geodynamics. Journal of Geodesy, 77, pp. 724-797.

    Modifications: 93-04-20. Phil Morgan.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    lat = np.asanyarray(lat)

    # Eqn p27.  UNESCO 1983.
    f = 2 * cte.OMEGA * np.sin( np.deg2rad(lat) )
    return f

def cndr(s, t, p):
    r"""
    Calculates conductivity ratio.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db]. The shape can be "broadcasted"

    Returns
    -------
    cndr : array_like
           conductivity ratio. R = C(s,t,p)/C(35,15(IPTS-68),0) [no units]

    See Also
    --------
    salds, sals, salrt

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested.

    Examples
    --------
    Data from UNESCO 1983 p9

    >>> import seawater.csiro as sw
    >>> t = sw.T90conv([0, 10, 0, 10, 10, 30])
    >>> p    = [0, 0, 1000, 1000, 0, 0]
    >>> s    = [25, 25, 25, 25, 40, 40]
    >>> sw.cndr(s, t, p)
    array([ 0.49800825,  0.65499015,  0.50624434,  0.66297496,  1.00007311,
            1.52996697])

    Notes
    -----
    TODO: Pressure broadcast feature need to be tested.

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap.
    in Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    Modifications: 93-04-21. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)

    T68 = T68conv(t)

    # DO A NEWTON-RAPHSON ITERATION FOR INVERSE INTERPOLATION OF Rt FROM S.
    Rx = np.sqrt( s / 35.0 ) # first guess at Rx = sqrt(Rt)
    SInc = sals( Rx**2, t ) # S Increment (guess) from Rx
    iloop = 0
    while True:
        Rx = Rx + ( s - SInc ) / salds( Rx, t - 15 )
        SInc = sals( Rx**2, t )
        iloop = iloop + 1
        dels = np.abs( SInc - s )
        if not (dels.all() > 1.0e-10) & (iloop < 100):
            break

    # ONCE Rt FOUND, CORRESPONDING TO EACH (S,T) EVALUATE R
    # eqn(4) p.8 UNESCO 1983
    d1 =  3.426e-2
    d2 =  4.464e-4
    d3 =  4.215e-1
    d4 = -3.107e-3

    e1 =  2.070e-5
    e2 = -6.370e-10
    e3 =  3.989e-15

    A  = ( d3 + d4 * T68 )
    B  = 1 + d1 * T68 + d2 * T68**2
    C  = p * ( e1 + e2 * p + e3 * p**2 )

    # eqn(6) p.9 UNESCO 1983.
    Rt    = Rx**2
    rt    = salrt(t)
    #Rtrt  = rt * Rt # NOTE: unused in the code, but present in the original
    D     = B - A * rt * Rt
    E     = rt * Rt * A * ( B + C )
    r     = np.sqrt( np.abs( D**2 + 4 * E ) ) - D
    r     = 0.5 * r/A
    return r


def sals(rt, t):
    r"""
    Salinity of sea water as a function of Rt and T.
    UNESCO 1983 polynomial.

    Parameters
    ----------
    rt : array_like
         :math:`rt(s,t) = \frac{C(s,t,0)}{C(35, t(\textrm{IPTS-68}), 0)}`
    t : array_like
        temperature [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    s : array_like
        salinity [psu (PSS-78)]

    See Also
    --------
    salt

    Notes
    -----
    TODO

    Examples
    --------
    Data from UNESCO 1983 p9

    >>> import seawater.csiro as sw
    >>> t = T90conv([15, 20, 5])
    >>> rt   = [  1, 1.0568875, 0.81705885]
    >>> sw.sals(rt, t)
    array([ 35.        ,  37.24562718,  27.99534701])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    Modifications: 93-04-17. Phil Morgan.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    rt, t = np.asanyarray(rt), np.asanyarray(t)

    # eqn (1) & (2) p6,7 unesco
    del_T68 = T68conv(t) - 15

    a0 =  0.0080
    a1 = -0.1692
    a2 = 25.3851
    a3 = 14.0941
    a4 = -7.0261
    a5 =  2.7081

    b0 =  0.0005
    b1 = -0.0056
    b2 = -0.0066
    b3 = -0.0375
    b4 =  0.0636
    b5 = -0.0144

    k  =  0.0162

    Rtx = (rt)**0.5
    del_S = ( ( del_T68 / ( 1 + k * del_T68 ) ) *
            ( b0 + ( b1 + ( b2 + ( b3 + ( b4 + b5 * Rtx )
            * Rtx ) * Rtx ) * Rtx ) * Rtx) )

    s = a0 + ( a1 + ( a2 + ( a3 + ( a4 + a5 * Rtx) *
               Rtx) * Rtx ) * Rtx ) * Rtx

    s = s + del_S

    return s

def salds(rtx, delt):
    r"""
    Calculates Salinity differential (:math:`\frac{dS}{d(\sqrt{Rt})}`) at
    constant temperature.

    Parameters
    ----------
    rtx : array_like
          :math:`\sqrt{rt}`
    delt : array_like
           t-15 [:math:`^\circ` C (IPTS-68)]

    Returns
    -------
    ds : array_like
         :math:`\frac{dS}{d rtx}`

    See Also
    --------
    cndr, salt

    Notes
    -----
    TODO

    Examples
    --------
    Data from UNESCO 1983 p9

    >>> import numpy as np
    >>> import seawater.csiro as sw
    >>> delt = T90conv([15, 20, 5])  - 15
    >>> rtx  = np.array([  1, 1.0568875, 0.81705885])**0.5
    >>> sw.salds(rtx, delt)
    array([ 78.31921607,  81.5689307 ,  68.19023687])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    Modifications: 93-04-21. Phil Morgan.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    rtx, delt = np.asanyarray(rtx), np.asanyarray(delt)

    #a0 =  0.0080 #TODO: unused in the code, but present in the original
    a1 = -0.1692
    a2 = 25.3851
    a3 = 14.0941
    a4 = -7.0261
    a5 =  2.7081

    #b0 =  0.0005 #TODO: unused in the code, but present in the original
    b1 = -0.0056
    b2 = -0.0066
    b3 = -0.0375
    b4 =  0.0636
    b5 = -0.0144

    k  =  0.0162

    ds = ( a1 + ( 2 * a2 + ( 3 * a3 + ( 4 * a4 + 5 * a5 * rtx ) * rtx )
          * rtx ) * rtx + ( delt / ( 1 + k * delt ) ) *
          ( b1 + ( 2 * b2 + ( 3 * b3 + ( 4 * b4 + 5 * b5 * rtx )
          * rtx ) * rtx) * rtx) )

    return ds

def salrt(t):
    r"""
    Equation for rt used in calculating salinity. UNESCO 1983 polynomial.

    .. math::
        rt(t) = \frac{C(35,t,0)}{C(35,15(\textrm{IPTS-68}), 0)}


    Parameters
    ----------
      t : array_like
          temperature [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    rt : array_like
    conductivity ratio  [no units]

    See Also
    --------
    salt

    Notes
    -----
    TODO

    Examples
    --------
    Data from UNESCO 1983 p9

    >>> import seawater.csiro as sw
    >>> t = T90conv([15, 20, 5])
    >>> sw.salrt(t)
    array([ 1.        ,  1.11649272,  0.77956585])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    Modifications: 93-04-17. Phil Morgan.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    t = np.asanyarray(t)

    #Eqn (3) p.7 UNESCO.
    T68 = T68conv(t)

    c0 =  0.6766097
    c1 =  2.00564e-2
    c2 =  1.104259e-4
    c3 = -6.9698e-7
    c4 =  1.0031e-9

    rt = c0 + ( c1 + ( c2 + ( c3 + c4 * T68) * T68) * T68 ) * T68
    return rt

def salt(r, t, p):
    r"""
    Calculates Salinity from conductivity ratio. UNESCO 1983 polynomial.

    Parameters
    ----------
    r : array_like
        conductivity ratio :math:`R = \frac{C(S,T,P)}{C(35,15(IPTS-68),0)}`
    t : array_like
        temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db]

    Returns
    -------
    s : array_like
        salinity [psu (PSS-78)]

    See Also
    --------
    sals, salrt, salrp

    Notes
    -----
    TODO

    Examples
    --------
    Data from UNESCO 1983 p9

    >>> import seawater.csiro as sw
    >>> r = [1, 1.2, 0.65]
    >>> t = sw.T90conv([15, 20, 5])
    >>> p = [0, 2000, 1500]
    >>> sw.salt(r, t, p)
    array([ 34.99999992,  37.24562765,  27.99534693])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    Modifications: 93-04-17. Phil Morgan.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    r, t, p = np.asanyarray(r), np.asanyarray(t), np.asanyarray(p)

    rt = salrt(t)
    rp = salrp(r, t, p )
    rt = r / ( rp * rt )
    s  = sals(rt, t)

    return s

def salrp(r, t, p):
    r"""
    Equation for Rp used in calculating salinity. UNESCO 1983 polynomial.

    .. math::
        Rp(S,T,P) = \frac{C(S,T,P)}{C(S,T,0)}


    Parameters
    ----------
    r : array_like
        conductivity ratio :math:`R = \frac{C(S,T,P)}{C(35,15(IPTS-68),0)}`
    t : array_like
        temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db]

    Returns
    -------
    rp : array_like
        conductivity ratio :math:`Rp(S,T,P) = \frac{C(S,T,P)}{C(S,T,0)}`

    See Also
    --------
    salt

    Notes
    -----
    TODO

    Examples
    --------

    >>> import seawater.csiro as sw
    >>> r = [1, 1.2, 0.65]
    >>> t = T90conv([15, 20, 5])
    >>> p = [0, 2000, 1500]
    >>> sw.salrp(r, t, p)
    array([ 1.        ,  1.01694294,  1.02048638])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    Modifications: 93-04-17. Phil Morgan.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    r, t, p = np.asanyarray(r), np.asanyarray(t), np.asanyarray(p)

    # eqn (4) p.8 unesco.
    T68 = T68conv(t)

    d1 =  3.426e-2
    d2 =  4.464e-4
    d3 =  4.215e-1
    d4 = -3.107e-3

    e1 =  2.070e-5
    e2 = -6.370e-10
    e3 =  3.989e-15

    rp = ( 1 + ( p * ( e1 + e2 * p + e3 * p**2 ) )
         / ( 1 + d1 * T68 + d2 * T68**2 + ( d3 + d4 * T68 ) * r) )

    return rp

def fp(s, p):
    r"""
    Freezing point of Sea Water using UNESCO 1983 polynomial.

    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    p : array_like
        pressure [db]

    Returns
    -------
    fp : array_like
        freezing point temperature [:math:`^\circ` C (ITS-90)]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    UNESCO DATA p.30

    >>> import seawater.csiro as sw
    >>> s = [[5, 10, 15, 20, 25, 30, 35, 40], [5, 10, 15, 20, 25, 30, 35, 40]]
    >>> p = [[ 0, 0, 0, 0, 0, 0, 0, 0], [500, 500, 500, 500, 500, 500, 500, 500]]
    >>> sw.fp(s, p)
    array([[-0.27369757, -0.54232831, -0.81142026, -1.0829461 , -1.35804594,
            -1.63748903, -1.9218401 , -2.2115367 ],
           [-0.65010724, -0.91873798, -1.18782992, -1.45935577, -1.73445561,
            -2.01389869, -2.29824976, -2.58794636]])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    Modifications: 93-04-20. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    s, p = np.asanyarray(s), np.asanyarray(p)


    #NOTE: P = P/10 # to convert db to Bar as used in UNESCO routines

    # eqn  p.29
    a0 = -0.0575
    a1 =  1.710523e-3
    a2 = -2.154996e-4
    b  = -7.53e-4

    fp = T90conv( a0 * s + a1 * s * (s)**0.5 + a2 * s**2 + b * p )

    return fp

def svel(s, t, p):
    r"""
    Sound Velocity in sea water using UNESCO 1983 polynomial.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    svel : array_like
           sound velocity  [m/s]


    See Also
    --------
    TODO

    Notes
    -----


    TODO: Add equation to docstring.

    Examples
    --------
    Data from Pond and Pickard Intro. Dynamical Oceanography 2nd ed. 1986

    >>> import seawater.csiro as sw
    >>> t = T90conv([[  0,  0,  0,  0,  0,  0], [ 10, 10, 10, 10, 10, 10], [ 20, 20, 20, 20, 20, 20], [ 30, 30, 30, 30, 30, 30], [ 40, 40, 40, 40, 40, 40]])
    >>> s = [[ 25, 25, 25, 35, 35, 35], [ 25, 25, 25, 35, 35, 35], [ 25, 25, 25, 35, 35, 35], [ 25, 25, 25, 35, 35, 35], [ 25, 25, 25, 35, 35, 35]]
    >>> p = [ 0, 5000, 10000, 0, 5000, 10000]
    >>> sw.svel(s, t, p)
    array([[ 1435.789875  ,  1520.358725  ,  1610.4074    ,  1449.13882813,
             1533.96863705,  1623.15007097],
           [ 1477.68316464,  1561.30635914,  1647.39267114,  1489.82233602,
             1573.40946928,  1658.99115504],
           [ 1510.31388348,  1593.59671798,  1676.80967748,  1521.4619731 ,
             1604.4762822 ,  1687.18305631],
           [ 1535.21434752,  1618.95631952,  1700.60547902,  1545.59485539,
             1628.97322783,  1710.06294277],
           [ 1553.44506636,  1638.02522336,  1719.15088536,  1563.20925247,
             1647.29949576,  1727.83176404]])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    Modifications: 93-04-20. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)

    # UNESCO 1983. eqn.33  p.46
    p = p/10  # convert db to bars as used in UNESCO routines
    T68 = T68conv(t)

    # eqn 34 p.46
    c00 = 1402.388
    c01 =    5.03711
    c02 =   -5.80852e-2
    c03 =    3.3420e-4
    c04 =   -1.47800e-6
    c05 =    3.1464e-9

    c10 =  0.153563
    c11 =  6.8982e-4
    c12 = -8.1788e-6
    c13 =  1.3621e-7
    c14 = -6.1185e-10

    c20 =  3.1260e-5
    c21 = -1.7107e-6
    c22 =  2.5974e-8
    c23 = -2.5335e-10
    c24 =  1.0405e-12

    c30 = -9.7729e-9
    c31 =  3.8504e-10
    c32 = -2.3643e-12

    Cw  =( ((((c32 * T68 + c31) * T68 + c30) * p +
          ((((c24 * T68 + c23) * T68 + c22) * T68 + c21) * T68 + c20)) * p +
          ((((c14 * T68 + c13) * T68 + c12) * T68 + c11) * T68 + c10)) * p +
          ((((c05 * T68 + c04) * T68 + c03) * T68 + c02) * T68 + c01)*T68+c00 )

    # eqn 35. p.47
    a00 =  1.389
    a01 = -1.262e-2
    a02 =  7.164e-5
    a03 =  2.006e-6
    a04 = -3.21e-8

    a10 =  9.4742e-5
    a11 = -1.2580e-5
    a12 = -6.4885e-8
    a13 =  1.0507e-8
    a14 = -2.0122e-10

    a20 = -3.9064e-7
    a21 =  9.1041e-9
    a22 = -1.6002e-10
    a23 =  7.988e-12

    a30 =  1.100e-10
    a31 =  6.649e-12
    a32 = -3.389e-13

    A = ( ((((a32 * T68 + a31) * T68 + a30) * p +
        (((a23 * T68 + a22) * T68 + a21) * T68 + a20)) * p +
        ((((a14 * T68 + a13) * T68 + a12) * T68 + a11) * T68 + a10)) * p +
        (((a04 * T68 + a03) * T68 + a02) * T68 + a01) * T68 + a00 )

    # eqn 36 p.47
    b00 = -1.922e-2
    b01 = -4.42e-5
    b10 =  7.3637e-5
    b11 =  1.7945e-7

    B = b00 + b01 * T68 + ( b10 + b11 * T68) * p

    # eqn 37 p.47
    d00 =  1.727e-3
    d10 = -7.9836e-6

    D = d00 + d10 * p

    # eqn 33 p.46
    svel = Cw + A * s + B * s * (s)**0.5 + D * s**2

    return svel

def pres(depth, lat):
    r"""
    Calculates pressure in dbars from depth in meters.

    Parameters
    ----------
    depth : array_like
            depth [meters]
    lat : array_like
          latitude in decimal degrees north [-90..+90]


    Returns
    -------
    p : array_like
           pressure [db]

    See Also
    --------
    pressure from depth (TODO)


    Notes
    -----


    Examples
    --------
    >>> import seawater.csiro as sw
    >>> depth, lat = 7321.45, 30
    >>> sw.pres(depth,lat)
    7500.0065130118019

    References
    ----------
    .. [1] Saunders, Peter M., 1981: Practical Conversion of Pressure to Depth.
    J. Phys. Oceanogr., 11, 573-574.
    doi: 10.1175/1520-0485(1981)011<0573:PCOPTD>2.0.CO;2

    Modifications: 93-06-25. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    depth, lat = np.asanyarray(depth), np.asanyarray(lat)

    X = np.sin( abs( np.deg2rad(lat) ) )
    C1 = 5.92E-3 + X**2 * 5.25E-3
    pres = ( (1-C1) - ( ( (1-C1)**2 ) - ( 8.84E-6 * depth ) )**0.5 ) / 4.42E-6
    return pres

def dist(lon, lat, units='km'):
    r"""
    Calculate distance between two positions on globe using the "Plane
    Sailing" method. Also uses simple geometry to calculate the bearing of
    the path between position pairs.

    Parameters
    ----------
    lon : array_like
          decimal degrees (+ve E, -ve W) [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [- 90.. +90]
    units : string, optional
            default kilometers

    Returns
    -------
    dist : array_like
           distance between positions in units
    phaseangle : array_like
                 angle of line between stations with x axis (East).
                 Range of values are -180..+180. (E=0, N=90, S=-90)

    See Also
    --------
    TODO

    Notes
    -----
    Usually used to create a distance vector to plot hydrographic data. However,
    pay attention to the phaseangle to avoid apples and oranges!

    Also not that the input order for the matlab version is lat,lon
    (alphabetic order), while this version is lon,lat (geometric order).

    Examples
    --------
    >>> import seawater.csiro as sw
    >>> lon = [35, 35]
    >>> lat = [41, 40]
    >>> sw.dist(lon, lat)
    (array([ 111.12]), array([-90.]))

    Create a distance vector

    >>> lon = np.arange(30,40,1)
    >>> lat = 35
    >>> np.cumsum(np.append(0, sw.dist(lon, lat, units='km')[0]))
    array([   0.        ,   91.02417516,  182.04835032,  273.07252548,
            364.09670065,  455.12087581,  546.14505097,  637.16922613,
            728.19340129,  819.21757645])

    References
    ----------
    .. [1] The PLANE SAILING method as described in "CELESTIAL NAVIGATION" 1989
    by Dr. P. Gormley. The Australian Antarctic Division.

    Modifications: 92-02-10. Phil Morgan and Steve Rintoul.
                   99-06-25. Lindsay Pender, name change distance to sw_dist.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    lon, lat = np.asanyarray(lon), np.asanyarray(lat)

    if lat.size == 1:
        lat = np.repeat(lat, lon.size)
    elif lon.size == 1:
        lon = np.repeat(lon, lat.size)

    npositions = max(lon.shape)

    ind = np.arange( 0, npositions-1, 1) # index to first of position pairs

    dlon = np.diff(lon, axis=0)
    if any( abs(dlon) > 180 ):
        flag = abs(dlon) > 180
        dlon[flag] = -np.sign( dlon[flag] ) * ( 360 - abs( dlon[flag] ) )

    latrad = abs( np.deg2rad(lat) )
    dep    = np.cos( ( latrad [ind+1] + latrad[ind] ) / 2 ) * dlon
    dlat   = np.diff( lat, axis=0 )
    dist   = cte.DEG2NM * ( dlat**2 + dep**2 )**0.5

    if units == 'km':
        dist = dist * cte.NM2KM

    # Calcualte angle to x axis
    phaseangle  = np.rad2deg( np.angle( dep + dlat * 1j ) )
    return dist, phaseangle

def satAr(s, t):
    r"""
    Solubility (saturation) of Argon (Ar) in sea water.

    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    satAr : array_like
            solubility of Ar [ml l :sup:`-1`]

    See Also
    --------
    satO2, satN2

    Notes
    -----
    TODO

    Examples
    --------
    Data from Weiss 1970

    >>> import seawater.csiro as sw
    >>> t = T90conv([[ -1, -1], [ 10, 10], [ 20, 20], [ 40, 40]])
    >>> s = [[ 20, 40], [ 20, 40], [ 20, 40], [ 20, 40]]
    >>> sw.satAr(s, t)
    array([[ 0.4455784 ,  0.38766011],
           [ 0.33970659,  0.29887756],
           [ 0.27660227,  0.24566428],
           [ 0.19861429,  0.17937698]])

    References
    ----------
    .. [1] Weiss, R. F. 1970. The Solubility of Nitrogen, Oxygen and Argon in
    Water and Seawater Deep-Sea Research Vol. 17, p. 721-735.
    doi:10.1016/0011-7471(70)90037-9

    Modifications: 97-11-05. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    s, t = np.asanyarray(s), np.asanyarray(t)

    # convert T to Kelvin
    t = cte.Kelvin + T68conv(t)

    # constants for Eqn (4) of Weiss 1970
    a1 = -173.5146
    a2 =  245.4510
    a3 =  141.8222
    a4 =  -21.8020
    b1 =   -0.034474
    b2 =    0.014934
    b3 =   -0.0017729

    # Eqn (4) of Weiss 1970
    lnC = ( a1 + a2 * ( 100/t ) + a3 * np.log( t/100 ) + a4 * ( t/100 ) +
          s * ( b1 + b2 * ( t/100 ) + b3 * ( ( t/100 )**2) ) )

    c = np.exp(lnC)

    return c

def satN2(s, t):
    r"""
    Solubility (saturation) of Nitrogen (N2) in sea water.

    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    satN2 : array_like
            solubility of N2  [ml l :sup:`-1`]

    See Also
    --------
    satO2, satAr

    Notes
    -----
    TODO

    Examples
    --------
    Data from Weiss 1970

    >>> import seawater.csiro as sw
    >>> t = T90conv([[ -1, -1], [ 10, 10], [ 20, 20], [ 40, 40]])
    >>> s = [[ 20, 40], [ 20, 40], [ 20, 40], [ 20, 40]]
    >>> sw.satN2(s, t)
    array([[ 16.27952432,  14.00784526],
           [ 12.64036196,  11.01277257],
           [ 10.46892822,   9.21126859],
           [  7.78163876,   6.95395099]])


    References
    ----------
    .. [1] Weiss, R. F. 1970. The Solubility of Nitrogen, Oxygen and Argon in
    Water and Seawater Deep-Sea Research Vol. 17, p. 721-735.
    doi:10.1016/0011-7471(70)90037-9

    Modifications: 97-11-05. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    s, t = np.asanyarray(s), np.asanyarray(t)

    # convert T to Kelvin
    t = cte.Kelvin + T68conv(t)

    # constants for Eqn (4) of Weiss 1970
    a1 = -172.4965
    a2 =  248.4262
    a3 =  143.0738
    a4 =  -21.7120
    b1 =   -0.049781
    b2 =    0.025018
    b3 =   -0.0034861

    # Eqn (4) of Weiss 1970
    lnC = ( a1 + a2 * ( 100/t ) + a3 * np.log( t/100 ) + a4 * ( t/100 ) +
            s * ( b1 + b2 * ( t/100 ) + b3 * ( ( t/100 )**2 ) ) )

    c = np.exp(lnC)
    return c

def satO2(s, t):
    r"""
    Solubility (saturation) of Oxygen (O2) in sea water.

    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [:math:`^\circ` C (ITS-68)]

    Returns
    -------
    satO2 : array_like
            solubility of O2  [ml l :sup:`-1` ]

    Notes
    -----
    TODO

    Examples
    --------
    Data from Weiss 1970

    >>> import seawater.csiro as sw
    >>> t = T90conv([[ -1, -1], [ 10, 10], [ 20, 20], [ 40, 40]])
    >>> s = [[ 20, 40], [ 20, 40], [ 20, 40], [ 20, 40]]
    >>> sw.satO2(s, t)
    array([[ 9.162056  ,  7.98404249],
           [ 6.95007741,  6.12101928],
           [ 5.64401453,  5.01531004],
           [ 4.0495115 ,  3.65575811]])

    References
    ----------
    .. [1] Weiss, R. F. 1970. The Solubility of Nitrogen, Oxygen and Argon in
    Water and Seawater Deep-Sea Research Vol. 17, p. 721-735.
    doi:10.1016/0011-7471(70)90037-9

    Modifications: 97-11-05. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    s, t = np.asanyarray(s), np.asanyarray(t)

    # convert T to Kelvin
    t = cte.Kelvin + T68conv(t)

    # constants for Eqn (4) of Weiss 1970
    a1 = -173.4292
    a2 =  249.6339
    a3 =  143.3483
    a4 =  -21.8492
    b1 =   -0.033096
    b2 =    0.014259
    b3 =   -0.0017000

    # Eqn (4) of Weiss 1970
    lnC = ( a1 + a2 * ( 100/t ) + a3 * np.log( t/100 ) + a4 * ( t/100 ) +
            s * ( b1 + b2 * ( t/100 ) + b3 * ( ( t/100 )**2 ) ) )

    c = np.exp(lnC)
    return c

def dens0(s, t):
    r"""
    Density of Sea Water at atmospheric pressure.

    Parameters
    ----------
    s(p=0) : array_like
             salinity [psu (PSS-78)]
    t(p=0) : array_like
             temperature [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    dens0(s, t) : array_like
                  density  [kg m :sup:`3`] of salt water with properties
                  (s, t, p=0) 0 db gauge pressure

    See Also
    --------
    svan, dens, smow, seck, pden

    Notes
    -----
    TODO

    Examples
    --------
    Data from UNESCO Tech. Paper in Marine Sci. No. 44, p22

    >>> import seawater.csiro as sw
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> sw.dens0(s, t)
    array([  999.842594  ,   999.842594  ,   995.65113374,   995.65113374,
            1028.10633141,  1028.10633141,  1021.72863949,  1021.72863949])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere equation
    of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
    doi:10.1016/0198-0149(81)90122-9

    Modifications: 92-11-05. Phil Morgan.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    s, t = np.asanyarray(s), np.asanyarray(t)

    T68 = T68conv(t)

    # UNESCO 1983 eqn(13) p17
    b0 =  8.24493e-1
    b1 = -4.0899e-3
    b2 =  7.6438e-5
    b3 = -8.2467e-7
    b4 =  5.3875e-9

    c0 = -5.72466e-3
    c1 =  1.0227e-4
    c2 = -1.6546e-6

    d0 = 4.8314e-4
    dens0 =( smow(t) + ( b0 + ( b1 + ( b2 + ( b3 + b4 * T68 ) * T68 ) * T68 )
    * T68 ) * s + ( c0 + ( c1 + c2 * T68 ) * T68 ) * s * (s)**0.5 + d0 * s**2 )
    return dens0

def smow(t):
    r"""
    Density of Standard Mean Ocean Water (Pure Water) using EOS 1980.

    Parameters
    ----------
    t : array_like
        temperature [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    dens(t) : array_like
              density  [kg m :sup:`3`]

    See Also
    --------
    svan, dens, dens0, seck, pden

    Notes
    -----
    Standard Mean Ocean Water (SMOW) is the water collected in the deep ocean
    used as a reference.

    Examples
    --------
    Data from UNESCO Tech. Paper in Marine Sci. No. 44, p22

    >>> import seawater.csiro as sw
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> sw.smow(t)
    array([ 999.842594  ,  999.842594  ,  995.65113374,  995.65113374,
            999.842594  ,  999.842594  ,  995.65113374,  995.65113374])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere equation
    of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
    doi:10.1016/0198-0149(81)90122-9

    Modifications: 92-11-05. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    t = np.asanyarray(t)

    a0 = 999.842594
    a1 =   6.793952e-2
    a2 =  -9.095290e-3
    a3 =   1.001685e-4
    a4 =  -1.120083e-6
    a5 =   6.536332e-9

    T68  = T68conv(t)
    dens = ( a0 + ( a1 + ( a2 + ( a3 + ( a4 + a5 * T68 ) * T68 ) * T68 )
             * T68 ) * T68 )
    return dens

def seck(s, t, p=0):
    r"""
    Secant Bulk Modulus (K) of Sea Water using Equation of state 1980. UNESCO
    polynomial implementation.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    k : array_like
        secant bulk modulus [bars]

    See Also
    --------
    dens

    Notes
    -----
    TODO

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22

    >>> import seawater.csiro as sw
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> p = [0, 10000, 0, 10000, 0, 10000, 0, 10000]
    >>> sw.seck(s, t, p)
    array([ 19652.21      ,  22977.2115    ,  22336.0044572 ,  25656.8196222 ,
            21582.27006823,  24991.99729129,  23924.21823158,  27318.32472464])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere equation
    of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
    doi:10.1016/0198-0149(81)90122-9

    Modifications: 92-11-05. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)

    # Compute compression terms
    p   = p/10.0 # convert from db to atmospheric pressure units
    T68 = T68conv(t)

    # Pure water terms of the secant bulk modulus at atmos pressure.
    # UNESCO eqn 19 p 18
    h3 = -5.77905E-7
    h2 =  1.16092E-4
    h1 =  1.43713E-3
    h0 =  3.239908   #[-0.1194975]

    AW = h0 + ( h1 + ( h2 + h3 * T68 ) * T68 ) * T68

    k2 =  5.2787E-8
    k1 = -6.12293E-6
    k0 =  8.50935E-5  #[+3.47718E-5]

    BW = k0 + ( k1 + k2 * T68 ) * T68

    e4 =    -5.155288E-5
    e3 =     1.360477E-2
    e2 =    -2.327105
    e1 =   148.4206
    e0 = 19652.21    #[-1930.06]

    KW  = e0 + ( e1 + ( e2 + ( e3 + e4 * T68 ) * T68 ) * T68 ) * T68 # eqn 19

    # Sea water terms of secant bulk modulus at atmos. pressure
    j0 = 1.91075E-4

    i2 = -1.6078E-6
    i1 = -1.0981E-5
    i0 =  2.2838E-3

    SR = (s)**0.5

    A  = AW + ( i0 + ( i1 + i2 * T68 ) * T68 + j0 * SR ) * s

    m2 =  9.1697E-10
    m1 =  2.0816E-8
    m0 = -9.9348E-7

    B = BW + ( m0 + ( m1 + m2 * T68 ) * T68 ) * s # eqn 18

    f3 =  -6.1670E-5
    f2 =   1.09987E-2
    f1 =  -0.603459
    f0 =  54.6746

    g2 = -5.3009E-4
    g1 =  1.6483E-2
    g0 =  7.944E-2

    K0 = ( KW + ( f0 + ( f1 + ( f2 + f3 * T68 ) * T68 ) * T68
            + ( g0 + ( g1 + g2 * T68 ) * T68 ) * SR ) * s ) # eqn 16

    K = K0 + ( A + B * p ) * p # eqn 15
    return K

def dens(s, t, p):
    r"""
    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    dens : array_like
           density  [kg m :sup:`3`]

    See Also
    --------
    dens0, seck

    Notes
    -----
    TODO

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22

    >>> import seawater.csiro as sw
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> p = [0, 10000, 0, 10000, 0, 10000, 0, 10000]
    >>> sw.dens(s, t, p)
    array([  999.842594  ,  1045.33710972,   995.65113374,  1036.03148891,
            1028.10633141,  1070.95838408,  1021.72863949,  1060.55058771])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J., Chen, C.T., Bradshaw, A., and Schleicher, K. A new
    high pressure equation of state for seawater. Deap-Sea Research., 1980,
    Vol27A, pp255-264. doi:10.1016/0198-0149(80)90016-3

    Modifications: 92-11-05. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)

    # UNESCO 1983. eqn.7  p.15
    densP0 = dens0(s, t)
    K      = seck(s, t, p)
    p      = p / 10.0 # convert from db to atm pressure units
    dens   = densP0 / ( 1-p / K )
    return dens

def pden(s, t, p, pr=0):
    r"""
    Calculates potential density of water mass relative to the specified
    reference pressure by pden = dens(S, ptmp, PR).

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].
    pr : number
         reference pressure [db], default = 0

    Returns
    -------
    pden : array_like
           potential density relative to the ref. pressure [kg m :sup:3]

    See Also
    --------
    ptmp, dens

    Notes
    -----
    The reference pressure is in "oceanographic" standards, so 0 db means at
    surface or 1 atm.

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22

    >>> import seawater.csiro as sw
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> p = [0, 10000, 0, 10000, 0, 10000, 0, 10000]
    >>> sw.pden(s, t, p)
    array([  999.842594  ,   999.79523994,   995.65113374,   996.36115932,
            1028.10633141,  1028.15738545,  1021.72863949,  1022.59634627])

    :math:`\sigma_{4}` (at 4000 db)

    >>> sw.pden(s, t, p, 4000) - 1000
    array([ 19.2895493 ,  19.33422519,  12.43271053,  13.27563816,
            46.30976432,  46.48818851,  37.76150878,  38.74500757])

    References
    ----------
    .. [1] A.E. Gill 1982. p.54  eqn 3.7.15 "Atmosphere-Ocean Dynamics"
    Academic Press: New York. ISBN: 0-12-283522-0

    Modifications: 92-04-06. Phil Morgan.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)
    pr = np.asanyarray(pr)

    pt   = ptmp(s, t, p, pr)
    pden = dens(s, pt, pr)
    return pden

def svan(s, t, p=0):
    r"""
    Specific Volume Anomaly calculated as svan = 1/dens(s,t,p) - 1/dens(35,0,p).

    Note that it is often quoted in literature as 1e8*units.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    svan : array_like
           specific volume anomaly  [m :sup:`3` kg :sup:`-1`]

    See Also
    --------
    dens

    Notes
    -----
    TODO

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22

    >>> import seawater.csiro as sw
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> p = ([0, 10000, 0, 10000, 0, 10000, 0, 10000])
    >>> sw.svan(s, t, p)
    array([  2.74953924e-05,   2.28860986e-05,   3.17058231e-05,
             3.14785290e-05,   0.00000000e+00,   0.00000000e+00,
             6.07141523e-06,   9.16336113e-06])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] S. Pond & G.Pickard 2nd Edition 1986 Introductory Dynamical
    Oceanography Pergamon Press Sydney. ISBN 0-08-028728-X

    Modifications: 92-11-05. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)

    svan = 1/dens( s, t, p ) - 1/dens( 35, 0.0, p )
    return svan

def gpan(s, t, p):
    r"""
    Geopotential Anomaly calculated as the integral of svan from the
    the sea surface to the bottom. THUS RELATIVE TO SEA SURFACE.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    gpan : array_like
           geopotential anomaly
           [m :sup:`3` kg :sup:`-1` Pa = m :sup:`2` s :sup:`-2` = J kg :sup:`-1`]

    See Also
    --------
    svan, gvel

    Notes
    -----
    Adapted method from Pond and Pickard (p76) to calculate gpan relative to
    sea surface whereas P&P calculated relative to the deepest common depth.
    Note that older literature may use units of "dynamic decimeter" for above.

    TODO: example with values that make some sense

    TODO: pass axis as argument

    Examples
    --------
    Data from Unesco Tech. Paper in Marine Sci. No. 44, p22

    >>> import numpy as np
    >>> import seawater.csiro as sw
    >>> s = [[0, 0, 0], [15, 15, 15], [30, 30, 30],[35,35,35]]
    >>> t = [[15]*3]*4
    >>> p = [[0], [250], [500], [1000]]
    >>> sw.gpan(s, t, p)
    array([[   0.        ,    0.        ,    0.        ],
           [  56.35465209,   56.35465209,   56.35465209],
           [  84.67266947,   84.67266947,   84.67266947],
           [ 104.95799186,  104.95799186,  104.95799186]])

    References
    ----------
    .. [1] S. Pond & G.Pickard 2nd Edition 1986 Introductory Dynamical
    Oceanography Pergamon Press Sydney. ISBN 0-08-028728-X

    Modifications: 92-11-05. Phil Morgan.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)
    s, t, p = np.broadcast_arrays(s, t, p)

    if (s.ndim != 2)  and (t.ndim != 2):
        raise ValueError('Arguments must be 2D arrays: n_depths, n_profiles')

    svn = svan(s, t, p)
    mean_svan  = 0.5*(svn[1:,:] + svn[0:-1,:])

    top = svn[0,:] * p[0,:] * cte.db2Pascal

    delta_ga = ( mean_svan * np.diff(p, axis=0) ) * cte.db2Pascal
    ga = np.cumsum( np.vstack( ( top, delta_ga ) ), axis=0 )

    return ga

def gvel(ga, distm, lat):
    r"""
    Calculates geostrophic velocity given the geopotential anomaly
    and position of each station.

    Parameters
    ----------
    ga : array_like
         geopotential anomaly relative to the sea surface.
    dist : array_like
           distance between stations [meters]
    lat : array_like
          lat to use for constant f determination

    Returns
    -------
    vel : array_like
           geostrophic velocity relative to the sea surface.
           dimension will be MxN-1 (N: stations)

    See Also
    --------
    svan, gpan

    Notes
    -----
    The original matlab version had gvel and gvel2 only, here the logic is
    "gvel2", where one must compute the distance first.

    TODO: dim(m, nstations-1) or pass axis?
          add example with a reference level.
          example with values that make some sense.

    Examples
    --------
    >>> import numpy as np
    >>> import seawater.csiro as sw
    >>> lon = [30, 30, 30]
    >>> lat = [30, 32, 35]
    >>> s = np.array([[0, 1, 2], [15, 16, 17], [30, 31, 32],[35,35,35]])
    >>> t = np.repeat(15, s.size).reshape(s.shape)
    >>> p = [[0], [250], [500], [1000]]
    >>> ga = sw.gpan(s,t,p)
    >>> distm = 1000.0 * sw.dist(lon, lat, units='km')[0]
    >>> sw.gvel(ga, distm, lat)
    array([[-0.        , -0.        ],
           [ 0.11385677,  0.07154215],
           [ 0.22436555,  0.14112761],
           [ 0.33366412,  0.20996272]])

    References
    ----------
    .. [1] S. Pond & G.Pickard 2nd Edition 1986 Introductory Dynamical
    Oceanography Pergamon Press Sydney. ISBN 0-08-028728-X

    Modifications: 92-03-26. Phil Morgan.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    ga, distm, lat = np.asanyarray(ga), np.asanyarray(distm), np.asanyarray(lat)

    f     = cor( ( lat[0:-1] + lat[1:] )/2 )
    lf    = f * distm
    vel   = -np.diff(ga, axis=1) / lf

    return vel

def cp(s, t, p):
    r"""
    Heat Capacity of Sea Water using UNESCO 1983 polynomial.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].

    Returns
    -------
    cp : array_like
         specific heat capacity [J kg :sup:`-1` C :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO

    Examples
    --------
    Data from Pond and Pickard Intro. Dynamical Oceanography 2nd ed. 1986

    >>> import seawater.csiro as sw
    >>> t = T90conv([[0, 0, 0, 0, 0, 0], [10, 10, 10, 10, 10, 10], [20, 20, 20, 20, 20, 20], [30, 30, 30, 30, 30, 30], [40, 40, 40, 40, 40, 40]])
    >>> s = [[25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35]]
    >>> p = [0, 5000, 10000, 0, 5000, 10000]
    >>> sw.cp(s, t, p)
    array([[ 4048.4405375 ,  3896.25585   ,  3807.7330375 ,  3986.53309476,
             3849.26094605,  3769.11791286],
           [ 4041.8276691 ,  3919.5550066 ,  3842.3111366 ,  3986.34061786,
             3874.72665865,  3804.415624  ],
           [ 4044.8438591 ,  3938.5978466 ,  3866.7400391 ,  3993.85441786,
             3894.99294519,  3828.29059113],
           [ 4049.0984351 ,  3952.0375476 ,  3882.9855526 ,  4000.68382238,
             3909.24271128,  3844.32151784],
           [ 4051.2244911 ,  3966.1132036 ,  3905.9162711 ,  4003.46192541,
             3923.89463092,  3868.28959814]])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    Modifications: Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)

    p = p/10. # to convert [db] to [bar] as used in UNESCO routines
    T68 = T68conv(t)

    # eqn 26 p.32
    c0 = 4217.4
    c1 =   -3.720283
    c2 =    0.1412855
    c3 =   -2.654387e-3
    c4 =    2.093236e-5

    a0 = -7.64357
    a1 =  0.1072763
    a2 = -1.38385e-3

    b0 =  0.1770383
    b1 = -4.07718e-3
    b2 =  5.148e-5

    Cpst0 = ( ( ( ( c4 * T68 + c3 ) * T68 + c2 ) * T68 + c1 ) * T68 + c0 +
            ( a0 + a1 * T68 + a2 * T68**2 ) * s +
            ( b0 + b1 * T68 + b2 * T68**2 ) * s *(s)**0.5 )

    # eqn 28 p.33
    a0 = -4.9592e-1
    a1 =  1.45747e-2
    a2 = -3.13885e-4
    a3 =  2.0357e-6
    a4 =  1.7168e-8

    b0 =  2.4931e-4
    b1 = -1.08645e-5
    b2 =  2.87533e-7
    b3 = -4.0027e-9
    b4 =  2.2956e-11

    c0 = -5.422e-8
    c1 =  2.6380e-9
    c2 = -6.5637e-11
    c3 =  6.136e-13

    del_Cp0t0 = ( (((((c3 * T68 + c2) * T68 + c1) * T68 + c0) * p +
                ((((b4 * T68 + b3) * T68 + b2) * T68 + b1) * T68 + b0)) * p +
                ((((a4 * T68 + a3) * T68 + a2) * T68 + a1) * T68 + a0)) * p )

    # eqn 29 p.34
    d0 =  4.9247e-3
    d1 = -1.28315e-4
    d2 =  9.802e-7
    d3 =  2.5941e-8
    d4 = -2.9179e-10

    e0 = -1.2331e-4
    e1 = -1.517e-6
    e2 =  3.122e-8

    f0 = -2.9558e-6
    f1 =  1.17054e-7
    f2 = -2.3905e-9
    f3 =  1.8448e-11

    g0 =  9.971e-8

    h0 =  5.540e-10
    h1 = -1.7682e-11
    h2 =  3.513e-13

    j1 = -1.4300e-12

    S3_2 = s * (s)**0.5

    del_Cpstp = ( (((((d4 * T68 + d3) * T68 + d2) * T68 + d1) * T68 + d0) * s +
                ((e2 * T68 + e1) * T68 + e0) * S3_2) * p +
                ((((f3 * T68 + f2) * T68 + f1) * T68 + f0) * s +
                 g0 * S3_2 ) * p**2 +
                (((h2 * T68 + h1) * T68 + h0) * s +
                j1 * T68 * S3_2) * p**3 )

    cp = Cpst0 + del_Cp0t0 + del_Cpstp

    return cp

def ptmp(s, t, p, pr=0):
    r"""
    Calculates potential temperature as per UNESCO 1983 report.

    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].
    pr : array_like
        reference pressure [db], default = 0

    Returns
    -------
    pt : array_like
         potential temperature relative to PR [:math:`^\circ` C (ITS-90)]

    See Also
    --------
    adtg, pden

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.csiro as sw
    >>> t = T90conv([[0, 0, 0, 0, 0, 0], [10, 10, 10, 10, 10, 10], [20, 20, 20, 20, 20, 20], [30, 30, 30, 30, 30, 30], [40, 40, 40, 40, 40, 40]])
    >>> s = [[25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35],  [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35], [25, 25, 25, 35, 35, 35]]
    >>> p = [0, 5000, 10000, 0, 5000, 10000]
    >>> sw.T68conv(sw.ptmp(s, t, p, pr=0))
    array([[  0.        ,  -0.30614418,  -0.96669485,   0.        ,
             -0.3855565 ,  -1.09741136],
           [ 10.        ,   9.35306331,   8.46840949,  10.        ,
              9.29063461,   8.36425752],
           [ 20.        ,  19.04376281,  17.94265   ,  20.        ,
             18.99845171,  17.86536441],
           [ 30.        ,  28.75124632,  27.43529911,  30.        ,
             28.72313484,  27.38506197],
           [ 40.        ,  38.46068173,  36.92544552,  40.        ,
             38.44979906,  36.90231661]])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Bryden, H. 1973. New Polynomials for thermal expansion, adiabatic
    temperature gradient and potential temperature of sea water.
    Deep-Sea Res. Vol20,401-408. doi:10.1016/0011-7471(73)90063-6

    Modifications: 92-04-06. Phil Morgan.
                   99-06-25. Lindsay Pender, Fixed transpose of row vectors.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    s, t, p = np.asanyarray(s), np.asanyarray(t), np.asanyarray(p)
    pr = np.asanyarray(pr)

    # theta1
    del_P  = pr - p
    del_th = del_P * adtg(s, t, p)
    th     = T68conv(t) + 0.5 * del_th
    q      = del_th

    # theta2
    del_th = del_P * adtg(s, T90conv(th), p + 0.5 * del_P )
    th     = th + ( 1 - 1/(2)**00.5 ) * ( del_th - q )
    q      = ( 2 - (2)**0.5 ) * del_th + ( -2 + 3/(2)**0.5 ) * q

    # theta3
    del_th = del_P * adtg( s, T90conv(th), p + 0.5 * del_P )
    th     = th + ( 1 + 1/(2)**0.5 ) * ( del_th - q )
    q      = ( 2 + (2)**0.5 ) * del_th + ( -2 -3/(2)**0.5 ) * q

    # theta4
    del_th = del_P * adtg( s, T90conv(th), p + del_P )
    pt     = T90conv( th + ( del_th - 2 * q ) / 6 )
    return pt

def temp(s, pt, p, pr=0):
    r"""
    Calculates temperature from potential temperature at the reference pressure
    PR and in situ pressure P.

    Parameters
    ----------
    s(p) : array_like
        salinity [psu (PSS-78)]
    pt(p) : array_like
         potential temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [db].
    pr : array_like
         reference pressure [db]

    Returns
    -------
    temp : array_like
           temperature [:math:`^\circ` C (ITS-90)]

    See Also
    --------
    ptmp

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.csiro as sw
    >>> s, t, p = 35, 15, 100
    >>> sw.temp(s, sw.ptmp(s, t, p), p)
    15.0

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Bryden, H. 1973. New Polynomials for thermal expansion, adiabatic
    temperature gradient and potential temperature of sea water. Deep-Sea Res.
    Vol20,401-408. doi:10.1016/0011-7471(73)90063-6

    Modifications: 92-04-06. Phil Morgan.
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-19. Filipe Fernandes, Reformulated docstring.
    """

    s, pt, p = np.asanyarray(s), np.asanyarray(pt), np.asanyarray(p)
    pr  = np.asanyarray(pr)

    # Carry out inverse calculation by swapping p0 & pr
    t = ptmp(s, pt, pr, p)

    return t

def swvel(length, depth):
    r"""
    Calculates surface wave velocity.

    length : array_like
            wave length
    depth : array_like
            water depth [meters]

    Returns
    -------
    speed : array_like
            surface wave speed [m s :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    TODO: add my wave function to extras

    Examples
    --------
    >>> import seawater.csiro as sw
    >>> sw.swvel(10, 100)
    3.9493270848342941

    Modifications: Lindsay Pender 2005
                   10-01-14. Filipe Fernandes, Python translation.
                   10-08-25. Filipe Fernandes, Reformulated docstring.
    """

    length, depth = np.asanyarray(length), np.asanyarray(depth)

    k = 2.0 * np.pi / length
    speed = ( cte.gdef * np.tanh(k * depth) / k )**0.5
    return speed

def test(fileout='matlab-test.txt'):
    r"""
    Copy of the matlab test.

    Execute test routines to verify SEAWATER Library routines for your
    platform. Prints output to file.

    Notes
    ------
    This is only to reproduce sw_test.m from the original. A better more
    complete test is performed via doctest.

    Modifications: Phil Morgan
                   03-12-12. Lindsay Pender, Converted to ITS-90.
                   10-01-14. Filipe Fernandes, Python translation.
    """
    f = open(fileout,'w')

    print >>f, '**********************************************'
    print >>f, '    TEST REPORT    '
    print >>f, ''
    print >>f, ' SEA WATER LIBRARY ', __version__
    print >>f, ''
    print >>f, ''
    # Show some info about this Python
    print >>f, 'python version:', sys.version
    print >>f, ' on ', uname()[0],uname()[-1], ' computer'
    print >>f, ''
    print >>f,  asctime( localtime() )
    print >>f, '**********************************************'
    print >>f, ''

    # test MAIN MODULE  ptmp
    module     = 'ptmp'
    submodules = 'adtg'

    print >>f, '*************************************'
    print >>f, '**  TESTING MODULE: ', module
    print >>f, '**  and SUB-MODULE: ', submodules
    print >>f, '*************************************'

    # test 1 - data from Unesco 1983 p45
    T = np.array([[ 0,  0,  0,  0,  0,  0],
            [10, 10, 10, 10, 10, 10],
            [20, 20, 20, 20, 20, 20],
            [30, 30, 30, 30, 30, 30],
            [40, 40, 40, 40, 40, 40]])

    T = T / 1.00024

    S = np.array([[25, 25, 25, 35, 35, 35],
            [25, 25, 25, 35, 35, 35],
            [25, 25, 25, 35, 35, 35],
            [25, 25, 25, 35, 35, 35],
            [25, 25, 25, 35, 35, 35]])

    P = np.array([[0, 5000, 10000, 0, 5000, 10000],
            [0, 5000, 10000, 0, 5000, 10000],
            [0, 5000, 10000, 0, 5000, 10000],
            [0, 5000, 10000, 0, 5000, 10000],
            [0, 5000, 10000, 0, 5000, 10000]])

    Pr = np.array([0, 0, 0, 0, 0, 0])

    UN_ptmp = np.array([[ 0, -0.3061, -0.9667,  0, -0.3856, -1.0974],
                    [10,  9.3531,  8.4684, 10,  9.2906,  8.3643],
                    [20, 19.0438, 17.9426, 20, 18.9985, 17.8654],
                    [30, 28.7512, 27.4353, 30, 28.7231, 27.3851],
                    [40, 38.4607, 36.9254, 40, 38.4498, 36.9023]])

    PT = ptmp(S, T, P, Pr)*1.00024

    # DISPLAY RESULTS
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, 'Comparison of accepted values from UNESCO 1983 '
    print >>f, ' (Unesco Tech. Paper in Marine Sci. No. 44, p45)'
    print >>f, '********************************************************'

    m,n = S.shape  # TODO: so many loops there must be a better way...
    for icol in range(0, n):
        print >>f, '   Sal  Temp  Press     PTMP       ptmp'
        print >>f, '  (psu)  (C)   (db)     (C)          (C)'
        result = np.vstack( ( S[:,icol], T[:,icol], P[:,icol],
        UN_ptmp[:,icol], PT[:,icol] ) )
        for iline in range(0, m):
            print >>f, " %4.0f  %4.0f   %5.0f   %8.4f  %11.5f" % tuple(result[:,iline])

        print >>f, ''

    # test MAIN MODULE  svan
    module     = 'svan'
    submodules = 'dens dens0 smow seck pden ptmp'

    print >>f, '*************************************'
    print >>f, '**  TESTING MODULE: ', module
    print >>f, '**  and SUB-MODULE: ', submodules
    print >>f, '*************************************'

    # test DATA FROM: Unesco Tech. Paper in Marine Sci. No. 44, p22
    s = np.array([0,     0,  0,     0, 35,    35, 35,   35])
    p = np.array([0, 10000,  0, 10000,  0, 10000,  0, 10000])
    t = np.array([0,     0, 30,    30,  0,     0, 30,    30]) / 1.00024

    UN_svan = np.array([2749.54, 2288.61, 3170.58, 3147.85,
                        0.0,    0.00,  607.14,  916.34])

    SVAN    = svan(s, t, p)

    # DISPLAY RESULTS
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, 'Comparison of accepted values from UNESCO 1983'
    print >>f, ' (Unesco Tech. Paper in Marine Sci. No. 44, p22)'
    print >>f, '********************************************************'
    print >>f, ''
    print >>f, '   Sal  Temp  Press        SVAN        svan'
    print >>f, '  (psu)  (C)   (db)    (1e-8*m3/kg)  (1e-8*m3/kg)'
    result = np.vstack([s, t, p, UN_svan, 1e+8*SVAN])
    for iline in range( 0, len(SVAN) ):
        print >>f,  " %4.0f  %4.0f   %5.0f   %11.2f    %11.3f" % tuple(result[:,iline])

    # test MAIN MODULE salt
    module     = 'salt'
    submodules = 'salrt salrp sals'
    print >>f, '*************************************'
    print >>f, '**  TESTING MODULE: ', module
    print >>f, '**  and SUB-MODULE: ', submodules
    print >>f, '*************************************'

    # test 1 - data from Unesco 1983 p9
    R    = np.array([  1,       1.2,       0.65]) # cndr = R
    T    = np.array([ 15,        20,          5]) / 1.00024
    P    = np.array([  0,      2000,       1500])
    #Rt   = np.array([  1, 1.0568875, 0.81705885])
    UN_S = np.array([35, 37.245628,  27.995347])

    S    = salt(R, T, P)

    # DISPLAY RESULTS
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, 'Comparison of accepted values from UNESCO 1983 '
    print >>f, '(Unesco Tech. Paper in Marine Sci. No. 44, p9)'
    print >>f, '********************************************************'
    print >>f, ''
    print >>f, '   Temp    Press       R              S           salt'
    print >>f, '   (C)     (db)    (no units)       (psu)          (psu) '
    table = np.vstack([T, P, R, UN_S, S])
    m,n = table.shape
    for iline in range( 0, n ):
        print >>f, " %4.0f       %4.0f  %8.2f      %11.6f  %14.7f" % tuple(table[:,iline])

    # test MAIN MODULE cndr
    module     = 'cndr'
    submodules = 'salds'
    print >>f, '*************************************'
    print >>f, '**  TESTING MODULE: ', module
    print >>f, '**  and SUB-MODULE: ', submodules
    print >>f, '*************************************'

    # test 1 - data from Unesco 1983 p9
    T    = np.array([  0, 10, 0, 10, 10, 30]) / 1.00024
    P    = np.array([  0,  0, 1000, 1000, 0, 0])
    S    = np.array([ 25, 25, 25, 25, 40, 40])
    UN_R = np.array([ 0.498088, 0.654990, 0.506244, 0.662975, 1.000073, 1.529967])
    R    = cndr(S, T, P)

    # DISPLAY RESULTS
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, 'Comparison of accepted values from UNESCO 1983 '
    print >>f, ' (Unesco Tech. Paper in Marine Sci. No. 44, p14)'
    print >>f, '********************************************************'
    print >>f, ''
    print >>f, '   Temp    Press       S            cndr         cndr'
    print >>f, '   (C)     (db)      (psu)        (no units)    (no units) '
    table = np.vstack([T, P, S, UN_R, R])
    m,n = table.shape
    for iline in range( 0, n ):
        print >>f, " %4.0f       %4.0f   %8.6f   %11.6f  %14.8f" % tuple(table[:,iline])

    # test MAIN MODULE depth
    module     = 'depth'
    print >>f, ''
    print >>f, '*************************************'
    print >>f, '**  TESTING MODULE: ', module
    print >>f, '*************************************'

    # test DATA - matrix "pressure", vector "lat"  Unesco 1983 data p30.
    lat = np.array([0, 30, 45, 90])
    P   = np.array([[  500,   500,   500,  500],
                 [ 5000,  5000,  5000, 5000],
                 [10000, 10000, 10000, 10000]])

    UN_dpth = np.array([[  496.65,  496.00,  495.34,  494.03],
                    [ 4915.04, 4908.56, 4902.08, 4889.13],
                    [ 9725.47, 9712.65, 9699.84, 9674.23]])

    dpth = depth(P, lat)

    # DISPLAY RESULTS
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, 'Comparison of accepted values from Unesco 1983 '
    print >>f, '(Unesco Tech. Paper in Marine Sci. No. 44, p28)'
    print >>f, '********************************************************'

    for irow in range(0, 3):
        print >>f, ''
        print >>f, '    Lat       Press     DPTH      dpth'
        print >>f, '  (degree)    (db)     (meter)    (meter)'
        table = np.vstack( [ lat, P[irow,:], UN_dpth[irow,:], dpth[irow,:] ] )
        m,n   = table.shape
        for iline in range(0, n):
            print >>f, "  %6.3f     %6.0f   %8.2f   %8.3f" % tuple(table[:,iline])

    # test MAIN MODULE fp
    module     = 'fp'
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, '**  TESTING MODULE: ', module
    print >>f, '********************************************************'

    # test 1 - UNESCO DATA p.30
    S     = np.array([ [5, 10, 15, 20, 25, 30, 35, 40],
                    [5, 10, 15, 20, 25, 30, 35, 40] ])

    P     = np.array([ [  0,   0,   0,   0,   0,   0,   0,   0],
                    [500, 500, 500, 500, 500, 500, 500, 500] ])


    UN_fp = np.array([[-0.274, -0.542, -0.812, -1.083, -1.358, -1.638, -1.922, -2.212],
                    [-0.650, -0.919, -1.188, -1.460, -1.735, -2.014, -2.299, -2.589] ])

    FP    = fp(S, P)

    # DISPLAY RESULTS
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, 'Comparison of accepted values from UNESCO 1983 '
    print >>f, ' (Unesco Tech. Paper in Marine Sci. No. 44, p30)'
    print >>f, '********************************************************'

    for irow in range(0, 2):
        print >>f, ''
        print >>f, '   Sal   Press      fp        fp'
        print >>f, '  (psu)   (db)      (C)        (C)'
        table = np.vstack( [ S[irow,:], P[irow,:], UN_fp[irow,:], FP[irow,:] ] )
        m,n   = table.shape
        for iline in range( 0, n ):
            print >>f, " %4.0f   %5.0f   %8.3f  %11.4f" % tuple(table[:,iline])

    # test MAIN MODULE cp
    module     = 'cp'
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, '**  TESTING MODULE: ', module
    print >>f, '********************************************************'

    # test 1 -
    # DATA FROM POND AND PICKARD INTRO. DYNAMICAL OCEANOGRAPHY 2ND ED. 1986
    T     = np.array([[ 0,  0,  0,  0,  0,  0],
                   [10, 10, 10, 10, 10, 10],
                   [20, 20, 20, 20, 20, 20],
                   [30, 30, 30, 30, 30, 30],
                   [40, 40, 40, 40, 40, 40]]) / 1.00024

    S     = np.array([[25, 25, 25, 35, 35, 35],
                   [25, 25, 25, 35, 35, 35],
                   [25, 25, 25, 35, 35, 35],
                   [25, 25, 25, 35, 35, 35],
                   [25, 25, 25, 35, 35, 35]])

    P     = np.array([[0, 5000, 10000, 0, 5000, 10000],
                   [0, 5000, 10000, 0, 5000, 10000],
                   [0, 5000, 10000, 0, 5000, 10000],
                   [0, 5000, 10000, 0, 5000, 10000],
                   [0, 5000, 10000, 0, 5000, 10000]])

    UN_cp = np.array([[4048.4,  3896.3,  3807.7,  3986.5,  3849.3,  3769.1],
                   [4041.8,  3919.6,  3842.3,  3986.3,  3874.7,  3804.4],
                   [4044.8,  3938.6,  3866.7,  3993.9,  3895.0,  3828.3],
                   [4049.1,  3952.0,  3883.0,  4000.7,  3909.2,  3844.3],
                   [4051.2,  3966.1,  3905.9,  4003.5,  3923.9,  3868.3]])

    CP    = cp(S, T, P)

    # DISPLAY RESULTS
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, 'Comparison of accepted values from UNESCO 1983 '
    print >>f, ' (Unesco Tech. Paper in Marine Sci. No. 44, p37)'
    print >>f, '********************************************************'

    m,n = S.shape
    for icol in range(0, n):
        print >>f, ''
        print >>f, '   Sal  Temp  Press      Cp        cp'
        print >>f, '  (psu)  (C)   (db)    (J/kg.C)   (J/kg.C)'
        result = np.vstack( [ S[:,icol], T[:,icol], P[:,icol],
                            UN_cp[:,icol], CP[:,icol] ] )
        for iline in range(0, m):
            print >>f, " %4.0f  %4.0f   %5.0f   %8.1f  %11.2f" % tuple(result[:,iline])

    # test MAIN MODULE svel
    module     = 'svel'
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, '**  TESTING MODULE: ', module
    print >>f, '********************************************************'

    # test 1 -
    # DATA FROM POND AND PICKARD INTRO. DYNAMICAL OCEANOGRAPHY 2ND ED. 1986
    T       = np.array([[  0,  0,  0,  0,  0,  0],
                     [ 10, 10, 10, 10, 10, 10],
                     [ 20, 20, 20, 20, 20, 20],
                     [ 30, 30, 30, 30, 30, 30],
                     [ 40, 40, 40, 40, 40, 40]]) / 1.00024

    S       = np.array([[ 25, 25, 25, 35, 35, 35],
                     [ 25, 25, 25, 35, 35, 35],
                     [ 25, 25, 25, 35, 35, 35],
                     [ 25, 25, 25, 35, 35, 35],
                     [ 25, 25, 25, 35, 35, 35]])

    P       = np.array([[ 0, 5000, 10000, 0, 5000, 10000],
                     [ 0, 5000, 10000, 0, 5000, 10000],
                     [ 0, 5000, 10000, 0, 5000, 10000],
                     [ 0, 5000, 10000, 0, 5000, 10000],
                     [ 0, 5000, 10000, 0, 5000, 10000]])

    UN_svel = np.array([[ 1435.8, 1520.4, 1610.4, 1449.1, 1534.0, 1623.2],
                     [ 1477.7, 1561.3, 1647.4, 1489.8, 1573.4, 1659.0],
                     [ 1510.3, 1593.6, 1676.8, 1521.5, 1604.5, 1687.2],
                     [ 1535.2, 1619.0, 1700.6, 1545.6, 1629.0, 1710.1],
                     [ 1553.4, 1638.0, 1719.2, 1563.2, 1647.3, 1727.8]])

    SVEL    = svel(S, T, P)

    # DISPLAY RESULTS
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, 'Comparison of accepted values from UNESCO 1983 '
    print >>f, ' (Unesco Tech. Paper in Marine Sci. No. 44, p50)'
    print >>f, '********************************************************'

    m,n = SVEL.shape
    for icol in range(0, n):
        print >>f, ''
        print >>f, '   Sal  Temp  Press     SVEL       svel'
        print >>f, '  (psu)  (C)   (db)     (m/s)       (m/s)'

        result = np.vstack( [ S[:,icol], T[:,icol], P[:,icol], UN_svel[:,icol], SVEL[:,icol] ] )
        for iline in range(0, m):
            print >>f, " %4.0f  %4.0f   %5.0f   %8.1f  %11.3f" % tuple(result[:,iline])

    # test SUBMODULES alpha beta aonb
    submodules     = 'alpha beta aonb'
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, '**  and SUB-MODULE: ', submodules
    print >>f, '********************************************************'

    # DATA FROM MCDOUOGALL 1987
    s    = 40
    PT   = 10
    p    = 4000
    beta_lit  = 0.72088e-03
    aonb_lit  = 0.34763
    alpha_lit = aonb_lit*beta_lit

    BETA  = beta( s, PT, p, pt=True)
    ALPHA = alpha(s, PT, p, pt=True)
    AONB  = aonb( s, PT, p, pt=True)

    # DISPLAY RESULTS
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, 'Comparison of accepted values from MCDOUGALL 1987 '
    print >>f, '********************************************************'

    print >>f, ''
    print >>f, '   Sal  Temp  Press     BETA       beta'
    print >>f, '  (psu)  (C)   (db)   (psu^-1)     (psu^-1)'
    table = np.hstack( [ s, PT, p, beta_lit, BETA ] )
    print >>f, " %4.0f  %4.0f   %5.0f   %11.4e  %11.5e" % tuple(table)

    print >>f, ''
    print >>f, '   Sal  Temp  Press     AONB       aonb'
    print >>f, '  (psu)  (C)   (db)   (psu C^-1)   (psu C^-1)'
    table = np.hstack( [s, PT, p, aonb_lit, AONB] )
    print >>f, " %4.0f  %4.0f   %5.0f   %8.5f  %11.6f" % tuple(table)

    print >>f, ''
    print >>f, '   Sal  Temp  Press     ALPHA       alpha'
    print >>f, '  (psu)  (C)   (db)    (psu^-1)     (psu^-1)'
    table = np.hstack( [ s, PT, p, alpha_lit, ALPHA ] )
    print >>f, " %4.0f  %4.0f   %5.0f   %11.4e  %11.4e" % tuple(table)

    # test MAIN MODULES  satO2 satN2 satAr
    module     = 'satO2 satN2 satAr'
    print >>f, ''
    print >>f, '********************************************************'
    print >>f, '**  TESTING MODULE: ', module
    print >>f, '********************************************************'

    # Data from Weiss 1970
    T      = np.array([[ -1, -1],
                    [ 10, 10],
                    [ 20, 20],
                    [ 40, 40]]) / 1.00024

    S      = np.array([[ 20, 40],
                    [ 20, 40],
                    [ 20, 40],
                    [ 20, 40]])

    lit_O2 = np.array([[ 9.162, 7.984],
                    [ 6.950, 6.121],
                    [ 5.644, 5.015],
                    [ 4.050, 3.656]])

    lit_N2 =  np.array([[ 16.28, 14.01],
                     [ 12.64, 11.01],
                     [ 10.47,  9.21],
                     [  7.78,  6.95]])

    lit_Ar =  np.array([[ 0.4456, 0.3877],
                     [ 0.3397, 0.2989],
                     [ 0.2766, 0.2457],
                     [ 0.1986, 0.1794]])

    O2     = satO2(S, T)
    N2     = satN2(S, T)
    Ar     = satAr(S, T)

    # DISPLAY RESULTS
    print >>f, ''
    print >>f, '************************************************************'
    print >>f, 'Comparison of accepted values from Weiss, R.F. 1979 '
    print >>f, '"The solubility of nitrogen, oxygen and argon in water'
    print >>f, ' and seawater." Deep-Sea Research., 1970, Vol 17, pp721-735.'
    print >>f, '************************************************************'

    m,n = S.shape
    for icol in range(0, n):
        print >>f, ''
        print >>f, '   Sal  Temp      O2         satO2'
        print >>f, '  (psu)  (C)      (ml/l)     (ml/l)'
        result = np.vstack( [ S[:,icol], T[:,icol],
                            lit_O2[:,icol], O2[:,icol] ] )
        for iline in range(0, m):
            print >>f, " %4.0f  %4.0f    %8.2f   %9.3f" % tuple(result[:,iline])

    for icol in range(0, n):
        print >>f, ''
        print >>f, '   Sal  Temp      N2         satN2'
        print >>f, '  (psu)  (C)      (ml/l)     (ml/l)'
        result = np.vstack( [ S[:,icol], T[:,icol],
                            lit_N2[:,icol], N2[:,icol] ] )
        for iline in range(0, m):
            print >>f, " %4.0f  %4.0f    %8.2f  %9.3f" % tuple(result[:,iline])

    for icol in range(0, n):
        print >>f, ''
        print >>f, '   Sal  Temp      Ar         satAr'
        print >>f, '  (psu)  (C)      (ml/l)     (ml/l)'
        result = np.vstack( [ S[:,icol], T[:,icol],
                              lit_Ar[:,icol], Ar[:,icol] ] )
        for iline in range(0, m):
            print >>f, " %4.0f  %4.0f     %8.4f  %9.4f" % tuple(result[:,iline])

if __name__=='__main__':
    import doctest
    doctest.testmod()
