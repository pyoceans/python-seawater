# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np
from seawater import constants as cte
import os

class Dict2Struc(object):
    r"""
    Open variables from a dictionary in a "matlab-like-structure"
    """
    def __init__(self, adict):
        self.__dict__.update(adict)

class Cache_npz(object):
    def __init__(self):
        self._cache = dict()
        self._default_path = os.path.join(os.path.dirname(__file__), 'data')

    def __call__(self, fname, datadir=None):
        if datadir is None:
            datadir = self._default_path
        fpath = os.path.join(datadir, fname)
        try:
            return self._cache[fpath]
        except KeyError:
            pass
        d = np.load(fpath)
        ret = Dict2Struc(d)
        self._cache[fpath] = ret
        return ret

_npz_cache = Cache_npz()

def read_data(fname, datadir=None):
    """
    Read variables from a numpy '.npz' file into a minimal
    class providing attribute access.

    A cache is used to avoid re-reading the same file.
    """
    return _npz_cache(fname, datadir=datadir)


def strip_mask(*args):
    """
    Process the standard arguments for efficient calculation.

    Return unmasked argments, plus a mask.

    The first argument, SA, is handled specially so that it can be

    This could be absorbed into a decorator, but it would
    require redefining functions to take the additional
    mask argument or kwarg.
    """
    mask = np.ma.getmaskarray(args[-1])
    SA = args[0]
    if SA.shape:
        SA = np.ma.asarray(SA)
        SA[ SA < 0] = np.ma.masked
        for a in args[:-1]:
            mask = np.ma.mask_or(mask, np.ma.getmask(a))
        newargs = [SA.filled(0)]
    elif SA < 0:
        SA = 0
        for a in args[1:-1]:
            mask = np.ma.mask_or(mask, np.ma.getmask(a))
        newargs = [SA]
    newargs.extend([np.ma.filled(a, 0) for a in args[1:]])
    newargs.append(mask)
    return newargs

def _gibbs(ns, nt, npr, SA, t, p):
    r"""
    Calculates specific Gibbs energy and its derivatives up to order 2 for
    seawater.

    The Gibbs function approach allows the calculation of internal energy,
    entropy, enthalpy, potential enthalpy and the chemical potentials of
    seawater as well as the freezing temperature, and the latent heats of
    freezing and of evaporation. These quantities were not available from
    EOS-80 but are essential for the accurate accounting of heat in the ocean
    and for the consistent and accurate treatment of air-sea and ice-sea heat
    fluxes.

    Parameters
    ----------
    ns : int
         order of SA derivative [0, 1 or 2 ]
    nt : int
         order of t derivative [0, 1 or 2 ]
    npr : int
          order of p derivative [0, 1 or 2 ]
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    gibbs : array_like
            Specific Gibbs energy or its derivatives.

            Gibbs energy (ns=nt=npr=0) has units of:
            [J kg :sup:`-1`]

            Absolute Salinity derivatives are output in units of:
            [(J kg :sup:`-1`) (g kg :sup:`-1`) :sup:`-ns`]

            Temperature derivatives are output in units of:
            [(J kg :sup:`-1`) K :sup:`-nt`]

            Pressure derivatives are output in units of:
            [(J kg :sup:`-1`) Pa :sup:`-npr`]

            The mixed derivatives are output in units of:
            [(J kg :sup:`-1`) (g kg :sup:`-1`) :sup:`-ns` K :sup:`-nt` Pa :sup:`-npr`]

    Notes
    -----
    The Gibbs function for seawater is that of TEOS-10 (IOC et al., 2010),
    being the sum of IAPWS-08 for the saline part and IAPWS-09 for the pure
    water part. These IAPWS releases are the officially blessed IAPWS
    descriptions of Feistel (2008) and the pure water part of Feistel (2003).
    Absolute Salinity, SA, in all of the GSW routines is expressed on the
    Reference-Composition Salinity Scale of 2008 (RCSS-08) of Millero et al.
    (2008).

    The derivatives are taken with respect to pressure in Pa, not withstanding
    that the pressure input into this routine is in dbar.

    References
    ----------
    .. [1] Feistel, R., 2003: A new extended Gibbs thermodynamic potential of
    seawater, Progr. Oceanogr., 58, 43-114.

    .. [2] Feistel, R., 2008: A Gibbs function for seawater thermodynamics
    for -6 to 80 :math:`^\circ` C and salinity up to 120 g kg :sup:`-1`,
    Deep-Sea Res. I, 55, 1639-1671.

    .. [3] IAPWS, 2008: Release on the IAPWS Formulation 2008 for the
    Thermodynamic Properties of Seawater. The International Association for the
    Properties of Water and Steam. Berlin, Germany, September 2008, available
    from http://www.iapws.org.  This Release is referred to as IAPWS-08.

    .. [4] IAPWS, 2009: Supplementary Release on a Computationally Efficient
    Thermodynamic Formulation for Liquid Water for Oceanographic Use. The
    International Association for the Properties of Water and Steam. Doorwerth,
    The Netherlands, September 2009, available from http://www.iapws.org.
    This Release is referred to as IAPWS-09.

    .. [5] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp. See section 2.6 and appendices A.6,  G and H.

    .. [6] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008:
    The composition of Standard Seawater and the definition of the
    Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.

    Modifications:
    2010-09-24. David Jackett, Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    # trick to use single number
    if SA.ndim == 0:
        SA = SA[np.newaxis]
    if (SA.size == 1 and SA[0] <= 0) or np.ma.all(SA <=0):
        nonzero_SA = False
    else:
        nonzero_SA = True
        SA = np.ma.masked_less_equal(SA, 0)
    _SA = SA
    _t = t
    _p = p

    SA = np.ma.filled(SA, 0)
    t = np.ma.filled(t, 20)
    p = np.ma.filled(p, 10)

    SA, t, p = np.broadcast_arrays(SA, t, p)

    gibbs = np.zeros(SA.shape, dtype=np.float) #use if all_masked becomes True
    all_masked = False

    # Ensure a full mask, so we can set elements if necessary.
    mask = np.ma.mask_or(np.ma.getmaskarray(_SA), np.ma.getmask(_t))
    mask = np.ma.mask_or(mask, np.ma.getmask(_p))
    mask = np.ma.mask_or(mask, SA<0)

    ipos = (SA > 0)
    inpos = ~ipos
    if np.all(ipos):
        ipos = slice(None) # more efficient for usual case


    x2 = cte.sfac * SA
    x = np.sqrt(x2)

    y = t * 0.025
    z = p * 1e-4 # The input pressure (p) is sea pressure in units of dbar.

    if (ns==0) & (nt==0) & (npr==0):
        g03 =( 101.342743139674 + z * ( 100015.695367145 +
        z * ( -2544.5765420363 + z * ( 284.517778446287 +
        z * ( -33.3146754253611 + ( 4.20263108803084 - 0.546428511471039 * z )
        * z ) ) ) ) +
        y * ( 5.90578347909402 + z * ( -270.983805184062 +
        z * ( 776.153611613101 + z * ( -196.51255088122 + ( 28.9796526294175 -
        2.13290083518327 * z ) * z ) ) ) +
        y * ( -12357.785933039 + z * ( 1455.0364540468 +
        z * ( -756.558385769359 + z * ( 273.479662323528 + z *
        ( -55.5604063817218 + 4.34420671917197 * z ) ) ) ) +
        y * ( 736.741204151612 + z * ( -672.50778314507 +
        z * ( 499.360390819152 + z * ( -239.545330654412 + ( 48.8012518593872 -
        1.66307106208905 * z ) * z ) ) ) +
        y * ( -148.185936433658 + z * ( 397.968445406972 +
        z * ( -301.815380621876 + ( 152.196371733841 - 26.3748377232802 * z ) *
        z ) ) +
        y * ( 58.0259125842571 + z * ( -194.618310617595 +
        z * ( 120.520654902025 + z * ( -55.2723052340152 +
        6.48190668077221 * z ) ) ) +
        y * ( -18.9843846514172 + y * ( 3.05081646487967 -
        9.63108119393062 * z ) +
        z * ( 63.5113936641785 + z * ( -22.2897317140459 +
        8.17060541818112 * z ) ) ) ) ) ) ) ) )

        if nonzero_SA:
            g08 = x2 * ( 1416.27648484197 + z * ( -3310.49154044839 +
            z * ( 384.794152978599 + z * ( -96.5324320107458 +
            ( 15.8408172766824 - 2.62480156590992 * z ) * z ) ) ) +
            x * ( -2432.14662381794 + x * ( 2025.80115603697 +
            y * ( 543.835333000098 + y * ( -68.5572509204491 +
            y * ( 49.3667694856254 + y * ( -17.1397577419788 +
            2.49697009569508 * y ) ) ) - 22.6683558512829 * z ) +
            x * ( -1091.66841042967 - 196.028306689776 * y +
            x * ( 374.60123787784 - 48.5891069025409 * x +
            36.7571622995805 * y ) + 36.0284195611086 * z ) +
            z * ( -54.7919133532887 + ( -4.08193978912261 -
            30.1755111971161 * z ) * z ) ) +
            z * ( 199.459603073901 + z * ( -52.2940909281335 +
            ( 68.0444942726459 - 3.41251932441282 * z ) * z ) ) +
            y * ( -493.407510141682 + z * ( -175.292041186547 +
            ( 83.1923927801819 - 29.483064349429 * z ) * z ) +
            y * ( -43.0664675978042 + z * ( 383.058066002476 + z *
            ( -54.1917262517112 + 25.6398487389914 * z ) ) +
            y * ( -10.0227370861875 - 460.319931801257 * z + y *
            ( 0.875600661808945 + 234.565187611355 * z ) ) ) ) )  +
            y * ( 168.072408311545 + z * ( 729.116529735046 +
            z * ( -343.956902961561 + z * ( 124.687671116248 + z *
            ( -31.656964386073 + 7.04658803315449 * z ) ) ) ) +
            y * ( 880.031352997204 + y * ( -225.267649263401 +
            y * ( 91.4260447751259 + y * ( -21.6603240875311 +
            2.13016970847183 * y ) +
            z * ( -297.728741987187 + ( 74.726141138756 -
            36.4872919001588 * z ) * z ) ) +
            z * ( 694.244814133268 + z * ( -204.889641964903 +
            ( 113.561697840594 - 11.1282734326413 * z ) * z ) ) ) +
            z * ( -860.764303783977 + z * ( 337.409530269367 +
            z * ( -178.314556207638 + ( 44.2040358308 -
            7.92001547211682 * z ) * z ) ) ) ) ) )

            g08[ipos] += x2[ipos] * ( 5812.81456626732 + 851.226734946706 *
            y[ipos] ) * np.log( x[ipos] )
        else:
            g08 = 0
        gibbs = g03 + g08


    elif (ns==1) & (nt==0) & (npr==0):
        if nonzero_SA:
            g08 =( 8645.36753595126 + z * ( -6620.98308089678 +
            z * ( 769.588305957198 + z * ( -193.0648640214916 +
            ( 31.6816345533648 - 5.24960313181984 * z ) * z ) ) ) +
            x * ( -7296.43987145382 + x * ( 8103.20462414788 +
            y * ( 2175.341332000392 + y * ( -274.2290036817964 +
            y * ( 197.4670779425016 + y * ( -68.5590309679152 +
            9.98788038278032 * y ) ) ) - 90.6734234051316 * z ) +
            x * ( -5458.34205214835 - 980.14153344888 * y +
            x * ( 2247.60742726704 - 340.1237483177863 * x +
            220.542973797483 * y ) + 180.142097805543 * z ) +
            z * ( -219.1676534131548 + ( -16.32775915649044 -
            120.7020447884644 * z ) * z ) ) +
            z * ( 598.378809221703 + z * ( -156.8822727844005 +
            ( 204.1334828179377 - 10.23755797323846 * z ) * z ) ) +
            y * ( -1480.222530425046 + z * ( -525.876123559641 +
            ( 249.57717834054571 - 88.449193048287 * z ) * z ) +
            y * ( -129.1994027934126 + z * ( 1149.174198007428 +
            z * ( -162.5751787551336 + 76.9195462169742 * z ) ) +
            y * ( -30.0682112585625 - 1380.9597954037708 * z + y *
            ( 2.626801985426835 + 703.695562834065 * z ) ) ) ) ) +
            y * ( 1187.3715515697959 + z * ( 1458.233059470092 +
            z * ( -687.913805923122 + z * ( 249.375342232496 + z *
            ( -63.313928772146 + 14.09317606630898 * z ) ) ) ) +
            y * ( 1760.062705994408 + y * ( -450.535298526802 +
            y * ( 182.8520895502518 + y * ( -43.3206481750622 +
            4.26033941694366 * y ) +
            z * ( -595.457483974374 + ( 149.452282277512 -
            72.9745838003176 * z ) * z ) ) +
            z * ( 1388.489628266536 + z * ( -409.779283929806 +
            ( 227.123395681188 - 22.2565468652826 * z ) * z ) ) ) +
            z * ( -1721.528607567954 + z * ( 674.819060538734 +
            z * ( -356.629112415276 + ( 88.4080716616 -
            15.84003094423364 * z ) * z ) ) ) ) ) )

            g08[ipos] = g08[ipos] + ( 11625.62913253464 + 1702.453469893412 *
            y[ipos] ) * np.log( x[ipos] )

            gibbs = 0.5 * cte.sfac * g08
        else:
            all_masked = True

    elif (ns==0) & (nt==1) & (npr==0):
        g03 = ( 5.90578347909402 + z * ( -270.983805184062 +
        z * ( 776.153611613101 + z * ( -196.51255088122 +
        ( 28.9796526294175 - 2.13290083518327 * z ) * z ) ) ) +
        y * ( -24715.571866078 + z * ( 2910.0729080936 +
        z * ( -1513.116771538718 + z * ( 546.959324647056 + z *
        ( -111.1208127634436 + 8.68841343834394 * z ) ) ) ) +
        y * ( 2210.2236124548363 + z * ( -2017.52334943521 +
        z * ( 1498.081172457456 + z * ( -718.6359919632359 +
        ( 146.4037555781616 - 4.9892131862671505 * z ) * z ) ) ) +
        y * ( -592.743745734632 + z * ( 1591.873781627888 +
        z * ( -1207.261522487504 + ( 608.785486935364 -
        105.4993508931208 * z ) * z ) ) +
        y * ( 290.12956292128547 + z * ( -973.091553087975 +
        z * ( 602.603274510125 + z * ( -276.361526170076 +
        32.40953340386105 * z ) ) ) +
        y * ( -113.90630790850321 + y * ( 21.35571525415769 -
        67.41756835751434 * z ) +
        z * ( 381.06836198507096 + z * ( -133.7383902842754 +
        49.023632509086724 * z ) ) ) ) ) ) ) )

        if nonzero_SA:
            g08 = x2 * ( 168.072408311545 + z * ( 729.116529735046 +
            z * ( -343.956902961561 + z * ( 124.687671116248 + z *
            ( -31.656964386073 + 7.04658803315449 * z ) ) ) ) +
            x * ( -493.407510141682 + x * ( 543.835333000098 + x *
            ( -196.028306689776 + 36.7571622995805 * x ) +
            y * ( -137.1145018408982 + y * ( 148.10030845687618 + y *
            ( -68.5590309679152 + 12.4848504784754 * y ) ) ) -
            22.6683558512829 * z ) + z * ( -175.292041186547 +
            ( 83.1923927801819 - 29.483064349429 * z ) * z ) +
            y * ( -86.1329351956084 + z * ( 766.116132004952 + z *
            ( -108.3834525034224 + 51.2796974779828 * z ) ) +
            y * ( -30.0682112585625 - 1380.9597954037708 * z + y *
            ( 3.50240264723578 + 938.26075044542 * z ) ) ) ) +
            y * ( 1760.062705994408 + y * ( -675.802947790203 +
            y * ( 365.7041791005036 + y * ( -108.30162043765552 +
            12.78101825083098 * y ) +
            z * ( -1190.914967948748 + ( 298.904564555024 -
            145.9491676006352 * z ) * z ) ) +
            z * ( 2082.7344423998043 + z * ( -614.668925894709 +
            ( 340.685093521782 - 33.3848202979239 * z ) * z ) ) ) +
            z * ( -1721.528607567954 + z * ( 674.819060538734 +
            z * ( -356.629112415276 + ( 88.4080716616 -
            15.84003094423364 * z ) * z ) ) ) ) )

            g08[ipos] += 851.226734946706 * x2[ipos] * np.log( x[ipos] )

            gibbs = (g03 + g08) * 0.025
        else:
            gibbs = g03

    elif (ns==0) & (nt==0) & (npr==1):
        g03 = ( 100015.695367145 + z * ( -5089.1530840726 +
        z * ( 853.5533353388611 + z * ( -133.2587017014444 +
        ( 21.0131554401542 - 3.278571068826234 * z ) * z ) ) ) +
        y * ( -270.983805184062 + z * ( 1552.307223226202 +
        z * ( -589.53765264366 + ( 115.91861051767 -
        10.664504175916349 * z ) * z ) ) +
        y * ( 1455.0364540468 + z * ( -1513.116771538718 +
        z * ( 820.438986970584 + z * ( -222.2416255268872 +
        21.72103359585985 * z ) ) ) +
        y * ( -672.50778314507 + z * ( 998.720781638304 +
        z * ( -718.6359919632359 + ( 195.2050074375488 -
        8.31535531044525 * z ) * z ) ) +
        y * ( 397.968445406972 + z * ( -603.630761243752 +
        ( 456.589115201523 - 105.4993508931208 * z ) * z ) +
        y * ( -194.618310617595 + y * ( 63.5113936641785 -
        9.63108119393062 * y +
        z * ( -44.5794634280918 + 24.511816254543362 * z ) ) +
        z * ( 241.04130980405 + z * ( -165.8169157020456 +
        25.92762672308884 * z ) ) ) ) ) ) ) )

        if nonzero_SA:
            g08 = x2 * ( -3310.49154044839 + z * ( 769.588305957198 +
            z * ( -289.5972960322374 + ( 63.3632691067296 -
            13.1240078295496 * z ) * z ) ) +
            x * ( 199.459603073901 + x * ( -54.7919133532887 +
            36.0284195611086 * x - 22.6683558512829 * y +
            ( -8.16387957824522 - 90.52653359134831 * z ) * z ) +
            z * ( -104.588181856267 + ( 204.1334828179377 -
            13.65007729765128 * z ) * z ) +
            y * ( -175.292041186547 + ( 166.3847855603638 -
            88.449193048287 * z ) * z +
            y * ( 383.058066002476 + y * ( -460.319931801257 +
            234.565187611355 * y ) +
            z * ( -108.3834525034224 + 76.9195462169742 * z ) ) ) ) +
            y * ( 729.116529735046 + z * ( -687.913805923122 +
            z * ( 374.063013348744 + z * ( -126.627857544292 +
            35.23294016577245 * z ) ) )  +
            y * ( -860.764303783977 + y * ( 694.244814133268 +
            y * ( -297.728741987187 + ( 149.452282277512 -
            109.46187570047641 * z ) * z ) +
            z * ( -409.779283929806 + ( 340.685093521782 -
            44.5130937305652 * z ) * z ) ) +
            z * ( 674.819060538734 + z * ( -534.943668622914 +
            ( 176.8161433232 - 39.600077360584095 * z ) * z ) ) ) ) )
        else:
            g08 = 0
        # Pressure derivative of the Gibbs function
        # in units of (J kg :sup:`-1`) (Pa :sup:`-1`) = m :sup:`3` kg :sup:`-1`
        gibbs = (g03 + g08) * 1e-8

    elif (ns==1) & (nt==1) & (npr==0):
        if nonzero_SA:
            g08 = ( 1187.3715515697959 + z * ( 1458.233059470092 +
            z * ( -687.913805923122 + z * ( 249.375342232496 + z *
            ( -63.313928772146 + 14.09317606630898 * z ) ) ) ) +
            x * ( -1480.222530425046 + x * ( 2175.341332000392 + x *
            ( -980.14153344888 + 220.542973797483 * x ) +
            y * ( -548.4580073635929 + y * ( 592.4012338275047 + y *
            ( -274.2361238716608 + 49.9394019139016 * y ) ) ) -
            90.6734234051316 * z ) + z * ( -525.876123559641 +
            ( 249.57717834054571 - 88.449193048287 * z ) * z ) +
            y * ( -258.3988055868252 + z * ( 2298.348396014856 +
            z * ( -325.1503575102672 + 153.8390924339484 * z ) ) +
            y * ( -90.2046337756875 - 4142.8793862113125 * z + y *
            ( 10.50720794170734 + 2814.78225133626 * z ) ) ) ) +
            y * ( 3520.125411988816 + y * ( -1351.605895580406 +
            y * ( 731.4083582010072 + y * ( -216.60324087531103 +
            25.56203650166196 * y ) +
            z * ( -2381.829935897496 + ( 597.809129110048 -
            291.8983352012704 * z ) * z ) ) +
            z * ( 4165.4688847996085 + z * ( -1229.337851789418 +
            ( 681.370187043564 - 66.7696405958478 * z ) * z ) ) ) +
            z * ( -3443.057215135908 + z * ( 1349.638121077468 +
            z * ( -713.258224830552 + ( 176.8161433232 -
            31.68006188846728 * z ) * z ) ) ) ) )

            g08[ipos] = g08[ipos] + 1702.453469893412 * np.log( x[ipos] )
            gibbs = 0.5 * cte.sfac * 0.025 * g08
            #mask[inpos] = True FIXME: commented by FF, g110 without nan didn't pass
        else:
            all_masked = True

    elif (ns==1) & (nt==0) & (npr==1):
        g08 = ( -6620.98308089678 + z * ( 1539.176611914396 +
        z * ( -579.1945920644748 + ( 126.7265382134592 -
        26.2480156590992 * z ) * z ) ) +
        x * ( 598.378809221703 + x * ( -219.1676534131548 +
        180.142097805543 * x - 90.6734234051316 * y +
        (-32.65551831298088 - 362.10613436539325 * z ) * z ) +
        z * ( -313.764545568801 + ( 612.4004484538132 -
        40.95023189295384 * z ) * z ) +
        y * ( -525.876123559641 + ( 499.15435668109143 -
        265.347579144861 * z ) * z +
        y * ( 1149.174198007428 + y * ( -1380.9597954037708 +
        703.695562834065 * y ) +
        z * ( -325.1503575102672 + 230.7586386509226 * z ) ) ) ) +
        y * ( 1458.233059470092 + z * ( -1375.827611846244 +
        z * ( 748.126026697488 + z * ( -253.255715088584 +
        70.4658803315449 * z ) ) )  +
        y * ( -1721.528607567954 + y * ( 1388.489628266536 +
        y * ( -595.457483974374 + ( 298.904564555024 -
        218.92375140095282 * z ) * z ) +
        z * ( -819.558567859612 + ( 681.370187043564 -
        89.0261874611304 * z ) * z ) ) +
        z * ( 1349.638121077468 + z * ( -1069.887337245828 +
        ( 353.6322866464 - 79.20015472116819 * z ) * z ) ) ) ) )

        # Derivative of the Gibbs function is in units of
        # (m :sup:`3` kg :sup:`-1`) / (g kg :sup:`-1`) = m :sup:`3` g :sup:`-1`
        # that is, it is the derivative of specific volume with respect to
        # Absolute Salinity measured in g kg :sup:`-1`

        gibbs = g08 * cte.sfac * 0.5e-8

    elif (ns==0) & (nt==1) & (npr==1):
        g03 = ( -270.983805184062 + z * ( 1552.307223226202 + z *
        ( -589.53765264366 +
        ( 115.91861051767 - 10.664504175916349 * z ) * z ) ) +
        y * ( 2910.0729080936 + z * ( -3026.233543077436 +
        z * ( 1640.877973941168 + z * ( -444.4832510537744 +
        43.4420671917197 * z ) ) ) +
        y * ( -2017.52334943521 + z * ( 2996.162344914912 +
        z * ( -2155.907975889708 + ( 585.6150223126464 -
        24.946065931335752 * z ) * z ) ) +
        y * ( 1591.873781627888 + z * ( -2414.523044975008 +
        ( 1826.356460806092 - 421.9974035724832 * z ) * z ) +
        y * ( -973.091553087975 + z * ( 1205.20654902025 + z *
        ( -829.084578510228 + 129.6381336154442 * z ) ) +
        y * ( 381.06836198507096 - 67.41756835751434 * y + z *
        ( -267.4767805685508 + 147.07089752726017 * z ) ) ) ) ) ) )

        if nonzero_SA:
            g08 = x2 * ( 729.116529735046 + z * ( -687.913805923122 +
            z * ( 374.063013348744 + z * ( -126.627857544292 +
            35.23294016577245 * z ) ) ) +
            x * ( -175.292041186547 - 22.6683558512829 * x +
            ( 166.3847855603638 - 88.449193048287 * z ) * z +
            y * ( 766.116132004952 + y * ( -1380.9597954037708 +
            938.26075044542 * y ) +
            z * ( -216.7669050068448 + 153.8390924339484 * z ) ) ) +
            y * ( -1721.528607567954 + y * ( 2082.7344423998043 +
            y * ( -1190.914967948748 + ( 597.809129110048 -
            437.84750280190565 * z ) * z ) +
            z * ( -1229.337851789418 + ( 1022.055280565346 -
            133.5392811916956 * z ) * z ) ) +
            z * ( 1349.638121077468 + z * ( -1069.887337245828 +
            ( 353.6322866464 - 79.20015472116819 * z ) * z ) ) ) )
        else:
            g08 = 0
        # Derivative of the Gibbs function is in units of (m :sup:`3` (K kg) )
        # that is, the pressure of the derivative in Pa.
        gibbs = (g03 + g08) * 2.5e-10

    elif (ns==2) & (nt==0) & (npr==0):
        g08 = 2.0 * ( 8103.20462414788 +
        y * ( 2175.341332000392 + y * ( -274.2290036817964 +
        y * ( 197.4670779425016 + y * ( -68.5590309679152 +
        9.98788038278032 * y ) ) ) - 90.6734234051316 * z ) +
        1.5 * x * ( -5458.34205214835 - 980.14153344888 * y +
        ( 4.0 / 3.0 ) * x * ( 2247.60742726704 - 340.1237483177863 * 1.25 *
        x + 220.542973797483 * y ) +
        180.142097805543 * z ) +
        z * ( -219.1676534131548 + ( -16.32775915649044 -
        120.7020447884644 * z ) * z ) )

        if nonzero_SA:
            tmp = ( ( -7296.43987145382 + z * ( 598.378809221703 +
            z * ( -156.8822727844005 + ( 204.1334828179377 -
            10.23755797323846 * z ) * z ) ) +
            y * ( -1480.222530425046 + z * ( -525.876123559641 +
            ( 249.57717834054571 - 88.449193048287 * z ) * z ) +
            y * ( -129.1994027934126 + z * ( 1149.174198007428 +
            z * ( -162.5751787551336 + 76.9195462169742 * z ) ) +
            y * ( -30.0682112585625 - 1380.9597954037708 * z +
            y * ( 2.626801985426835 + 703.695562834065 * z ) ) ) ) ) / x +
            ( 11625.62913253464 + 1702.453469893412 * y ) )
            g08[ipos] += tmp[ipos] / x2[ipos]

        gibbs = 0.25 * cte.sfac**2 * g08

    elif (ns==0) & (nt==2) & (npr==0):
        g03 = ( -24715.571866078 + z * ( 2910.0729080936 + z *
        ( -1513.116771538718 + z * ( 546.959324647056 + z *
        ( -111.1208127634436 + 8.68841343834394 * z ) ) ) ) +
        y * ( 4420.4472249096725 + z * ( -4035.04669887042 +
        z * ( 2996.162344914912 + z * ( -1437.2719839264719 +
        ( 292.8075111563232 - 9.978426372534301 * z ) * z ) ) ) +
        y * ( -1778.231237203896 + z * ( 4775.621344883664 +
        z * ( -3621.784567462512 + ( 1826.356460806092 -
        316.49805267936244 * z ) * z ) ) +
        y * ( 1160.5182516851419 + z * ( -3892.3662123519 +
        z * ( 2410.4130980405 + z * ( -1105.446104680304 +
        129.6381336154442 * z ) ) ) +
        y * ( -569.531539542516 + y * ( 128.13429152494615 -
        404.50541014508605 * z ) +
        z * ( 1905.341809925355 + z * ( -668.691951421377 +
        245.11816254543362 * z ) ) ) ) ) ) )

        if nonzero_SA:
            g08 = x2 * ( 1760.062705994408 + x * ( -86.1329351956084 +
            x * ( -137.1145018408982 + y * ( 296.20061691375236 +
            y * ( -205.67709290374563 + 49.9394019139016 * y ) ) )  +
            z * ( 766.116132004952 + z * ( -108.3834525034224 +
            51.2796974779828 * z ) ) +
            y * ( -60.136422517125 - 2761.9195908075417 * z +
            y * ( 10.50720794170734 + 2814.78225133626 * z ) ) ) +
            y * ( -1351.605895580406 + y * ( 1097.1125373015109 +
            y * ( -433.20648175062206 + 63.905091254154904 * y ) +
            z * ( -3572.7449038462437 + ( 896.713693665072 -
            437.84750280190565 * z ) * z ) ) +
            z * ( 4165.4688847996085 + z * ( -1229.337851789418 +
            ( 681.370187043564 - 66.7696405958478 * z ) * z ) ) ) +
            z * ( -1721.528607567954 + z * ( 674.819060538734 +
            z * ( -356.629112415276 + ( 88.4080716616 -
            15.84003094423364 * z ) * z ) ) ) )
        else:
            g08 = 0
        gibbs = (g03 + g08) * 0.000625

    elif (ns==0) & (nt==0) & (npr==2):
        g03 = ( -5089.1530840726 + z * ( 1707.1066706777221 +
        z * ( -399.7761051043332 + ( 84.0526217606168 -
        16.39285534413117 * z ) * z ) ) +
        y * ( 1552.307223226202 + z * ( -1179.07530528732 +
        ( 347.75583155301 - 42.658016703665396 * z ) * z ) +
        y * ( -1513.116771538718 + z * ( 1640.877973941168 +
        z * ( -666.7248765806615 + 86.8841343834394 * z ) ) +
        y * ( 998.720781638304 + z * ( -1437.2719839264719 +
        ( 585.6150223126464 - 33.261421241781 * z ) * z ) +
        y * ( -603.630761243752 + ( 913.178230403046 -
        316.49805267936244 * z ) * z +
        y * ( 241.04130980405 + y * ( -44.5794634280918 +
        49.023632509086724 * z ) +
        z * ( -331.6338314040912 + 77.78288016926652 * z ) ) ) ) ) ) )

        if nonzero_SA:
            g08 = x2 * ( 769.588305957198 + z * ( -579.1945920644748 +
            ( 190.08980732018878 - 52.4960313181984 * z ) * z ) +
            x * ( -104.588181856267 + x * ( -8.16387957824522 -
            181.05306718269662 * z ) +
            ( 408.2669656358754 - 40.95023189295384 * z ) * z +
            y * ( 166.3847855603638 - 176.898386096574 * z + y *
            ( -108.3834525034224 + 153.8390924339484 * z ) ) ) +
            y * ( -687.913805923122 + z * ( 748.126026697488 +
            z * ( -379.883572632876 + 140.9317606630898 * z ) ) +
            y * ( 674.819060538734 + z * ( -1069.887337245828 +
            ( 530.4484299696 - 158.40030944233638 * z ) * z ) +
            y * ( -409.779283929806 + y * ( 149.452282277512 -
            218.92375140095282 * z ) +
            ( 681.370187043564 - 133.5392811916956 * z ) * z ) ) ) )
        else:
            g08 = 0
        # Second derivative of the Gibbs function with respect to pressure,
        # measured in Pa; units of (J kg :sup:`-1`) (Pa :sup:`-2`).
        gibbs = (g03 + g08) * 1e-16
    else:
        raise NameError('Wrong Combination of order/variables')

    gibbs = np.ma.array(gibbs, mask=mask, copy=False)
    if all_masked:
        gibbs[:] = np.ma.masked
    return gibbs

def _entropy_part(SA, t, p):
    r"""
    Calculates entropy, except that it does not evaluate any terms that are
    functions of Absolute Salinity alone.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    entropy_part : array_like
                   entropy minus the terms that due to SA alone
                   [J kg :sup:`-1` K :sup:`-1`]

    Notes
    -----
    By not calculating these terms, which are a function only of Absolute
    Salinity, several unnecessary computations are avoided (including saving
    the computation of a natural logarithm). These terms are a necessary part
    of entropy, but are not needed when calculating potential temperature from
    in situ temperature.

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p, mask = strip_mask(SA, t, p)

    x2 = cte.sfac * SA
    x = np.sqrt(x2)
    y = t * 0.025
    z = p * 1e-4

    g03 = ( z * ( -270.983805184062 +
    z * ( 776.153611613101 + z * ( -196.51255088122 + ( 28.9796526294175 -
    2.13290083518327 * z ) * z ) ) ) +
    y * ( -24715.571866078 + z * ( 2910.0729080936 +
    z * ( -1513.116771538718 + z * ( 546.959324647056 + z *
    ( -111.1208127634436 + 8.68841343834394 * z ) ) ) ) +
    y * ( 2210.2236124548363 + z * ( -2017.52334943521 +
    z * ( 1498.081172457456 + z * ( -718.6359919632359 +
    ( 146.4037555781616 - 4.9892131862671505 * z ) * z ) ) ) +
    y * ( -592.743745734632 + z * ( 1591.873781627888 +
    z * ( -1207.261522487504 + ( 608.785486935364 -
    105.4993508931208 * z ) * z ) ) +
    y * ( 290.12956292128547 + z * ( -973.091553087975 +
    z * ( 602.603274510125 + z * ( -276.361526170076 +
    32.40953340386105 * z ) ) ) +
    y * ( -113.90630790850321 + y *
    ( 21.35571525415769 - 67.41756835751434 * z ) +
    z * ( 381.06836198507096 + z * ( -133.7383902842754 +
    49.023632509086724 * z ) ) ) ) ) ) ) )

    # TODO? short-circuit this if SA is zero
    g08 = x2 * ( z * ( 729.116529735046 +
    z * ( -343.956902961561 + z * ( 124.687671116248 + z * ( -31.656964386073 +
    7.04658803315449 * z ) ) ) ) +
    x * ( x * ( y * ( -137.1145018408982 + y * ( 148.10030845687618 +
    y * ( -68.5590309679152 + 12.4848504784754 * y ) ) ) -
    22.6683558512829 * z ) + z * ( -175.292041186547 +
    ( 83.1923927801819 - 29.483064349429 * z ) * z ) +
    y * ( -86.1329351956084 + z * ( 766.116132004952 +
    z * ( -108.3834525034224 + 51.2796974779828 * z ) ) +
    y * ( -30.0682112585625 - 1380.9597954037708 * z +
    y * ( 3.50240264723578 + 938.26075044542 * z ) ) ) ) +
    y * ( 1760.062705994408 + y * ( -675.802947790203 +
    y * ( 365.7041791005036 + y * ( -108.30162043765552 +
    12.78101825083098 * y ) +
    z * ( -1190.914967948748 + ( 298.904564555024 -
    145.9491676006352 * z ) * z ) ) +
    z * ( 2082.7344423998043 + z * ( -614.668925894709 +
    ( 340.685093521782 - 33.3848202979239 * z ) * z ) ) ) +
    z * ( -1721.528607567954 + z * ( 674.819060538734 +
    z * ( -356.629112415276 + ( 88.4080716616 -
    15.84003094423364 * z ) * z ) ) ) ) )

    entropy_part = -( g03 + g08 )  * 0.025

    return np.ma.array(entropy_part, mask=mask, copy=False)

def _gibbs_pt0_pt0(SA, pt0):
    r"""
    Calculates the second derivative of the specific Gibbs function with
    respect to temperature at zero sea pressure or _gibbs(0,2,0,SA,t,0).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt0 : array_like
          potential temperature relative to 0 dbar [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    gibbs_pt0_pt0 : array_like
                    TODO: write the eq. for the second derivative of the
                    specific Gibbs function. FIXME: [units]

    Notes
    -----
    This library function is called by both "pt_from_CT(SA,CT)"
    and "pt0_from_t(SA,t,p)".

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, pt0, mask = strip_mask(SA, pt0)

    x2 = cte.sfac * SA
    x = np.sqrt(x2)
    y = pt0 * 0.025

    g03 = ( -24715.571866078 +
    y * ( 4420.4472249096725 +
    y * ( -1778.231237203896 +
    y * ( 1160.5182516851419 +
    y * ( -569.531539542516 + y * 128.13429152494615) ) ) ) )

    g08 = x2 * ( 1760.062705994408 + x * ( -86.1329351956084 +
    x * ( -137.1145018408982 + y * ( 296.20061691375236 +
    y * ( -205.67709290374563 + 49.9394019139016 * y ) ) ) +
    y * ( -60.136422517125 + y * 10.50720794170734 ) ) +
    y * ( -1351.605895580406 + y * ( 1097.1125373015109 +
    y * ( -433.20648175062206 + 63.905091254154904 * y ) ) ) )

    gibbs_pt0_pt0 = ( g03 + g08 ) * 0.000625

    return np.ma.array(gibbs_pt0_pt0, mask=mask, copy=False)

def _entropy_part_zerop(SA, pt0):
    r"""
    Calculates entropy at a sea surface (p = 0 dbar), except that it does not
    evaluate any terms that are functions of Absolute Salinity alone.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt0 : array_like
          potential temperature relative to 0 dbar [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    entropy_part_zerop : array_like
                         [J kg :sup:`-1` K :sup:`-1`]

    Notes
    -----
    By not calculating these terms, which are a function only of Absolute
    Salinity, several unnecessary computations are avoided (including saving
    the computation of a natural logarithm). These terms are a necessary part
    of entropy, but are not needed when calculating potential temperature from
    in situ temperature.

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, pt0, mask = strip_mask(SA, pt0)

    x2 = cte.sfac * SA
    x = np.sqrt(x2)
    y = pt0 * 0.025

    g03 = y * ( -24715.571866078 + y * ( 2210.2236124548363 +
    y * ( -592.743745734632 + y * ( 290.12956292128547 +
    y * ( -113.90630790850321 + y * 21.35571525415769) ) ) ) )

    g08 = x2 * ( x * ( x * ( y * ( -137.1145018408982 + y *
    ( 148.10030845687618 +
    y * ( -68.5590309679152 + 12.4848504784754 * y ) ) ) ) +
    y * ( -86.1329351956084 + y * ( -30.0682112585625 + y *
    3.50240264723578 ) ) ) +
    y * ( 1760.062705994408 + y * ( -675.802947790203 +
    y * ( 365.7041791005036 + y * ( -108.30162043765552 +
    12.78101825083098 * y ) ) ) ) )

    entropy_part_zerop = -( g03 + g08 ) * 0.025

    return np.ma.array(entropy_part_zerop, mask=mask, copy=False)

def _enthalpy_SSO_0_CT25(p):
    r"""
    Calculates enthalpy at the Standard Ocean Salinity (SSO) and at a
    Conservative Temperature of zero degrees C (CT=0), as a function of
    pressure (p [dbar]) or enthalpy_CT25(35.16504,0,p).

    Parameters
    ----------
    p : array_like
        pressure [dbar]

    Returns
    -------
    enthalpy_CT25 : array_like
                    enthalpy_CT25 at (SSO, CT = 0, p), 25-term equation.
                    [J kg :sup:`-1`]

    Notes
    -----
    Uses a streamlined version of the 25-term CT version of the Gibbs function,
    that is, a streamlined version of the code "enthalpy_CT25(SA,CT,p)"

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    p = np.asanyarray(p)
    mask = np.ma.getmask(p)
    p = np.ma.filled(p, 0)

    SSO = cte.SSO

    a0 = 1 + SSO * (2.0777716085618458e-3 + np.sqrt(SSO) *
    3.4688210757917340e-6)
    a1 = 6.8314629554123324e-6
    b0 = 9.9984380290708214e2 + SSO * (2.8925731541277653e0 + SSO *
    1.9457531751183059e-3)
    b1 = 0.5 * (1.1930681818531748e-2 + SSO * 5.9355685925035653e-6)
    b2 = -2.5943389807429039e-8
    A = b1 - np.sqrt(b1**2 - b0 * b2)
    B = b1 + np.sqrt(b1**2 - b0 * b2)

    part = ( a0 * b2 - a1 * b1) / (b2 * (B - A) )

    enthalpy_SSO_0_CT25 = cte.db2Pascal * ( ( a1 / (2*b2) ) *
    np.log( 1 + p * ( 2 * b1 + b2 * p ) / b0 ) + part *
    np.log( 1 + ( b2 * p * (B - A) ) / (A * (B + b2 * p ) ) ) )

    return np.ma.array(enthalpy_SSO_0_CT25, mask=mask, copy=False)

def _specvol_SSO_0_CT25(p):
    r"""
    Calculates specific volume at the Standard Ocean Salinity (SSO) and
    Conservative Temperature of zero degrees C (CT=0), as a function of
    pressure (p [dbar]) or spec_vol_CT25(35.16504,0,p).

    Parameters
    ----------
    p : array_like
        pressure [dbar]

    Returns
    -------
    specvol_SSO_0_CT25 : array_like
                         Specific volume at (SSO, CT=0, p), 25-term equation.
                         [m :sup:`3` kg :sup:`-1`]

    Notes
    -----
    It uses a streamlined version of the 25-term CT version of specific volume
    that is, a streamlined version of the code "rho_alpha_beta_CT25(SA,CT,p)"

    Modifications
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    p = np.asanyarray(p)
    # No need to strip mask and replace it here; the calculation is simple.

    SSO = cte.SSO
    specvol_SSO_0_CT25 = ( (1.00000000e+00 + SSO * ( 2.0777716085618458e-003 +
    np.sqrt(SSO) * 3.4688210757917340e-006) + p * 6.8314629554123324e-006) /
    (9.9984380290708214e+002 + SSO * ( 2.8925731541277653e+000 + SSO *
    1.9457531751183059e-003) + p * ( 1.1930681818531748e-002 + SSO *
    5.9355685925035653e-006 + p * -2.5943389807429039e-008) ) )

    return specvol_SSO_0_CT25

"""
Salinity lib functions
"""
def _SP_from_SA_Baltic(SA, lon, lat):
    r"""
    Calculates Practical Salinity (SP) for the Baltic Sea, from a value
    computed analytically from Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
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
    SP_from_SA, SP_from_Sstar

    Notes
    -----
    This program will only produce Practical Salinity values for the Baltic Sea.

    Examples
    --------
    >>> import seawater.library as lib
    >>> SA = [6.6699, 6.7738, 6.9130, 7.3661, 7.5862, 10.3895]
    >>> lon, lat = 20, 59
    >>> lat = 59
    >>> lib._SP_from_SA_Baltic(SA, lon, lat)
    masked_array(data = [6.56825466873 6.67192351682 6.8108138311 7.26290579519 7.4825161269
     10.2795794748],
                 mask = [False False False False False False],
           fill_value = 1e+20)
    <BLANKLINE>

    References
    ----------
    .. [1] Feistel, R., S. Weinreben, H. Wolf, S. Seitz, P. Spitzer, B. Adel,
    G. Nausch, B. Schneider and D. G. Wright, 2010c: Density and Absolute
    Salinity of the Baltic Sea 2006-2009.  Ocean Science, 6, 3-24.
    http://www.ocean-sci.net/6/3/2010/os-6-3-2010.pdf

    .. [2] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    .. [3] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm
    for estimating Absolute Salinity in the global ocean. Submitted to Ocean
    Science. A preliminary version is available at Ocean Sci. Discuss.,
    6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """
    SA, lon, lat = map(np.ma.masked_invalid, (SA, lon, lat))
    lon, lat, SA = np.broadcast_arrays(lon, lat, SA)

    xb1, xb2, xb3 = 12.6, 7., 26.
    xb1a, xb3a = 45., 26.
    yb1, yb2, yb3 = 50., 59., 69.

    inds_baltic = (xb2 < lon) & (lon < xb1a) & (yb1 < lat) & (lat < yb3)
    if not inds_baltic.sum():
        return None

    SP_baltic = np.ma.masked_all(SA.shape, dtype=np.float)

    xx_left = np.interp( lat[inds_baltic], [yb1,yb2,yb3], [xb1,xb2,xb3])
    xx_right = np.interp( lat[inds_baltic], [yb1,yb3], [xb1a,xb3a] )

    inds_baltic1 = (   (xx_left <= lon[inds_baltic])
                     & (lon[inds_baltic] <= xx_right))

    if not inds_baltic1.sum():
        return None

    SP_baltic[inds_baltic[inds_baltic1]] = ( ( 35 / ( cte.SSO - 0.087 ) )
                             * ( SA[inds_baltic[inds_baltic1]] - 0.087) )

    return SP_baltic

def _SA_from_SP_Baltic(SP, lon, lat):
    r"""
    Calculates Absolute Salinity in the Baltic Sea, from Practical Salinity.
    Since SP is non-negative by definition, this function changes any negative
    input values of SP to be zero.

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
    SA_baltic : masked array, or None if there are no Baltic positions
                Absolute salinity [g kg :sup:`-1`]

    See Also
    --------
    SA_from_SP, Sstar_from_SP, SA_Sstar_from_SP

    Notes
    -----
    This program will only produce Absolute Salinity values for the Baltic Sea.

    Examples
    --------
    >>> import seawater.library as lib
    >>> SP = [6.5683, 6.6719, 6.8108, 7.2629, 7.4825, 10.2796]
    >>> lon, lat = 20, 59
    >>> lib._SA_from_SP_Baltic(SP, lon, lat)
    masked_array(data = [6.66994543234 6.77377643074 6.91298613806 7.36609419189 7.58618383714
     10.389520571],
                 mask = [False False False False False False],
           fill_value = 1e+20)
    <BLANKLINE>

    References
    ----------
    .. [1] Feistel, R., S. Weinreben, H. Wolf, S. Seitz, P. Spitzer, B. Adel,
    G. Nausch, B. Schneider and D. G. Wright, 2010c: Density and Absolute
    Salinity of the Baltic Sea 2006-2009.  Ocean Science, 6, 3-24.
    http://www.ocean-sci.net/6/3/2010/os-6-3-2010.pdf

    .. [2] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    .. [3] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm
    for estimating Absolute Salinity in the global ocean. Submitted to Ocean
    Science. A preliminary version is available at Ocean Sci. Discuss.,
    6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SP, lon, lat = map(np.ma.masked_invalid, (SP, lon, lat))
    lon, lat, SP = np.broadcast_arrays(lon, lat, SP)

    xb1, xb2, xb3 = 12.6, 7.0, 26.0
    xb1a, xb3a = 45.0, 26.0
    yb1, yb2, yb3 = 50.0, 59.0, 69.0

    inds_baltic = (xb2 < lon) & (lon < xb1a) & (yb1 < lat) & (lat < yb3)
    if not inds_baltic.sum():
        return None

    SA_baltic = np.ma.masked_all(SP.shape, dtype=np.float)

    xx_left = np.interp( lat[inds_baltic], [yb1,yb2,yb3], [xb1,xb2,xb3])
    xx_right = np.interp( lat[inds_baltic], [yb1,yb3], [xb1a,xb3a] )

    inds_baltic1 = (   (xx_left <= lon[inds_baltic])
                     & (lon[inds_baltic] <= xx_right))
    if not inds_baltic1.sum():
        return None

    SA_baltic[inds_baltic[inds_baltic1]] = ((( cte.SSO - 0.087 ) / 35 )
                                    * SP[inds_baltic[inds_baltic1]] + 0.087)

    return SA_baltic


def _delta_SA(p, lon, lat):
    r"""
    Calculates the Absolute Salinity anomaly, SA - SR, in the open ocean by
    spatially interpolating the global reference data set of delta_SA to the
    location of the seawater sample.

    Parameters
    ----------
    p : array_like, maximum 1D
        pressure [dbar]
    lon : array_like, maximum 1D
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like, maximum 1D
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    delta_SA : masked array; masked where no nearby ocean is found in data
               Absolute Salinity anomaly [g kg :sup:`-1`]

    See Also
    --------
    _dsa_add_barrier, _dsa_add_mean

    Notes
    -----
    The Absolute Salinity Anomaly in the Baltic Sea is evaluated separately,
    since it is a function of Practical Salinity, not of space. The present
    function returns a delta_SA of zero for data in the Baltic Sea. The correct
    way of calculating Absolute Salinity in the Baltic Sea is by calling
    SA_from_SP.

    The mask is only set when the observation is well and truly on dry
    land; often the warning flag is not set until one is several hundred
    kilometers inland from the coast.

    Examples
    --------
    >>> import seawater.library as lib
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat = 188, 4
    >>> lib._delta_SA(p, lon, lat)
    masked_array(data = [0.000167785807437 0.00026867590804 0.000665539507353 0.00269430342286
     0.00562666390947 0.00939665321653],
                 mask = [False False False False False False],
           fill_value = 1e+20)
    <BLANKLINE>

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
    of seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm
    for estimating Absolute Salinity in the global ocean.  Submitted to Ocean
    Science. A preliminary version is available at Ocean Sci. Discuss.,
    6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    ????-??-??. David Jackett.
    2010-07-23. Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    2011-02-10. Bjørn Ådlandsvik, Several bug fixes. Arguments must be scalars or 1D arrays.
    """

    # Input argument handling
    # -----------------------

    # For now, convert masked to unmasked, and use nan internally.
    # Maybe reverse the strategy later, or, more likely, use separate
    # data and mask internally.
    p, lon, lat = [np.ma.filled(var, np.nan) for var in (p, lon, lat)]
    p   = np.atleast_1d(p).astype(np.float) # must be float for interpolation
    # Make all arrays of same shape, raise ValueError if not compatible
    p, lon, lat = np.broadcast_arrays(p, lon, lat)

    # hack to force 1-D (FIXME)
    shape = p.shape
    p, lon, lat = p.flatten(), lon.flatten(), lat.flatten()

    # Check 1D
    if p.ndim != 1:
        raise ValueError, 'Arguments must be scalars or 1D arrays'

    # Put longitudes in 0-360
    lon %= 360

    # Read data file
    # --------------
    data = read_data("gsw_data_v2_0.npz")

    delta_SA_ref = np.ma.masked_invalid(data.delta_SA_ref)
    lats_ref = data.lats_ref
    longs_ref = data.longs_ref
    p_ref = data.p_ref                # Depth levels
    # Local number of depth levels
    ndepth_ref = np.ma.masked_invalid(data.ndepth_ref).astype(np.int8)
    ndepth_ref -= 1 # change to 0-based index

    # Grid resolution
    dlongs_ref = longs_ref[1] - longs_ref[0]
    dlats_ref = lats_ref[1] - lats_ref[0]

    # Find horizontal indices bracketing the position
    # -----------------------------------------------

    # Find indsx0 such that
    #   lons_ref[indsx0] <= lon < lons_ref[indsx0+1]
    # Border cases:
    #   indsx0 = lons_ref.size - 2 for
    indsx0 = (lon-longs_ref[0]) / dlongs_ref
    indsx0 = indsx0.astype(np.int)
    indsx0 = np.clip(indsx0, 0, longs_ref.size-2)

    # Find indsy0 such that
    #   lats_ref[indsy0] <= lat < lats_ref[indsy0+1]
    # Border cases:
    #   indsy0 = 0                 for lat < -86 = lats_refs[0]
    #   indsy0 = lats_ref.size - 2 for lat = 90 = lats_refs[-1]
    indsy0 = (lat-lats_ref[0]) / dlats_ref
    indsy0 = indsy0.astype(np.int)
    indsy0 = np.clip(indsy0, 0, lats_ref.size-2)

    indsz0 = np.searchsorted(p_ref, p, side='right') - 1

    nmax = np.c_[ ndepth_ref[indsy0, indsx0],
                  ndepth_ref[indsy0, indsx0+1],
                  ndepth_ref[indsy0+1, indsx0+1],
                  ndepth_ref[indsy0+1, indsx0] ]
    nmax = nmax.max(axis=1)

    deepmask = indsz0 > nmax
    p[deepmask] = p_ref[nmax[deepmask]]

    indsz0 = np.clip(indsz0, 0, p_ref.size-2)

    inds0 = (indsz0 + indsy0 * delta_SA_ref.shape[0]
                    + indsx0 * delta_SA_ref.shape[0] * delta_SA_ref.shape[1])

    data_indices = np.c_[indsx0, indsy0, indsz0, inds0]
    data_inds = data_indices[:,2]

    r1 = ( lon - longs_ref[indsx0] ) / dlongs_ref
    s1 = ( lat - lats_ref[indsy0] )  / dlats_ref
    t1 = ( p - p_ref[indsz0] ) / ( p_ref[indsz0+1] - p_ref[indsz0] )

    nksum = 0
    no_levels_missing = 0

    sa_upper = np.nan * ( np.ones(data_inds.shape) )
    sa_lower = np.nan * ( np.ones(data_inds.shape) )
    delta_SA = np.nan * ( np.ones(data_inds.shape) )

    for k in range(indsz0.min(), indsz0.max()+1):
        inds_k = np.where(indsz0 == k)[0]
        nk = inds_k.size

        if nk > 0:
            nksum = nksum + nk
            indsx = indsx0[inds_k]
            indsy = indsy0[inds_k]
            indsz = k * np.ones( indsx.shape, dtype='int64' )
            inds_di = (data_inds == k) # level k interpolation
            dsa = np.nan * np.ones( (4, p.size) )

            dsa[0, inds_k] = delta_SA_ref[indsz, indsy, indsx]
            dsa[1, inds_k] = delta_SA_ref[indsz, indsy, indsx+1]
            dsa[2, inds_k] = delta_SA_ref[indsz, indsy+1, indsx+1]
            dsa[3, inds_k] = delta_SA_ref[indsz, indsy+1, indsx]

            inds = ( (260. <= lon) & (lon <= 295.217) & (0. <= lat) & (lat <=
                      19.55) & (indsz0 == k) )
            """ TODO: describe add_barrier """
            if inds.any():
                dsa[:,inds] = _dsa_add_barrier( dsa[:,inds], lon[inds],
                lat[inds], longs_ref[indsx0[inds]], lats_ref[indsy0[inds]],
                dlongs_ref, dlats_ref)

            inds = np.where( (np.isnan(np.sum(dsa, axis=0))) & (indsz0==k) )[0]
            """ TODO: describe add_mean """
            if inds.size !=0:
                dsa[:,inds] = _dsa_add_mean(dsa[:,inds])

            # level k+1 interpolation
            sa_upper[inds_di] = ( ( 1 - s1[inds_di] ) * ( dsa[0, inds_k] +
            r1[inds_di] * ( dsa[1, inds_k] - dsa[0, inds_k] ) ) +
            s1[inds_di] * ( dsa[3, inds_k] +
            r1[inds_di] * ( dsa[2, inds_k] - dsa[3,inds_k] ) ) )

            dsa = np.nan * np.ones( (4, p.size) )
            dsa[0, inds_k] = delta_SA_ref[indsz+1, indsy, indsx]
            dsa[1, inds_k] = delta_SA_ref[indsz+1, indsy, indsx+1]
            dsa[2, inds_k] = delta_SA_ref[indsz+1, indsy+1, indsx+1]
            dsa[3, inds_k] = delta_SA_ref[indsz+1, indsy+1, indsx]

            inds = ( (260. <= lon) & (lon <= 295.217) & (0 <= lat) & (lat <=
            19.55) & (indsz0 == k) )
            """ TODO: describe add_barrier """
            if inds.any():
                dsa[:,inds] = _dsa_add_barrier( dsa[:,inds], lon[inds],
                lat[inds], longs_ref[indsx0[inds]], lats_ref[indsy0[inds]],
                dlongs_ref, dlats_ref)

            inds = ( np.isnan( np.sum(dsa, axis=0) ) ) & (indsz0==k)
            """ TODO: describe add_mean """
            dsa[:,inds] = _dsa_add_mean(dsa[:,inds])

            sa_lower[inds_di] = ( ( 1 - s1[inds_di] ) * ( dsa[0, inds_k] +
            r1[inds_di] * ( dsa[1, inds_k] - dsa[0,inds_k] ) ) +
            s1[inds_di] * ( dsa[3, inds_k] +
            r1[inds_di] * ( dsa[2, inds_k] - dsa[3, inds_k] ) ) )

            inds_diff = np.where( np.isfinite( sa_upper[inds_di] ) &
                                       np.isnan( sa_lower[inds_di] ) )[0]
            inds_di = np.where(inds_di)[0] #TODO: terrible solution, but works
            if inds_diff.size != 0:
                sa_lower[inds_di[inds_diff]] = sa_upper[inds_di[inds_diff]]

            delta_SA[inds_di] = ( sa_upper[inds_di] + t1[inds_di] *
                                ( sa_lower[inds_di] - sa_upper[inds_di] ) )
        else:
            no_levels_missing = no_levels_missing + 1

    delta_SA = np.ma.masked_invalid(delta_SA)
    # hack to force 1-D (FIXME)
    return np.reshape(delta_SA, shape)

def _dsa_add_barrier(dsa, lon, lat, longs_ref, lats_ref, dlongs_ref, dlats_ref):
    r"""
    Adds a barrier through Central America (Panama) and then averages over the
    appropriate side of the barrier.

    Parameters
    ----------
    dsa : array_like
          Absolute Salinity anomaly of the 4 adjacent neighbors [g kg :sup:`-1`]
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
          Absolute Salinity anomaly of the 4 adjacent neighbors  [g kg :sup:`-1`]

    Notes
    -----
    originally inside "_delta_SA"

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """
    longs_pan = np.array([260.0000, 272.5900, 276.5000, 278.6500,
                                              280.7300, 295.2170])
    lats_pan = np.array([19.5500, 13.9700, 9.6000, 8.1000, 9.3300, 0])

    lats_lines0 = np.interp( lon, longs_pan, lats_pan )
    lats_lines1 = np.interp( longs_ref, lats_pan, longs_pan )
    lats_lines2 = np.interp( (longs_ref+dlongs_ref), lats_pan, longs_pan )

    above_line = np.bool_( np.ones(4) )
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
            above_line[3] = False

        if lats_lines2[k0] <= lats_ref[k0]:
            above_line[1] = True
        else:
            above_line[1] = False

        if lats_lines2[k0] <= (lats_ref[k0] + dlats_ref):
            above_line[2] = True
        else:
            above_line[2] = False

        # indices of different sides of CA line
        inds = ( above_line != above_line0 )
        dsa[inds,k0] = np.nan

    dsa_mean = dsa.mean()
    inds_nan = np.isnan( dsa_mean )

    if inds_nan.any():
        no_nan = len(np.where(inds_nan))

        for kk in range(0,no_nan):
            col = inds_nan[kk]
            inds_kk = np.where( np.isnan( dsa[:,col] ) )[0]
            Inn = np.where( ~np.isnan( dsa[:,col] ) )[0]

            if Inn.size == 0:
                dsa[inds_kk,col] = dsa[Inn,col].mean()


    delta_SA = dsa
    return delta_SA

def _dsa_add_mean(dsa):
    r"""
    Replaces NaN's with nanmean of the 4 adjacent neighbors

    Parameters
    ----------
    dsa : array_like
          Absolute Salinity anomaly of the 4 adjacent neighbors [g kg :sup:`-1`]

    Returns
    -------
    delta_SA : array_like
               Absolute Salinity anomaly of the 4 adjacent neighbours
               [g kg :sup:`-1`]

    Notes
    -----
    originally inside "_delta_SA"

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    dsa_mean = dsa.mean(axis = 0)
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

if __name__ == '__main__':
    import doctest
    doctest.testmod()
