#TODO: create a CTD class that will take Conductive, in situ t and Pressure (lon, lat) and output SA
#TODO: check_dim for p in all "p" functions
#FIXME: some function return values even with NaN in the input, check this behavior (also present in the original). print doctest to check which one

from __future__ import division

import numpy as np
from seawater import constants as cte
import os

""" The global data set of Absolute Salinity Anomaly:
data/gsw_data_v2_0.pkl
TODO: Add a option for a local Atlas
"""
try:
    import cPickle as pickle
except:
    import pickle

datadir = os.sep.join([os.path.dirname(__file__), 'data/'])

"""
Section A: Library functions
"""

def _gibbs(ns, nt, npr, SA, t, p):
    r"""
    Calculates specific Gibbs energy and its derivatives up to order 2 for seawater.

    The Gibbs function approach allows the calculation of internal energy, entropy, enthalpy, potential enthalpy and the chemical potentials of seawater as well as the freezing temperature, and the latent heats of freezing and of evaporation. These quantities were not available from EOS-80 but are essential for the accurate accounting of heat in the ocean and for the consistent and accurate treatment of air-sea and ice-sea heat fluxes.

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
            The Gibbs energy (when ns = nt = npr = 0) has units of: [ J kg :sup:`-1` ]
            The Absolute Salinity derivatives are output in units of: [ (J kg :sup:`-1`) (g kg :sup:`-1`) :sup:`-ns` ]
            The temperature derivatives are output in units of: [ (J kg :sup:`-1`) K :sup:`-nt` ]
            The pressure derivatives are output in units of: [ (J kg :sup:`-1`) Pa :sup:`-npr` ]
            The mixed derivatives are output in units of: [ (J kg :sup:`-1`) (g kg :sup:`-1`) :sup:`-ns` K :sup:`-nt` Pa :sup:`-npr` ]

    Notes
    -----
    The Gibbs function for seawater is that of TEOS-10 (IOC et al., 2010), being the sum of IAPWS-08 for the saline part and IAPWS-09 for the pure water part. These IAPWS releases are the officially blessed IAPWS descriptions of Feistel (2008) and the pure water part of Feistel (2003). Absolute Salinity, SA, in all of the GSW routines is expressed on the Reference-Composition Salinity Scale of 2008 (RCSS-08) of Millero et al. (2008).

    The derivatives are taken with respect to pressure in Pa, not withstanding that the pressure input into this routine is in dbar.

    References
    ----------
    .. [1] Feistel, R., 2003: A new extended Gibbs thermodynamic potential of seawater, Progr. Oceanogr., 58, 43-114.

    .. [2] Feistel, R., 2008: A Gibbs function for seawater thermodynamics for -6 to 80 :math:`^\circ` C and salinity up to 120 g kg :sup:`-1`, Deep-Sea Res. I, 55, 1639-1671.

    .. [3] IAPWS, 2008: Release on the IAPWS Formulation 2008 for the Thermodynamic Properties of Seawater. The International Association for the Properties of Water and Steam. Berlin, Germany, September 2008, available from http://www.iapws.org.  This Release is referred to as IAPWS-08.

    .. [4] IAPWS, 2009: Supplementary Release on a Computationally Efficient Thermodynamic Formulation for Liquid Water for Oceanographic Use. The International Association for the Properties of Water and Steam. Doorwerth, The Netherlands, September 2009, available from http://www.iapws.org.  This Release is referred to as IAPWS-09.

    .. [5] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org See section 2.6 and appendices A.6,  G and H.

    .. [6] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: The composition of Standard Seawater and the definition of the Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.

    Modifications:
    2010-09-24. David Jackett, Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    # trick to use single number
    if SA.ndim == 0:
        SA = SA[np.newaxis]

    SA[SA < 0] = 0 # ensure that SA is non-negative.

    x2 = cte.sfac * SA
    x = np.sqrt(x2)
    y = t * 0.025
    z = p * 1e-4 # The input pressure (p) is sea pressure in units of dbar.

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

        gibbs = 0.5 * cte.sfac * g08

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

        gibbs = (g03 + g08) * 0.025

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

        # Pressure derivative of the Gibbs function in units of (J kg :sup:`-1`) (Pa :sup:`-1`) = m :sup:`3` kg :sup:`-1`
        gibbs = (g03 + g08) * 1e-8

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
        gibbs = 0.5 * cte.sfac * 0.025 * g08

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

        # Derivative of the Gibbs function is in units of (m :sup:`3` kg :sup:`-1`) / (g kg :sup:`-1`) = m :sup:`3` g :sup:`-1` that is, it is the derivative of specific volume with respect to Absolute Salinity measured in g kg :sup:`-1`.
        gibbs = g08 * cte.sfac * 0.5e-8

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

        # Derivative of the Gibbs function is in units of (m :sup:`3` (K kg) ) that is, the pressure of the derivative in Pa.
        gibbs = (g03 + g08) * 2.5e-10

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

        gibbs = 0.25 * cte.sfac**2 * g08

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

        # Second derivative of the Gibbs function with respect to pressure, measured in Pa; units of (J kg :sup:`-1`) (Pa :sup:`-2`).
        gibbs = (g03 + g08) * 1e-16
    else:
        raise NameError('Wrong Combination of order/variables')

    return gibbs

def _entropy_part(SA, t, p):
    r"""
    Calculates entropy, except that it does not evaluate any terms that are functions of Absolute Salinity alone.

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
                   entropy minus the terms that due to SA alone [J kg :sup:`-1` K :sup:`-1`]

    Notes
    -----
    By not calculating these terms, which are a function only of Absolute Salinity, several unnecessary computations are avoided (including saving the computation of a natural logarithm). These terms are a necessary part of entropy, but are not needed when calculating potential temperature from in situ temperature.

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    if SA.shape:
        SA[SA < 0] = 0
    elif SA < 0:
        SA = 0

    x2 = cte.sfac * SA
    x = np.sqrt(x2)
    y = t * 0.025
    z = p * 1e-4

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

    entropy_part = -( g03 + g08 )  * 0.025

    return entropy_part

def _gibbs_pt0_pt0(SA, pt0):
    r"""
    Calculates the second derivative of the specific Gibbs function with respect to temperature at zero sea pressure or _gibbs(0,2,0,SA,t,0).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt0 : array_like
          potential temperature relative to 0 dbar [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    gibbs_pt0_pt0 : array_like
                    TODO: write the eq. for the second derivative of the specific Gibbs function.
                    FIXME: [units]

    Notes
    -----
    This library function is called by both "pt_from_CT(SA,CT)" and "pt0_from_t(SA,t,p)".

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, pt0 = np.asanyarray(SA), np.asanyarray(pt0)

    # Ensure that SA is non-negative. Allows for array and single number
    if SA.shape:
        SA[SA < 0] = 0
    elif SA < 0:
        SA = 0

    x2 = cte.sfac * SA
    x = np.sqrt(x2)
    y = pt0 * 0.025

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

def _entropy_part_zerop(SA, pt0):
    r"""
    Calculates entropy at a sea surface (p = 0 dbar), except that it does not evaluate any terms that are functions of Absolute Salinity alone.

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
    By not calculating these terms, which are a function only of Absolute Salinity, several unnecessary computations are avoided (including saving the computation of a natural logarithm). These terms are a necessary part of entropy, but are not needed when calculating potential temperature from in situ temperature.

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, pt0 = np.asanyarray(SA), np.asanyarray(pt0)

    # Ensure that SA is non-negative
    if SA.shape:
        SA[SA < 0] = 0
    elif SA < 0:
        SA = 0

    x2 = cte.sfac * SA
    x = np.sqrt(x2)
    y = pt0 * 0.025

    g03 = y * ( -24715.571866078 + y * ( 2210.2236124548363 + \
    y * ( -592.743745734632 + y * ( 290.12956292128547 + \
    y * ( -113.90630790850321 + y * 21.35571525415769) ) ) ) )

    g08 = x2 * ( x * ( x * ( y * ( -137.1145018408982 + y * ( 148.10030845687618 + \
    y * ( -68.5590309679152 + 12.4848504784754 * y ) ) ) ) + \
    y * ( -86.1329351956084 + y * ( -30.0682112585625 + y * 3.50240264723578 ) ) ) + \
    y * ( 1760.062705994408 + y * ( -675.802947790203 + \
    y * ( 365.7041791005036 + y * ( -108.30162043765552 + 12.78101825083098 * y ) ) ) ) )

    entropy_part_zerop = -( g03 + g08 ) * 0.025

    return entropy_part_zerop

def _enthalpy_SSO_0_CT25(p):
    r"""
     Calculates enthalpy at the Standard Ocean Salinity (SSO) and at a Conservative Temperature of zero degrees C (CT=0), as a function of pressure (p [dbar]) or enthalpy_CT25(35.16504,0,p).

    Parameters
    ----------
    p : array_like
        pressure [dbar]

    Returns
    -------
    enthalpy_CT25 : array_like
                    enthalpy_CT25 at (SSO, CT = 0, p), 25-term equation. [J kg :sup:`-1`]

    Notes
    -----
    It Uses a streamlined version of the 25-term CT version of the Gibbs function, that is, a streamlined version of the code "enthalpy_CT25(SA,CT,p)"

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    p = np.asanyarray(p)

    SSO = cte.SSO * np.ones( p.shape )

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

def _specvol_SSO_0_CT25(p):
    r"""
    Calculates specific volume at the Standard Ocean Salinity (SSO) and Conservative Temperature of zero degrees C (CT=0), as a function of pressure (p [dbar]) or spec_vol_CT25(35.16504,0,p).

    Parameters
    ----------
    p : array_like
        pressure [dbar]

    Returns
    -------
    specvol_SSO_0_CT25 : array_like
                         Specific volume at (SSO, CT=0, p), 25-term equation. [m :sup:`3` kg :sup:`-1`]

    Notes
    -----
    It uses a streamlined version of the 25-term CT version of specific volume that is, a streamlined version of the code "rho_alpha_beta_CT25(SA,CT,p)"

    Modifications
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    p = np.asanyarray(p)

    SSO = cte.SSO * np.ones( p.shape )
    specvol_SSO_0_CT25 = (1.00000000e+00 + SSO * ( 2.0777716085618458e-003 +np.sqrt(SSO) * 3.4688210757917340e-006) + p * 6.8314629554123324e-006) / (9.9984380290708214e+002 + SSO * ( 2.8925731541277653e+000 + SSO * 1.9457531751183059e-003) + p * ( 1.1930681818531748e-002 + SSO * 5.9355685925035653e-006 + p * -2.5943389807429039e-008) )

    return specvol_SSO_0_CT25

def _check_dim(prop1, prop2):
    r"""
    Broadcast prop1 to the shape of prop2.
    Prop1 can be scalar, row equal or column equal to prop2.
    """

    # Comic book guy would say: "Worst function ever!"
    if prop1.ndim == 1:
        prop1 = prop1.flatten()

    if (prop1.ndim == 1) & (prop1.size == 1):
        prop1 = prop1 * np.ones( prop2.shape )
    elif (prop1.ndim == 1) & (prop2.ndim != 1):
        if prop1.size == prop2.shape[1]:
            prop1 = prop1 * np.ones(prop2.shape)
        elif prop1.size == prop2.shape[0]:
            prop1 = prop1[:,np.newaxis] * np.ones(prop2.shape)
        else:
            raise NameError('Blahrg')

    if prop1.ndim == 0:
        prop1 = prop1 * np.ones(prop2.shape)

    return prop1

"""
Salinity lib functions
"""
def _SP_from_SA_Baltic(SA, lon, lat):
    r"""
    Calculates Practical Salinity (SP) for the Baltic Sea, from a value computed analytically from Absolute Salinity.

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
    >>> import seawater.gibbs as gsw
    >>> SA = [6.6699, 6.7738, 6.9130, 7.3661, 7.5862, 10.3895]
    >>> lon = 20
    >>> lat = 59
    >>> gsw._SP_from_SA_Baltic(SA, lon, lat)
    array([  6.56825467,   6.67192352,   6.81081383,   7.2629058 ,
             7.48251613,  10.27957947])

    References
    ----------
    .. [1] Feistel, R., S. Weinreben, H. Wolf, S. Seitz, P. Spitzer, B. Adel, G. Nausch, B. Schneider and D. G. Wright, 2010c: Density and Absolute Salinity of the Baltic Sea 2006-2009.  Ocean Science, 6, 3-24.
    http://www.ocean-sci.net/6/3/2010/os-6-3-2010.pdf

    .. [2] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.
    Available from http://www.TEOS-10.org

    .. [3] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, lon, lat = np.asanyarray(SA), np.asanyarray(lon), np.asanyarray(lat)
    lon, lat = _check_dim(lon, SA), _check_dim(lat, SA)

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
            SP_baltic[inds[inds1]] = ( 35 / ( cte.SSO - 0.087 ) ) * ( SA[inds[inds1]] - 0.087)

        SP_baltic = np.reshape( SP_baltic, lon.shape )

    return SP_baltic

def _SA_from_SP_Baltic(SP, lon, lat):
    r"""
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
                Absolute salinity [g kg :sup:`-1`]

    See Also
    --------
    SA_from_SP, Sstar_from_SP, SA_Sstar_from_SP

    Notes
    -----
    This program will only produce Absolute Salinity values for the Baltic Sea.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SP = [6.5683, 6.6719, 6.8108, 7.2629, 7.4825, 10.2796]
    >>> lon = 20
    >>> lat = 59
    >>> gsw._SA_from_SP_Baltic(SP, lon, lat)
    array([  6.66994543,   6.77377643,   6.91298614,   7.36609419,
             7.58618384,  10.38952057])

    References
    ----------
    .. [1] Feistel, R., S. Weinreben, H. Wolf, S. Seitz, P. Spitzer, B. Adel, G. Nausch, B. Schneider and D. G. Wright, 2010c: Density and Absolute Salinity of the Baltic Sea 2006-2009.  Ocean Science, 6, 3-24.
    http://www.ocean-sci.net/6/3/2010/os-6-3-2010.pdf

    .. [2] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    .. [3] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SP, lon, lat = np.asanyarray(SP), np.asanyarray(lon), np.asanyarray(lat)
    lon, lat = _check_dim(lon, SP), _check_dim(lat, SP)

    xb1, xb2, xb3 = 12.6, 7., 26.
    xb1a, xb3a = 45., 26.
    yb1, yb2, yb3 = 50., 59., 69.

    inds_baltic = (xb2 < lon) & (lon < xb1a) & (yb1 < lat) & (lat < yb3)
    SA_baltic = np.ones( SP.shape ) * np.nan

    if list(inds_baltic):
        xx_left = np.interp( lat[inds_baltic], [yb1,yb2,yb3], [xb1,xb2,xb3])
        xx_right = np.interp( lat[inds_baltic], [yb1,yb3], [xb1a,xb3a] )
        inds_baltic1 = (xx_left <= lon[inds_baltic]) & (lon[inds_baltic] <= xx_right)
        SA_baltic[inds_baltic[inds_baltic1]] = ( ( cte.SSO - 0.087 ) / 35 ) * SP[inds_baltic[inds_baltic1]] + 0.087

    return SA_baltic

def _delta_SA(p, lon, lat):
    r"""
    Calculates the Absolute Salinity anomaly, SA - SR, in the open ocean by spatially interpolating the global reference data set of delta_SA to the location of the seawater sample.

    Parameters
    ----------
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    delta_SA : array_like
               Absolute Salinity anomaly [g kg :sup:`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    _dsa_add_barrier, _dsa_add_mean

    Notes
    -----
    The Absolute Salinity Anomaly in the Baltic Sea is evaluated separately, since it is a function of Practical Salinity, not of space. The present function returns a delta_SA of zero for data in the Baltic Sea. The correct way of calculating Absolute Salinity in the Baltic Sea is by calling _SA_from_SP.

    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometers inland from the coast.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat = 188, 4
    >>> gsw._delta_SA(p, lon, lat)
    (array([ 0.00016779,  0.00026868,  0.00066554,  0.0026943 ,  0.00562666,
            0.00939665]), array([ True,  True,  True,  True,  True,  True], dtype=bool))
    >>> gsw._delta_SA(p, -80, 45)[1]
    array([ True,  True,  True,  True,  True, True], dtype=bool)

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean.  Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    ????-??-??. David Jackett.
    2010-07-23. Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    p, lon, lat = np.asanyarray(p), np.asanyarray(lon), np.asanyarray(lat)
    lon, lat = _check_dim(lon, p), _check_dim(lat, p)

    if (lon < 0).any():
        lon[lon < 0] = lon[lon < 0] + 360.

    data = pickle.load( open(os.path.join(datadir + 'gsw_data_v2_0.pkl'), 'rb') )

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
    indsz0 = np.sum( (P_REF >= P), axis=0 ) -1 #FIXME: -1(?)

    nmax = np.c_[ ndepth_ref[indsy0, indsx0], \
                  ndepth_ref[indsy0, indsx0+1], \
                  ndepth_ref[indsy0+1, indsx0+1], \
                  ndepth_ref[indsy0+1, indsx0] ]
    nmax = np.int64( np.nanmax(nmax, axis=1) )

    inds1 = np.where(indsz0 > nmax)[0] # casts deeper than GK maximum
    if inds1.size != 0:
        # have reset p here to reset indsz0
        p[inds1] = p_ref[nmax[inds1]-1]

    #FIXME: Ugly matlab matrix dot multiplication, there must be a better way...
    P_REF = np.dot( np.ones(p_ref.size)[:,np.newaxis], p[np.newaxis,:] )
    P = np.dot( p_ref[:,np.newaxis], np.ones(p.size)[np.newaxis,:] )
    indsz0 = np.sum( (P_REF >= P), axis=0 ) -1

    inds = (indsz0 == p_ref.size-1)
    indsz0[inds] = p_ref.size - 2

    inds0 = indsz0 + indsy0 * delta_SA_ref.shape[0] + indsx0 * delta_SA_ref.shape[0] * delta_SA_ref.shape[1]

    data_indices = np.c_[indsx0, indsy0, indsz0, inds0]
    data_inds = data_indices[:,2]

    r1 = ( lon - longs_ref[indsx0] ) / ( longs_ref[indsx0+1] - longs_ref[indsx0] )
    s1 = ( lat - lats_ref[indsy0] ) / ( lats_ref[indsy0+1] - lats_ref[indsy0] )
    t1 = np.float64( ( p - p_ref[indsz0] ) ) / np.float64( ( p_ref[indsz0+1] - p_ref[indsz0] ) )

    nksum = 0
    no_levels_missing = 0

    sa_upper = np.nan * ( np.ones(data_inds.shape) )
    sa_lower = np.nan * ( np.ones(data_inds.shape) )
    delta_SA = np.nan * ( np.ones(data_inds.shape) )

    #for k in range(0, p_ref.size-1):
    for k in range(0, p_ref.size):
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
            dsa[1, inds_k] = delta_SA_ref[indsz, indsy, indsx+1]   # inds0 + ny*nz
            dsa[2, inds_k] = delta_SA_ref[indsz, indsy+1, indsx+1] # inds0 + ny*nz + nz
            dsa[3, inds_k] = delta_SA_ref[indsz, indsy+1, indsx] #  inds0 + nz

            #inds = np.where( (260. <= lon) & (lon <= 295.217) & (0. <= lat) & (lat <= 19.55) & (indsz0 == k) )[0]
            inds = ( (260. <= lon) & (lon <= 295.217) & (0. <= lat) & (lat <= 19.55) & (indsz0 == k) )
            """ TODO: describe add_barrier """
            if inds.any():
                dsa[:,inds] = _dsa_add_barrier( dsa[:,inds], lon[inds], lat[inds], longs_ref[indsx0[inds]], lats_ref[indsy0[inds]], dlongs_ref, dlats_ref)
                #print "add_barrier %s" % dsa

            inds = np.where( ( np.isnan( np.sum(dsa, axis=0) ) ) & (indsz0==k))[0]
            #print "inds = %s" % inds
            #raw_input("Press Enter to continue...")
            """ TODO: describe add_mean """
            if inds.size !=0:
                #print "dsa = %s" % dsa[:,inds]
                dsa[:,inds] = _dsa_add_mean(dsa[:,inds])
                #print "add_mean %s" % dsa

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
            if inds.any():
                dsa[:,inds] = _dsa_add_barrier( dsa[:,inds], lon[inds], lat[inds], longs_ref[indsx0[inds]], lats_ref[indsy0[inds]], dlongs_ref, dlats_ref)
                #print "add_barrier %s" % dsa

            inds = ( np.isnan( np.sum(dsa, axis=0) ) ) & (indsz0==k)

            """ TODO: describe add_mean """
            dsa[:,inds] = _dsa_add_mean(dsa[:,inds])
            #print "add_mean %s" % dsa

            sa_lower[inds_di] = ( 1 - s1[inds_di] ) * ( dsa[0, inds_k] + \
            r1[inds_di] * ( dsa[1, inds_k] - dsa[0,inds_k] ) ) + \
            s1[inds_di] * ( dsa[3, inds_k] + \
            r1[inds_di] * ( dsa[2, inds_k] - dsa[3, inds_k] ) )

            inds_different = np.isfinite(sa_upper[inds_di]) & np.isnan(sa_lower[inds_di])
            sa_lower[inds_di[inds_different]] = sa_upper[inds_di[inds_different]]

            delta_SA[inds_di] = sa_upper[inds_di] + t1[inds_di] * ( sa_lower[inds_di] - sa_upper[inds_di] )
        else:
            no_levels_missing = no_levels_missing + 1

    in_ocean = ~np.isfinite(delta_SA)
    delta_SA[in_ocean] = 0
    in_ocean = (in_ocean == False) # reverse True <-> False
    return delta_SA, in_ocean

def _dsa_add_barrier(dsa, lon, lat, longs_ref, lats_ref, dlongs_ref, dlats_ref):
    r"""
    Adds a barrier through Central America (Panama) and then averages over the appropriate side of the barrier.

    Parameters
    ----------
    dsa : array_like
          Absolute Salinity anomaly of the 4 adjacent neighbors  [g kg :sup:`-1`]
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
    dsa = np.asanyarray(dsa)
    lon, lat = np.asanyarray(lon), np.asanyarray(lat)
    longs_ref, lats_ref = np.asanyarray(longs_ref), np.asanyarray(lats_ref)
    dlongs_ref, dlats_ref = np.asanyarray(dlongs_ref), np.asanyarray(dlats_ref)

    longs_pan = np.array([260.0000, 272.5900, 276.5000, 278.6500, 280.7300, 295.2170])
    lats_pan = np.array([19.5500, 13.9700, 9.6000, 8.1000, 9.3300, 0])

    lats_lines0 = np.interp(lon, longs_pan, lats_pan)

    lats_lines1 = np.interp( longs_ref, lats_pan, longs_pan)
    lats_lines2 = np.interp( (longs_ref+dlongs_ref), lats_pan, longs_pan)

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
          Absolute Salinity anomaly of the 4 adjacent neighbors  [g kg :sup:`-1`]

    Returns
    -------
    delta_SA : array_like
          Absolute Salinity anomaly of the 4 adjacent neighbours  [g kg :sup:`-1`]

    Notes
    -----
    originally inside "_delta_SA"

    Modifications:
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    dsa = np.asanyarray(dsa)

    #FIXME: should be nanmean here in the original
    # the best approach should be a complete re-write using masked array and elimination in_ocean
    #print "dsa%s" % dsa
    dsa_mean = dsa.mean(axis = 0)
    #print "dsa_mean %s" % dsa_mean
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

"""
Section B: functions
"""

def entropy(SA, t, p):
    r"""
    Calculates specific entropy of seawater.

    The specific entropy of seawater :math:`\eta` is given by:

    .. math::
        \eta(SA, t, p) = -g_T = \frac{\partial g}{\partial T}\Big|_{SA,p}

    When taking derivatives with respect to *in situ* temperature, the symbol :math:`T` will be used for temperature in order that these derivatives not be confused with time derivatives.

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
    entropy : array_like
              specific entropy [J kg :sup:`-1` K :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.entropy(SA, t, p)
    array([ 400.38942528,  395.43817843,  319.8664982 ,  146.79088159,
             98.64734087,   62.79150873])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1 = 0, 1

    entropy = -_gibbs(n0, n1, n0, SA, t, p)

    return entropy

def rho(SA, t, p):
    r"""
    Calculates in situ density of seawater from Absolute Salinity and in situ temperature.

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
    rho : array_like
          in situ density [kg m :sup:`-3`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.rho(SA, t, p)
    array([ 1021.84017319,  1022.26268993,  1024.42771594,  1027.79020181,
            1029.83771473,  1032.00240412])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.8.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1 = 0, 1

    rho = 1. / _gibbs(n0, n0, n1, SA, t, p)

    return rho

def cp(SA, t, p):
    r"""
    Calculates the isobaric heat capacity of seawater.

    SA : array_like
        Absolute salinity [g kg :sup:`-1`]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]

    Returns
    -------
    cp : array_like
         heat capacity of seawater [J kg :sup:`-1` K :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.cp(SA, t, p)
    array([ 4002.88800396,  4000.98028393,  3995.54646889,  3985.07676902,
            3973.59384348,  3960.18408479])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n2 = 0, 2

    cp = -( t + cte.Kelvin ) * _gibbs(n0, n2, n0, SA, t, p)

    return cp

def helmholtz_energy(SA, t, p):
    r"""
    Calculates the Helmholtz energy of seawater.


    The specific Helmholtz energy of seawater :math:`f` is given by:

    .. math::
        f(SA, t, p) = g - (p + P_0) \nu = g - (p + P_0) \frac{\partial g}{\partial P}\Big|_{SA,T}

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
    Helmholtz_energy : array_like
                       Helmholtz energy [J kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.helmholtz_energy(SA, t, p)
    array([-5985.58288209, -5830.81845224, -3806.96617841,  -877.66369421,
            -462.17033905,  -245.50407205])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.13.

    Modifications:
    2010-08-26. Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1 = 0, 1

    helmholtz_energy = _gibbs(n0, n0, n0, SA, t, p) - \
                    ( cte.db2Pascal * p + 101325 ) * _gibbs(n0, n0, n1, SA, t, p)

    return helmholtz_energy

def internal_energy(SA, t, p):
    r"""
    Calculates the Helmholtz energy of seawater.

    The specific internal energy of seawater :math:`u` is given by:

    .. math::
        u(SA, t, p) = g + (T_0 + t)\eta - (p + P_0)\nu = g - (T_0 + t)\frac{\partial g}{\partial T}\Big|_{SA,p} - (p + P_0)\frac{\partial g}{\partial P}\Big|_{SA,T}

    where :math:`T_0` is the Celsius zero point, 273.15 K and :math:`P_0` = 101 325 Pa is the standard atmosphere pressure.

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
    internal_energy (u) : array_like
                          specific internal energy [J kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.internal_energy(SA, t, p)
    array([ 114906.23847309,  113426.57417062,   90860.81858842,
             40724.34005719,   27162.66600185,   17182.50522667])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.11.1)

    Modifications:
    2010-08-22. Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1 = 0, 1

    internal_energy = _gibbs(n0, n0, n0, SA, t, p) - \
                    (cte.Kelvin + t) * _gibbs(n0, n1, n0, SA, t, p) - \
                    (cte.db2Pascal * p + 101325) * _gibbs(n0, n0, n1, SA, t, p)

    return internal_energy

def sound_speed(SA, t, p):
    r"""
    Calculates the speed of sound in seawater.


    The speed of sound in seawater :math:`c` is given by:

    .. math::
        c(SA, t, p) = \sqrt{ \partial P  / \partial \rho |_{SA,\eta}} = \sqrt{(\rho\kappa)^{-1}} = g_P \sqrt{g_{TT}/(g^2_{TP} - g_{TT}g_{PP})}

    Note that in these expressions, since sound speed is in m s :sup`-1` and density has units of kg m :sup:`-3` it follows that the pressure of the partial derivatives must be in Pa and the isentropic compressibility :math:`kappa` must have units of Pa :sup:`-1`. The sound speed c produced by both the SIA and the GSW software libraries (appendices M and N) has units of m s :sup:`-1`.

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
    sound_speed : array_like
                  speed of sound in seawater [m s :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.sound_speed(SA, t, p)
    array([ 1542.61580359,  1542.70353407,  1530.84497914,  1494.40999692,
            1487.37710252,  1483.93460908])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.17.1)

    Modifications:
    2010-07-23. David Jackett, Paul Barker and Trevor McDougall.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1, n2 = 0, 1, 2

    g_tt = _gibbs(n0, n2, n0, SA, t, p)
    g_tp = _gibbs(n0, n1, n1, SA, t, p)

    sound_speed = _gibbs(n0, n0, n1, SA, t, p) * \
    np.sqrt( g_tt / ( g_tp * g_tp - g_tt * _gibbs(n0, n0, n2, SA, t, p ) ) )

    return sound_speed

def adiabatic_lapse_rate(SA, t, p):
    r"""
    Calculates the adiabatic lapse rate of sea water.

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
    adiabatic_lapse_rate : array_like
                           Adiabatic lapse rate [K Pa :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    The output is in unit of degrees Celsius per Pa, (or equivalently K/Pa) not in units of K/dbar

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.adiabatic_lapse_rate(SA, t, p)
    array([  2.40350282e-08,   2.38496700e-08,   2.03479880e-08,
             1.19586543e-08,   9.96170718e-09,   8.71747270e-09])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.22.1).

    Modifications:
    2010-08-26. Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1, n2 = 0, 1, 2

    adiabatic_lapse_rate = - _gibbs(n0, n1, n1, SA, t, p) / ( _gibbs(n0, n2, n0, SA, t, p ) )

    return adiabatic_lapse_rate

def chem_potential_relative(SA, t, p):
    r"""
    Calculates the adiabatic lapse rate of sea water.

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
    chem_potential_relative : array_like
                              relative chemical potential [J kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.chem_potential_relative(SA, t, p)
    array([ 79.4254481 ,  79.25989214,  74.69154859,  65.64063719,
            61.22685656,  57.21298557])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-08-26. Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1 = 0, 1

    chem_potential_relative = _gibbs(n1, n0, n0, SA, t, p)

    return chem_potential_relative

def specvol(SA, t, p):
    r"""
    Calculates the specific volume of seawater.

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
    specvol : array_like
              specific volume [m :sup:`3` kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.specvol(SA, t, p)
    array([ 0.00097863,  0.00097822,  0.00097615,  0.00097296,  0.00097103,
            0.00096899])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.7.

    Modifications:
    2010-08-26. David Jackett & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1 = 0, 1

    specvol = _gibbs(n0, n0, n1, SA, t, p)

    return specvol

def conservative_t(SA, t, p):
    r"""
    Calculates Conservative Temperature of seawater from in situ temperature.

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
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.conservative_t(SA, t, p)
    array([ 28.80991983,  28.43922782,  22.78617689,  10.22618927,
             6.82721363,   4.32357575])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 3.3.

    Modifications:
    2010-08-26. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    #pt0 = potential_t(SA, t, p) #NOTE: original call a faster version (pt0_from_t) instead of potential_t
    pt0 = pt0_from_t(SA, t, p)

    CT = CT_from_pt(SA, pt0)

    return CT

def potential_t(SA, t, p, pr=0):
    r"""
    Calculates potential temperature with the general reference pressure, pr, from in situ temperature.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]
    pr : int, float, optional
         reference pressure, default = 0

    Returns
    -------
    pt : array_like
         potential temperature [:math:`^\circ` C (ITS-90)]

    See Also
    --------
    _entropy_part

    Notes
    -----
    This function calls "entropy_part" which evaluates entropy except for the parts which are a function of Absolute Salinity alone. A faster routine exists pt0_from_t(SA,t,p) if pr is indeed zero dbar.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.potential_t(SA, t, p)
    array([ 28.78319682,  28.42098334,  22.7849304 ,  10.23052366,
             6.82923022,   4.32451057])
    >>> gsw.potential_t(SA, t, p, pr = 1000)
    array([ 29.02665528,  28.662375  ,  22.99149634,  10.35341725,
             6.92732954,   4.4036    ])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 3.1.

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010: A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater. To be submitted to Ocean Science Discussions.

    Modifications:
    2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)
    pr = np.asanyarray(pr)

    n0, n2 = 0, 2

    SA[SA < 0] = 0

    s1 = SA * 35. / cte.SSO

    pt = t + ( p - pr ) * ( 8.65483913395442e-6  - \
    s1 * 1.41636299744881e-6 - \
    ( p + pr ) * 7.38286467135737e-9 + \
    t * ( -8.38241357039698e-6 + \
    s1 * 2.83933368585534e-8 + \
    t * 1.77803965218656e-8 + \
    ( p + pr ) * 1.71155619208233e-10 ) )

    dentropy_dt = cte.cp0 / ( (cte.Kelvin + pt) * ( 1 - 0.05 * ( 1 - SA / cte.SSO ) ) )

    true_entropy_part = _entropy_part(SA, t, p)

    for Number_of_iterations in range(0,2,1):
        pt_old = pt
        dentropy = _entropy_part(SA, pt_old, pr) - true_entropy_part
        pt = pt_old - dentropy / dentropy_dt # this is half way through the modified method
        ptm = 0.5 * (pt + pt_old)
        dentropy_dt = -_gibbs(n0, n2, n0, SA, ptm, pr)
        pt = pt_old - dentropy / dentropy_dt

    # maximum error of 6.3x10^-9 degrees C for one iteration.
    # maximum error is 1.8x10^-14 degrees C for two iterations (two iterations is the default, "for Number_of_iterations = 1:2).
    # These errors are over the full "oceanographic funnel" of McDougall et al. (2010), which reaches down to p = 8000 dbar.

    return pt

def enthalpy(SA, t, p):
    r"""
    Calculates the specific enthalpy of seawater.


    The specific enthalpy of seawater :math:`h` is given by:

    .. math::
        h(SA, t, p) = g + (T_0 + t)\eta = g - (T_0 + t) \frac{\partial g}{\partial T}\Big|_{SA,p}

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
    enthalpy : array_like
               specific enthalpy [J kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.enthalpy(SA, t, p)
    array([ 115103.26047838,  114014.8036012 ,   92179.9209311 ,
             43255.32838089,   33087.21597002,   26970.5880448 ])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See appendix A.11.

    Modifications:
    2010-08-26. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1 = 0, 1

    enthalpy = _gibbs(n0, n0, n0, SA, t, p) - ( t + cte.Kelvin ) * _gibbs(n0, n1, n0, SA, t, p)

    return enthalpy

def chem_potential_water(SA, t, p):
    r"""
    Calculates the chemical potential of water in seawater.

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
    chem_potential_water : array_like
                           chemical potential of water in seawater [J kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.chem_potential_water(SA, t, p)
    array([-8545.56114628, -8008.08554834, -5103.98013987,  -634.06778275,
            3335.56680347,  7555.43444597])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-09-28. Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    SA[ SA < 0] = 0 # ensure that SA is non-negative

    x2 = cte.sfac * SA

    x = np.sqrt(x2)
    y = t * 0.025
    z = p * 1e-4 # Note that the input pressure (p) is sea pressure in units of dbar

    #TODO: Check why this doesn't call _gibbs
    g03_g = 101.342743139674 + z * ( 100015.695367145 + \
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

    g08_g = x2 * ( 1416.27648484197 + \
        x * ( -2432.14662381794 + x * ( 2025.80115603697 + \
        y * ( 543.835333000098 + y * ( -68.5572509204491 + \
        y * ( 49.3667694856254 + y * ( -17.1397577419788 + 2.49697009569508 * y ) ) ) - 22.6683558512829 * z ) + \
        x * ( -1091.66841042967 - 196.028306689776 * y + \
        x * ( 374.60123787784 - 48.5891069025409 * x + 36.7571622995805 * y ) + 36.0284195611086 * z ) + \
        z * ( -54.7919133532887 + ( -4.08193978912261 - 30.1755111971161 * z ) * z ) ) + \
        z * ( 199.459603073901 + z * ( -52.2940909281335 + ( 68.0444942726459 - 3.41251932441282 * z ) * z ) ) + \
        y * ( -493.407510141682 + z * ( -175.292041186547 + ( 83.1923927801819 - 29.483064349429 * z ) * z ) + \
        y * ( -43.0664675978042 + z * ( 383.058066002476 + z * ( -54.1917262517112 + 25.6398487389914 * z ) ) + \
        y * ( -10.0227370861875 - 460.319931801257 * z + y * ( 0.875600661808945 + 234.565187611355 * z ) ) ) ) ) + \
        y * ( 168.072408311545 ) )

    g_SA_part = 8645.36753595126 + \
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
        y * ( 1187.3715515697959)

    chem_potential_water =  g03_g + g08_g  - 0.5 * cte.sfac * SA * g_SA_part

    return chem_potential_water

def chem_potential_salt(SA, t, p):
    r"""
    Calculates the chemical potential of salt in seawater.

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
    chem_potential_salt : array_like
                          chemical potential of salt in seawater [J kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.chem_potential_salt(SA, t, p)
    array([-8466.13569818, -7928.8256562 , -5029.28859129,  -568.42714556,
            3396.79366004,  7612.64743154])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.9.

    Modifications:
    2010-09-28. Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    chem_potential_salt = chem_potential_relative(SA, t, p) + \
                          chem_potential_water(SA, t, p)

    return chem_potential_salt

def isochoric_heat_cap(SA, t, p):
    r"""
    Calculates the isochoric heat capacity of seawater.

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
    isochoric_heat_cap : array_like
                         isochoric heat capacity [J kg :sup:`-1` K :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.isochoric_heat_cap(SA, t, p)
    array([ 3928.13708702,  3927.27381633,  3941.36418525,  3966.26126146,
            3960.50903222,  3950.13901342])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.21.

    Modifications:
    2010-08-26. Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1, n2 = 0, 1, 2

    g_tt = _gibbs(n0, n2, n0, SA, t, p)
    g_tp = _gibbs(n0, n1, n1, SA, t, p)
    g_pp = _gibbs(n0, n0, n2, SA, t, p)

    isochoric_heat_cap = -(cte.Kelvin + t) * (g_tt - g_tp * g_tp / g_pp)

    return isochoric_heat_cap

def kappa(SA, t, p):
    r"""
    Calculates the isentropic compressibility of seawater.

    When the entropy and Absolute Salinity are held constant while the pressure is changed, the isentropic and isohaline compressibility :math:`kappa` is obtained:

    .. math::
        \kappa(SA, t, p) = \rho^{-1}\frac{\partial \rho}{\partial P}\Big|_{SA,\eta} = -\nu^{-1}\frac{\partial \nu}{\partial P}\Big|_{SA,\eta} = \rho^{-1}\frac{\partial \rho}{\partial P}\Big|_{SA,\theta} = -\nu^{-1}\frac{\partial \nu}{\partial P}\Big|_{SA,\theta} =
        -\frac{ (g_{TP}^2 - g_{TT} g_{PP} ) }{g_P g_{TT}}

    The isentropic and isohaline compressibility is sometimes called simply the isentropic compressibility (or sometimes the "adiabatic compressibility"), on the unstated understanding that there is also no transfer of salt during the isentropic or adiabatic change in pressure.

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
    kappa : array_like
            Isentropic compressibility [Pa :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    The output is Pascal and not dbar.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.kappa(SA, t, p)
    array([  4.11245799e-10,   4.11029072e-10,   4.16539558e-10,
             4.35668338e-10,   4.38923693e-10,   4.40037576e-10])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqns. (2.16.1) and the row for kappa in Table P.1 of appendix P

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1, n2 = 0, 1, 2

    g_tt = _gibbs(n0, n2, n0, SA, t, p)
    g_tp = _gibbs(n0, n1, n1, SA, t, p)

    kappa = ( g_tp * g_tp - g_tt * _gibbs(n0, n0, n2, SA, t, p) ) / ( _gibbs(n0, n0, n1, SA, t, p ) * g_tt)

    return kappa

def kappa_const_t(SA, t, p):
    r"""
    Calculates isothermal compressibility of seawater at constant in situ temperature.

    .. math::
        \kappa^t(SA, t, p) = \rho^{-1}\frac{\partial \rho}{\partial P}\Big|_{SA,T} = -\nu^{-1}\frac{\partial \nu}{\partial P}\Big|_{SA,T} = -\frac{g_{PP}}{g_P}

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
    kappa : array_like
            Isothermal compressibility [Pa :sup:`-1`]

    See Also
    --------
    TODO

    Notes
    -----
    This is the compressibility of seawater at constant in situ temperature.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.kappa_const_t(SA, t, p)
    array([  4.19071646e-10,   4.18743202e-10,   4.22265764e-10,
             4.37735100e-10,   4.40373818e-10,   4.41156577e-10])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.15.1)

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1, n2 = 0, 1, 2

    kappa = - _gibbs(n0, n0, n2, SA, t, p) / _gibbs(n0, n0, n1, SA, t, p)

    return kappa

def pot_rho(SA, t, p, pr=0):
    r"""
    Calculates potential density of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]
    pr : int, float, optional
         reference pressure, default = 0

    Returns
    -------
    pot_rho : array_like
              potential density  [kg m :sup:`-3`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.pot_rho(SA, t, p)
    array([ 1021.79814581,  1022.05248442,  1023.89358365,  1026.66762112,
            1027.10723087,  1027.40963126])
    >>> gsw.pot_rho(SA, t, p, pr=1000)
    array([ 1025.95554512,  1026.21306986,  1028.12563226,  1031.1204547 ,
            1031.63768355,  1032.00240412])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 3.4.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    pt = potential_t(SA, t, p, pr=pr)

    pot_rho = rho(SA, pt, pr)

    return pot_rho

def specvol_anom(SA, t, p):
    r"""
    Calculates specific volume anomaly from Absolute Salinity, in situ temperature and pressure, using the full TEOS-10 Gibbs function.

    The reference value of Absolute Salinity is SSO and the reference value of Conservative Temperature is equal to 0 :math:`^\circ` C.

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
    specvol_anom : array_like
                   specific volume anomaly  [m :sup:`3` kg :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.specvol_anom(SA, t, p)
    array([  6.01044463e-06,   5.78602432e-06,   4.05564999e-06,
             1.42198662e-06,   1.04351837e-06,   7.63964850e-07])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (3.7.3)

    Modifications:
    2010-08-26. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1 = 0, 1

    SSO = cte.SSO * np.ones( SA.shape )

    CT0 = np.zeros( SSO.shape )

    pt_zero = pt_from_CT(SSO, CT0)

    t_zero = potential_t(SSO, pt_zero, 0, p)

    specvol_anom = _gibbs(n0, n0, n1, SA, t, p) - \
    _gibbs(n0, n0, n1, SSO, t_zero, p)

    return specvol_anom

def alpha_wrt_t(SA, t, p):
    r"""
    Calculates the thermal expansion coefficient of seawater with respect to in situ temperature.

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
    alpha_wrt_t : array_like
                  thermal expansion coefficient [K :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.alpha_wrt_t(SA, t, p)
    array([ 0.0003256 ,  0.00032345,  0.00028141,  0.00017283,  0.00014557,
            0.00012836])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.18.1)

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1 = 0, 1

    alpha_wrt_t = _gibbs(n0, n1, n1, SA, t, p) / _gibbs(n0, n0, n1, SA, t, p)

    return alpha_wrt_t

def alpha_wrt_CT(SA, t, p):
    r"""
    Calculates the thermal expansion coefficient of seawater with respect to Conservative Temperature.

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
    alpha_wrt_CT : array_like
                   thermal expansion coefficient [K :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.alpha_wrt_CT(SA, t, p)
    array([ 0.00032471,  0.00032272,  0.00028118,  0.00017314,  0.00014627,
            0.00012943])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.18.3).

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1, n2 = 0, 1, 2

    #pt0 = potential_t(SA, t, p) #NOTE: original call a faster version "pt0_from_t" instead of potential_t
    pt0 = pt0_from_t(SA, t, p)

    factor = -cte.cp0 / ( (cte.Kelvin + pt0) * _gibbs(n0, n2, n0, SA, t, p ) )

    alpha_wrt_CT = factor  * ( _gibbs(n0, n1, n1, SA, t, p) / _gibbs(n0, n0, n1, SA, t, p ) )

    return alpha_wrt_CT

def alpha_wrt_pt(SA, t, p):
    r"""
    Calculates the thermal expansion coefficient of seawater with respect to potential temperature, with a reference pressure of zero.

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
    alpha_wrt_pt : array_like
                   thermal expansion coefficient [K :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.alpha_wrt_pt(SA, t, p)
    array([ 0.00032562,  0.00032355,  0.00028164,  0.00017314,  0.00014623,
            0.00012936])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.18.2).

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1, n2 = 0, 1, 2

    #pt0 = potential_t(SA, t, p) #NOTE: original call a faster version (pt0_from_t) instead of potential_t
    pt0 = pt0_from_t(SA, t, p)

    factor = _gibbs(n0, n2, n0, SA, pt0, 0) / _gibbs(n0, n2, n0, SA, t, p)

    alpha_wrt_pt = factor * ( _gibbs(n0, n1, n1, SA, t, p) / _gibbs( n0, n0, n1, SA, t, p ) )

    return alpha_wrt_pt

def beta_const_t(SA, t, p):
    r"""
    Calculates the saline (i.e. haline) contraction coefficient of seawater at constant in situ temperature.

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
    beta_const_t : array_like
                   saline contraction coefficient [kg g :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.beta_const_t(SA, t, p)
    array([ 0.00073112,  0.00073107,  0.00073602,  0.00075381,  0.00075726,
            0.00075865])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.19.1)

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1 = 0, 1

    beta_const_t = -_gibbs(n1, n0, n1, SA, t, p) / _gibbs(n0, n0, n1, SA, t, p)

    return beta_const_t

def beta_const_CT(SA, t, p):
    r"""
    Calculates the saline (i.e. haline) contraction coefficient of seawater at constant Conservative Temperature.

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
    beta_const_CT : array_like
                    saline contraction coefficient [kg g :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.beta_const_CT(SA, t, p)
    array([ 0.00071749,  0.00071765,  0.00072622,  0.00075051,  0.00075506,
            0.00075707])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.19.3)

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1, n2 = 0, 1, 2

    #pt0 = potential_t(SA, t, p) #NOTE: original call a faster version (pt0_from_t) instead of potential_t
    pt0 = pt0_from_t(SA, t, p)

    gp = _gibbs(n0, n0, n1, SA, t, p)

    factora = _gibbs(n1, n1, n0, SA, t, p) - \
    _gibbs(n1, n0, n0, SA, pt0, 0) / (cte.Kelvin + pt0)

    factor = factora / ( gp * _gibbs(n0, n2, n0, SA, t, p) )

    beta_const_CT = _gibbs(n0, n1, n1, SA, t, p) * factor - \
    _gibbs(n1, n0, n1, SA, t, p) / gp

    return beta_const_CT

def beta_const_pt(SA, t, p):
    r"""
    Calculates the saline (i.e. haline) contraction coefficient of seawater at constant potential temperature with a reference pressure of 0 dbar.

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
    beta_const_pt : array_like
                    saline contraction coefficient [kg g :sup:`-1`]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.beta_const_pt(SA, t, p)
    array([ 0.00073112,  0.00073106,  0.00073599,  0.00075375,  0.00075712,
            0.00075843])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Eqn. (2.19.2)

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0, n1, n2 = 0, 1, 2

    #pt0 = potential_t(SA, t, p) #NOTE: original call a faster version "pt0_from_t" instead of potential_t
    pt0 = pt0_from_t(SA, t, p)

    gp = _gibbs(n0, n0, n1, SA, t, p)

    factora = _gibbs(n1, n1, n0, SA, t, p) - \
    _gibbs(n1, n1, n0, SA, pt0, 0)

    factor = factora / ( gp * _gibbs(n0, n2, n0, SA, t, p) )

    beta_const_pt = _gibbs(n0, n1, n1, SA, t, p) * factor - \
    _gibbs(n1, n0, n1, SA, t, p) / gp

    return beta_const_pt

def osmotic_coefficient(SA, t, p):
    r"""
    Calculates the osmotic coefficient of seawater.

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
    osmotic_coefficient : array_like
                          osmotic coefficient of seawater [unitless]

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
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.osmotic_coefficient(SA,t , p)
    array([ 0.90284718,  0.90298624,  0.90238866,  0.89880927,  0.89801054,
            0.89767912])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-09-28. Trevor McDougall and Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    n0 = 0

    molal = molality(SA) # molality of seawater in mol kg :sup:`-1`
    part = molal * cte.R * ( cte.Kelvin + t )

    #FIXME: SAzero: ValueError: shape mismatch: objects cannot be broadcast to a single shape
    SAzero = np.zeros( SA.shape )

    osmotic_coefficient = ( _gibbs(n0, n0, n0, SAzero, t, p) - chem_potential_water(SA, t, p) ) / part

    return osmotic_coefficient

def molality(SA):
    r"""
    Calculates the molality of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]

    Returns
    -------
    molal : array_like
            seawater molality [mol kg :sup:`-1`]

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
    >>> gsw.molality(SA)
    array([ 1.14508476,  1.15122708,  1.15581223,  1.14971265,  1.14593231,
            1.14578877])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-09-28. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA = np.asanyarray(SA)

    # only >= than zero
    Isalty = (SA >= 0)

    molal = np.nan * np.zeros( SA.shape )
    # molality of seawater in mol kg :sup:`-1`
    molal[Isalty] = SA[Isalty] / (cte.M_S * ( 1000 - SA[Isalty] ) )

    return molal

def ionic_strength(SA):
    r"""
    Calculates the ionic strength of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]

    Returns
    -------
    ionic_strength : array_like
                     ionic strength of seawater [mol kg :sup:`-1`]

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
    >>> gsw.ionic_strength(SA)
    array([ 0.71298118,  0.71680567,  0.71966059,  0.71586272,  0.71350891,
            0.71341953])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See Table L.1.

    .. [2] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: The composition of Standard Seawater and the definition of the Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72. see Eqns. 5.9 and 5.12.

    Modifications:
    2010-09-28. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA = np.asanyarray(SA)

    Z_2 = 1.2452898 # the valence factor of sea salt

    molal = molality(SA) # molality of seawater in mol kg :sup:`-1`

    ionic_strength = 0.5 * Z_2 * molal

    return ionic_strength

"""
Section C: extra functions for Temperature
TODO: study the possibility of transforming this into a DataSet like class with "properties" defining the different temperatures.
NOTE: section B calls: CT_from_pt(SA, pt0), pt_from_CT(SSO, CT0), pt0_from_t(SA,t, p)
pt_from_t(SA, t, p, pr=0) was renamed to potential_t(SA, t, p, pr=0) in Section B
"""

def CT_from_pt(SA, pt): #NOTE: used in conservative_t(SA, t, p)
    r"""
    Calculates Conservative Temperature of seawater from potential temperature (whose reference sea pressure is zero dbar).

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

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
    >>> gsw.CT_from_pt(SA, pt)
    array([ 28.80992302,  28.43914426,  22.78624661,  10.22616561,
             6.82718342,   4.32356518])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 3.3.

    Modifications:
    2010-08-05. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, pt = np.asanyarray(SA), np.asanyarray(pt)

    SA[SA < 0] = 0

    x2 = cte.sfac * SA
    x = np.sqrt(x2)
    y = pt * 0.025 # normalize for F03 and F08

    pot_enthalpy =  61.01362420681071 + y * ( 168776.46138048015 + \
    y * ( -2735.2785605119625 + y * ( 2574.2164453821433 + \
    y * ( -1536.6644434977543 + y * ( 545.7340497931629 + \
    ( -50.91091728474331 - 18.30489878927802 * y ) * y ) ) ) ) ) + \
    x2 * ( 268.5520265845071 + y * ( -12019.028203559312 + \
    y * ( 3734.858026725145 + y * ( -2046.7671145057618 + \
    y * ( 465.28655623826234 + ( -0.6370820302376359 - \
    10.650848542359153 * y ) * y ) ) ) ) + \
    x * ( 937.2099110620707 + y * ( 588.1802812170108 + \
    y * ( 248.39476522971285 + ( -3.871557904936333 - \
    2.6268019854268356 * y ) * y ) ) + \
    x * ( -1687.914374187449 + x * ( 246.9598888781377 + \
    x * ( 123.59576582457964 - 48.5891069025409 * x ) ) + \
    y * ( 936.3206544460336 + \
    y * ( -942.7827304544439 + y * ( 369.4389437509002 + \
    ( -33.83664947895248 - 9.987880382780322 * y ) * y ) ) ) ) ) )

    """
    The above polynomial for pot_enthalpy is the full expression for potential enthalpy in terms of SA and pt, obtained from the Gibbs function as below.

    It has simply collected like powers of x and y so that it is computationally faster than calling the Gibbs function twice as is done in the commented code below. When this code below is run, the results are identical to calculating pot_enthalpy as above, to machine precision.

    n0, n1 = 0, 1
    pot_enthalpy = _gibbs(n0, n0, n0, SA, pt, 0) - (cte.Kelvin + pt) * _gibbs(n0, n1, n0, SA, pt, 0)

    #----------------This is the end of the alternative code------------------
    #%timeit gsw.CT_from_pt(SA, pt)
    #1000 loops, best of 3: 1.34 ms per loop <- calling gibbs
    #1000 loops, best of 3: 254 us per loop <- standard
    """

    CT = pot_enthalpy / cte.cp0

    return CT

def pt_from_CT(SA, CT): #NOTE: used in specvol_anom(SA,t, p) inside gibbs.py
    r"""
    Calculates potential temperature (with a reference sea pressure of zero dbar) from Conservative Temperature.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    pt : array_like
         potential temperature referenced to a sea pressure of zero dbar [:math:`^\circ` C (ITS-90)]

    See Also
    --------
    TODO

    Notes
    -----
    This function uses 1.5 iterations through a modified Newton-Raphson (N-R) iterative solution procedure, starting from a rational-function-based initial condition for both pt and dCT_dpt.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.pt_from_CT(SA, CT)
    array([ 28.78317705,  28.4209556 ,  22.78495347,  10.23053439,
             6.82921659,   4.32453484])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See sections 3.1 and 3.3.

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater.  To be submitted to Ocean Science Discussions.

    Modifications:
    2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, CT = np.asanyarray(SA), np.asanyarray(CT)

    SA[SA < 0] = 0

    s1 = SA * 35. / cte.SSO

    a0 = -1.446013646344788e-2
    a1 = -3.305308995852924e-3
    a2 =  1.062415929128982e-4
    a3 =  9.477566673794488e-1
    a4 =  2.166591947736613e-3
    a5 =  3.828842955039902e-3

    b0 =  1.000000000000000e+0
    b1 =  6.506097115635800e-4
    b2 =  3.830289486850898e-3
    b3 =  1.247811760368034e-6

    a5CT = a5 * CT
    b3CT = b3 * CT
    CT_factor = ( a3 + a4 * s1 + a5CT )
    pt_num = a0 + s1 * ( a1 + a2 * s1 ) + CT * CT_factor
    pt_den = b0 + b1 * s1 + CT * ( b2 + b3CT )
    pt = pt_num / pt_den

    dCT_dpt = pt_den / (CT_factor + a5CT - (b2 + b3CT + b3CT) * pt)

    # start the 1.5 iterations through the modified Newton-Rapshon iterative method
    CT_diff = CT_from_pt(SA, pt) - CT
    pt_old = pt
    pt = pt_old - CT_diff /dCT_dpt # 1/2-way through the 1st modified N-R loop
    ptm = 0.5 * (pt + pt_old)

    #This routine calls gibbs_pt0_pt0(SA,pt0) to get the second derivative of the Gibbs function with respect to temperature at zero sea pressure.

    dCT_dpt = -(ptm + cte.Kelvin) * _gibbs_pt0_pt0(SA, ptm) / cte.cp0
    pt = pt_old - CT_diff / dCT_dpt # end of 1st full modified N-R iteration
    CT_diff = CT_from_pt(SA, pt) - CT
    pt_old = pt
    pt = pt_old - CT_diff / dCT_dpt # 1.5 iterations of the modified N-R method

    return pt

def t_from_CT(SA, CT, p):
    r"""
    Calculates in situ temperature from Conservative Temperature of seawater.

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
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]

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
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.t_from_CT(SA, CT, p)
    array([ 28.78558023,  28.43287225,  22.81032309,  10.26001075,
             6.8862863 ,   4.40362445])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See sections 3.1 and 3.3.

    Modifications:
    2010-08-26. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, CT, p = np.asanyarray(SA), np.asanyarray(CT), np.asanyarray(p)

    pt0 = pt_from_CT(SA, CT)

    t = potential_t(SA, pt0, 0, p)

    return t

def pt0_from_t(SA, t, p): #NOTE: potential_t does has the same result (only slower)
    r"""
    Calculates potential temperature with reference pressure, pr = 0 dbar. The present routine is computationally faster than the more general function "potential_t(SA, t, p, pr)" which can be used for any reference pressure value.

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
    pt0 : array_like
          potential temperature relative to 0 dbar [:math:`^\circ` C (ITS-90)]

    See Also
    --------
    _entropy_part, _gibbs_pt0_pt0, _entropy_part_zerop

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> t = [28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> gsw.pt0_from_t(SA, t, p)
    array([ 28.78319682,  28.42098334,  22.7849304 ,  10.23052366,
             6.82923022,   4.32451057])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 3.1.

    .. [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater.  To be submitted to Ocean Science Discussions.

    Modifications:
    2010-08-26. Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t, p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    s1 = SA * (35. / cte.SSO)

    pt0 = t + p * ( 8.65483913395442e-6 - \
             s1 *   1.41636299744881e-6 - \
              p *   7.38286467135737e-9 + \
              t * (-8.38241357039698e-6 + \
             s1 *   2.83933368585534e-8 + \
              t *   1.77803965218656e-8 + \
              p *   1.71155619208233e-10 ) )

    dentropy_dt = cte.cp0 / ( (cte.Kelvin + pt0) * ( 1 - 0.05 * ( 1 - SA / cte.SSO ) ) )

    true_entropy_part = _entropy_part(SA, t, p)

    for Number_of_iterations in range(0,2,1):
        pt0_old = pt0
        dentropy = _entropy_part_zerop(SA, pt0_old) - true_entropy_part
        pt0 = pt0_old - dentropy / dentropy_dt # this is half way through the modified method
        pt0m = 0.5 * (pt0 + pt0_old);
        dentropy_dt = -_gibbs_pt0_pt0(SA, pt0m)
        pt0 = pt0_old - dentropy / dentropy_dt

    """
    maximum error of 6.3x10^-9 degrees C for one iteration.
    maximum error is 1.8x10^-14 degrees C for two iterations (two iterations is the default, "for Number_of_iterations = 1:2")

    These errors are over the full "oceanographic funnel" of McDougall et al. (2010), which reaches down to p = 8000 dbar.
    """

    return pt0

def t_from_entropy(SA, entropy, t_type='pt'):
    r"""
    Calculates potential temperature with reference pressure pr = 0 dbar or Conservative temperature from entropy.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    entropy : array_like
              specific entropy [J kg :sup:`-1` K :sup:`-1`]
    t_type : str, optional
             'pt' for potential temperature [:math:`^\circ` C (ITS-90)]
             'CT' for Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    t : array_like
        potential temperature [:math:`^\circ` C (ITS-90)] (Default) or Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    See Also
    --------
    _gibbs_pt0_pt0

    Notes
    -----
    TODO

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> entropy = [400.3892, 395.4378, 319.8668, 146.7910, 98.6471, 62.7919]
    >>> gsw.t_from_entropy(SA, entropy)
    array([ 28.78317983,  28.42095483,  22.78495274,  10.23053207,
             6.82921333,   4.32453778])
    >>> gsw.t_from_entropy(SA, entropy, t_type='CT')
    array([ 28.80990279,  28.43919923,  22.78619927,  10.22619767,
             6.82719674,   4.32360295])


    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See appendix  A.10.

    Modifications:
    2010-10-13. Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, entropy = np.asanyarray(SA), np.asanyarray(entropy)

    SA[SA < 0] = 0

    n0 = 0
    n1 = 1

    part1 = 1 - SA / cte.SSO
    part2 = 1 - 0.05 * part1
    ent_SA = (cte.cp0 / cte.Kelvin) * part1 * ( 1 - 1.01 * part1)
    c = (entropy - ent_SA) * part2 / cte.cp0
    pt = cte.Kelvin * (np.exp(c) - 1)
    dentropy_dt = cte.cp0 / ( (cte.Kelvin + pt) * part2) #this is the initial value of dentropy_dt

    for Number_of_iterations in range(0,3):
        pt_old = pt
        dentropy = entropy_from_t(SA, pt_old, t_type='pt') - entropy
        pt = pt_old - dentropy / dentropy_dt #this is half way through the modified method
        ptm = 0.5 * (pt + pt_old)
        dentropy_dt = -_gibbs_pt0_pt0(SA, ptm)
        t = pt_old - dentropy / dentropy_dt

    if t_type == 'CT':
        t = CT_from_pt(SA, t)

    return t

def entropy_from_t(SA, t, t_type='pt'):
    r"""
    Calculates specific entropy of seawater.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
        temperature [:math:`^\circ` C]
    t_type : str, optional
             'pt' for potential temperature [:math:`^\circ` C (ITS-90)]
             'CT' for Conservative Temperature [:math:`^\circ` C (TEOS-10)]

    Returns
    -------
    entropy : array_like
              specific entropy [J kg :sup:`-1` K :sup:`-1`]

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
    >>> pt = [28.7832, 28.4210, 22.7850, 10.2305, 6.8292, 4.3245]
    >>> gsw.entropy_from_t(SA, pt)
    array([ 400.38946744,  395.43839949,  319.86743859,  146.79054828,
             98.64691006,   62.79135672])
    >>> CT = [28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]
    >>> gsw.entropy_from_t(SA, CT, 'CT')
    array([ 400.38916315,  395.43781023,  319.86680989,  146.79103279,
             98.64714648,   62.79185763])

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See appendix A.10.

    Modifications:
    2010-10-13. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, t = np.asanyarray(SA), np.asanyarray(t)

    SA[SA < 0] = 0

    n0, n1 = 0, 1

    if t_type == 'pt':
        entropy = -_gibbs(n0, n1, n0, SA, t, 0)
    elif t_type == 'CT':
        pt0 = pt_from_CT(SA, t)
        entropy = -_gibbs(n0, n1, n0, SA, pt0, 0)

    return entropy

"""
Section D: extra functions for Depth, Pressure and Distance
TODO: create a class for Depth pressure conversions
"""

def  z_from_p(p, lat):
    r"""
    Calculates height from sea pressure using the computationally-efficient 25-term expression for density in terms of SA, CT and p.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    p : array_like
        pressure [dbar]

    Returns
    -------
    z : array_like
        height [m]

    See Also
    --------
    _enthalpy_SSO_0_CT25


    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lat = 4
    >>> gsw.z_from_p(p, lat)
    array([  -9.94460074,  -49.71817465, -124.2728275 , -248.47044828,
           -595.82618014, -992.0931748 ])

    Notes
    -----
    At sea level z = 0, and since z (HEIGHT) is defined to be positive upwards, it follows that while z is positive in the atmosphere, it is NEGATIVE in the ocean.

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    ,, [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater.

    .. [3] Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.

    Modifications:
    2010-08-26. Trevor McDougall, Claire Roberts-Thomson & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    p, lat = np.asanyarray(p), np.asanyarray(lat)

    X     = np.sin( np.deg2rad(lat) )
    sin2  = X**2
    B     = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2 )
    A     = -0.5 * cte.gamma * B
    C     = _enthalpy_SSO_0_CT25(p)
    z     = -2 * C / ( B + np.sqrt( B**2 - 4 * A * C ) )

    return z

def  p_from_z(z, lat):
    r"""
    Calculates sea pressure from height using computationally-efficient 25-term expression for density, in terms of SA, CT and p.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    z : array_like
        height [m]

    Returns
    -------
    p : array_like
        pressure [dbar]

    See Also
    --------
    _specvol_SSO_0_CT25, _enthalpy_SSO_0_CT25

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> z = [10, 50, 125, 250, 600, 1000]
    >>> lat = 4.
    >>> gsw.p_from_z(z, lat)
    array([  -10.05521794,   -50.2711751 ,  -125.6548857 ,  -251.23284504,
            -602.44050752, -1003.07609807])
    >>> z = [-9.94460074, -49.71817465, -124.2728275, -248.47044828, -595.82618014, -992.0931748]
    >>> gsw.p_from_z(z, lat)
    array([   10.,    50.,   125.,   250.,   600.,  1000.])

    Notes
    -----
    Height (z) is NEGATIVE in the ocean.  Depth is -z. Depth is not used in the gibbs library.

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    ,, [2] McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term expression for the density of seawater in terms of Conservative Temperature, and related properties of seawater.

    .. [3] Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.

    .. [4] Saunders, P. M., 1981: Practical conversion of pressure to depth. Journal of Physical Oceanography, 11, 573-574.

    Modifications:
    2010-08-26. Trevor McDougall, Claire Roberts-Thomson and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    z, lat = np.asanyarray(z), np.asanyarray(lat)

    X     = np.sin( np.deg2rad(lat) )
    sin2  = X**2
    gs    = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2 )
    # get the first estimate of p from Saunders (1981)
    c1 =  5.25e-3 * sin2 + 5.92e-3
    p  = -2 * z / ( (1-c1) + np.sqrt( (1-c1) * (1-c1) + 8.84e-6 * z ) )
    # end of the first estimate from Saunders (1981)
    df_dp = cte.db2Pascal * _specvol_SSO_0_CT25(p) # initial value of the derivative of f
    f     = _enthalpy_SSO_0_CT25(p) + gs * ( z - 0.5 * cte.gamma * ( z**2 ) )
    p_old = p
    p     = p_old - f / df_dp
    pm    = 0.5 * (p + p_old)
    df_dp = cte.db2Pascal * _specvol_SSO_0_CT25(pm)
    p     = p_old - f / df_dp

    return p

def grav(lat, p=0):
    r"""
    Calculates acceleration due to gravity as a function of latitude and as a function of pressure in the ocean.

    Parameters
    ----------
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    p : number or array_like. Default p = 0
        pressure [dbar]

    Returns
    -------
    g : array_like
        gravity [m s :sup:`2`]

    See Also
    --------
    TODO

    Notes
    -----
    In the ocean z is negative.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> lat = [-90, -60, -30, 0]
    >>> p = 0
    >>> gsw.grav(lat, p)
    array([ 9.83218621,  9.81917886,  9.79324926,  9.780327  ])
    >>> gsw.grav(45)
    9.8061998770458008

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    .. [2] Moritz (2000) Goedetic reference system 1980. J. Geodesy, 74, 128-133.

    .. [3] Saunders, P.M., and N.P. Fofonoff (1976) Conversion of pressure to depth in the ocean. Deep-Sea Res.,pp. 109 - 111.

    Modifications:
    2010-07-23. Trevor McDougall & Paul Barker
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    lat, p, = np.asanyarray(lat), np.asanyarray(p)

    X = np.sin( np.deg2rad(lat) )
    sin2 = X**2
    gs = 9.780327 * ( 1.0 + ( 5.2792e-3 + ( 2.32e-5 * sin2 ) ) * sin2)
    z = z_from_p(p, lat)
    grav = gs * (1 - cte.gamma * z) # z is the height corresponding to p
    return grav

def distance(lon, lat, p=None):
    r"""
    Calculates the distance in met res between successive points in the vectors lon and lat, computed using the Haversine formula on a spherical earth of radius 6,371 km, being the radius of a sphere having the same volume as Earth. For a spherical Earth of radius 6,371,000 m, one nautical mile is 1,853.2488 m, thus one degree of latitude is 111,194.93 m.

    Haversine formula:
        R = earth's radius (mean radius = 6,371 km)

    .. math::
        a = \sin^2(\delta \text{lat}/2) + \cos(\text{lat}_1) \cos(\text{lat}_2) \sin^2(\delta \text{lon}/2)

        c = 2 \times \text{atan2}(\sqrt{a}, \sqrt{(1-a)})

        d = R \times c

    Parameters
    ----------
    lon : array_like
          decimal degrees east [0..+360] or [-180 ... +180]
    lat : array_like
          latitude in decimal degrees north [-90..+90]
    p : number or array_like. Default p = 0
        pressure [dbar]

    Returns
    -------
    dist: array_like
          distance between points on a spherical Earth at pressure (p) [m]

    See Also
    --------
    TODO

    Notes
    -----
    Distances are probably good to better than 1\% of the "true" distance on the ellipsoidal earth.
    The check value below differ from the original online docs at "http://www.teos-10.org/pubs/gsw/html/gsw_distance.html" but agree with the result.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> lon = [159, 220]
    >>> lat = [-35, 35]
    >>> gsw.distance(lon, lat)
    array([[ 10030974.652916]])
    >>> p = [200, 1000]
    >>> gsw.distance(lon, lat, p)
    array([[ 10030661.63878009]])

    References
    ----------
    .. [1] http://www.eos.ubc.ca/~rich/map.html

    Modifications:
    2000-11-06. Rich Pawlowicz
    2010-07-28. Paul Barker and Trevor McDougall
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    lon, lat = np.asanyarray(lon), np.asanyarray(lat)

    if (lon.size == 1) & (lat.size == 1):
        raise NameError('more than one point is needed to compute distance')
    elif lon.ndim != lat.ndim:
        raise NameError('lon, lat must have the same dimension')

    if (lon.size == 1) & (lat.size != 1): # fill lon with lat size
        lon = _check_dim(lon, lat)
    elif (lat.size == 1) & (lon.size != 1): # fill lat with lon size
        lat = _check_dim(lat, lon)

    if (lon.ndim == 1) & (lat.ndim == 1): #NOTE: Ugly way to "matlabsize it"
        lon = lon[np.newaxis,:]
        lat = lat[np.newaxis,:]

    # check for lon/lat and p dimensions
    if p == None:
        p = np.zeros( lon.shape )
    else:
        p = np.asanyarray(p)

    if p.ndim > lat.ndim:
        lon = _check_dim(lon, p)
        lat = _check_dim(lat, p)
    elif p.ndim == 1:
        p = _check_dim(p, lon)

    dlon = np.deg2rad( np.diff(lon, axis=1) )
    dlat = np.deg2rad( np.diff(lat, axis=1) )

    a = ( np.sin(dlat/2.) )**2 + np.cos( np.deg2rad( lat[:,:-1] ) ) * np.cos( np.deg2rad( lat[:,1:] ) ) * ( np.sin(dlon/2.) )**2

    angles = 2. * np.arctan2( np.sqrt(a), np.sqrt(1-a) )

    p_mid = 0.5 * (   p[:,0:-1] +   p[:,0:-1] )
    lat_mid = 0.5 * ( lat[:,:-1] + lat[:,1:] )

    z = z_from_p(p_mid, lat_mid) #NOTE: z is height and is negative in the ocean

    distance = (cte.a + z) * angles

    return distance

"""
Section E: extra functions for Salinity
TODO: study the possibility of transforming this into a DataSet like class with "properties" defining the different salinities.
NOTE: ---
"""

""" cndr_from_SP == sw.cndr """
from  seawater.csiro import cndr as cndr_from_SP
""" SP_from_cndr == sw.salt """
from  seawater.csiro import salt as SP_from_cndr

def SA_from_SP(SP, p, lon, lat):
    r"""
    Calculates Absolute Salinity from Practical Salinity.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    _delta_SA, _SA_from_SP_Baltic

    Notes
    -----
    Since SP is non-negative by definition, this function changes any negative input values of SP to be zero.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SP = [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat = 188, 4
    >>> gsw.SA_from_SP(SP, p, lon, lat)
    (array([ 34.71177971,  34.89152372,  35.02554774,  34.84723008,
            34.7366296 ,  34.73236186]), array([ True,  True,  True,  True,  True,  True], dtype=bool))

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SP, p, lon, lat = np.asanyarray(SP), np.asanyarray(p), np.asanyarray(lon), np.asanyarray(lat)

    p = _check_dim(p, SP)
    lon, lat = _check_dim(lon, SP), _check_dim(lat, SP)

    SP[SP < 0] = 0

    if (lon < 0).any():
        lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(SP)

    SA = np.NaN * np.zeros( SP.shape )
    dSA = np.NaN * np.zeros( SP.shape )

    in_ocean = np.bool_( np.ones( SP.shape ) )

    dSA[inds], in_ocean[inds] = _delta_SA( p[inds], lon[inds], lat[inds] )
    SA[inds] = ( cte.SSO / 35 ) * SP[inds] + dSA[inds]
    SA_baltic = _SA_from_SP_Baltic( SP, lon, lat )

    indsbaltic = ~np.isnan(SA_baltic)

    SA[indsbaltic] = SA_baltic[indsbaltic]

    return SA, in_ocean

def SA_from_Sstar(Sstar, p, lon, lat):
    r"""
    Calculates Absolute Salinity from Preformed Salinity.

    Parameters
    ----------
    Sstar : array_like
            Preformed Salinity [g kg :sup:`-1`]
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    _delta_SA

    Notes
    -----
    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometers inland from the coast.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> Sstar = [34.7115, 34.8912, 35.0247, 34.8436, 34.7291, 34.7197]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat = 188, 4
    >>> gsw.SA_from_Sstar(Sstar, p, lon, lat)
    (array([ 34.71172651,  34.89156271,  35.02559848,  34.84723731,
            34.736696  ,  34.73238548]), array([ True,  True,  True,  True,  True,  True], dtype=bool))

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    Sstar, p, lon, lat = np.asanyarray(Sstar), np.asanyarray(p), np.asanyarray(lon), np.asanyarray(lat)

    p = _check_dim(p, Sstar)
    lon = _check_dim(lon, Sstar)
    lat = _check_dim(lat, Sstar)

    if (lon < 0).any():
        lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(Sstar)
    SA = np.NaN * np.zeros( Sstar.shape )
    dSA = np.NaN * np.zeros( Sstar.shape )

    in_ocean = np.bool_( np.ones( Sstar.shape ) )

    dSA[inds], in_ocean[inds] = _delta_SA( p[inds], lon[inds], lat[inds] )

    SA[inds] = Sstar[inds] + ( 1 + cte.r1 ) * dSA[inds]

    #NOTE: In the Baltic Sea, SA = Sstar, and note that _delta_SA returns zero for dSA in the Baltic.

    return SA, in_ocean

def SP_from_SA(SA, p, lon, lat):
    r"""
    Calculates Practical Salinity from Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    _delta_SA, _SP_from_SA_Baltic

    Notes
    -----
    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometers inland from the coast.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat =  188, 4
    >>> gsw.SP_from_SA(SA, p, lon, lat)
    (array([ 34.54872019,  34.72747639,  34.86055202,  34.68097006,
            34.56797054,  34.56003796]), array([ True,  True,  True,  True,  True,  True], dtype=bool))

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, p, lon, lat = np.asanyarray(SA), np.asanyarray(p), np.asanyarray(lon), np.asanyarray(lat)

    p = _check_dim(p, SA)
    lat = _check_dim(lat, SA)
    lon = _check_dim(lon, SA)

    if (lon < 0).any():
        lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(SA)
    SP = np.NaN * np.zeros( SA.shape )
    dSA = np.NaN * np.zeros( SA.shape )

    in_ocean = np.bool_( np.ones( SA.shape ) )

    dSA[inds], in_ocean[inds] = _delta_SA( p[inds], lon[inds], lat[inds] )

    SP[inds] = (35./cte.SSO) * ( SA[inds] - dSA[inds] )

    SP_baltic = _SP_from_SA_Baltic( SA, lon, lat )

    indsbaltic = ~np.isnan(SP_baltic)

    SP[indsbaltic] = SP_baltic[indsbaltic]

    return SP, in_ocean

def Sstar_from_SA(SA, p, lon, lat):
    r"""
    Converts Preformed Salinity from Absolute Salinity.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    Sstar : array_like
            Preformed Salinity [g kg :sup:`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    _delta_SA

    Notes
    -----
    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometers inland from the coast.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat =  188, 4
    >>> gsw.Sstar_from_SA(SA, p, lon, lat)
    (array([ 34.71157349,  34.89113729,  35.02470152,  34.84356269,
            34.729004  ,  34.71971452]), array([ True,  True,  True,  True,  True,  True], dtype=bool))

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SA, p, lon, lat = np.asanyarray(SA), np.asanyarray(p), np.asanyarray(lon), np.asanyarray(lat)

    p = _check_dim(p, SA)
    lon = _check_dim(lon, SA)
    lat = _check_dim(lat, SA)

    if (lon < 0).any():
        lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(SA)
    Sstar = np.NaN * np.zeros( SA.shape )
    dSA = np.NaN * np.zeros( SA.shape )

    in_ocean = np.bool_( np.ones( SA.shape ) )

    dSA[inds], in_ocean[inds] = _delta_SA( p[inds], lon[inds], lat[inds] )

    Sstar[inds] =  SA[inds] - ( 1 + cte.r1 ) * dSA[inds]
    #NOTE: In the Baltic Sea, SA = Sstar, and note that _delta_SA returns zero for dSA in the Baltic.

    return Sstar, in_ocean

def SP_from_Sstar(Sstar, p, lon, lat):
    r"""
    Calculates Practical Salinity from Preformed Salinity.

    Parameters
    ----------
    Sstar : array_like
            Preformed Salinity [g kg :sup:`-1`]
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    _delta_SA, _SP_from_SA_Baltic

    Notes
    -----
    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometers inland from the coast.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> Sstar = [34.7115, 34.8912, 35.0247, 34.8436, 34.7291, 34.7197]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat =  188, 4
    >>> gsw.SP_from_Sstar(Sstar, p, lon, lat)
    (array([ 34.54864705,  34.72753881,  34.8605505 ,  34.68100719,
            34.56806609,  34.56002351]), array([ True,  True,  True,  True,  True,  True], dtype=bool))

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    Sstar, p, lon, lat = np.asanyarray(Sstar), np.asanyarray(p), np.asanyarray(lon), np.asanyarray(lat)

    p = _check_dim(p, Sstar)
    lon = _check_dim(lon, Sstar)
    lat = _check_dim(lat, Sstar)

    if (lon < 0).any():
        lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(Sstar)
    SP = np.NaN * np.zeros( Sstar.shape )
    dSA = np.NaN * np.zeros( Sstar.shape )

    in_ocean = np.bool_( np.ones( Sstar.shape ) )

    dSA[inds], in_ocean[inds] = _delta_SA( p[inds], lon[inds], lat[inds] )
    SP[inds] = (35/cte.SSO) * ( Sstar[inds] + cte.r1 * dSA[inds] )

    # In the Baltic Sea, SA = Sstar.
    SP_baltic = _SP_from_SA_Baltic( Sstar, lon, lat )

    indsbaltic = ~np.isnan(SP_baltic)

    SP[indsbaltic] = SP_baltic[indsbaltic]

    return SP, in_ocean

def Sstar_from_SP(SP, p, lon, lat):
    r"""
    Calculates Preformed Salinity from Absolute Salinity.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    Sstar : array_like
            Preformed Salinity [g kg :sup:`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    _delta_SA, _SA_from_SP_Baltic

    Notes
    -----
    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometers inland from the coast.

    Since SP is non-negative by definition, this function changes any negative input values of SP to be zero.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SP = [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat =  188, 4
    >>> gsw.Sstar_from_SP(SP, p, lon, lat)
    (array([ 34.7115532 ,  34.89116101,  35.02464926,  34.84359277,
            34.7290336 ,  34.71967638]), array([ True,  True,  True,  True,  True,  True], dtype=bool))

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    ,, [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall and Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SP, p, lon, lat = np.asanyarray(SP), np.asanyarray(p), np.asanyarray(lon), np.asanyarray(lat)

    p = _check_dim(p, SP)
    lon = _check_dim(lon, SP)
    lat = _check_dim(lat, SP)

    SP[SP < 0] = 0

    if (lon < 0).any():
        lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(SP)
    Sstar = np.NaN * np.zeros( SP.shape )
    dSA = np.NaN * np.zeros( SP.shape )

    in_ocean = np.bool_( np.ones( SP.shape ) )

    dSA[inds], in_ocean[inds] = _delta_SA( p[inds], lon[inds], lat[inds] )
    Sstar[inds] = (cte.SSO/35.) * SP[inds] - cte.r1 * dSA[inds]

    # In the Baltic Sea, SA == Sstar.
    Sstar_baltic = _SA_from_SP_Baltic( SP, lon, lat )

    indsbaltic = ~np.isnan(Sstar_baltic)

    Sstar[indsbaltic] = Sstar_baltic[indsbaltic]

    return Sstar, in_ocean

def SA_Sstar_from_SP(SP, p, lon, lat):
    r"""
    Calculates Absolute Salinity and Preformed Salinity from Practical Salinity.

    Parameters
    ----------
    SP : array_like
         salinity [psu (PSS-78)], unitless
    p : array_like
        pressure [dbar]
    lon : array_like
          decimal degrees east [0..+360] or [-180..+180]
    lat : array_like
          decimal degrees (+ve N, -ve S) [-90..+90]

    Returns
    -------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    Sstar : array_like
            Preformed Salinity [g kg :sup:`-1`]
    in_ocean : False, if [lon, lat] are a long way from the ocean
               True, [lon, lat] are in the ocean.

    See Also
    --------
    _delta_SA, _SA_from_SP_Baltic

    Notes
    -----
    The in_ocean flag is only set when the observation is well and truly on dry land; often the warning flag is not set until one is several hundred kilometers inland from the coast.

    Since SP is non-negative by definition, this function changes any negative input values of SP to be zero.

    Examples
    --------
    >>> import seawater.gibbs as gsw
    >>> SP = [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lon, lat =  188, 4
    >>> gsw.SA_Sstar_from_SP(SP, p, lon, lat)
    (array([ 34.71177971,  34.89152372,  35.02554774,  34.84723008,
            34.7366296 ,  34.73236186]), array([ 34.7115532 ,  34.89116101,  35.02464926,  34.84359277,
            34.7290336 ,  34.71967638]), array([ True,  True,  True,  True,  True,  True], dtype=bool))

    References
    ----------
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.5 and appendices A.4 and A.5.

    .. [2] McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm for estimating Absolute Salinity in the global ocean. Submitted to Ocean Science. A preliminary version is available at Ocean Sci. Discuss., 6, 215-242.
    http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf

    Modifications:
    2010-07-23. David Jackett, Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    SP, p, lon, lat = np.asanyarray(SP), np.asanyarray(p), np.asanyarray(lon), np.asanyarray(lat)

    p = _check_dim(p, SP)
    lat = _check_dim(lat, SP)
    lon = _check_dim(lon, SP)

    SP[SP < 0] = 0

    if (lon < 0).any():
        lon[lon < 0] = lon[lon < 0] + 360.

    inds = np.isfinite(SP)

    SA = np.NaN * np.zeros( SP.shape )
    Sstar = np.NaN * np.zeros( SP.shape )
    dSA = np.NaN * np.zeros( SP.shape )

    in_ocean = np.bool_( np.ones( SP.shape ) )

    dSA[inds], in_ocean[inds] = _delta_SA( p[inds], lon[inds], lat[inds] )

    SA[inds] = ( cte.SSO / 35 ) * SP[inds] + dSA[inds]
    Sstar[inds] = ( cte.SSO / 35 ) * SP[inds] - cte.r1 * dSA[inds]

    SA_baltic = _SA_from_SP_Baltic( SP, lon, lat )

    indsbaltic = ~np.isnan(SA_baltic)

    SA[indsbaltic] = SA_baltic[indsbaltic]
    Sstar[indsbaltic] = SA_baltic[indsbaltic]
    #NOTE: In the Baltic Sea, Sstar == SA.

    return SA, Sstar, in_ocean

def SA_from_rho(rho, t, p):
    r"""
    Calculates the Absolute Salinity of a seawater sample, for given values of its density, in situ temperature and sea pressure (in dbar).

    One use for this function is in the laboratory where a measured value of the in situ density :math:`\rho` of a seawater sample may have been made at the laboratory temperature :math:`t` and at atmospheric pressure :math:`p`. The present function will return the Absolute Salinity SA of this seawater sample.

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
    This is expressed on the Reference-Composition Salinity Scale of Millero et al. (2008).

    After two iterations of a modified Newton-Raphson iteration, the error in SA is typically no larger than 2 :math:`^\times` 10 :sup:`-13` [g kg :sup:`-1`]

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
    .. [1] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010: Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp. See section 2.5.

    .. [2] Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: The composition of Standard Seawater and the definition of the Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.

    Modifications:
    2010-08-23. Trevor McDougall & Paul Barker.
    2010-12-09. Filipe Fernandes, Python translation from gsw toolbox.
    """

    rho, t, p = np.asanyarray(rho), np.asanyarray(t), np.asanyarray(p)

    n0, n1 = 0, 1
    v_lab = np.ones( rho.shape ) / rho
    v_0 = _gibbs(n0, n0, n1, 0, t, p)
    v_120 = _gibbs(n0, n0, n1, 120, t, p)
    SA = 120 * ( v_lab - v_0 ) / ( v_120 - v_0 ) # initial estimate of SA
    Ior = (SA < 0) | (SA > 120)
    v_SA = ( v_120 - v_0 ) / 120 # initial estimate of v_SA, the SA derivative of v

    for iter in range(0,2):
        SA_old = SA
        delta_v = _gibbs(n0, n0, n1, SA_old, t, p) - v_lab
        SA = SA_old - delta_v / v_SA # this is half way through the modified N-R method
        SA_mean = 0.5 * ( SA + SA_old )
        v_SA = _gibbs(n1, n0, n1, SA_mean, t, p)
        SA = SA_old - delta_v / v_SA

    SA[Ior] = np.NaN

    return SA

"""
Section F: Classes
"""

class Gibbs:
    r"""
    Class that aggregate all SA, t, p functions.

    Parameters
    ----------
    SA : array_like
         Absolute salinity [g kg :sup:`-1`]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        pressure [dbar]

    """
    def __init__(self, SA=None, t=None, p=None):
        self.SA, self.t, self.p = np.asanyarray(SA), np.asanyarray(t), np.asanyarray(p)

    def entropy(self):
        """
        See entropy docstring
        """
        return entropy(self.SA, self.t, self.p)

    def rho(self):
        """
        See rho docstring
        """
        return rho(self.SA, self.t, self.p)

    def cp(self):
        """
        See cp docstring
        """
        return cp(self.SA, self.t, self.p)

    def helmholtz_energy(self):
        """
        See helmholtz_energy docstring
        """
        return helmholtz_energy(self.SA, self.t, self.p)

    def internal_energy(self):
        """
        See internal_energy docstring
        """
        return internal_energy(self.SA, self.t, self.p)

    def sound_speed(self):
        """
        See sound_speed docstring
        """
        return sound_speed(self.SA, self.t, self.p)

    def adiabatic_lapse_rate(self):
        """
        See adiabatic_lapse_rate docstring
        """
        return adiabatic_lapse_rate(self.SA, self.t, self.p)

    def chem_potential_relative(self):
        """
        See chem_potential_relative docstring
        """
        return chem_potential_relative(self.SA, self.t, self.p)

    def specvol(self):
        """
        See specvol docstring
        """
        return specvol(self.SA, self.t, self.p)

    def conservative_t(self):
        """
        See conservative_t docstring
        """
        return conservative_t(self.SA, self.t, self.p)

    def potential_t(self, pr=0):
        """
        See potential_t docstring
        """
        return potential_t(self.SA, self.t, self.p)

    def enthalpy(self):
        """
        See enthalpy docstring
        """
        return enthalpy(self.SA, self.t, self.p)

    def chem_potential_water(self):
        """
        See chem_potential_water docstring
        """
        return chem_potential_water(self.SA, self.t, self.p)

    def chem_potential_salt(self):
        """
        See chem_potential_salt docstring
        """
        return chem_potential_salt(self.SA, self.t, self.p)

    def isochoric_heat_cap(self):
        """
        See isochoric_heat_cap docstring
        """
        return isochoric_heat_cap(self.SA, self.t, self.p)

    def kappa(self):
        """
        See kappa docstring
        """
        return kappa(self.SA, self.t, self.p)

    def kappa_const_t(self):
        """
        See kappa_const_t docstring
        """
        return kappa_const_t(self.SA, self.t, self.p)

    def pot_rho(self, pr=0):
        """
        See pot_rho docstring
        """
        return pot_rho(self.SA, self.t, self.p)

    def specvol_anom(self):
        """
        See specvol_anom docstring
        """
        return specvol_anom(self.SA, self.t, self.p)

    def alpha_wrt_t(self):
        """
        See alpha_wrt_t docstring
        """
        return alpha_wrt_t(self.SA, self.t, self.p)

    def alpha_wrt_CT(self):
        """
        See alpha_wrt_CT docstring
        """
        return alpha_wrt_CT(self.SA, self.t, self.p)

    def alpha_wrt_pt(self):
        """
        See alpha_wrt_pt docstring
        """
        return alpha_wrt_pt(self.SA, self.t, self.p)

    def beta_const_t(self):
        """
        See beta_const_t docstring
        """
        return beta_const_t(self.SA, self.t, self.p)

    def beta_const_CT(self):
        """
        See beta_const_CT docstring
        """
        return beta_const_CT(self.SA, self.t, self.p)

    def beta_const_pt(self):
        """
        See beta_const_pt docstring
        """
        return beta_const_pt(self.SA, self.t, self.p)

    def osmotic_coefficient(self):
        """
        See osmotic_coefficient docstring
        """
        return osmotic_coefficient(self.SA, self.t, self.p)

    def molality(self):
        """
        See molality docstring
        """
        return molality(self.SA)

    def ionic_strength(self):
        """
        See ionic_strength docstring
        """
        return ionic_strength(self.SA)

class Dict2Struc(object):
    r"""
    Open variables from a dictionary in a "matlab-like-structure"
    """
    def __init__(self, adict):
        self.__dict__.update(adict)

if __name__=='__main__':
    r"""
    This test only the Gibbs class
    """
    def test_print(STP, method, comp_value=None):
        """
        Run a function test mimicking the original logic. This is done to allow for a direct comparison of the result from the Matlab to the python package.
        """

        if comp_value is None:
            comp_value = method

        # test for floating differences with: computed - check value >= defined precision
        exec( "unequal = (gsw_cv." +comp_value+ " - STP." +method+ "() ) >= gsw_cv." +comp_value+ "_ca")

        width = 23
        if unequal.any():
            print "%s: Failed" % method.rjust(width)
        else:
            # test if check value is identical to computed value
            if eval( "( gsw_cv." +comp_value+ "[~np.isnan(gsw_cv."+comp_value+")] == STP." +method+ "()[~np.isnan(STP."+method+"())] ).all()" ):
                print "%s: Passed" % method.rjust(width)
            else:
                # test for differences in case their aren't equal. This is an attempt to place all tests together (i.e. term25 and small float differences that will appear)
                exec("nmax = ( gsw_cv."+comp_value+" - STP."+method+"() )[~np.isnan(gsw_cv."+comp_value+")].max()")
                exec("nmin = ( gsw_cv."+comp_value+" - STP."+method+"() )[~np.isnan(gsw_cv."+comp_value+")].min()")
                print "%s: Passed, but small diff ranging from: %s to %s" % ( method.rjust(width), nmax, nmin)

    """ test data: data/gsw_cv.pkl
    (1) a reference cast (for the McD_Klocker streamfunction)
    (2) three vertical profiles of (SP, t, p) at known long & lat
    (3) the outputs of all the GSW functions for these 3 profiles
    (4) the required accuracy of all these outputs.
    """

    data = pickle.load( open(os.path.join(datadir + 'gsw_cv.pkl'),'rb') )
    gsw_cv = Dict2Struc(data) # then type data.<tab> to navigate through your variables

    STP = Gibbs(gsw_cv.SA_from_SP, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)

    test_print(STP, "entropy")
    test_print(STP, "rho")
    test_print(STP, "cp")
    test_print(STP, "helmholtz_energy", "Helmholtz_energy")
    test_print(STP, "internal_energy")
    test_print(STP, "sound_speed")
    test_print(STP, "adiabatic_lapse_rate")
    test_print(STP, "chem_potential_relative", "chem_potential")
    test_print(STP, "specvol")
    test_print(STP, "molality")
    test_print(STP, "ionic_strength")
    test_print(STP, "potential_t", "pt_from_t")
    test_print(STP, "conservative_t", "CT_from_t")
    test_print(STP, "enthalpy")
    test_print(STP, "alpha_wrt_t")
    test_print(STP, "alpha_wrt_CT")
    test_print(STP, "alpha_wrt_pt")
    test_print(STP, "beta_const_t")
    test_print(STP, "beta_const_CT")
    test_print(STP, "beta_const_pt")
    test_print(STP, "chem_potential_water")
    test_print(STP, "chem_potential_salt")
    test_print(STP, "isochoric_heat_cap")
    test_print(STP, "kappa")
    test_print(STP, "kappa_const_t")
    test_print(STP, "osmotic_coefficient")
    test_print(STP, "pot_rho")
    test_print(STP, "specvol_anom")
