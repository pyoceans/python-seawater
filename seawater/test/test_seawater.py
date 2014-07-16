# -*- coding: utf-8 -*-
#
# test_seawater.py
#
# purpose:  Test to mimic matlab sw_test.m function.
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.github.io/
# created:  05-Aug-2013
# modified: Wed 16 Jul 2014 10:15:04 AM BRT
#
# obs:  Check the unittest for a more thorough code testing.
#


from __future__ import division

import sys
from platform import uname
from time import asctime, localtime

import numpy as np
import seawater as sw


def test(fileout='python-test.txt'):
    r"""Copy of the Matlab test.

    Modifications: Phil Morgan
                   03-12-12. Lindsay Pender, Converted to ITS-90.
    """
    f = open(fileout, 'w')
    asterisks = '*' * 76

    f.write(asterisks)
    f.write('\n    TEST REPORT    ')
    f.write('\n')
    f.write('\n SEA WATER LIBRARY %s' % sw.__version__)
    f.write('\n')
    # Show some info about this Python.
    f.write('\npython version: %s' % sys.version)
    f.write('\n on %s computer %s' % (uname()[0], uname()[-1]))
    f.write('\n')
    f.write('\n')
    f.write(asctime(localtime()))
    f.write('\n')
    f.write('\n')

    # Test main module  ptmp.
    f.write('\n%s' % asterisks)
    f.write('\n**  TESTING MODULE: ptmp')
    f.write('\n**  and SUB-MODULE: adtg')
    f.write('\n%s' % asterisks)
    f.write('\n')
    f.write('\n')

    # Test 1 - data from Unesco 1983 p45.
    T = np.array([[0,  0,  0,  0,  0,  0],
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

    UN_ptmp = np.array([[0, -0.3061, -0.9667,  0, -0.3856, -1.0974],
                        [10,  9.3531,  8.4684, 10,  9.2906,  8.3643],
                        [20, 19.0438, 17.9426, 20, 18.9985, 17.8654],
                        [30, 28.7512, 27.4353, 30, 28.7231, 27.3851],
                        [40, 38.4607, 36.9254, 40, 38.4498, 36.9023]])

    PT = sw.ptmp(S, T, P, Pr) * 1.00024

    # Display results.
    f.write('\n')
    f.write('\n%s' % asterisks)
    f.write('\nComparison of accepted values from UNESCO 1983 ')
    f.write('\n (Unesco Tech. Paper in Marine Sci. No. 44, p45)')
    f.write('\n%s' % asterisks)
    f.write('\n')

    m, n = S.shape  # TODO: so many loops there must be a better way.
    for icol in range(0, n):
        f.write('\n   Sal  Temp  Press     PTMP       ptmp')
        f.write('\n  (psu)  (C)   (db)     (C)          (C)\n')
        result = np.vstack((S[:, icol], T[:, icol], P[:, icol],
                            UN_ptmp[:, icol], PT[:, icol]))
        for iline in range(0, m):
            f.write(" %4.0f  %4.0f   %5.0f   %8.4f  %11.5f\n" %
                    tuple(result[:, iline]))

    # Test main module svan.
    f.write('\n%s' % asterisks)
    f.write('\n**  TESTING MODULE: svan')
    f.write('\n**  and SUB-MODULE: dens dens0 smow seck pden ptmp')
    f.write('\n%s' % asterisks)

    # Test data FROM: Unesco Tech. Paper in Marine Sci. No. 44, p22.
    s = np.array([0,     0,  0,     0, 35,    35, 35,   35])
    p = np.array([0, 10000,  0, 10000,  0, 10000,  0, 10000])
    t = np.array([0,     0, 30,    30,  0,     0, 30,    30]) / 1.00024

    UN_svan = np.array([2749.54, 2288.61, 3170.58, 3147.85,
                        0.0,    0.00,  607.14,  916.34])

    SVAN = sw.svan(s, t, p)

    # DISPLAY RESULTS
    f.write('\n')
    f.write('\n%s' % asterisks)
    f.write('\n')
    f.write('\nComparison of accepted values from UNESCO 1983')
    f.write('\n (Unesco Tech. Paper in Marine Sci. No. 44, p22)')
    f.write('\n%s' % asterisks)
    f.write('\n')
    f.write('\n   Sal  Temp  Press        SVAN        svan')
    f.write('\n  (psu)  (C)   (db)    (1e-8*m3/kg)  (1e-8*m3/kg)\n')
    result = np.vstack([s, t, p, UN_svan, 1e+8 * SVAN])
    for iline in range(0, len(SVAN)):
        f.write(" %4.0f  %4.0f   %5.0f   %11.2f    %11.3f\n" %
                tuple(result[:, iline]))

    # Test main module salt.
    f.write('\n%s' % asterisks)
    f.write('\n**  TESTING MODULE: salt')
    f.write('\n**  and SUB-MODULE: salrt salrp sals')
    f.write('\n%s' % asterisks)
    f.write('\n')

    # Test 1 - data from Unesco 1983 p9.
    R = np.array([1, 1.2, 0.65])  # cndr = R.
    T = np.array([15, 20, 5]) / 1.00024
    P = np.array([0, 2000, 1500])
    #Rt   = np.array([  1, 1.0568875, 0.81705885])
    UN_S = np.array([35, 37.245628,  27.995347])

    S = sw.salt(R, T, P)

    # Display results.
    f.write('\n')
    f.write('\n%s' % asterisks)
    f.write('\nComparison of accepted values from UNESCO 1983 ')
    f.write('\n(Unesco Tech. Paper in Marine Sci. No. 44, p9)')
    f.write('\n%s' % asterisks)
    f.write('\n')

    f.write('\n   Temp    Press       R              S           salt')
    f.write('\n   (C)     (db)    (no units)       (psu)          (psu)\n')
    table = np.vstack([T, P, R, UN_S, S])
    m, n = table.shape
    for iline in range(0, n):
        f.write(" %4.0f       %4.0f  %8.2f      %11.6f  %14.7f\n" %
                tuple(table[:, iline]))

    # Test main module cndr.
    f.write('\n%s' % asterisks)
    f.write('\n**  TESTING MODULE: cndr')
    f.write('\n**  and SUB-MODULE: salds')
    f.write('\n%s' % asterisks)

    # Test 1 - data from Unesco 1983 p9.
    T = np.array([0, 10, 0, 10, 10, 30]) / 1.00024
    P = np.array([0,  0, 1000, 1000, 0, 0])
    S = np.array([25, 25, 25, 25, 40, 40])
    UN_R = np.array([0.498088, 0.654990, 0.506244, 0.662975, 1.000073,
                     1.529967])
    R = sw.cndr(S, T, P)

    # Display results.
    f.write('\n')
    f.write('\n%s' % asterisks)
    f.write('\nComparison of accepted values from UNESCO 1983 ')
    f.write('\n (Unesco Tech. Paper in Marine Sci. No. 44, p14)')
    f.write('\n%s' % asterisks)
    f.write('\n')
    f.write('\n')

    f.write('\n   Temp    Press       S            cndr         cndr')
    f.write('\n   (C)     (db)      (psu)        (no units)    (no units)\n')
    table = np.vstack([T, P, S, UN_R, R])
    m, n = table.shape
    for iline in range(0, n):
        f.write(" %4.0f       %4.0f   %8.6f   %11.6f  %14.8f\n" %
                tuple(table[:, iline]))

    # Test main module depth.
    f.write('\n%s' % asterisks)
    f.write('\n**  TESTING MODULE: depth')
    f.write('\n%s' % asterisks)

    # Test data - matrix "pressure", vector "lat"  Unesco 1983 data p30.
    lat = np.array([0, 30, 45, 90])
    P = np.array([[500,   500,   500,  500],
                  [5000,  5000,  5000, 5000],
                  [10000, 10000, 10000, 10000]])

    UN_dpth = np.array([[496.65,  496.00,  495.34,  494.03],
                        [4915.04, 4908.56, 4902.08, 4889.13],
                        [9725.47, 9712.65, 9699.84, 9674.23]])

    dpth = sw.dpth(P, lat)

    # Display results.
    f.write('\n')
    f.write('\n%s' % asterisks)
    f.write('\nComparison of accepted values from Unesco 1983 ')
    f.write('\n(Unesco Tech. Paper in Marine Sci. No. 44, p28)')
    f.write('\n%s' % asterisks)
    f.write('\n')

    f.write('\n')
    for irow in range(0, 3):
        f.write('\n    Lat       Press     DPTH      dpth')
        f.write('\n  (degree)    (db)     (meter)    (meter)\n')
        table = np.vstack([lat, P[irow, :], UN_dpth[irow, :], dpth[irow, :]])
        m, n = table.shape
        for iline in range(0, n):
            f.write("  %6.3f     %6.0f   %8.2f   %8.3f\n" %
                    tuple(table[:, iline]))

    # Test main module fp.
    f.write('\n%s' % asterisks)
    f.write('\n**  TESTING MODULE: fp')
    f.write('\n%s' % asterisks)

    # Test 1 - UNESCO data p.30.
    S = np.array([[5, 10, 15, 20, 25, 30, 35, 40],
                  [5, 10, 15, 20, 25, 30, 35, 40]])

    P = np.array([[0,   0,   0,   0,   0,   0,   0,   0],
                  [500, 500, 500, 500, 500, 500, 500, 500]])

    UN_fp = np.array([[-0.274, -0.542, -0.812, -1.083, -1.358, -1.638, -1.922,
                       -2.212], [-0.650, -0.919, -1.188, -1.460, -1.735,
                                 -2.014, -2.299, -2.589]])

    FP = sw.fp(S, P)

    # Display results.
    f.write('\n')
    f.write('\n%s' % asterisks)
    f.write('\nComparison of accepted values from UNESCO 1983 ')
    f.write('\n (Unesco Tech. Paper in Marine Sci. No. 44, p30)')
    f.write('\n%s' % asterisks)
    f.write('\n')

    f.write('\n')
    for irow in range(0, 2):
        f.write('\n   Sal   Press      fp        fp')
        f.write('\n  (psu)   (db)      (C)        (C)\n')
        table = np.vstack([S[irow, :], P[irow, :], UN_fp[irow, :],
                           FP[irow, :]])
        m, n = table.shape
        for iline in range(0, n):
            f.write(" %4.0f   %5.0f   %8.3f  %11.4f\n" %
                    tuple(table[:, iline]))

    # Test main module cp.
    f.write('\n%s' % asterisks)
    f.write('\n**  TESTING MODULE: cp')
    f.write('\n%s' % asterisks)

    # Test 1.
    # Data from Pond and Pickard Intro. Dynamical Oceanography 2nd ed. 1986
    T = np.array([[0,  0,  0,  0,  0,  0],
                  [10, 10, 10, 10, 10, 10],
                  [20, 20, 20, 20, 20, 20],
                  [30, 30, 30, 30, 30, 30],
                  [40, 40, 40, 40, 40, 40]]) / 1.00024

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

    UN_cp = np.array([[4048.4,  3896.3,  3807.7,  3986.5,  3849.3,  3769.1],
                      [4041.8,  3919.6,  3842.3,  3986.3,  3874.7,  3804.4],
                      [4044.8,  3938.6,  3866.7,  3993.9,  3895.0,  3828.3],
                      [4049.1,  3952.0,  3883.0,  4000.7,  3909.2,  3844.3],
                      [4051.2,  3966.1,  3905.9,  4003.5,  3923.9,  3868.3]])

    CP = sw.cp(S, T, P)

    # Display results.
    f.write('\n')
    f.write('\n%s' % asterisks)
    f.write('\nComparison of accepted values from UNESCO 1983 ')
    f.write('\n (Unesco Tech. Paper in Marine Sci. No. 44, p37)')
    f.write('\n%s' % asterisks)
    f.write('\n')

    m, n = S.shape
    f.write('\n')
    for icol in range(0, n):
        f.write('\n   Sal  Temp  Press      Cp        cp')
        f.write('\n  (psu)  (C)   (db)    (J/kg.C)   (J/kg.C)\n')
        result = np.vstack([S[:, icol], T[:, icol], P[:, icol],
                            UN_cp[:, icol], CP[:, icol]])
        for iline in range(0, m):
            f.write(" %4.0f  %4.0f   %5.0f   %8.1f  %11.2f\n" %
                    tuple(result[:, iline]))

    # Test main module svel.
    f.write('\n%s' % asterisks)
    f.write('\n**  TESTING MODULE: svel')
    f.write('\n%s' % asterisks)

    # Test 1.
    # Data from Pond and Pickard Intro. Dynamical Oceanography 2nd ed. 1986
    T = np.array([[0,  0,  0,  0,  0,  0],
                  [10, 10, 10, 10, 10, 10],
                  [20, 20, 20, 20, 20, 20],
                  [30, 30, 30, 30, 30, 30],
                  [40, 40, 40, 40, 40, 40]]) / 1.00024

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

    UN_svel = np.array([[1435.8, 1520.4, 1610.4, 1449.1, 1534.0, 1623.2],
                        [1477.7, 1561.3, 1647.4, 1489.8, 1573.4, 1659.0],
                        [1510.3, 1593.6, 1676.8, 1521.5, 1604.5, 1687.2],
                        [1535.2, 1619.0, 1700.6, 1545.6, 1629.0, 1710.1],
                        [1553.4, 1638.0, 1719.2, 1563.2, 1647.3, 1727.8]])

    SVEL = sw.svel(S, T, P)

    # Display results.
    f.write('\n')
    f.write('\n%s' % asterisks)
    f.write('\nComparison of accepted values from UNESCO 1983 ')
    f.write('\n (Unesco Tech. Paper in Marine Sci. No. 44, p50)')
    f.write('\n%s' % asterisks)
    f.write('\n')

    m, n = SVEL.shape
    f.write('\n')
    for icol in range(0, n):
        f.write('\n   Sal  Temp  Press     SVEL       svel')
        f.write('\n  (psu)  (C)   (db)     (m/s)       (m/s)\n')

        result = np.vstack([S[:, icol], T[:, icol], P[:, icol],
                            UN_svel[:, icol], SVEL[:, icol]])
        for iline in range(0, m):
            f.write(" %4.0f  %4.0f   %5.0f   %8.1f  %11.3f\n" %
                    tuple(result[:, iline]))

    # Test submodules alpha beta aonb.
    f.write('\n%s' % asterisks)
    f.write('\n**  and SUB-MODULE: alpha beta aonb')
    f.write('\n%s' % asterisks)

    # Data from McDougall 1987.
    s = 40
    PT = 10
    p = 4000
    beta_lit = 0.72088e-03
    aonb_lit = 0.34763
    alpha_lit = aonb_lit * beta_lit

    BETA = sw.beta(s, PT, p, pt=True)
    ALPHA = sw.alpha(s, PT, p, pt=True)
    AONB = sw.aonb(s, PT, p, pt=True)

    # Display results.
    f.write('\n')
    f.write('\n%s' % asterisks)
    f.write('\nComparison of accepted values from MCDOUGALL 1987 ')
    f.write('\n%s' % asterisks)
    f.write('\n')
    f.write('\n')

    f.write('\n   Sal  Temp  Press     BETA       beta')
    f.write('\n  (psu)  (C)   (db)   (psu^-1)     (psu^-1)\n')
    table = np.hstack([s, PT, p, beta_lit, BETA])
    f.write(" %4.0f  %4.0f   %5.0f   %11.4e  %11.5e\n" %
            tuple(table))

    f.write('\n   Sal  Temp  Press     AONB       aonb')
    f.write('\n  (psu)  (C)   (db)   (psu C^-1)   (psu C^-1)\n')
    table = np.hstack([s, PT, p, aonb_lit, AONB])
    f.write(" %4.0f  %4.0f   %5.0f   %8.5f  %11.6f\n" %
            tuple(table))

    f.write('\n   Sal  Temp  Press     ALPHA       alpha')
    f.write('\n  (psu)  (C)   (db)    (psu^-1)     (psu^-1)\n')
    table = np.hstack([s, PT, p, alpha_lit, ALPHA])
    f.write(" %4.0f  %4.0f   %5.0f   %11.4e  %11.4e\n" %
            tuple(table))

    # Test main moduleS  satO2 satN2 satAr.
    f.write('\n%s' % asterisks)
    f.write('\n**  TESTING MODULE: satO2 satN2 satAr')
    f.write('\n%s' % asterisks)
    f.write('\n')

    # Data from Weiss 1970.
    T = np.array([[-1, -1],
                  [10, 10],
                  [20, 20],
                  [40, 40]]) / 1.00024

    S = np.array([[20, 40],
                  [20, 40],
                  [20, 40],
                  [20, 40]])

    lit_O2 = np.array([[9.162, 7.984],
                       [6.950, 6.121],
                       [5.644, 5.015],
                       [4.050, 3.656]])

    lit_N2 = np.array([[16.28, 14.01],
                       [12.64, 11.01],
                       [10.47,  9.21],
                       [7.78,  6.95]])

    lit_Ar = np.array([[0.4456, 0.3877],
                       [0.3397, 0.2989],
                       [0.2766, 0.2457],
                       [0.1986, 0.1794]])

    O2 = sw.satO2(S, T)
    N2 = sw.satN2(S, T)
    Ar = sw.satAr(S, T)

    # Display results.
    f.write('\n')
    f.write('\n%s' % asterisks)
    f.write('\nComparison of accepted values from Weiss, R.F. 1979 ')
    f.write('\n"The solubility of nitrogen, oxygen and argon in water')
    f.write('\n and seawater." Deep-Sea Research., 1970, Vol 17, pp721-735.')
    f.write('\n%s' % asterisks)
    f.write('\n')

    m, n = S.shape
    f.write('\n')
    for icol in range(0, n):
        f.write('\n   Sal  Temp      O2         satO2')
        f.write('\n  (psu)  (C)      (ml/l)     (ml/l)\n')
        result = np.vstack([S[:, icol], T[:, icol],
                            lit_O2[:, icol], O2[:, icol]])
        for iline in range(0, m):
            f.write(" %4.0f  %4.0f    %8.2f   %9.3f\n" %
                    tuple(result[:, iline]))

    for icol in range(0, n):
        f.write('\n   Sal  Temp      N2         satN2')
        f.write('\n  (psu)  (C)      (ml/l)     (ml/l)\n')
        result = np.vstack([S[:, icol], T[:, icol],
                            lit_N2[:, icol], N2[:, icol]])
        for iline in range(0, m):
            f.write(" %4.0f  %4.0f    %8.2f  %9.3f\n" %
                    tuple(result[:, iline]))

    for icol in range(0, n):
        f.write('\n   Sal  Temp      Ar         satAr')
        f.write('\n  (psu)  (C)      (ml/l)     (ml/l)\n')
        result = np.vstack([S[:, icol], T[:, icol],
                            lit_Ar[:, icol], Ar[:, icol]])
        for iline in range(0, m):
            f.write(" %4.0f  %4.0f     %8.4f  %9.4f\n" %
                    tuple(result[:, iline]))


if __name__ == '__main__':
    test()
