# -*- coding: utf-8 -*-
#
# test_result_comparison.py
#
# purpose:  Compare a dataset from python vs matlab (run on octave).
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.github.io/
# created:  05-Aug-2013
# modified: Wed 21 Aug 2013 05:43:57 PM BRT
#
# obs:
#


from __future__ import division

import os
import unittest

import numpy as np
from oct2py import octave
#from pandas import read_csv

import seawater as sw
from seawater.constants import c3515

path = os.path.join(os.path.dirname(__file__), 'seawater_v3_3')
_ = octave.addpath(octave.genpath(path))


def compare_results(name, function, args, values):
    args = [values.get(arg) for arg in args]
    res = function(*args)  # Python call.
    # FIXME: Test only the first output when multiple outputs are present.
    if isinstance(res, tuple):
        nout = len(res)
        res = res[0]
    else:
        nout = 1
    val = octave.call('sw_%s' % name, *args, nout=nout)  # Octave call.
    if nout > 1:
        val = val[0]
    val, res = val.squeeze(), res.squeeze()
    np.testing.assert_allclose(val, res)


class OctaveResultComparison(unittest.TestCase):
    def setUp(self):
        # TODO: More tests with station data.
        #kw = dict(comment='#', header=5, index_col=0)
        #st61 = read_csv('Endeavor_Cruise-88_Station-61.csv', **kw)
        #st64 = read_csv('Endeavor_Cruise-88_Station-64.csv', **kw)
        kw = dict(comments='#', skiprows=6, delimiter=',')
        st61 = np.loadtxt('Endeavor_Cruise-88_Station-61.csv', **kw)
        st64 = np.loadtxt('Endeavor_Cruise-88_Station-64.csv', **kw)

        latst = 36. + 40.03 / 60., 37. + 39.93 / 60.
        lonst = -(70. + 59.59 / 60.), -71.
        Sal = np.c_[st61[:, 2], st64[:, 2]]
        Temp = np.c_[st61[:, 1], st64[:, 1]]
        Pres = np.c_[st61[:, 0], st61[:, 0]]
        Gpan = sw.gpan(Sal, Temp, Pres)

        self.values = dict(r=np.array([56.4125, 56.3161, 50.6703, 38.1345,
                                       35.0565, 32.9865]) / c3515,
                           s=np.array([34.5487, 34.7275, 34.8605, 34.6810,
                                       34.5680, 34.5600]),
                           t=np.array([28.7856, 28.4329, 22.8103, 10.2600,
                                       6.8863, 4.4036]),
                           p=np.array([10., 50., 125., 250., 600., 1000.]),
                           pref=np.array([0.0]),
                           pt=np.array([28.8099, 28.4392, 22.7862, 10.2262,
                                        6.8272, 4.3236]),
                           rtx=np.array([0.99353194]),
                           delt=np.array([13.7856]),
                           rt=np.array([1.32968079, 1.32094651, 1.18368907,
                                        0.89332541, 0.81977076, 0.76703445]),
                           units='km',
                           lon=np.array([-30.] * 6),
                           lat=np.linspace(-22., -21., 6.),
                           length=np.array([100.] * 6),
                           latst=latst, lonst=lonst, Sal=Sal, Temp=Temp,
                           Pres=Pres, Gpan=Gpan)

    def tearDown(self):
        unittest.TestCase.tearDown(self)

    # library.py  # NOTE: Matlab does not have version of T68conv and T90conv.
    def test_cndr(self):
        name = 'cndr'
        function, args = (sw.cndr, ('s', 't', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_salds(self):
        name = 'salds'
        function, args = (sw.salds, ('rtx', 'delt'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_salrp(self):
        name = 'salrp'
        function, args = (sw.salrp, ('r', 't', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_salrt(self):
        name = 'salrt'
        function, args = (sw.salrt, ('t',))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_seck(self):
        name = 'seck'
        function, args = (sw.seck, ('s', 't', 'pref'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_sals(self):
        name = 'sals'
        function, args = (sw.sals, ('rt', 't'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    # extras.py
    def test_dist(self):
        name = 'dist'
        function, args = (sw.dist, ('lat', 'lon', 'units'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_f(self):
        name = 'f'
        function, args = (sw.f, ('lat',))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_satAr(self):
        name = 'satAr'
        function, args = (sw.satAr, ('s', 't'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_satN2(self):
        name = 'satN2'
        function, args = (sw.satN2, ('s', 't'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_satO2(self):
        name = 'satO2'
        function, args = (sw.satO2, ('s', 't'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_swvel(self):
        name = 'swvel'
        function, args = (sw.swvel, ('length', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    # eos80.py
    def test_adtg(self):
        name = 'adtg'
        function, args = (sw.adtg, ('s', 't', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_alpha(self):
        name = 'alpha'
        function, args = (sw.alpha, ('s', 't', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_aonb(self):
        name = 'aonb'
        function, args = (sw.aonb, ('s', 't', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_beta(self):
        name = 'beta'
        function, args = (sw.beta, ('s', 't', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_bfrq1d(self):
        name = 'bfrq'
        function, args = (sw.bfrq, ('s', 't', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_dpth(self):
        name = 'dpth'
        function, args = (sw.dpth, ('p', 'lat'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_g(self):
        name = 'g'
        function, args = (sw.g, ('lat',))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_salt(self):
        name = 'salt'
        function, args = (sw.salt, ('r', 't', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_fp(self):
        name = 'fp'
        function, args = (sw.fp, ('s', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_svel(self):
        name = 'svel'
        function, args = (sw.svel, ('s', 't', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_pres(self):
        name = 'pres'
        function, args = (sw.pres, ('p', 'lat'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_dens0(self):
        name = 'dens0'
        function, args = (sw.dens0, ('s', 't'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_smow(self):
        name = 'smow'
        function, args = (sw.smow, ('t'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_dens(self):
        name = 'dens'
        function, args = (sw.dens, ('s', 't', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_pden(self):
        name = 'pden'
        function, args = (sw.pden, ('s', 't', 'p', 'pref'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_cp(self):
        name = 'cp'
        function, args = (sw.cp, ('s', 't', 'p'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_ptmp(self):
        name = 'ptmp'
        function, args = (sw.ptmp, ('s', 't', 'p', 'pref'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_temp(self):
        name = 'temp'
        function, args = (sw.temp, ('s', 'pt', 'p', 'pref'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    # geostrophic.py
    def test_svan(self):
        name = 'svan'
        function, args = (sw.svan, ('s', 't', 'pref'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_bfrq(self):
        name = 'bfrq'
        function, args = (sw.bfrq, ('Sal', 'Temp', 'Pres'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_gpan(self):
        name = 'gpan'
        function, args = (sw.gpan, ('Sal', 'Temp', 'Pres'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)

    def test_gvel(self):
        name = 'gvel'
        function, args = (sw.gvel, ('Gpan', 'latst', 'lonst'))
        compare_results(name=name, function=function, args=args,
                        values=self.values)


if __name__ == '__main__':
    unittest.main()
