# -*- coding: utf-8 -*-
#
# test_result_comparison.py
#
# purpose:  Compare a dataset from python vs matlab (run on octave).
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.github.io/
# created:  05-Aug-2013
# modified: Mon 19 Aug 2013 08:39:17 PM BRT
#
# obs:
#


from __future__ import division

import os
import unittest

import numpy as np
from oct2py import octave
from pandas import read_csv

import seawater as sw
from seawater.constants import c3515

path = os.path.join(os.path.dirname(__file__), 'seawater_v3_3')
path = './seawater_v3_3'
_ = octave.addpath(octave.genpath(path))

kw = dict(comment='#', header=5, index_col=0)
st61 = read_csv('Endeavor_Cruise-88_Station-61.csv', **kw)
st64 = read_csv('Endeavor_Cruise-88_Station-64.csv', **kw)
latst = 36. + 40.03 / 60., 37. + 39.93 / 60.
lonst = -(70. + 59.59 / 60.), -71.
Sal = np.c_[st61['S'].values, st64['S'].values]
Temp = np.c_[st61['t'].values, st64['t'].values]
Pres = np.c_[st61.index.values.astype(float),
                st64.index.values.astype(float)]
Gpan = sw.gpan(Sal, Temp, Pres)

values = dict(r=np.array([56.4125, 56.3161, 50.6703, 38.1345, 35.0565,
                          32.9865]) / c3515,
              s=np.array([34.5487, 34.7275, 34.8605, 34.6810, 34.5680,
                          34.5600]),
              t=np.array([28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]),
              p=np.array([10., 50., 125., 250., 600., 1000.]),
              pref=np.array([0.0]),
              pt=np.array([28.8099, 28.4392, 22.7862, 10.2262, 6.8272, 4.3236]),
              rtx=np.array([0.99353194]),
              delt=np.array([13.7856]),
              rt=np.array([1.32968079, 1.32094651, 1.18368907, 0.89332541,
                           0.81977076, 0.76703445]),
              units='km',
              lon=np.array([-30.] * 6),
              lat=np.linspace(-22., -21., 6.),
              length=np.array([100.] * 6),
              latst=latst,
              lonst=lonst,
              Sal=Sal,
              Temp=Temp,
              Pres=Pres,
              Gpan=Gpan)

def compare_results(name, function, args):
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
    np.testing.assert_equal(val, res)


# All seawater functions to test.
functions = dict({
    # library.py  # NOTE: Matlab does not have version of T68conv and T90conv.
    'cndr': (sw.cndr, ('s', 't', 'p')),  # FIXME: Diff. at 5th decimal place.
    'salds' : (sw.salds, ('rtx', 'delt')),  # Matlab version not vectorized.
    'salrp': (sw.salrp, ('r', 't', 'p')),
    'salrt': (sw.salrt, ('t',)),
    'seck': (sw.seck, ('s', 't', 'pref')),
    'sals' : (sw.sals, ('rt', 't')),
    # extras.py
    'dist': (sw.dist, ('lat', 'lon', 'units')),
    'f': (sw.f, ('lat',)),
    'satAr': (sw.satAr, ('s', 't')),
    'satN2': (sw.satN2, ('s', 't')),
    'satO2': (sw.satO2, ('s', 't')),
    'swvel': (sw.swvel, ('length', 'p')),
    'adtg': (sw.adtg, ('s', 't', 'p')),
    'alpha': (sw.alpha, ('s', 't', 'p')),
    'aonb': (sw.aonb, ('s', 't', 'p')),
    'beta': (sw.beta, ('s', 't', 'p')),
    'bfrq': (sw.bfrq, ('s', 't', 'p')),
    'dpth': (sw.dpth, ('p', 'lat')),
    'g': (sw.g, ('lat',)),
    'salt': (sw.salt, ('r', 't', 'p')),
    'fp': (sw.fp, ('s', 'p')),
    'svel': (sw.svel, ('s', 't', 'p')),
    'pres': (sw.pres, ('p', 'lat')),
    'dens0': (sw.dens0, ('s', 't')),
    'smow': (sw.smow, ('t')),
    'dens': (sw.dens, ('s', 't', 'p')),
    'pden': (sw.pden, ('s', 't', 'p', 'pref')),
    'cp': (sw.cp, ('s', 't', 'p')),
    'ptmp': (sw.ptmp, ('s', 't', 'p', 'pref')),
    'temp': (sw.temp, ('s', 'pt', 'p', 'pref')),
    # geostrophic.py
    'svan': (sw.svan, ('s', 't', 'pref')),
    'bfrq': (sw.bfrq, ('Sal', 'Temp', 'Pres')),
    'gpan': (sw.gpan, ('Sal', 'Temp', 'Pres')),
    'gvel': (sw.gvel, ('Gpan', 'latst', 'lonst')),
    })

class OctaveResultComparison(unittest.TestCase):
    def setUp(self):
        self.values = values
        self.functions = functions

    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def test_all_functions(self):
        for name, (function, args) in self.functions.iteritems():
            compare_results(name=name, function=function, args=args)

if __name__ == '__main__':
    unittest.main()
