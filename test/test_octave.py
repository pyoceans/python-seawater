# -*- coding: utf-8 -*-
#
# test_octave.py
#
# purpose:  Compare a dataset from python vs matlab (run on octave).
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.github.io/
# created:  05-Aug-2013
# modified: Mon 05 Aug 2013 09:57:15 AM BRT
#
# obs:
#


from __future__ import division

import os
import sys
from collections import OrderedDict

import numpy as np
from oct2py import octave
from pandas import read_csv

import seawater as sw
from seawater.constants import c3515

try:
    path = sys.argv[1]
except IndexError:
    path = "./seawater_v3_3"

if not os.path.exists(path):
    raise ValueError("matlab seawater path %s not found" % path)

_ = octave.addpath(octave.genpath(path))

kw = dict(comment='#', header=5, index_col=0)
st61 = read_csv('Endeavor_Cruise-88_Station-61.csv', **kw)
st64 = read_csv('Endeavor_Cruise-88_Station-64.csv', **kw)
latst = 36. + 40.03 / 60., 37. + 39.93 / 60.
lonst = -(70. + 59.59 / 60.), -71.

Sal=np.c_[st61['S'].values, st64['S'].values]
Temp=np.c_[st61['t'].values, st64['t'].values]
Pres=np.c_[st61.index.values.astype(float), st64.index.values.astype(float)]
Gpan = sw.gpan(Sal, Temp, Pres)

def compare_results(name, function, args):
    args = [values.get(arg) for arg in args]

    try:  # Python.
        res = function(*args)
    except:
        print('%s: python runtime error' % name)
        raise
        return 'no_python'

    # FIXME: Testing only the first output when multiple outputs are present.
    nout = 1
    if isinstance(res, tuple):
        nout = len(res)
        res = res[0]

    try:  # Octave.
        val = octave.call('sw_%s' % name, *args, verbose=True, nout=nout)
        if nout > 1:
            val = val[0]
    except Oct2PyError:
        print('%s: Octave runtime error' % name)
        print("python:\n%s" % res)
        return 'no_octave'

    val = val.squeeze()
    res = res.squeeze()

    try:
        perfect = (val == res).all()
    except:
        print('%s: Comparison failed' % name)
        print("octave:\n%s" % val)
        print("python:\n%s" % res)
        return 'no_comparison'
    if np.allclose(val, res, rtol=1e-15, atol=0):
        print('%s: Passed' % name)
        return 'passed'
    else:
        print('%s: Failed' % name)
        print("octave:\n%s" % val)
        print("python:\n%s" % res)
        return 'failed'
    print('')


values = dict(r=np.array([56.4125, 56.3161, 50.6703, 38.1345, 35.0565,
                          32.9865]) / c3515,
              s=np.array([34.5487, 34.7275, 34.8605, 34.6810, 34.5680,
                           34.5600]),
              t=np.array([28.7856, 28.4329, 22.8103, 10.2600, 6.8863, 4.4036]),
              p=np.array([10., 50., 125., 250., 600., 1000.]),
              pref=np.array([0.0]),
              pt=np.array([28.8099, 28.4392, 22.7862, 10.2262, 6.8272,
                           4.3236]),
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
              Gpan=Gpan
              )

# Functions.
library = OrderedDict({
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

if __name__ == '__main__':
    outcomes = ['passed', 'no_octave', 'no_python', 'failed', 'no_comparison']
    results = dict([(k, list()) for k in outcomes])

    for name, (function, args) in library.iteritems():
        ret = compare_results(name=name, function=function, args=args)
        results[ret].append(name)

    print('\nSummary:')
    print('passed:\n  %s' % '\n  '.join(results['passed']))
    print('octave call failed:\n  %s' % '\n  '.join(results['no_octave']))
    print('python call failed:\n  %s' % '\n  '.join(results['no_python']))
    print('results did not match:\n  %s' % '\n  '.join(results['failed']))
    print('comparison failed:\n  %s' % '\n  '.join(results['no_comparison']))
    print('')
