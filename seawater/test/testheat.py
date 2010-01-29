#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Python test script for the heat functions in the seawater module
#
# This script tries to recreate the heat tables from pp 37, 41, and 45
# in the UNESCO 1983 report. The functions tested are heatcap, adtgrad,
# and temppot.

# Reference:
# UNESCO 1983
# N.P. Fofonoff and R.C. Millard Jr.,
# Algorithms for computation of fundamental properties of seawater,
# Unesco technical papers in marine science 44.

from numpy import array, arange
from seawater import cp, adtg, ptmp

T = 10.0*arange(5)
P = 1000.0*arange(11)
S = array([0., 30., 35., 40.])

deg = 'deg'  # degree symbol in iso-latin1 encoding ?

def printheatcap():
    S = array([0., 30., 35., 40.])
    entry = "      %6.1f"
    print 
    print "                     SPECIFIC HEAT SEAWATER Cp [ J/(Kg K) ]"
    print "PRESSURE                     TEMPERATURE",
    print deg+"C  IPTS-68                 SALINTY:", "%2.0f" % S[0]
    print "DECIBARS %8d %11d %11d %11d %11d" % tuple(T)
    for p in P:
        print "%6.0f" % p, 5*entry % tuple(cp(S[0],T,p))
    for s in S[1:]:
        print 68*"-" + "SALINTIY:", "%2.0f" % s
        for p in P:
            print "%6.0f" % p, 5*entry % tuple(cp(s,T,p))

def printadtgrad():
    S = array([25., 30., 35., 40.])
    entry = "      %6.4f"
    print 
    print "                     ADIABATIC LAPSE RATE  [" \
          + deg + "C/1000 Decibars]"
    print "PRESSURE                     TEMPERATURE",
    print deg+"C  IPTS-68                 SALINTY:", "%2.0f" % S[0]
    print "DECIBARS %8d %11d %11d %11d %11d" % tuple(T)
    for p in P:
        print "%6.0f" % p, 5*entry % tuple(1000*adtg(S[0],T,p))
    for s in S[1:]:
        print 68*"-" + "SALINTIY:", "%2.0f" % s
        for p in P:
            print "%6.0f" % p, 5*entry % tuple(1000*adtg(s,T,p))

def printtemppot():
    S = array([25., 30., 35., 40.])
    entry = "     %7.4f"
    print 
    print "        POTENTIAL TEMPERATURE " + deg + \
          "C (Ref. Pres. = 0.0)"
    print "PRESSURE                     TEMPERATURE",
    print deg+"C  IPTS-68                 SALINTY:", "%2.0f" % S[0]
    print "DECIBARS %8d %11d %11d %11d %11d" % tuple(T)
    for p in P:
        print "%6.0f" % p, 5*entry % tuple(ptmp(S[0],T,p))
    for s in S[1:]:
        print 68*"-" + "SALINTIY:", "%2.0f" % s
        for p in P:
            print "%6.0f" % p, 5*entry % tuple(ptmp(s,T,p))

printheatcap()
raw_input("press a key to continue:")
printadtgrad()
raw_input("press a key to continue:")
printtemppot()