#! /usr/bin/env python
# -*- coding: utf-8 -*-
# --- encoding: iso-8859-1 ---

# Python test script for the miscellaneous functions in the seawater module
#
# This script tries to recreate miscellaneous tables from pp 28, 30, and 50
# in the UNESCO 1983 report. The functions tested are depth, freezept,
# and soundvel

# Reference:
# UNESCO 1983
# N.P. Fofonoff and R.C. Millard Jr.,
# Algorithms for computation of fundamental properties of seawater,
# Unesco technical papers in marine science 44.

from seawater import depth, fp, svel
from numpy import array, arange

deg = 'deg'  # degree symbol in iso-latin1 encoding ?

def printdepth():
    LAT = array([0., 30., 45., 60., 90.])
    P   = arange(0, 10001, 1000)
    P[0] = 500

    print 
    print "                TABLE OF DEPTH (METERS)"
    print
    print "                   LATITUDE (DEGREES)"
    print "PRESSURE"
    print "DECIBARS", ("%7.0f" + 4*"%11.0f") % tuple(LAT)
    print 
    entry = "    %7.2f"
    for p in P:
        print "%6.0f" % p, 5*entry % tuple(depth(p, LAT))

def printfreezept():
    S = arange(5, 41, 5)
    P = arange(0, 501, 100)
    
    print 
    print "          FREEZING POINT TEMPERATURE " + deg + "C"
    print
    print "                   SALINITY PSS-78"
    print
    print "PRESSURE", 8*"%8.0f" % tuple(S)
    print "DECIBARS"
    entry =  "  %6.3f"
    for p in P:
        print "%9.0f" % p, 8*entry % tuple(fp(S, p))

def printsoundvel():
    T = 10.0*arange(5)
    P = 1000.0*arange(11)
    S = array([25., 30., 35., 40.])
    entry = "      %6.1f"
    print 
    print "                               SOUND SPEED IN SEAWATER U [m/s]"
    print "PRESSURE                     TEMPERATURE",
    print deg+"C  IPTS-68                 SALINTY:", "%2.0f" % S[0]
    print "DECIBARS %8d %11d %11d %11d %11d" % tuple(T)
    for p in P:
        print "%6.0f" % p, 5*entry % tuple(svel(S[0],T,p))
    for s in S[1:]:
        print 68*"-" + "SALINTIY:", "%2.0f" % s
        for p in P:
            print "%6.0f" % p, 5*entry % tuple(svel(s,T,p))

printdepth()
raw_input("press any key to continue:")
printfreezept()
raw_input("press any key to continue:")
printsoundvel()