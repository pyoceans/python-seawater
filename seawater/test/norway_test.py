#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Simple test script for the seawater package

Test the functions against standard values in the
UNESCO 1983 report.

UNESCO 1983 
N.P. Fofonoff and R.C. Millard Jr.,
Algorithms for computation of fundamental properties of seawater,
Unesco technical papers in marine science, 44.
"""

from seawater import *

format1 = "Computed: %25s = "
format2 = "Check value                         = "

print
# Check value from UNESCO 1983, p. 20
print "Checking svan"
print
print "S = 40, T = 40 degC, P = 10000 dbar"
print format1 % "svan(40, 40, 10000)", svan(40, 40, 10000)
print format2, "981.30210E-8"

print
# Check value from UNESCO 1983, p. 20
print "Checking dens"
print
print "S = 40, T = 40 degC, P = 10000 dbar"
print format1 % "dens(40, 40, 10000)", dens(40, 40, 10000)
print format2, 1059.82037

print
# Check value from UNESCO 1983, p. 11
print "Checking salt"
print
print "Salinity = 40.0000"
print "cond = 1.888091, T = 40 degC, P = 10000 dbar"
print format1 % "salt(1.888091, 40, 10000)", salt(1.888091, 40, 10000)
print format2, 40.0000

print
# Check value from UNESCO 1983, p. 11
print "Checking cndr"
print
print "S = 40, T = 40 degC, P = 10000 dbar"
print format1 % "cndr(40, 40, 10000)", cndr(40, 40, 10000)
print format2, 1.888091

print
# Check value from UNESCO 1983, p. 35
print "Checking cp"
print
print "S = 40, T = 40 degC, P = 10000 dbar"
print format1 % "cp(40, 40, 10000)", cp(40, 40, 10000)
print format2, "3849.500" 

print
# Check value from UNESCO 1983, p. 36
print "Checking adtg"
print
print "S = 40, T = 40 degC, P = 10000 dbar"
print format1 % "adtg(40, 40, 10000)", adtg(40, 40, 10000)
print format2, "3.255976E-4"

print
# Check value from UNESCO 1983, p. 44
print "Checking ptmp"
print 
print "S = 40, T = 40 degC, P = 10000 dbar, Pref = 0"
print format1 % "ptmp(40, 40, 10000)", ptmp(40, 40, 10000)
print format2, 36.89073

print
# Check value from UNESCO 1983, p. 30
print "Checking fp"
print
print "S = 40, p = 500 dbar"
print format1 % "fp(40, 500)", fp(40, 500)
print format2, -2.588567

print
# Check value from UNESCO 1983, p. 49
print "Checking svel"
print 
print "S = 40, T = 40 degC,  P = 10000 dbar"
print format1 % "svel(40, 40, 10000)", svel(40, 40, 10000)
print format2, 1731.995

print
# Check value from UNESCO 1983, p. 28
print "Checking depth"
print
print "P = 10000 dbar, latitude = 30 degrees"
print format1 % "depth(10000, 30)", depth(10000, 30)
print format2, 9712.653