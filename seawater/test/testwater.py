#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Python test script for most functions in the water module
#
# The script tries to recreate table A3.1 in Gill's book
# Most of the salinity deratives are missing
# (since they are not implemented in the water module yet).
#
# The functions fp and beta are not tested by this script,
# while some functions are only tested indirectly as they are
# used by other functions.
#
# The results for alpha, the termal expansion coefficient are
# not totally identical to the table.
#
# Reference:
# Adrian E. Gill, Atmosphere-Ocean Dynamics
# Academic Press, 1982

from seawater import dens, alpha, cp, ptmp, svel #drhods

def values(S, T, P):
    '''Compute a row in the table'''
    #return (P/10, S, T, dens(S,T,P)-1000,
            #drhods(S,T,P), 10**7*alpha(S,T,P),  cp(S,T,P),
            #1000*ptmp( S, T, P ), svel( S, T, P ) )
    return ( P/10, S, T, dens( S, T, P ) - 1000,
            10**7 * alpha( S, T, P ),  cp( S, T, P ),
            1000 * ptmp( S, T, P ), svel( S, T, P ) )

# OUTPUT format for the rows
format = '%4.0f %4.0f %4.0f %7.3f %6.0f %6.0f %6.0f %8.1f'

S = 35

P = 0
for T in [-2, 0, 2, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31]:
    print format % values( S, T, P )

P = 1000
for T in [-2, 0, 2, 4, 7, 10, 13, 16, 19]:
    print format % values( S, T, P )

P = 2000
for T in [-2, 0, 2, 4, 7]:
    print format % values( S, T, P )

P = 3000
for T in [-2, 0, 2, 4]:
    print format % values( S, T, P )

P = 4000
for T in [-2, 0, 2, 4]:
    print format % values( S, T, P )

P = 5000
for T in [-2, 0, 2]:
    print format % values( S, T, P )

P = 6000
for T in [-2, 0, 2]:
    print format % values( S, T, P )
    