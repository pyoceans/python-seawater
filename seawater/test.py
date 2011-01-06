#TODO: organize in the same order as original version
#TODO: incorporate into gibbs

import seawater.gibbs as gsw
import seawater.csiro as sw
import numpy as np

try:
    import cPickle as pickle
except:
    import pickle

""" load test data """
class Dict2Struc(object):
    """ all the variables from a dict in a "matlab-like-structure" """
    def __init__(self, adict):
        self.__dict__.update(adict)

data = pickle.load( open('gsw_cv.pkl','rb') )
gsw_cv = Dict2Struc(data) # then type dat.<tab> to navigate through your variables

""" SA_from_SP """
SA_chck_cast = gsw.SA_from_SP(gsw_cv.SP_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0]
ISA_from_SP = (gsw_cv.SA_from_SP - SA_chck_cast) >= gsw_cv.SA_from_SP_ca

if ISA_from_SP.any():
    print "SA_from_SP:   Failed. Note that this will cause many other programmes in the GSW toolbox to fail"
else:
    print "SA_from_SP:   Passed"

""" z_from_p """
z_from_p = gsw.z_from_p(gsw_cv.p_chck_cast, gsw_cv.lat_chck_cast)
Iz_from_p = (gsw_cv.z_from_p - z_from_p) >= gsw_cv.z_from_p_ca

if Iz_from_p.any():
    print "z_from_p:   Failed"
else:
    print "z_from_p:   Passed"

""" grav """
grav = gsw.grav(gsw_cv.lat_chck_cast, gsw_cv.p_chck_cast )
Igrav = (gsw_cv.grav - grav) >= gsw_cv.grav_ca

if Igrav.any():
    print "grav:   Failed"
else:
    print "grav:   Passed"

""" gsw/sw f """
f = sw.cor(gsw_cv.lat_chck_cast)

If = (gsw_cv.f - f) >= gsw_cv.f_ca

if If.any():
    print "f:   Failed"
else:
    print "f:   Passed"