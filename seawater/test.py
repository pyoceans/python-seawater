try:
    import cPickle as pickle
except:
    import pickle

import seawater.gibbs as gsw

""" FIXME: temporary solution until SA is OK"""
import scipy.io as sio
SA_chck_cast = sio.loadmat('SA_chck_cast.mat', squeeze_me=True)['SA_chck_cast']

""" load test data """
class Dict2Struc(object):
    """ all the variables from a dict in a "matlab-like-structure" """
    def __init__(self, adict):
        self.__dict__.update(adict)

data = pickle.load( open('gsw_cv.pkl','rb') )
gsw_cv = Dict2Struc(data) # then type dat.<tab> to navigate through your variables

""" z_from_p """
z_from_p = gsw.z_from_p(gsw_cv.p_chck_cast, gsw_cv.lat_chck_cast)
Iz_from_p = ( (gsw_cv.z_from_p - z_from_p) >= gsw_cv.z_from_p_ca ).nonzero()

if Iz_from_p[0].size != 0:
    print "gsw.z_from_p:   Failed\n"
else:
    print "gsw.z_from_p:   Passed\n"

""" grav """
grav = gsw.grav(gsw_cv.lat_chck_cast, gsw_cv.p_chck_cast )
Igrav = ( (gsw_cv.grav - grav) >= gsw_cv.grav_ca ).nonzero()

if Igrav[0].size != 0:
    print "gsw.grav:   Failed\n"
else:
    print "gsw.grav:   Passed\n"

""" molality """
molality = gsw.molality(SA_chck_cast)
Imolality = ( (gsw_cv.molality - molality) >= gsw_cv.molality_ca ).nonzero()

if Imolality[0].size != 0:
    print "gsw_molality:   Failed\n"
else:
    print "gsw_molality:   Passed\n"


""" ionic_strength """
ionic_strength = gsw.ionic_strength(SA_chck_cast)
Iionic_strength = ( (gsw_cv.ionic_strength - ionic_strength) >= gsw_cv.ionic_strength_ca ).nonzero()

if Iionic_strength[0].size != 0:
    print "gsw_ionic_strength:   Failed\n"
else:
    print "gsw_ionic_strength:   Passed\n"