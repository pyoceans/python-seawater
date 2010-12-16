try:
    import cPickle as pickle
except:
    import pickle

import seawater.gibbs as gsw
import seawater.csiro as sw

""" FIXME: temporary solution until SA is OK"""
import scipy.io as sio
#SA_chck_cast = sio.loadmat('SA_chck_cast.mat', squeeze_me=True)['SA_chck_cast']

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
    print "z_from_p:   Failed"
else:
    print "z_from_p:   Passed"

""" grav """
grav = gsw.grav(gsw_cv.lat_chck_cast, gsw_cv.p_chck_cast )
Igrav = ( (gsw_cv.grav - grav) >= gsw_cv.grav_ca ).nonzero()

if Igrav[0].size != 0:
    print "grav:   Failed"
else:
    print "grav:   Passed"


""" SA_from_SP """
#TODO:
SA_chck_cast = gsw.SA_from_SP(gsw_cv.SP_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)
ISA_from_SP = ( (gsw_cv.SA_from_SP - SA_chck_cast) >= gsw_cv.SA_from_SP_ca).nonzero()

if ISA_from_SP[0].size != 0:
    print "SA_from_SP:   Failed. Note that this will cause many other programmes in the GSW toolbox to fail"
else:
    print "SA_from_SP:   Passed"

""" molality """
molality = gsw.molality(SA_chck_cast)
Imolality = ( (gsw_cv.molality - molality) >= gsw_cv.molality_ca ).nonzero()

if Imolality[0].size != 0:
    print "molality:   Failed"
else:
    print "molality:   Passed"


""" ionic_strength """
ionic_strength = gsw.ionic_strength(SA_chck_cast)
Iionic_strength = ( (gsw_cv.ionic_strength - ionic_strength) >= gsw_cv.ionic_strength_ca ).nonzero()

if Iionic_strength[0].size != 0:
    print "ionic_strength:   Failed"
else:
    print "ionic_strength:   Passed"


""" gsw/sw f """
f = sw.cor(gsw_cv.lat_chck_cast)

If = ( (gsw_cv.f - f) >= gsw_cv.f_ca ).nonzero()

if If[0].size != 0:
    print "f:   Failed"
else:
    print "f:   Passed"


""" CT_from_t """
CT_chck_cast = gsw.CT_from_t(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
ICT_from_t = ( (gsw_cv.CT_from_t - CT_chck_cast) >= gsw_cv.CT_from_t_ca ).nonzero()
if ICT_from_t[0].size != 0:
    print "CT_from_t:   Failed. Note that this will cause many other programmes in the GSW toolbox to fail."
else:
    print "CT_from_t:   Passed"



""" pt_from_t """
pt = gsw.pt_from_t(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast, gsw_cv.pr)
Ipt_from_t = ( (gsw_cv.pt_from_t - pt) >= gsw_cv.pt_from_t_ca ).nonzero()

if Ipt_from_t[0].size != 0:
    print "pt_from_t:   Failed"
else:
    print "pt_from_t:   Passed"


""" entropy_from_t 'pt' """
entropy_from_pt =  gsw.entropy_from_t(SA_chck_cast, pt)
Ientropy_from_pt = ( (gsw_cv.entropy_from_pt - entropy_from_pt) >= gsw_cv.entropy_from_pt_ca ).nonzero()

if Ientropy_from_pt[0].size != 0:
    print "entropy_from_pt:   Failed"
else:
    print "entropy_from_pt:   Passed"


""" entropy_from_t 'CT'"""
entropy_from_CT =  gsw.entropy_from_t(SA_chck_cast, CT_chck_cast, t_type='CT')
Ientropy_from_CT = ( (gsw_cv.entropy_from_CT - entropy_from_CT) >= gsw_cv.entropy_from_CT_ca).nonzero()

if Ientropy_from_CT[0].size != 0:
    print "entropy_from_CT:   Failed"
else:
    print "entropy_from_CT:   Passed"


""" CT_from_pt """
CT_from_pt = gsw.CT_from_pt(SA_chck_cast, pt)
ICT_from_pt = ( (gsw_cv.CT_from_pt - CT_from_pt) >= gsw_cv.CT_from_pt_ca ).nonzero()

if ICT_from_pt[0].size != 0:
    print "CT_from_pt:   Failed"
else:
    print "CT_from_pt:   Passed"


""" pt0_from_t """
pt0 = gsw.pt0_from_t(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
Ipt0 = ( (gsw_cv.pt0 - pt0) >= gsw_cv.pt0_ca ).nonzero()

if Ipt0[0].size != 0:
    print "pt0_from_t:   Failed"
else:
    print "pt0_from_t:   Passed"



""" pt_from_CT """
pt_from_CT = gsw.pt_from_CT(SA_chck_cast, CT_chck_cast)
Ipt_from_CT = ( (gsw_cv.pt - pt_from_CT) >= gsw_cv.pt_ca ).nonzero()

if Ipt_from_CT[0].size != 0:
    print "pt_from_CT:   Failed"
else:
    print "pt_from_CT:   Passed"



""" entropy """
entropy = gsw.entropy(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
Ientropy = ( (gsw_cv.entropy - entropy) >= gsw_cv.entropy_ca ).nonzero()

if Ientropy[0].size != 0:
    print "entropy:   Failed"
else:
    print "entropy:   Passed"

""" pt_from_entropy """
pt_from_entropy =  gsw.t_from_entropy(SA_chck_cast, entropy)
Ipt_from_entropy = ( (gsw_cv.pt_from_entropy - pt_from_entropy) >= gsw_cv.pt_from_entropy_ca ).nonzero()

if Ipt_from_entropy[0].size != 0:
    print "pt_from_entropy:   Failed"
else:
    print "pt_from_entropy:   Passed"


""" CT_from_entropy """
CT_from_entropy =  gsw.t_from_entropy(SA_chck_cast, entropy, 'CT')
ICT_from_entropy = ( (gsw_cv.CT_from_entropy - CT_from_entropy) >= gsw_cv.CT_from_entropy_ca ).nonzero()

if ICT_from_entropy[0].size != 0:
    print "CT_from_entropy:   Failed"
else:
    print "CT_from_entropy:   Passed"