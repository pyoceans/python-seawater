#TODO: organize in the same order as original version

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


""" entropy_from_t 'pt' """ #FIXME: pass, but small float are detected, investigate further
entropy_from_pt =  gsw.entropy_from_t(SA_chck_cast, pt)
Ientropy_from_pt = (gsw_cv.entropy_from_pt - entropy_from_pt) >= gsw_cv.entropy_from_pt_ca

if Ientropy_from_pt.any():
    print "entropy_from_pt:   Failed"
else:
    print "entropy_from_pt:   Passed"


""" entropy_from_t 'CT'"""
entropy_from_CT =  gsw.entropy_from_t(SA_chck_cast, CT_chck_cast, t_type='CT')
Ientropy_from_CT = (gsw_cv.entropy_from_CT - entropy_from_CT) >= gsw_cv.entropy_from_CT_ca

if Ientropy_from_CT.any():
    print "entropy_from_CT:   Failed"
else:
    print "entropy_from_CT:   Passed"


""" CT_from_pt """
CT_from_pt = gsw.CT_from_pt(SA_chck_cast, pt) #FIXME: pass, but small float are detected, investigate further
ICT_from_pt = (gsw_cv.CT_from_pt - CT_from_pt) >= gsw_cv.CT_from_pt_ca

if ICT_from_pt.any():
    print "CT_from_pt:   Failed"
else:
    print "CT_from_pt:   Passed"


""" pt_from_CT """
pt_from_CT = gsw.pt_from_CT(SA_chck_cast, CT_chck_cast)
Ipt_from_CT = (gsw_cv.pt - pt_from_CT) >= gsw_cv.pt_ca

if Ipt_from_CT.any():
    print "pt_from_CT:   Failed"
else:
    print "pt_from_CT:   Passed"


""" CT_from_entropy """ #FIXME: pass, but small float are detected, investigate further
CT_from_entropy =  gsw.t_from_entropy(SA_chck_cast, entropy, 'CT')
ICT_from_entropy = (gsw_cv.CT_from_entropy - CT_from_entropy) >= gsw_cv.CT_from_entropy_ca

if ICT_from_entropy.any():
    print "CT_from_entropy:   Failed"
else:
    print "CT_from_entropy:   Passed"


""" sigma_CT """
sigma0_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast)
Isigma_CT = (gsw_cv.sigma0_CT - sigma0_CT) >= gsw_cv.sigma0_CT_ca

if Isigma_CT.any():
    print "sigma_CT at 0 db:   Failed"
else:
    print "sigma_CT at 0 db:   Passed"

sigma1_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast, 1000)
Isigma_CT = (gsw_cv.sigma1_CT - sigma1_CT) >= gsw_cv.sigma1_CT_ca

if Isigma_CT.any():
    print "sigma_CT at 1000 db:   Failed"
else:
    print "sigma_CT at 1000 db:   Passed"

sigma2_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast, 2000)
Isigma_CT = (gsw_cv.sigma2_CT - sigma2_CT) >= gsw_cv.sigma2_CT_ca

if Isigma_CT.any():
    print "sigma_CT at 2000 db:   Failed"
else:
    print "sigma_CT at 2000 db:   Passed"

sigma3_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast, 3000)
Isigma_CT = (gsw_cv.sigma3_CT - sigma3_CT) >= gsw_cv.sigma3_CT_ca

if Isigma_CT.any():
    print "sigma_CT at 3000 db:   Failed"
else:
    print "sigma_CT at 3000 db:   Passed"

sigma4_CT = gsw.sigma_CT(SA_chck_cast, CT_chck_cast, 4000)
Isigma_CT = (gsw_cv.sigma4_CT - sigma4_CT) >= gsw_cv.sigma4_CT_ca

if Isigma_CT.any():
    print "sigma_CT at 4000 db:   Failed"
else:
    print "sigma_CT at 4000 db:   Passed"

""" enthalpy """
enthalpy = gsw.enthalpy(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
Ienthalpy = (gsw_cv.enthalpy - enthalpy) >= gsw_cv.enthalpy_ca

if Ienthalpy.any():
    print "enthalpy:   Failed"
else:
    print "enthalpy:   Passed"

""" t_from_CT """ #FIXME: pass, but small float are detected, investigate further
t_from_CT =  gsw.t_from_CT(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
It_from_CT = (gsw_cv.t_chck_cast - t_from_CT) >= gsw_cv.t_from_CT_ca

if It_from_CT.any():
    print "t_from_CT:   Failed"
else:
    print "t_from_CT:   Passed"

""" enthalpy_CT """
enthalpy_CT =  gsw.enthalpy(SA_chck_cast,CT_chck_cast, gsw_cv.p_chck_cast, t_type='CT')
Ienthalpy_CT = (gsw_cv.enthalpy_CT - enthalpy_CT) >= gsw_cv.enthalpy_CT_ca

if Ienthalpy_CT.any():
    print "enthalpy_CT:   Failed"
else:
    print "enthalpy_CT:   Passed"

""" enthalpy_CT25 """ #FIXME: pass, but small float are detected, investigate further
enthalpy_CT25 =  gsw.enthalpy(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast, t_type='CT', term25=True)[0]
Ienthalpy_CT25 = (gsw_cv.enthalpy_CT25 - enthalpy_CT25) >= gsw_cv.enthalpy_CT25_ca

if Ienthalpy_CT25.any():
    print "enthalpy_CT25:   Failed"
else:
    print "enthalpy_CT25:   Passed"