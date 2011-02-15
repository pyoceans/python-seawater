""" test data: data/gsw_cv.npz
(1) a reference cast (for the McD_Klocker streamfunction)
(2) three vertical profiles of (SP, t, p) at known long & lat
(3) the outputs of all the GSW functions for these 3 profiles
(4) the required accuracy of all these outputs.
"""

import numpy as np
import seawater.gibbs as gsw
import os

class Dict2Struc(object):
    """
    Open variables from a dictionary in a "matlab-like-structure"
    """
    def __init__(self, adict):
        self.__dict__.update(adict)

def read_data(fname):
    datadir = os.path.join(os.path.dirname(__file__), '../data')
    d = np.load(os.path.join(datadir, fname))
    return Dict2Struc(d)

def test_print(func, comp=None):
    """
    Run a function test mimicking the original logic. This is done to allow for
    a direct comparison of the result from the Matlab to the python package.
    """
    width = 23
    if comp is None:
        comp = func

    # test if check value is identical to computed value
    if eval( "( gsw_cv." +comp+ "[~np.isnan(gsw_cv."+comp+")] == " +func+
                                        "[~np.isnan("+func+")] ).all()" ):
        print "%s: Passed" % func.rjust(width)
    else:
        # floating diffs with: computed - check value >= defined precision
        try:
            exec( "unequal = (gsw_cv." +comp+ " - " +func+
                              " ) >= gsw_cv." +comp+ "_ca" )
        except:
            exec( "unequal = (gsw_cv." +comp+ " - " +func+
                              " ) >= gsw_cv." +func+ "_ca" )
        if unequal.any():
            print "%s: Failed" % func.rjust(width)
        else:
            # for term25 and small float differences that will appear
            exec("nmax = ( (gsw_cv."+comp+" - "+func+
                            ")[~np.isnan(gsw_cv."+comp+")]).max()")
            exec("nmin = ( (gsw_cv."+comp+" - "+func+
                            ")[~np.isnan(gsw_cv."+comp+")]).min()")
            print "%s: Passed, diff range from: %s to %s" % ( func.rjust(width), nmax, nmin)

# load test data
gsw_cv = read_data("gsw_cv.npz")

# longitude is stored as uint8
gsw_cv.long_chck_cast = np.float64(gsw_cv.long_chck_cast)

# Test Section
#FIXME: _delta_SA should accept 2-D or make the flatten() inside so the user can use 2D

""" Absolute Salinity (SA) and Preformed Salinity (Sstar) """
SP_chck_cast = gsw_cv.SP_chck_cast.flatten()
p_chck_cast = gsw_cv.p_chck_cast.flatten()
long_chck_cast = lon = [ 142.,  183.,   20.]*45
lat_chck_cast = lat = [ 11. ,   9.5,  59. ]*45

SA_chck_cast = gsw.SA_from_SP(SP_chck_cast, p_chck_cast, long_chck_cast, lat_chck_cast)
#SA_chck_cast = gsw.SA_from_SP(gsw_cv.SP_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)
test_print("SA_chck_cast", "SA_from_SP")

#Sstar_from_SP = gsw.Sstar_from_SP(gsw_cv.SP_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)
Sstar_from_SP = gsw.Sstar_from_SP(SP_chck_cast, p_chck_cast, long_chck_cast, lat_chck_cast)
test_print("Sstar_from_SP")

#SA_SA_Sstar_from_SP, Sstar_SA_Sstar_from_SP = gsw.SA_Sstar_from_SP(gsw_cv.SP_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)[0:2]
SA_SA_Sstar_from_SP, Sstar_SA_Sstar_from_SP = gsw.SA_Sstar_from_SP(SP_chck_cast, p_chck_cast, long_chck_cast, lat_chck_cast)[0:2]
test_print("SA_SA_Sstar_from_SP")
test_print("Sstar_SA_Sstar_from_SP")

#FIXME: fails with nans, passes without
#FIXME: match_args_return ?
""" Conservative Temperature (CT) """
#CT_chck_cast = gsw.conservative_t(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
#test_print("CT_chck_cast", "CT_from_t")

""" other conversions between temperatures, salinities, pressure and height """
#t_from_CT =  gsw.t_from_CT(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("t_from_CT", "t_chck_cast") #NOTE: diffs are also found in the original

#pt = gsw.potential_t(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast,gsw_cv.pr)
#test_print("pt", "pt_from_t")

#CT_from_pt = gsw.CT_from_pt(SA_chck_cast, pt)
#test_print("CT_from_pt") #NOTE: diffs are also found in the original

#pot_enthalpy = gsw.pot_enthalpy_from_pt(SA_chck_cast, pt)
#test_print("pot_enthalpy")

#pt0 = gsw.pt0_from_t(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
#test_print("pt0")

#pt_from_CT = gsw.pt_from_CT(SA_chck_cast,CT_chck_cast)
#test_print("pt_from_CT", "pt")

""" More on salinity """
Sstar_from_SA = gsw.Sstar_from_SA(SA_chck_cast, p_chck_cast, long_chck_cast, lat_chck_cast)
#Sstar_from_SA = gsw.Sstar_from_SA(SA_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)
test_print("Sstar_from_SA")

#SA_from_Sstar = gsw.SA_from_Sstar(SA_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)
SA_from_Sstar = gsw.SA_from_Sstar(SA_chck_cast, p_chck_cast, long_chck_cast, lat_chck_cast)
SA_from_Sstar = np.reshape( SA_from_Sstar, gsw_cv.SP_chck_cast.shape )
test_print("SA_from_Sstar")

#SP_from_SA = gsw.SP_from_SA(SA_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)
SP_from_SA = gsw.SP_from_SA(SA_chck_cast, p_chck_cast, long_chck_cast, lat_chck_cast)
SP_from_SA = np.reshape( SP_from_SA, gsw_cv.SP_chck_cast.shape )
test_print("SP_from_SA", "SP_chck_cast") #NOTE: diffs are also found in the original

#SP_from_Sstar = gsw.SP_from_SA(SA_chck_cast, gsw_cv.p_chck_cast, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast)
SP_from_Sstar = gsw.SP_from_SA(SA_chck_cast, p_chck_cast, long_chck_cast, lat_chck_cast)
test_print("SP_from_Sstar")

""" pressure /depth """
z_from_p = gsw.z_from_p(gsw_cv.p_chck_cast, gsw_cv.lat_chck_cast)
test_print("z_from_p")

#FIXME: output is 1D when input is 2D
#p_from_z = gsw.p_from_z(z_from_p, gsw_cv.lat_chck_cast)
p_from_z = np.reshape( gsw.p_from_z(z_from_p, gsw_cv.lat_chck_cast), gsw_cv.p_from_z.shape )
test_print("p_from_z") #FIXME: diffs are not found in the original

""" more temperatures conversions """
from seawater.csiro import T90conv
t90_from_t68 = T90conv(gsw_cv.t_chck_cast, t_type='T68')
test_print("t90_from_t68")

t90_from_t48 = T90conv(gsw_cv.t_chck_cast, t_type='T48')
test_print("t90_from_t48")

""" density and enthalpy,  based on the 25-term expression for density """
#rho_CT25 = gsw.rho_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("rho_CT25")

#rho_CT25rab, alpha_CT25rab, beta_CT25rab = gsw.rho_alpha_beta_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("rho_CT25rab")
#test_print("alpha_CT25rab")
#test_print("beta_CT25rab")

#specvol_CT25 = gsw.specvol_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("specvol_CT25")

#specvol_anom_CT25 = gsw.specvol_anom_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("specvol_anom_CT25")

#enthalpy_CT25 =  gsw.enthalpy_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("enthalpy_CT25")

#enthalpy_diff_CT25 =  gsw.enthalpy_diff_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast_shallow, gsw_cv.p_chck_cast_deep)
#test_print("enthalpy_diff_CT25")

""" water column properties,  based on the 25-term expression for density """
#n2, p_mid_n2 = gsw.Nsquared_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast, gsw_cv.lat_chck_cast)
#test_print("n2")
#test_print("p_mid_n2")

#Tu, Rsubrho, p_mid_TuRsr = gsw.Turner_Rsubrho_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("Tu")
#test_print("Rsubrho")
#test_print("p_mid_TuRsr")

#IPVfN2, p_mid_IPVfN2 = gsw.IPV_vs_fNsquared_ratio_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast, gsw_cv.pr)
#test_print("IPVfN2")
#test_print("p_mid_IPVfN2")

""" geostrophic streamfunctions,  based on the 25-term expression for density """
#geo_strf_dyn_height = gsw.geo_strf_dyn_height(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("geo_strf_dyn_height")

#geo_strf_dyn_height_pc, dh_pmid = gsw.geo_strf_dyn_height_pc(SA_chck_cast, CT_chck_cast, gsw_cv.delta_p_chck_cast)
#test_print("geo_strf_dyn_height_pc")
#test_print("dh_pmid")

#geo_strf_McD_Klocker = gsw.geo_strf_McD_Klocker(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast, gsw_cv.Neutral_Density, gsw_cv.p_Neutral_Density)
#test_print("geo_strf_McD_Klocker")

#geo_strf_McD_Klocker_pc, mk_p_mid = gsw.geo_strf_McD_Klocker_pc(SA_chck_cast, CT_chck_cast, gsw_cv.delta_p_chck_cast, gsw_cv.Neutral_Density[0], 3)
#test_print("geo_strf_McD_Klocker_pc")
#test_print("mk_p_mid")

#geo_strf_Montgomery = gsw.geo_strf_Montgomery(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("geo_strf_Montgomery")

#geo_strf_Cunningham = gsw.geo_strf_Cunningham(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("geo_strf_Cunningham")

#geo_str_velocity, gv_mid_lat, gv_mid_long = gsw.geostrophic_velocity(geo_strf_dyn_height, gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast, gsw_cv.p_chck_cast)
#test_print("geo_str_velocity")
#test_print("gv_mid_lat")
#test_print("gv_mid_long")

""" neutral and non-linear properties,  based on the 25-term expression for density """
#cabbeling_CT25 = gsw.cabbeling_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("cabbeling_CT25")

#thermobaric_CT25 = gsw.thermobaric_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("thermobaric_CT25")

#isopycnal_slope_ratio_CT25 = gsw.isopycnal_slope_ratio_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast, gsw_cv.pr)
#test_print("isopycnal_slope_ratio_CT25")

#G_CT_CT25, p_mid_G_CT_CT25 = gsw.isopycnal_vs_ntp_CT_ratio_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast, gsw_cv.pr)
#test_print("G_CT_CT25")
#test_print("p_mid_G_CT_CT25")

#ntpptCT_CT25 = gsw.ntp_pt_vs_CT_ratio_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("ntpptCT_CT25")

""" basic thermodynamic properties """
SA_chck_cast = np.reshape( SA_chck_cast, gsw_cv.z_from_p.shape ) #FIXME: won't be need once _delta_SA is fixed
gsw_cv.p_chck_cast = np.reshape( gsw_cv.p_chck_cast, gsw_cv.z_from_p.shape )

rho = gsw.rho(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("rho")

pot_rho = gsw.pot_rho(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast, gsw_cv.pr)
test_print("pot_rho")

specvol = gsw.specvol(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("specvol")

specvol_anom = gsw.specvol_anom(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("specvol_anom")

alpha_wrt_CT = gsw.alpha_wrt_CT(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("alpha_wrt_CT")

alpha_wrt_pt = gsw.alpha_wrt_pt(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("alpha_wrt_pt")

alpha_wrt_t = gsw.alpha_wrt_t(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("alpha_wrt_t")

beta_const_CT = gsw.beta_const_CT(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("beta_const_CT")

beta_const_pt = gsw.beta_const_pt(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("beta_const_pt")

beta_const_t = gsw.beta_const_t(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("beta_const_t")

entropy = gsw.entropy(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("entropy")

internal_energy = gsw.internal_energy(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("internal_energy")

enthalpy = gsw.enthalpy(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("enthalpy")

cp = gsw.cp(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("cp")

isochoric_heat_cap = gsw.isochoric_heat_cap(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("isochoric_heat_cap")

chem_potential =  gsw.chem_potential_relative(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("chem_potential")

chem_potential_water =  gsw.chem_potential_water(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("chem_potential_water")

chem_potential_salt =  gsw.chem_potential_salt(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("chem_potential_salt")

Helmholtz_energy = gsw.helmholtz_energy(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("Helmholtz_energy")

sound_speed = gsw.sound_speed(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("sound_speed")

kappa = gsw.kappa(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("kappa")

kappa_const_t = gsw.kappa_const_t(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("kappa_const_t")

adiabatic_lapse_rate = gsw.adiabatic_lapse_rate(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("adiabatic_lapse_rate")

molality = gsw.molality(SA_chck_cast)
test_print("molality")

ionic_strength = gsw.ionic_strength(SA_chck_cast)
test_print("ionic_strength")

osmotic_coefficient = gsw.osmotic_coefficient(SA_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("osmotic_coefficient")

#t_maxden, pt_maxden, CT_maxden = gsw.temps_maxdensity(SA_chck_cast, gsw_cv.p_chck_cast)
#test_print("t_maxden")
#test_print("pt_maxden")
#test_print("CT_maxden")

""" basic thermodynamic properties in terms of CT and pt """
#rho_CT = gsw.rho_CT(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("rho_CT")

#rho_CTrab, alpha_CTrab, beta_CTrab = gsw.rho_alpha_beta_CT(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("rho_CTrab")
#test_print("alpha_CTrab")
#test_print("beta_CTrab")

#specvol_CT = gsw.specvol_CT25(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("specvol_CT")

#specvol_anom_CT = gsw.specvol_anom_CT(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("specvol_anom_CT")

#sigma0_CT = gsw.sigma0_CT(SA_chck_cast, CT_chck_cast)
#test_print("sigma0_CT")

#sigma1_CT = gsw.sigma1_CT(SA_chck_cast, CT_chck_cast)
#test_print("sigma1_CT")

#sigma2_CT = gsw.sigma2_CT(SA_chck_cast, CT_chck_cast)
#test_print("sigma2_CT")

#sigma3_CT = gsw.sigma3_CT(SA_chck_cast, CT_chck_cast)
#test_print("sigma3_CT")

#sigma4_CT = gsw.sigma4_CT(SA_chck_cast, CT_chck_cast)
#test_print("sigma4_CT")

#enthalpy_CT =  gsw.enthalpy_CT(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("enthalpy_CT")

#enthalpy_diff_CT =  gsw.enthalpy_diff_CT(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast_shallow, gsw_cv.p_chck_cast_deep)
#test_print("enthalpy_diff_CT")

#FIXME: pt is not passing
#entropy_from_pt =  gsw.entropy_from_t(SA_chck_cast, pt, t_type='pt')
#test_print("entropy_from_pt") #NOTE: diffs are also found in the original

#FIXME: CT_chck_cast is not passing
#entropy_from_CT =  gsw.entropy_from_t(SA_chck_cast, CT_chck_cast, t_type='CT')
#test_print("entropy_from_CT")

CT_from_entropy =  gsw.t_from_entropy(SA_chck_cast, entropy, t_type='CT')
test_print("CT_from_entropy") #FIXME: diffs are not found in the original

pt_from_entropy =  gsw.t_from_entropy(SA_chck_cast, entropy, t_type='pt')
test_print("pt_from_entropy") #FIXME: diffs are not found in the original

""" derivatives of enthalpy,  entropy,  CT and pt """
#[CT_SA,  CT_pt] = gsw.CT_first_derivatives(SA_chck_cast, pt)
#test_print("CT_SA")
#test_print("CT_pt")

#CT_SA_SA, CT_SA_pt, CT_pt_pt = gsw.CT_second_derivatives(SA_chck_cast, pt)
#test_print("CT_SA_SA")
#test_print("CT_SA_pt")
#test_print("CT_pt_pt")

#h_SA, h_CT, h_P = gsw.enthalpy_first_derivatives(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("h_SA")
#test_print("h_CT")
#test_print("h_P")

#h_SA_SA, h_SA_CT, h_CT_CT = gsw.enthalpy_second_derivatives(SA_chck_cast, CT_chck_cast, gsw_cv.p_chck_cast)
#test_print("h_SA_SA")
#test_print("h_SA_CT")
#test_print("h_CT_CT")

#eta_SA, eta_CT = gsw.entropy_first_derivatives(SA_chck_cast, CT_chck_cast)
#test_print("eta_SA")
#test_print("eta_CT")

#eta_SA_SA, eta_SA_CT, eta_CT_CT = gsw.entropy_second_derivatives(SA_chck_cast, CT_chck_cast)
#test_print("eta_SA_SA")
#test_print("eta_SA_CT")
#test_print("eta_CT_CT")

#pt_SA,  pt_CT = gsw.pt_first_derivatives(SA_chck_cast, CT_chck_cast)
#test_print("pt_SA")
#test_print("pt_CT")

#pt_SA_SA, pt_SA_CT, pt_CT_CT = gsw.pt_second_derivatives(SA_chck_cast, CT_chck_cast)
#test_print("pt_SA_SA")
#test_print("pt_SA_CT")
#test_print("pt_CT_CT")

""" planet earth properties """
from seawater.csiro import cor
f = cor(gsw_cv.lat_chck_cast)
test_print("f")

grav = gsw.grav(gsw_cv.lat_chck_cast, gsw_cv.p_chck_cast)
test_print("grav")

distance = gsw.distance(gsw_cv.long_chck_cast, gsw_cv.lat_chck_cast, gsw_cv.p_chck_cast)
test_print("distance") #FIXME: diffs are not found in the original

""" Absolute Salinity from direct density measurements:- a laboratory function """
SA_from_rho = gsw.SA_from_rho(rho, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("SA_from_rho")

#sigma0_pt = gsw.sigma0_pt(SA_chck_cast, pt0)
#test_print("sigma0_pt")

""" Practical Salinity (SP):- PSS-78 """
cndr = gsw.cndr_from_SP(gsw_cv.SP_chck_cast, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("cndr") #FIXME: diffs are not found in the original

SP_from_cndr = gsw.SP_from_cndr(cndr, gsw_cv.t_chck_cast, gsw_cv.p_chck_cast)
test_print("SP_from_cndr") #FIXME: diffs are not found in the original
