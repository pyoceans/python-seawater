# -*- coding: utf-8 -*-

"""Unit check for standard profiles for the 
   Gibbs Sea Water python package

   Functionality similar to matlab_test.py
   without the nice output
"""

# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-03-03

import os
import unittest
import functools   # require python 2.5
import numpy as np
import seawater.gibbs as gsw
from seawater.library import Dict2Struc

# ------------------------------------
# Table of functions with arguments
# ------------------------------------

# Codes for non-tested functions
# NI: Not Implemented
# NG: In seawater.library not in seawater.gibbs
# NV: No test value
# NA: Don't fit the testing scheme
# ERR: Gives error

# Could perhaps be auto-generated
function_arguments = {

#NA    'CT_first_derivatives' : ('SA', 'pt')
    'CT_from_entropy' : ('SA', 'entropy'),
    'CT_from_pt'      : ('SA', 'pt'),
    'CT_from_t'       : ('SA', 't', 'p'),
#NA    'CT_second_derivatives' : ('SA', 'pt'),

    'Helmholtz_energy' : ('SA', 't', 'p'),
    
#NI    'IPV_vs_fNsquared_ratio_CT25 :

#NI    'Nsquared_CT25' :

#NA    'SA_Sstar_from_SP' : ('SP', 'p', 'long', 'lat'),
    'SA_from_SP' : ('SP', 'p', 'long', 'lat'),
#NG    'SA_from_SP_Baltic' : ('SP', 'long', 'lat'),
    'SA_from_Sstar' : ('Sstar', 'p', 'long', 'lat'),
    'SA_from_rho' : ('rho', 't', 'p'),
    'SP_from_SA'  : ('SA', 'p', 'long', 'lat'),
#NG    'SP_from_SA_Baltic' : ('SA', 'long', 'lat'),
    'SP_from_Sstar' : ('Sstar', 'p', 'long', 'lat'),
    'SP_from_cndr'  : ('R', 't', 'p'),
    'Sstar_from_SA' : ('SA', 'p', 'long', 'lat'),
    'Sstar_from_SP' : ('SP', 'p', 'long', 'lat'),

#NI    'Turner_Rsubrho_CT25' : ('SA', 'CT', 'p'),

    'adiabatic_lapse_rate' : ('SA', 't', 'p'),
    'alpha_wrt_CT'         : ('SA', 't', 'p'),
    'alpha_wrt_pt'         : ('SA', 't', 'p'),
    'alpha_wrt_t'          : ('SA', 't', 'p'),
    
    'beta_const_CT' : ('SA', 't', 'p'),
    'beta_const_pt' : ('SA', 't', 'p'),
    'beta_const_t'  : ('SA', 't', 'p'),
    
#NI    'cabelling_CT25'          : ('SA', 't', 'p'),
#NV?    'chem_potential_relative' : ('SA', 't', 'p'),
    'chem_potential_salt'     : ('SA', 't', 'p'),
    'chem_potential_water'    : ('SA', 't', 'p'),
    'cp'                      : ('SA', 't', 'p'),
    
#NV   'delta_SA' : ('p', 'long', 'lat'),
    'distance' : ('long', 'lat', 'p'),

    'enthalpy' : ('SA', 't', 'p'),
#NI    'enthalpy_CT' : ('SA', 'CT', 'p'),
    'enthalpy_CT25' : ('SA', 'CT', 'p'),
#NI    'enthalpy_SSO_0_CT25' : ('p',),
#NA    'enthalpy_diff_CT' : ('SA', 'CT', 'p', 'p'),
#NA    'enthalpy_diff_CT25' : ('SA', 'CT', 'p', 'p),
#NA    'enthalpy_first_derivatives' : ('SA', 'CT', 'p'),
#NA    'enthalpy_second_derivatives' : ('SA', 'CT', 'p'),
    'entropy' : ('SA', 't', 'p'),
#NA    'entropy_first_derivatives' : ('SA', 'CT'),
    'entropy_from_CT' : ('SA', 'CT'),
    'entropy_from_pt' : ('SA', 'pt'),
#?    'entropy_part' : ?
#?    'entropy_part_zerop' : ?
#NA    'entropy_second_derivatives' : ('SA', 'pt'),
    
#?    'f' : ('lat'),

#NA    'geo_strf_Cunningham' :
#NA    'geo_strf_McD_Klocker' :
#NA    'geo_strf_McD_Klocker_pc' :
#NA    'geo_strf_Montgomery' :
#NA    'geo_strf_dyn_height' :
#NA    'geo_strf_dyn_height_pc' :
#NA?    'geostrophic_velocity' :
#NA    'gibbs' :
#NV    'gibbs_pt0_pt0' :
    'grav' : ('lat', 'p'),

    'internal_energy'    : ('SA', 't', 'p'),
#NA    'interp_McD_Klocker' :
#NA    'interp_SA_CT' :
    'ionic_strength'     : ('SA',),
    'isochoric_heat_cap' : ('SA', 't', 'p'),
#NI    'isopycnal_slope_ratio_CT25' :
#NI    'isopycnal_vs_ntp_CT_ratio_CT25' :

    'kappa'         : ('SA', 't', 'p'),
    'kappa_const_t' : ('SA', 't', 'p'),
    
    'molality'      : ('SA',),

#NI    'ntp_pt_vs_CT_ratio_CT25' : ('SA', 't', 'p'),

    'osmotic_coefficient' : ('SA', 't', 'p'),

    'p_from_z'             : ('z', 'lat'),
    'pot_enthalpy_from_pt' : ('SA', 'pt'),
#NA    'pot_rho' : ('SA', 't', 'p', 'pr'),
    'pt0_from_t'      : ('SA', 't', 'p'),
#NA    'pt_first_derivatives' :
    'pt_from_CT'      : ('SA', 'CT'),
    'pt_from_entropy' : ('SA', 'entropy'),
#NA    'pt_from_t' : ('SA', 't', 'p', 'pr'),
#NA    'pt_second_derivatives' :

    'rho' : ('SA', 't', 'p'),
#NI    'rho_CT' : ('SA', 'CT', 'p'),
    'rho_CT25' : ('SA', 'CT', 'p'),
#NA    'rho_alpha_beta_CT' : ('SA', 'CT', 'p'),
#NA    'rho_alpha_beta_CT25' : ('SA', 'CT', 'p'),




    
            }


# ------------------------------------------
# Read data file with check value profiles
# ------------------------------------------

datadir = os.path.join(os.path.dirname(gsw.__file__), 'data')
fname = 'gsw_cv.npz'
cv = Dict2Struc(np.load(os.path.join(datadir, fname)))

# ---------------------------
# Aliases in the cv struct
# ---------------------------

# Make aliases for some values to be used as arguments

# Note: there seems to be a bug in the dataset
# The arrays cv.SA_from_Sstar is wrong
# The error is also in mat-file
#
# cv.SA_from_Sstar != cv.SA_from_SP
# cv.SA_from_Sstar != gsw.SA_from_Sstar(Sstar, p, lon, lat)
# cv.SA_from_Sstar == gsw.SA_from_Sstar(SA, p, lon, lat)
#
# Bug work-around
cv.SA_from_Sstar = cv.SA_from_SP


cv.SA_chck_cast      = cv.SA_from_SP
cv.CT_chck_cast      = cv.CT_from_t
cv.pt_chck_cast      = cv.pt_from_t
cv.Sstar_chck_cast   = cv.Sstar_from_SA 
##cv.Sstar_chck_cast   = cv.SA_chck_cast    
cv.entropy_chck_cast = cv.entropy
cv.z_chck_cast       = cv.z_from_p
cv.rho_chck_cast     = cv.rho
cv.R_chck_cast       = cv.cndr

# Make aliases for check values whose names
# does not match the function
not_match = {'pt0_from_t' : 'pt0',
             'pt_from_CT' : 'pt',
             'pot_enthalpy_from_pt' : 'pot_enthalpy',
            }
for f in not_match:
    setattr(cv, f, getattr(cv, not_match[f]))
    setattr(cv, f+'_ca', getattr(cv, not_match[f]+'_ca'))

# ---------------------
# Generic test method
# ---------------------

def generic_test(self, func=None, argnames=None):
    """Generic test function, to be spesialized by functools.partial"""
    # Transform argument names to name convention in cv dataset
    args = [getattr(cv, a+'_chck_cast') for a in argnames]
    # Perform the function call
    out = getattr(gsw, func)(*args)
    # Check that the maximal error is less than the given tolerance
    maxdiff = np.nanmax(abs(out - getattr(cv, func)))
    #print func, maxdiff, getattr(cv, func+'_ca')
    self.assertTrue(maxdiff < getattr(cv, func+'_ca'))


# --------------------------------------------------------
# Dictionary of functions with corresponding test methods
# --------------------------------------------------------

function_test = {}
for f in function_arguments:
    function_test[f] = functools.partial(generic_test, 
                      func=f, argnames=function_arguments[f])

# ---------------------------
# Auto-generated TestCase 
# ---------------------------

class Test_profiles(unittest.TestCase):

    for f in function_test:
        method_def = ( "test_" + f + 
            " = lambda self: function_test['" + f + "'](self)" )
        exec(method_def)


if __name__ == '__main__':
    # A more verbose output
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_profiles)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()

