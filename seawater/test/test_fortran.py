# -*- coding: utf-8 -*-

"""
Testing values from the Fortran GSW library, version 1

Values from check_values.f90

These values are independent of values
calculated by python package or matlab toolbox

"""

# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-03-16

import sys
import unittest
import functools
import numpy.testing as npt
import seawater.gibbs as gsw
import seawater.gibbs25 as gsw25
import seawater.library as gswl

# ------------------------------------------
# Check values from check_values.f90
# in Fortran GSW library, version 1
# ------------------------------------------

# Argument values
arg_values = {
       'SP'  :   35.52764437773386,
       't'   :   25.5,
       'CT'  :   25.4805463842239,
       'p'   : 1023.0,
       'pr'  :    0.0,
       'SA'  :   35.7,
       'lon' :  201.0,
       'lat' :  -21.0,
       'pt'  :   25.2720983155409,
             }

# Check values from check_values.f90
# Fortran_name : (check_value, decimals)
#   decimals found by trial and error
check_values =        {
    'Asal'            : (35.7,                   4),
    'PSal_from_ASal'  : (35.5276443777339,       4),
    'alpha_t'         : (0.0003098378393192645, 18),
    'beta_t'          : (0.0007257297978386655, 32),
    'cp'              : (3974.42541259729,      32),
    'ctmp'            : (25.4805463842239,      13),
    'ptmp0_from_ctmp' : (25.5,                  13),
    'density'         : (1027.95249315662,      12),
    'enthalpy'        : (110776.712408975,       9),
    'entropy'         : (352.81879771528,       12),
    'kappa'           : (4.033862685464779e-10, 24), # corrected to dbar 
    'kappa_t'         : (4.104037946151349e-10, 24), # corrected to dbar
    'pden'            : (1023.66254941185,      11),
    'ptmp'            : (25.2720983155409,      13),
    'ptmp_inverse'    : (25.5,                  13),
    'specvol'         : (0.0009728076021579713, 32),
    'svel'            : (1552.93372863425,      11),
                      }

# Translation between fortran and python functions + arguments
python_function =     {
    'Asal'            : (gsw.SA_from_SP,    ('SP', 'p', 'lon', 'lat')),
    'PSal_from_ASal'  : (gsw.SP_from_SA,    ('SA', 'p', 'lon', 'lat')),
    'alpha_t'         : (gsw.alpha_wrt_t,   ('SA', 't', 'p')),
    'beta_t'          : (gsw.beta_const_t,  ('SA', 't', 'p')),
    'cp'              : (gsw.cp,            ('SA', 't', 'p')),
    'ctmp'            : (gsw.CT_from_pt,    ('SA', 't')),
    'ptmp0_from_ctmp' : (gsw.pt_from_CT,    ('SA', 'CT')), 
    'density'         : (gsw.rho,           ('SA', 't', 'p')),
    'enthalpy'        : (gsw.enthalpy,      ('SA', 't', 'p')),
    'entropy'         : (gsw.entropy,       ('SA', 't', 'p')),
    'kappa'           : (gsw.kappa,         ('SA', 't', 'p')),
    'kappa_t'         : (gsw.kappa_const_t, ('SA', 't', 'p')),
    'pden'            : (gsw.pot_rho,       ('SA', 't', 'p', 'pr')),
    'ptmp'            : (gsw.pt_from_t,     ('SA', 't', 'p', 'pr')),
    'ptmp_inverse'    : (gsw.pt_from_t,     ('SA', 'pt', 'pr', 'p')),
    'specvol'         : (gsw.specvol,       ('SA', 't', 'p')),
    'svel'            : (gsw.sound_speed,   ('SA', 't', 'p')),
                      }

# Generic test function to be specialised by functools.partial
def generic_test(self, fname):
    func, arg_names = python_function[fname]
    out_check, decimal = check_values[fname]
    args = [arg_values[s] for s in arg_names]
    out = func(*args)
    #print fname, abs(out-out_check)
    npt.assert_almost_equal(out, out_check, decimal=decimal)

# Make a dictionary of test methods
test_func = {}
for f in python_function:
    test_func[f] = functools.partial(generic_test, fname=f)


# The unit test case
class Test_Fortran(unittest.TestCase):

    # Get test methods into the test class
    for f in test_func:
        method_def = ( "test_" + f + 
            " = lambda self: test_func['" + f + "'](self)" )
        exec(method_def)

# --------------------------------------------------------
# Test Gibbs function 
#
# Check values from GSW_Library_5.F90 in the SIA library
# --------------------------------------------------------

# arguments to gibbs function : check value, number of correct decimals
gibbs_check =         {
   (0,0,0,35,26.85,0) : (-5113.70064124,          8),  
   (1,0,0,35,26.85,0) : (78.5928261339,          10),
   (0,1,0,35,26.85,0) : (-374.452000830,          9),
   (0,0,1,35,26.85,0) : (0.977858058750E-03,     15),
   (2,0,0,35,26.85,0) : (2.24755137017,          11),
   (1,1,0,35,26.85,0) : (0.789935187192,         12),
   (1,0,1,35,26.85,0) : (-0.716680931996E-06,    18),
   (0,2,0,35,26.85,0) : (-13.3358337534,         10),
   (0,1,1,35,26.85,0) : (0.304607508052E-06,     18),
   (0,0,2,35,26.85,0) : (-0.410939723950E-12,    24),
                      }
    
class Test_Gibbs(unittest.TestCase):
    """Testing the Gibbs function"""
 
    def test_gibbs(self):
        """Testing Gibbs function with derivatives"""
        for args in gibbs_check:
            out = gswl._gibbs(*args)
            out_check, decimal = gibbs_check[args]
            #print k, abs(out-out_check)
            npt.assert_almost_equal(out, out_check, decimal)
        
# ------------------------------------------------
# Other check values 
# from GSW_Library_F.F90 in the SIA library
# also in Values_GSW.F90
# ------------------------------------------------

SA, t, p = 35.0, 20.0, 1000.0
arg_values2 = {
        'SA' : SA,
        't'  : t,
        'p'  : p,
        'CT' : gsw.CT_from_t(SA, t, p)
             }

check_values2 = {
    'alpha_ct'  : (2.69418609861E-04,     8),
    'alpha_pt0' : (2.69753733317E-04,     8),
    'beta_ct'   : (7.23213672954E-04,     9),
    'beta_pt0'  : (7.31582583383E-04,    11),
    'cabb_ct'   : (0.896907383083E-05,    6),
    'cabb_pt0'  : (0.875963154048E-05,   32),
    'thrmb_ct'  : (0.172708365652E-11,   13),  # corrected for dbar
    'thrmb_pt0' : (0.170945045984E-07,   32), 
                }

python_function2 = {
      'alpha_ct'   : (gsw.alpha_wrt_CT,   ('SA', 't', 'p')),
      'alpha_pt0'  : (gsw.alpha_wrt_pt,   ('SA', 't', 'p')),
      'beta_ct'    : (gsw.beta_const_CT,  ('SA', 't', 'p')),
      'beta_pt0'   : (gsw.beta_const_pt,  ('SA', 't', 'p')),
      # Cabbeling, large relative error, correct??
      'cabb_ct'    : (gsw25.cabbeling_CT25, ('SA', 'CT', 'p')),
     #'cabb_pt0'   : No python function?
      'thrmb_ct'   : (gsw25.thermobaric_CT25, ('SA', 'CT', 'p')),
     #'thrmb_pt0'  : No python function?
                   }

# Generic test function to be specialised by functools.partial
def generic_test2(self, fname):
    func, arg_names = python_function2[fname]
    out_check, decimal = check_values2[fname]
    args = [arg_values2[s] for s in arg_names]
    out = func(*args)
    #print fname, abs(out-out_check)
    npt.assert_almost_equal(out, out_check, decimal=decimal)

# Make a dictionary of test methods
test_func2 = {}
for f in python_function2:
    test_func2[f] = functools.partial(generic_test2, fname=f)

class Test_Fortran2(unittest.TestCase):

    # Get test methods into the test class
    for f in test_func2:
        method_def = ( "test_" + f + 
            " = lambda self: test_func2['" + f + "'](self)" )
        exec(method_def)

if __name__ == '__main__':
    suite1 = unittest.TestLoader().loadTestsFromTestCase(Test_Fortran)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(Test_Gibbs)
    suite3 = unittest.TestLoader().loadTestsFromTestCase(Test_Fortran2)
    suite = unittest.TestSuite([suite1, suite2, suite3])
    a = unittest.TextTestRunner(verbosity=2).run(suite)
    if a.errors or a.failures: sys.exit(256)


           

