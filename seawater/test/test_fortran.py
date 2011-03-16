# -*- coding: utf-8 -*-

"""
Testing values from the Fortran GSW library, version 1

Values from check_values.f90

These values are independent of values
calculated by python package or matlab toolbox

"""

# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-03-16

import unittest
import functools
import numpy.testing as npt
import seawater.gibbs as gsw

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
#   decimals found by testing
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
    #print abs(out-out_check)
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
    

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Fortran)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()

           

