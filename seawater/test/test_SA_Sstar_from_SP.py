# -*- coding: utf-8 -*-

"""Unit test for SA_Sstar_from_SP gibbs seawater module"""

# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-02-27

import unittest
import numpy as np
import seawater.gibbs as gsw

class Test_SA_Sstar_from_SP(unittest.TestCase):
    """Test interface to the SA_Sstar_from_SP function"""

    def test_scalar(self):
        """Should return a two-tuple of scalars for scalar input"""
        SP = 35.0
        p = 1000.0
        lon, lat = 2, 66
        output = gsw.SA_Sstar_from_SP(SP, p, lon, lat)
        self.assertEqual(type(output), tuple)
        self.assertEqual(len(output), 2)
        self.assertTrue(np.isscalar(output[0]))
    
    def test_standard_values(self):
        "Check some standard values"
        SP =   [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
        p  =   [     10,      50,     125,     250,     600,    1000]
        lon, lat = 188, 4
        SA_standard = np.array((34.711779712893147,
                                34.891523721622328,
                                35.025547737221643,
                                34.847230081708567,
                                34.736629598766619,
                                34.732361864645107))
        Sstar_standard = np.array((34.711553202053111,
                                   34.891161009146472,
                                   35.024649258886718,
                                   34.843592772087710,
                                   34.729033602488833,
                                   34.719676382802788))
        
        SA, Sstar = gsw.SA_Sstar_from_SP(SP, p, lon, lat)
        self.assertTrue(np.all(abs(SA - SA_standard) < 1.0e-18))
        self.assertTrue(np.all(abs(Sstar - Sstar_standard) < 1.0e-18))


if __name__ == '__main__':
    unittest.main()







