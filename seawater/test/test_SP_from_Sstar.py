# -*- coding: utf-8 -*-

"""Unit code test for the SA_from_SP function"""

# Gjør dette til generell test-baltic funksjon

# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-02-13

import unittest
import numpy as np
import seawater.gibbs as gsw

class Test_SP_from_Sstar(unittest.TestCase):

    def test_scalar_input(self):
        """Return a scalar for scalar arguments"""
        Sstar = 35.0
        p = 0.0
        lon, lat = 2.0, 66.0
        SP = gsw.SP_from_Sstar(Sstar, p, lon, lat)
        self.assertTrue(np.isscalar(SP))


    def test_in_baltic(self):
        """Test that the function works in the baltic"""

        Sstar = 35.0
        p = 20.0
        lon, lat = 20.0, 59.0
        SP = gsw.SP_from_Sstar(Sstar, p, lon, lat)
        self.assertTrue(np.isscalar(SP))


    def test_in_baltic2(self):
        """Test that the function works in the baltic"""

        Sstar = 35.0
        p = 20.0
        lon = [-10.0, 20.0]
        lat = 59.0
        SP = gsw.SP_from_Sstar(Sstar, p, lon, lat)
        self.assertTrue(np.shape(SP) == (2,))




if __name__ == '__main__':
    unittest.main()
