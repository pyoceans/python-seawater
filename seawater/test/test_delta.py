# -*- coding: utf-8 -*-

"""Unit test for _delta_SA"""


# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-02-07

import unittest
import numpy as np
import seawater.gibbs as gsw

class Test_delta_SA(unittest.TestCase):
    """Test the interface to the SA_from_SP function"""

    def test_scalar_input(self):
        """Accept scalar input, return sigle element array"""
        # May be changed to return scalar
        p = 0.0
        lon, lat = 2.0, 66.0
        dsa = gsw._delta_SA(p, lon, lat)[0]
        self.assertTrue(dsa.shape == (1,))

    def test_single_array_input(self):
        """Accept and return single element arrays"""
        p = np.array([0.0])
        lon, lat = np.array([2.0]), np.array([66.0])
        dsa = gsw._delta_SA(p, lon, lat)[0]
        self.assertTrue(dsa.shape == (1,))

    def test_not_1D(self):
        """Raise ValueError for nD, n > 1"""
        p = np.array([[10, 100], [20, 200]])
        lon, lat = 2.0, 66.0
        self.assertRaises(ValueError, gsw._delta_SA, p, lon, lat)

    def test_non_compatible(self):
        """Raise ValueError if input not broadcastable to same shape"""
        p = np.array([10, 100, 1000])
        lon = [0.0, 2.0]
        lat = 66.0
        self.assertRaises(ValueError, gsw._delta_SA, p, lon, lat)
        
    def test_check_values(self):
        """Check some standard values"""
        p = np.array([10, 50, 125, 250, 600, 1000])
        lon, lat = 188, 4
        gsa_fasit = np.array([ 0.00016779, 0.00026868, 0.00066554,
                               0.0026943 , 0.00562666, 0.00939665])
        gsa = gsw._delta_SA(p, lon, lat)[0]
        self.assertTrue(np.all(np.abs(gsa - gsa_fasit) < 1.0E-8))
                             
    def test_extreme_positions(self):
        """Test that extreme positions works"""
        # May include test with value error for illegal positions
        p    = 0.0
        lon1 = [ -200, -180.0, 0.0, 180.0, 359.9, 360.0, 400.0]
        lat1 = 66.0
        gsa1 = gsw._delta_SA(p, lon1, lat1)[0]
        lon2 = -10.0
        lat2 = [ -110, -90.0, -89.0, 0.0, 89.9, 90.0, 130.0]
        gsa2 = gsw._delta_SA(p, lon2, lat2)[0]

    def test_callable_(self):
        """test that call from SA_from_SP works"""
        SP = [33, 34, 35]
        p = 100
        lon, lat = 2, 66
        SA = gsw.SA_from_SP(SP, p, lon, lat)

    


if __name__ == '__main__':
    unittest.main()

  




  
