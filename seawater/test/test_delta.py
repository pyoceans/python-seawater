# -*- coding: utf-8 -*-

"""Unit test for _delta_SA"""


# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-02-07

import unittest
import numpy as np
import seawater.gibbs.library as gswl

class Test_delta_SA(unittest.TestCase):
    """Test the interface to the SA_from_SP function"""

    def test_scalar_input(self):
        """Accept scalar input, return scalar output"""
        p = 0.0
        lon, lat = 2.0, 66.0
        dsa = gswl._delta_SA(p, lon, lat)
        self.assertFalse(hasattr(dsa, '__iter__'))

    def test_single_array_input(self):
        """Accept and return single element arrays"""
        p = np.array([0.0])
        lon, lat = np.array([2.0]), np.array([66.0])
        dsa = gswl._delta_SA(p, lon, lat)
        self.assertTrue(dsa.shape == (1,))

    def test_non_compatible(self):
        """Raise ValueError if input not broadcastable to same shape"""
        p = np.array([10, 100, 1000])
        lon = [0.0, 2.0]
        lat = 66.0
        self.assertRaises(ValueError, gswl._delta_SA, p, lon, lat)

    def test_nan(self):
        """Accept and return not-a-number, do not mess up finite values"""
        p = np.array([10.0, np.nan, 100.0])
        lon, lat = 2, 66
        dsa  = gswl._delta_SA(p, lon, lat)
        dsa0 = gswl._delta_SA(p[0], lon, lat)
        # Should return NaN for NaN
        self.assertTrue(np.isnan(dsa[1]))
        # Should return correct value for not NaN
        self.assertEqual(dsa[0], dsa0)

    def test_masked(self):
        """Accept and return correctly masked arrays"""
        p = np.array([10, 50, -99, 100.0])  # One pressure value missing
        p = np.ma.masked_less(p, 0)
        lon, lat = 2, 66
        dsa  = gswl._delta_SA(p, lon, lat)
        dsa0 = gswl._delta_SA(p[0], lon, lat)
        # Return array should have the same mask
        self.assertTrue(np.all(dsa.mask == p.mask))
        # Correct value for non-masked entries
        self.assertEqual(dsa[0], dsa0)

    def test_extreme_positions(self):
        """Test that extreme positions works"""
        # May include test with value error for illegal positions
        p    = 0.0
        lon1 = [-200, -180.0, 0.0, 180.0, 359.9, 360.0, 400.0]
        lat1 = 66.0
        gsa1 = gswl._delta_SA(p, lon1, lat1)
        lon2 = -10.0
        lat2 = [-110, -90.0, -89.0, 0.0, 89.9, 90.0, 130.0]
        gsa2 = gswl._delta_SA(p, lon2, lat2)


    def test_check_standard(self):
        r"""Check some standard values

        Standard values from
            http://www.teos-10.org/pubs/gsw/html/gsw_delta_SA.html
        """
        
        p =  [10, 50, 125, 250, 600, 1000]
        lon, lat = 188, 4
        dsa_standard = np.array((0.000167785807437,
                                 0.000268675908040,
                                 0.000665539507353,
                                 0.002694303422857,
                                 0.005626663909471,
                                 0.009396653216531))
        dsa = gswl._delta_SA(p, lon, lat)
        self.assertTrue(np.all(np.abs(dsa - dsa_standard) < 1.0E-15))


if __name__ == '__main__':
    unittest.main()







