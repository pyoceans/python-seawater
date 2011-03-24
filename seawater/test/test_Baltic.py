# -*- coding: utf-8 -*-

"""Test seawater functions for correct behaviour in the Baltic"""

# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-03-04

import unittest
import numpy as np
import seawater.gibbs as gsw
import seawater.gibbs.library as gswl

class Test_in_Baltic(unittest.TestCase):

    def test_scalar_inside(self):
        """Should return np.array([True]) for point inside Baltic"""
        lon, lat = 20, 59
        inside = gswl._in_Baltic(lon, lat)
        self.assertEqual(inside, np.array([True]))

    def test_scalar_outside(self):
        """Should return np.array([False]) for point outside"""
        # longitude outside
        lon, lat = -10, 59       # North Atlantic  
        inside = gswl._in_Baltic(lon, lat)
        self.assertEqual(inside, np.array([False]))

        # latitude outside
        lon, lat = 20, 77        # Barents Sea
        inside = gswl._in_Baltic(lon, lat)
        self.assertEqual(inside, np.array([False]))

        # in rectangle, outside Baltic
        lon, lat = 10.0, 66.0    # Norwegian Sea
        inside = gswl._in_Baltic(lon, lat)
        self.assertEqual(inside, np.array([False]))

        # outside rectangle, on land inside Baltic polygon
        lon, lat = 40.0, 52.0    # Russia
        inside = gswl._in_Baltic(lon, lat)
        self.assertEqual(inside, np.array([False]))

    def test_array_input(self):
        """Accept arrays and sequences as input, respect shape"""
        lon = [-10, 20]   # outside, inside
        lat = 59
        inside = gswl._in_Baltic(lon, lat)
        self.assertTrue(np.all(inside == np.array([False, True])))

        # test 2D array
        lon = np.array([[-10, 20], [21, 22]])
        lat = 59
        inside = gswl._in_Baltic(lon, lat)
        self.assertTrue(np.all(
              inside == np.array([[False, True],[True, True]])))

    def test_nan_input(self):
        """Accept not-a-number, return False"""
        lon = [np.nan, 20, -10]
        lat = 59
        inside = gswl._in_Baltic(lon, lat)
        self.assertTrue(np.all(inside == np.array([False, True, False])))

    def test_masked_input(self):
        """Accept masked arrays, return False for masked entries"""
        # Point inside Baltic
        lon = np.ma.MaskedArray([20, 20, 20], mask=[True, False, False])
        lat = np.ma.MaskedArray([59, 59, 59], mask=[True, False, True])
        inside = gswl._in_Baltic(lon, lat)
        self.assertTrue(np.all(inside == np.array([False, True, False])))

        # Point outside Baltic
        lon = np.ma.MaskedArray([-10, -10], mask=[True, False])
        lat = 59
        inside = gswl._in_Baltic(lon, lat)
        self.assertTrue(np.all(inside == np.array([False, False])))

# -----------------------------------------
        
class Test_SA_from_SP_Baltic(unittest.TestCase):

    def test_standard(self):
        """Test standard values"""
        SP  = [6.5683, 6.6719, 6.8108, 7.2629, 7.4825, 10.2796]
        lon, lat = 20, 59

        SA = gswl._SA_from_SP_Baltic(SP,lon,lat)
        SA_standard = np.array([6.6699,  6.7738,
                                6.9130,  7.3661,
                                7.5862, 10.3895])
        self.assertTrue(np.all(abs(SA - SA_standard) < 1.0e-4))

    def test_scalar(self):
        """Return a 1-element 1D masked array for a scalar"""

        # position in Baltic, not masked
        SP = 6.5683
        lon, lat = 20, 59
        SA = gswl._SA_from_SP_Baltic(SP, lon, lat)
        self.assertEqual(SA.shape, (1,))
        self.assertFalse(SA.mask)

        # position outside Baltic, masked
        lon = -10
        SA = gswl._SA_from_SP_Baltic(SP, lon, lat)
        self.assertEqual(SA.shape, (1,))
        self.assertTrue(SA.mask)

    def test_in_out(self):
        """Test array with positions inside and outside"""
        SP = [6.5683, 35.0]
        lon = [20, -10]
        lat = 59
        SA = gswl._SA_from_SP_Baltic(SP, lon, lat)
        #print SA, SA.mask
        self.assertEqual(SA.shape, (2,))
        self.assertFalse(SA.mask[0])          # Inside Baltic
        self.assertTrue(SA.mask[1])           # Outside
        
    def test_masked_array(self):
        """Handle masked arrays,
           mask output where not in Baltic or masked input
        """
        SP = np.ma.MaskedArray([6.5683, 35.0, 10.0, 10.0],
                               mask=[False, False, True, False])
        lon = np.ma.MaskedArray([20, -10, 20, 20],
                               mask=[False, False, False, True])
        lat = 59
        SA = gswl._SA_from_SP_Baltic(SP, lon, lat)
        self.assertTrue(np.all(SA.mask == [False, True, True, True]))
        self.assertTrue(np.isfinite(SA[0]))  # ordinary value, inside

    def test_nans(self):
        """Handle Not-a-Numbers,
           masked output for point outside or nan in input
        """
        SP = [6.5683, 35.0, np.nan, 10.0]
        lon = [20, -10, 20, np.nan]
        lat = 59
        SA = gswl._SA_from_SP_Baltic(SP, lon, lat)
        self.assertTrue(np.all(SA.mask == [False, True, True, True]))
        self.assertTrue(np.isfinite(SA[0]))  # ordinary value, inside

# -----------------------------------------------------

class test_SP_from_SA_Baltic(unittest.TestCase):

    def test_standard(self):
        """Test standard values"""
        SA  = [6.6699, 6.7738, 6.9130, 7.3661, 7.5862, 10.3895]
        lon, lat = 20, 59
        SP = gswl._SP_from_SA_Baltic(SA, lon, lat)
        SP_standard  = np.array([6.5683,  6.6719,
                                 6.8108,  7.2629,
                                 7.4825, 10.2796])
        self.assertTrue(np.all(abs(SP - SP_standard) < 1.0e-4))






if __name__ == '__main__':
    unittest.main()


        
