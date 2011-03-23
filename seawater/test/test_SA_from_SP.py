# -*- coding: utf-8 -*-

"""Unit code test for the SA_from_SP function"""

# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-02-13

import unittest
import numpy as np
import seawater.gibbs as gsw

class Test_SA_from_SP(unittest.TestCase):
    """Test the interface to the SA_from_SP function"""

    def test_scalar_input(self):
        """Return a scalar for scalar arguments"""
        SP = 35.0
        p = 0.0
        lon, lat = 2.0, 66.0
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        self.assertTrue(np.isscalar(SA))

    def test_single_array_input(self):
        """Accept and return single element arrays"""
        SP = np.array([35.0])
        p = np.array([0.0])
        lon, lat = np.array([2.0]), np.array([66.0])
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        self.assertTrue(SA.shape == (1,))

    def test_list_input(self):
        """Accept list input"""
        SP = [30.0, 35.0]
        p  = [10.0, 10.0]
        lon, lat = [2.0, 2.0], [65.0, 66.0]
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        self.assertTrue(SA.shape == (2,))

    def test_array_shape(self):
        """Respect the shape of input arrays"""
        SP = np.array([[30.0, 31.0, 32.0], [33.0, 34.0, 35.0]])
        p  = np.array([[10.0, 20.0, 30.0], [10.0, 20.0, 30.0]])
        lon, lat = 2.0, 66.0
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        self.assertTrue(SA.shape == SP.shape)

    def test_broadcast(self):
        """Follow usual broadcast rules"""
        SP = np.array([[30.0, 31.0, 32.0], [33.0, 34.0, 35.0]])
        p  = [10.0, 20.0, 30.0]
        lon, lat = 2.0, 66.0
        common_shape = np.broadcast(SP,p,lon,lat).shape
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        self.assertTrue(SA.shape == common_shape)
        
    def test_incompatible(self):
        """Raise ValueError for incompatible shapes"""
        SP  = np.array([[30.0, 31.0, 32.0], [33.0, 34.0, 35.0]])
        p1  = np.array([10.0, 20.0])
        p2  = np.array([10.0, 20.0, 30.0, 40.0])
        lon, lat = 2.0, 66.0
        self.assertRaises(ValueError, gsw.SA_from_SP, SP, p1, lon, lat)
        self.assertRaises(ValueError, gsw.SA_from_SP, SP, p2, lon, lat)


    def test_nan(self):
        """Accept and return not-a-number, do not mess up finite values"""
        SP = np.array([34.0, np.nan, 35.0])
        p = 100.0
        lon, lat = 2, 66
        SA  = gsw.SA_from_SP(SP, p, lon, lat)
        SA0 = gsw.SA_from_SP(SP[0], p, lon, lat)
        # Should return NaN for NaN
        self.assertTrue(np.isnan(SA[1]))
        # Should return correct value for not NaN
        self.assertEqual(SA[0], SA0)

    def test_masked(self):
        """Accept and return correctly masked arrays"""
        SP = np.array([33, 34, -99, 35])  # One salinity value missing
        p  = [0, 10, 50, 100]
        lon, lat = 2, 66
        SP = np.ma.masked_where(SP < 0, SP) 
        SA  = gsw.SA_from_SP(SP, p, lon, lat)
        SA0 = gsw.SA_from_SP(SP[0], p[0], lon, lat)
        # Return array should have the same mask
        self.assertTrue(np.all(SA.mask == SP.mask))
        # Correct value for non-masked entries
        self.assertEqual(SA[0], SA0)



    def test_correct(self):
        r"""Test correct answer for some standard values

        Standard values from 
        http://www.teos-10.org/pubs/gsw/html/gsw_SA_from_SP.html
        """

        SP = [34.5487, 34.7275, 34.8605, 34.6810, 34.5680, 34.5600]
        p = [10, 50, 125, 250, 600, 1000]
        lon, lat = 188, 4
        SA_official = np.array((34.711779712893147,
                                34.891523721622328,
                                35.025547737221643,
                                34.847230081708567,
                                34.736629598766619,
                                34.732361864645107))
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        self.assertTrue(np.all(SA == SA_official))

    def test_in_baltic(self):
        """Test that the function works in the baltic"""

        SP = [6.5683, 6.6719, 6.8108, 7.2629, 7.4825, 10.2796]
        lon, lat = 20, 59
        SA_correct = np.array([6.66994543, 6.77377643,  6.91298614,
                               7.36609419, 7.58618384, 10.38952057])
        p = 20  # Arbitrary value
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        self.assertTrue(np.all(abs(SA - SA_correct) < 1.0E-8))

if __name__ == '__main__':
    unittest.main()

  




  
