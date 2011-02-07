# -*- coding: utf-8 -*-

"""Unit code test for interfaces to GSW functions"""


# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-02-07

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
        #self.assertTrue(np.isscalar(SA))
        self.assertTrue(np.isscalar(SA[0]))

    def test_single_array_input(self):
        """Accept and return single element arrays"""
        SP = np.array([35.0])
        p = np.array([0.0])
        lon, lat = np.array([2.0]), np.array([66.0])
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        #self.assertTrue(SA.shape == (1,))
        self.assertTrue(SA[0].shape == (1,))

    def test_list_input(self):
        """Accept list input"""
        SP = [30.0, 35.0]
        p  = [10.0, 10.0]
        lon, lat = [2.0, 2.0], [65.0, 66.0]
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        #self.assertTrue(SA.shape == (2,))
        self.assertTrue(SA[0].shape == (2,))[0]

    def test_array_shape(self):
        """Respect the shape of input arrays"""
        SP = np.array([[30.0, 31.0, 32.0], [33.0, 34.0, 35.0]])
        p  = np.array([[10.0, 20.0, 30.0], [10.0, 20.0, 30.0]])
        lon, lat = 2.0, 66.0
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        #self.assertTrue(SA.shape == SP.shape)
        self.assertTrue(SA[0].shape == SP.shape)

    def test_broadcast(self):
        """Follow usual broadcast rules"""
        SP = np.array([[30.0, 31.0, 32.0], [33.0, 34.0, 35.0]])
        p  = [10.0, 20.0, 30.0]
        lon, lat = 2.0, 66.0
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        self.assertTrue(SA.shape == (SP+p).shape)
        
    def test_incompatible(self):
        """Raise ValueError for incompatible shapes"""
        SP  = np.array([[30.0, 31.0, 32.0], [33.0, 34.0, 35.0]])
        p1  = np.array([10.0, 20.0])
        p2  = np.array([10.0, 20.0, 30.0, 40.0])
        lon, lat = 2.0, 66.0
        #SA = gsw.SA_from_SP(SP, p, lon, lat)
        #SA = gsw.SA_from_SP(SP, p, lon, lat)[0]
        self.assertRaises(ValueError, gsw.SA_from_SP, SP, p1, lon, lat)
        self.assertRaises(ValueError, gsw.SA_from_SP, SP, p2, lon, lat)

        
        



if __name__ == '__main__':
    unittest.main()

  




  
