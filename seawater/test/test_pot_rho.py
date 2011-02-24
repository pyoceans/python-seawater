# -*- coding: utf-8 -*-

"""Unit test for potential density function in gibbs seawater module"""

# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-02-07

import unittest
import numpy as np
import seawater.gibbs as gsw

class Test_pot_rho(unittest.TestCase):
    """Test the rho function"""

    def test_scalar_input(self):
        """Accept scalar input, return a scalar"""
        SA = 35.0
        t  = 10.0
        p  = 1000.0
        pr = 0
        rho = gsw.pot_rho(SA, t, p, pr)
        self.assertTrue(np.isscalar(rho))

    def test_array_input(self):
        """Accept array input, return broadcasted shape"""

        # Test 1D
        SA = np.array([30,32,34,36])
        t  = 10.0
        p  = np.array([0, 100, 1000, 4000])
        rho = gsw.pot_rho(SA, t, p)
        self.assertTrue(rho.shape == np.broadcast(SA,t,p).shape)
        
        # Test 2D
        SA = np.array([30,32,34,36])
        t = np.array([[2,4,6,8], [12,14,16,18]])
        p = np.array([100.0])
        rho = gsw.pot_rho(SA, t, p)
        self.assertTrue(rho.shape == np.broadcast(SA,t,p).shape)

        # Test 3D
        SA = 30 + np.linspace(0,1,24).reshape((2,3,4))
        t = 18.0
        p = 0
        pr = 4000
        rho = gsw.pot_rho(SA, t, p, pr)
        self.assertTrue(rho.shape == np.broadcast(SA,t,p).shape)
    
    def test_list_input(self):
        """Lists may be used instead of arrays on input"""

        # Test 1D
        SA = [30, 32, 34, 36]
        t  = 10.0
        p  = [0, 100, 1000, 4000]
        pr = 1000.0
        rho = gsw.pot_rho(SA, t, p, pr)
        self.assertTrue(rho.shape == (4,))
        
        # Test 2D
        SA = [30,32,34,36]
        t  = [[2,4,6,8], [12,14,16,18]]
        p  = 100.0
        rho = gsw.pot_rho(SA, t, p)
        self.assertTrue(rho.shape == (2,4))

    def test_non_compatible(self):
        """Raise ValueError if input not broadcastable to same shape"""
        SA = np.array([34, 35])
        t = 10.0
        p = np.array([10, 100, 1000])
        self.assertRaises(ValueError, gsw.rho, SA, t, p)

    def test_nan(self):
        """Accept and return not-a-number, do not mess up finite values"""
        SA = np.array([34.0, np.nan, 35.0])
        t = 10.0
        p = 0
        pr = 1000
        rho = gsw.pot_rho(SA, t, p, pr)
        rho0 = gsw.pot_rho(SA[0], t, p, pr)
        # Should return NaN for NaN
        self.assertTrue(np.isnan(rho[1]))
        # Should return correct value for not NaN
        self.assertEqual(rho[0], rho0)

    def test_masked(self):
        """Accept and return correctly masked arrays"""
        SA = np.array([33, 34, -99, 35])  # One salinity value missing
        t  = [18, 17,  12, 8]
        p  = [0, 10, 50, 100]
        SA = np.ma.masked_where(SA < 0, SA) # Make masked array
        rho = gsw.pot_rho(SA, t, p)
        rho0 = gsw.pot_rho(SA[0], t[0], p[0])
        # Return array should have the same mask
        self.assertTrue(np.all(rho.mask == SA.mask))
        # Correct value for non-masked entries
        self.assertEqual(rho[0], rho0)

    def test_standard(self):
        r"""Test some standard values
        
        Standard values from
            http://www.teos-10.org/pubs/gsw/html/gsw_rho.html
        """
        # 
        SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
        t  = [28.7856, 28.4329, 22.8103, 10.2600,  6.8863,  4.4036]
        p  = [     10,      50,     125,     250,     600,    1000]
        pr = 0
        
        rho = gsw.pot_rho(SA, t, p, pr)
        rho_standard = np.array((1021.798145811089,
                                 1022.052484416980,
                                 1023.893583651958,
                                 1026.667621124443,
                                 1027.107230868492,
                                 1027.409631264134))

        self.assertTrue(np.all(abs(rho-rho_standard) < 1.0e-12))


if __name__ == '__main__':
    unittest.main()







