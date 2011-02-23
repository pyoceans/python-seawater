# -*- coding: utf-8 -*-

"""Unit test for interface to the _gibbs function
   in the Gibbs seawater library
"""

import unittest
import numpy as np
from seawater.library import _gibbs

class Test_gibbs(unittest.TestCase):
    """Test the _gibbs function"""

    def test_legal_derivatives(self):
        """Test that derivatives up to second order are accepted"""
        SA, t, p = 35, 10, 100
        # Zeroth order
        g000 = _gibbs(0, 0, 0, SA, t, p)
        # First order
        g100 = _gibbs(1, 0, 0, SA, t, p)
        g010 = _gibbs(0, 1, 0, SA, t, p)
        g001 = _gibbs(0, 0, 1, SA, t, p)
        # Second order
        g110 = _gibbs(1, 1, 0, SA, t, p)
        g110 = _gibbs(1, 0, 1, SA, t, p)
        g110 = _gibbs(0, 1, 1, SA, t, p)
        g110 = _gibbs(2, 0, 0, SA, t, p)
        g110 = _gibbs(0, 2, 0, SA, t, p)
        g110 = _gibbs(0, 0, 2, SA, t, p)

    def test_third_order_derivatives(self):
        """Test some illegal derivatives of third order"""
        SA, t, p = 35, 10, 100
        self.assertRaises(ValueError, _gibbs, 1, 1, 1, SA, t, p)
        self.assertRaises(ValueError, _gibbs, 0, 2, 1, SA, t, p)
        self.assertRaises(ValueError, _gibbs, 3, 0, 0, SA, t, p)

    def test_non_integer_derivatives(self):
        """Non-integer derivatives are not allowed"""
        SA, t, p = 35, 10, 100
        ns, nt, npr = 0.5, 1.0, 0.5
        self.assertRaises(ValueError, _gibbs, ns, nt, npr, SA, t, p)

    def test_return_masked(self):
        """Return values should be masked arrays of rank >= 1"""
        # Scalar input
        SA, t, p = 35, 10, 100
        g101 = _gibbs(1, 0, 1, SA, t, p)
        self.assertTrue(np.ma.isMaskedArray(g101))
        self.assertTrue(g101.ndim > 0)

        # 0D arrays
        SA = np.array(35.0)
        t = np.array(10)
        p = np.array(100)
        g101 = _gibbs(1, 0, 1, SA, t, p)
        self.assertTrue(np.ma.isMaskedArray(g101))
        self.assertTrue(g101.ndim > 0)

        # Array and list input
        SA = np.array([31, 33, 35, 35.2])
        t = [12, 11, 8, 7.7]
        p = 100
        g002 = _gibbs(0, 0, 2, SA, t, p)
        self.assertTrue(np.ma.isMaskedArray(g002))

    def test_broadcasting(self):
        """Handle 2D array and broadcasting correctly"""
        SA = np.array([[31,33],[34, 34.5]])
        t  = [10, 15]
        p = 100
        g100 = _gibbs(1, 0, 0, SA, t, p)
        self.assertEqual(g100.shape, (2,2))
        
    def test_non_compatible(self):
        """Return ValueError if not broadcastable input"""
        SA = np.array([[31,33],[34, 34.5]])
        t  = [10, 15, 20]
        p = np.array([[0, 1000]])
        self.assertRaises(ValueError,  _gibbs, 0, 0, 0, SA, t, p)       

    def test_nan(self):
        """Accept and return not-a-number, do not mess up finite values"""
        SA = np.array([34.0, np.nan, 35.0])
        t = 10.0
        p = 0
        g101 = _gibbs(1, 0, 1, SA, t, p)
        # Should return NaN for NaN
        self.assertTrue(np.isnan(g101[1]))
        # Should return correct value for not NaN
        self.assertEqual(g101[0], _gibbs(1, 0, 1, SA[0], t, p))

    def test_masked(self):
        """Accept and return correctly masked arrays"""
        SA = np.array([33, 34, -99, 35])  # One salinity value missing
        t  = [18, 17,  12, 8]
        p  = [0, 10, 50, 100]
        SA = np.ma.masked_less(SA, 0) # Make masked array
        g020 = _gibbs(0, 2, 0, SA, t, p)
        g101 = _gibbs(1, 0, 1, SA, t, p)
        # Return array should have the same mask
        self.assertTrue(np.all(g020.mask == SA.mask))
        self.assertTrue(np.all(g101.mask == SA.mask))
        # Correct value for non-masked entries
        self.assertEqual(g020[0], _gibbs(0, 2, 0, SA[0], t[0], p[0]))

    def test_all_masked(self):
        """Handle totally masked arrays"""
        SA = np.ma.MaskedArray([-99, 33], mask=[True, True])
        t = 10
        p = 100
        g001 = _gibbs(0, 0, 1, SA, t, p)
        g110 = _gibbs(1, 1, 0, SA, t, p)
        self.assertTrue(np.all(g001.mask))
        self.assertTrue(np.all(g110.mask))

    def test_zero_salinity(self):
        """Handle zero salinity correctly"""

        # If ns > 0, SA = 0 should be masked
        #   at least matlab explicitly puts a NaN in this case
        SA, t, p = [30, 0], 10, 100
        g110 = _gibbs(1, 1, 0, SA, t, p)
        self.assertTrue(g110.mask[1])    
        
        # If ns = 0, SA = 0 should return an unmasked value
        #   according to matlab
        SA, t, p = 0, 10, 100
        g011 = _gibbs(0, 1, 1, SA, t, p)
        self.assertFalse(np.all(g011.mask))   
        self.assertTrue(np.isfinite(g011[0])) 
                
if __name__ == '__main__':
    unittest.main()







