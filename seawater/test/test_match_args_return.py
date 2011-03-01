# -*- coding: utf-8 -*-

"""Unit test for the match_args_return decorator"""

# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-03-01

import unittest
import numpy as np
from seawater.library import match_args_return

class Test_match_args_return(unittest.TestCase):

    def setUp(self):
        # Define a simple function to work with
        # and decorate it explicitly
        def sum3(x, y, z):
            return x+y+z
        self.Sum3 = match_args_return(sum3)

    def test_scalar_input(self):
        """Scalar input should give float scalar output"""

        out = self.Sum3(2, 3, 4)

        self.assertTrue(np.isscalar(out))
        # The scalar has type numpy.float64 not
        # python's ordinary float
        # This behaviour should perhaps be changed
        #self.assertEqual(type(out), float)  # Gives failure
        self.assertEqual(type(out), np.float64)

    def test_list_input(self):
        """Accept lists and return array"""

        out = self.Sum3([1,1,1], [2,2,2], [3,3,3])

        self.assertEqual(type(out), np.ndarray)
        self.assertEqual(out.dtype, float)
        self.assertEqual(out.shape, (3,))

    def test_array_input(self):
        """Accept and return arrays"""

        x = np.array([[2.0, 3.14], [1.1, 2.2]])
        y = 3.0
        z = np.array([1, 1])
        out = self.Sum3(x, y, z)

        self.assertEqual(type(out), np.ndarray)
        self.assertEqual(out.dtype, float)
        self.assertEqual(out.shape, (2,2))
        
    def test_masked_input(self):
        """Accept and return masked arrays,
           without messing up non-masked values"""

        x = np.ma.MaskedArray([1, 2, 3, 4],
                              mask=[False, True, False, False])
        y = [1, 2, 3, 4]
        z = np.ma.MaskedArray([1, 2, 3, 4],
                              mask=[False, False, False, True])
        out = self.Sum3(x, y, z)

        # Return a masked array
        self.assertTrue(np.ma.isMaskedArray(out))
        # Output mask is correct
        self.assertTrue(np.all(out.mask == [False,True,False,True]))
        # Unmasked values OK
        self.assertEqual(out[0], 3)
        self.assertEqual(out[2], 9)

    def test_nan_input(self):
        """Propagate nan's without compromizing finite values"""
        x = np.array([1, np.nan, 3, 4])
        y = [1, 2, 3, 4]
        z = [1, 2, 3, np.nan]
        out = self.Sum3(x, y, z)

        # Return value is an array
        self.assertEqual(type(out), np.ndarray)
        # nan-s in correct places
        self.assertTrue(np.isnan(out[1]))
        self.assertTrue(np.isnan(out[3]))
        # Uncontaminated values should work
        self.assertEqual(out[0], 3)
        self.assertEqual(out[2], 9)

    def test_projection(self):
        """Test propagation of non-used nans and masks"""
        def f(x, y, z):  # Project onto first argument
            return x
        Pr0 = match_args_return(f)
        x = [1, 2, 3, 4]
        y = np.array([1, np.nan, 3, 4])
        z = np.ma.MaskedArray([1, 2, 3, 4],
                              mask=[False, False, False, True])
        out = Pr0(x, y, z)

        # Returns a masked array instead of ndarray
        # But the mask is correct (False everywhere)
        # self.assertEqual(type(out), np.ndarray) # gives failure
        self.assertTrue(np.ma.isMaskedArray(out))
        self.assertFalse(np.any(out.mask))

        # Check correct values, in particular no nans.
        self.assertTrue(np.all(x == out))

        
if __name__ == '__main__':
    unittest.main()
    
