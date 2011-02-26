# -*- coding: utf-8 -*-

"""Unit test for entropy_from_pt in gibbs seawater module"""

import unittest
import numpy as np
import seawater.gibbs as gsw


class Test_entropy_from_pt(unittest.TestCase):
    """Test the entropy_from_pt function"""

    def test_masked(self):
        SA = np.array([35, -99])
        SA = np.ma.masked_where(SA < 0, SA)
        pt = 10.0
        s = gsw.entropy_from_pt(SA, pt)
        self.assertTrue(np.all(s.mask == SA.mask))
        
    def test_standard(self):
        SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
        pt = [28.7832, 28.4210, 22.7850, 10.2305, 6.8292, 4.3245]
        s = gsw.entropy_from_pt(SA, pt)
        s_correct = np.array([400.38946744, 395.43839949, 319.86743859,
                              146.79054828,  98.64691006,  62.79135672])
        self.assertTrue(np.all(abs(s - s_correct) < 1.0e-8))




if __name__ == '__main__':
    unittest.main()







