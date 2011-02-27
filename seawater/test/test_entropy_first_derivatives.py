# -*- coding: utf-8 -*-

"""Unit test for entropy_first_derivatives in the gibbs seawater module"""

# Bjørn Ådlandsvik <bjorn@imr.no>
# 2011-02-26

import unittest
import numpy as np
import seawater.gibbs as gsw

class Test_entropy_first_derivatives(unittest.TestCase):
    """Test the entropy_first_derivatives function"""

    def test_scalar(self):
        """Should return a two-tuple of scalars for scalar input"""
        SA = 35.0
        CT = 10.0
        output = gsw.entropy_first_derivatives(SA, CT)
        self.assertEqual(type(output), tuple)
        self.assertEqual(len(output), 2)
        self.assertTrue(np.isscalar(output[0]))
    
    def test_standard_values(self):
        """Test some standard values"""
        SA = [34.7118, 34.8915, 35.0256, 34.8472, 34.7366, 34.7324]
        CT = [28.8099, 28.4392, 22.7862, 10.2262,  6.8272,  4.3236]

        eta_SA_standard = np.array((-0.263286800711655,
                                    -0.263977276574528,
                                    -0.255367497912925,
                                    -0.238066586439561,
                                    -0.234438260606436,
                                    -0.232820684341694))
        eta_CT_standard = np.array((13.221031210083824,
                                    13.236911191313675,
                                    13.489004628681361,
                                    14.086599016583795,
                                    14.257729576432077,
                                    14.386429945649411))
        eta_SA, eta_CT = gsw.entropy_first_derivatives(SA, CT)
        self.assertTrue(np.all(abs(eta_SA - eta_SA_standard) < 1.0e-15))
        self.assertTrue(np.all(abs(eta_CT - eta_CT_standard) < 1.0e-16))



if __name__ == '__main__':
    unittest.main()













