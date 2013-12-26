# -*- coding: utf-8 -*-
#
# test_input_shapes.py
#
# purpose:  Test gpan with various input shapes.
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.github.io/
# created:  06-Aug-2013
# modified: Tue 06 Aug 2013 09:54:30 AM BRT
#
# obs:
#

import unittest

import numpy as np
import seawater as sw
from scipy.io import loadmat

# Original data.
mat = loadmat('130624_stdAndNoENSO_smoothed_analysis.mat')

cp = mat.get('cp')
dens = mat.get('dens')
t_mean = mat.get('t_mean')
s_mean_comp = mat.get('s_mean_comp')
pt_mean_comp = mat.get('pt_mean_comp')
steric_height = mat.get('steric_height')

lon = mat.get('lon').astype(float).squeeze()
lat = mat.get('lat').astype(float).squeeze()
depth = mat.get('depths').astype(float).squeeze()

# Using first dimension from the data-set.
s_mean, pt_mean = s_mean_comp[0], pt_mean_comp[0]


def test_array(arr1, arr2):
    try:
        np.testing.assert_equal(arr1, arr2)
    except AssertionError:
        return False
    return True

if __name__ == '__main__':
    # 1D.
    steric_height_new = sw.gpan(s_mean[:, 0, 0], t_mean[:, 0, 0], depth)
    if test_array(steric_height[:, 0, 0], steric_height_new):
        print("1D sw.gpan test passed.")
    else:
        print("1D sw.gpan test failed.")

    # 2D.
    steric_height_new = sw.gpan(s_mean[..., 1], t_mean[..., 1],
                                depth[..., None])
    if test_array(steric_height[..., 1], steric_height_new):
        print("2D sw.gpan test passed.")
    else:
        print("2D sw.gpan test failed.")

    # 3D.
    steric_height_new = sw.gpan(s_mean, t_mean, depth[..., None, None])

    if test_array(steric_height, steric_height_new):
        print("3D sw.gpan test passed.")
    else:
        print("3D sw.gpan test failed.")
