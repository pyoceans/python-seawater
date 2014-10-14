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

import numpy as np
import seawater as sw

# Load data.
with np.load('shapes.npz') as shapes:
    cp = shapes['cp']
    lon = shapes['lon']
    lat = shapes['lat']
    dens = shapes['dens']
    depth = shapes['depth']
    s_mean = shapes['s_mean']
    t_mean = shapes['t_mean']
    pt_mean = shapes['pt_mean']
    steric_height = shapes['steric_height']


# 1D.
steric_height_new = sw.gpan(s_mean[:, 0, 0], t_mean[:, 0, 0], depth)
np.testing.assert_array_almost_equal(steric_height[:, 0, 0], steric_height_new)

# 2D.
steric_height_new = sw.gpan(s_mean[..., 1], t_mean[..., 1],
                            depth[..., None])
np.testing.assert_array_almost_equal(steric_height[..., 1], steric_height_new)

# 3D.
steric_height_new = sw.gpan(s_mean, t_mean, depth[..., None, None])
np.testing.assert_array_almost_equal(steric_height, steric_height_new)
