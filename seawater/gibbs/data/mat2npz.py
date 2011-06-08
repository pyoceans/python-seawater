#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# mat2npz.py
#
# purpose:  Convert matlab file from TEOS-10 group to a npz file
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  06-Jun-2011
# modified: Mon 06 Jun 2011 07:07:06 PM EDT
#
# obs:
#

import scipy.io as sio
import numpy as np

gsw_data_v3_0 = sio.loadmat('gsw_data_v3_0.mat', squeeze_me=True) 

# save compare values in a separate file
gsw_cv = gsw_data_v3_0['gsw_cv']
del gsw_data_v3_0['gsw_cv']

# save files
def dict2npz(adict, npzfile):
    """Save a dict as a npz file."""
    for k in adict:
        #eval(k+'='+adict[k])
        print(k+'='+adict[k])


"""
 # gsw_data_v3_0       gsw_data_v2_0
 sigma_2_ref_cast
      SA_ref_cast
         lats_ref
     lat_ref_cast
      deltaSA_ref <--> delta_SA_ref
       ndepth_ref
            p_ref
        ocean_ref
        longs_ref
    long_ref_cast
 gamma_n_ref_cast
       p_ref_cast
      CT_ref_cast

     # new value:
         SAAR_ref
           SR_ref
     version_date
   version_number
    gsw_demo_data

   # check value:
           gsw_cv
"""
