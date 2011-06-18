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
# modified: Sat 18 Jun 2011 02:15:27 PM EDT
#
# obs:
#

import scipy.io as sio
import numpy as np

data_ver = 'v3_0'
gsw_data = sio.loadmat('gsw_data_'+data_ver+'.mat', squeeze_me=True)

""" save compare values in a separate file
"""
gsw_cv = gsw_data['gsw_cv']
del gsw_data['gsw_cv']

cv_vars = []
cv_exact = []
for name in gsw_cv.dtype.names:
    cmd = ('%s = np.atleast_1d(gsw_cv[\'%s\'])[0]' % (name,name))
    exec(cmd)
    if 'exact' in name:
        cv_exact.append('%s=%s' % (name,name))
    else:
        cv_vars.append('%s=%s' % (name,name))

# check values
var_list = ', '.join(cv_vars)
exec('np.savez("gsw_cv_%s", %s)' % (data_ver, var_list))

# exact values (separated because np.savez does not accept > 255 args)
var_list = ', '.join(cv_exact)
exec('np.savez("gsw_exact_%s", %s)' % (data_ver, var_list))

# delta SA Atlas
ref_table = []
for k in gsw_data:
    if '__' not in k:
        if 'deltaSA_ref' in k:
             name =  'delta_SA_ref'
        else:
            name = k
        
        cmd = ('%s = np.atleast_1d(gsw_data[\'%s\'])' % (name, k))
        ref_table.append('%s=%s' % (name, name))
        exec(cmd)

var_list = ', '.join(ref_table)
exec('np.savez("gsw_data_%s", %s)' % (data_ver, var_list))

"""
 #  GSW_DATA_V3_0      GSW_DATA_V2_0
 sigma_2_ref_cast <--> same
      SA_ref_cast <--> same
         lats_ref <--> same
     lat_ref_cast <--> same
      deltaSA_ref <--> delta_SA_ref
       ndepth_ref <--> same
            p_ref <--> same
        ocean_ref <--> same
        longs_ref <--> same
    long_ref_cast <--> same
 gamma_n_ref_cast <--> same
       p_ref_cast <--> same
      CT_ref_cast <--> same

  # NEW VARIABLES:
         SAAR_ref <--> NA
           SR_ref <--> NA
     version_date <--> NA
   version_number <--> NA
    gsw_demo_data <--> NA

   # CHECK VALUE:
           gsw_cv <--> TODO
"""
