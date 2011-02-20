# -*- coding: utf-8 -*-

"""Perform all unit tests in the present dictionary"""

import os
import glob

# Find all test scripts
test_scripts = glob.glob('test*.py')
test_scripts.remove('test_all.py')  # test_all.py is not a unit-test

# Perform the unit-tests
for script in test_scripts:
    print 'running %s' % script
    os.system('python %s' % script)
    


