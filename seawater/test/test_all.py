#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Perform all unit tests in the present dictionary"""

import os
import glob

# Find all test scripts
test_scripts = glob.glob('test*.py')
test_scripts.remove('test_all.py')  # test_all.py is not a unit-test


errors = []
# Perform the unit-tests
for script in test_scripts:
    print 'running %s' % script
    a = os.system('python %s' % script)
    if a != 0:
        errors.append(script)

print "\n --- test_all.py finished ---"

#
print "  Ran %d test scripts" % len(test_scripts)
if not errors:
    print "  no problems"
else:
    print "  problems in :"
    print errors
    
    


