# -*- coding: utf-8 -*-

"""Perform all unit tests in the present dictionary"""

# Last modified 2011-04-06 
# This script should now work both with python2 and python3

import sys
import os
import glob

# Set executables for python2 and python3
python2 = 'python'   # python 2.x
python3 = 'python3'  # python 3.x

# Choose the correct executable
if sys.version[0] == 3:   # Version 3.x
    python = python3
else:
    python = python2

# Find all test scripts
test_scripts = glob.glob('test*.py')
test_scripts.remove('test_all.py')  # test_all.py is not a unit-test

# Perform the unit-tests and accumulate errors
errors = []
for script in test_scripts:
    sys.stdout.write('running %s\n' % script)
    a = os.system('%s %s' % (python, script))
    if a != 0:
        errors.append(script)

# Display a summary
sys.stdout.write("\n --- test_all.py finished ---\n")
sys.stdout.write("  Ran %d test scripts\n" % len(test_scripts))
if not errors:
    sys.stdout.write("  no problems\n")
else:
    sys.stdout.write("  problems in :\n")
    sys.stdout.write(str(errors))
