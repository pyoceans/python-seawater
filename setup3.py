#!/usr/bin/env python
# -*- coding: utf-8 -*-

# setup file for the seawater package,
# works with both python2 and python3
# Based on distutils
#
# To install under python3 do:
# python3 setup3.py build
# sudo python3 setup3.py install
#


from distutils.core import setup
try:                 # Python 3
  from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:  # Python 2
  from distutils.command.build_py import build_py


classifiers = """\
Development Status :: 5 - Production/Stable
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
License :: OSI Approved :: MIT License
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Education
Topic :: Software Development :: Libraries :: Python Modules
"""

setup(name             = 'seawater',
      version          = '2.0.0',
      author           = 'Filipe Fernandes, Eric Firing, Bjørn Ådlandsvik',
      author_email     = 'ocefpaf@gmail.com',
      maintainer       = 'Filipe Fernandes',
      maintainer_email = 'ocefpaf@gmail.com',
      url              = 'http://pypi.python.org/pypi/seawater/',
      description      = 'Seawater Libray for Python',
      long_description = """\
This python package contains a python translation for two Matlab user
contributed toolboxes. (1) the
`seawater <http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm>`_ (EOS-80)
and (2) the `gibbs <http://www.teos-10.org/software.htm>`_ seawater (TEOS-10).
""",
      download_url     = 'http://pypi.python.org/packages/source/s/seawater/',
      #classifiers      = filter(None, classifiers.split("\n")),
      platforms        = 'any',
      cmdclass = {'build_py': build_py},
      #packages = find_packages(),
      packages = ['seawater', 'seawater.gibbs', 'seawater.csiro'],
      package_data     = {'':['data/*.npz']},
      zip_safe = False,
      license          = 'MIT',
      keywords         = 'oceanography seawater',
    )
