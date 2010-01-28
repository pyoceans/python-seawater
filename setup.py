#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

classifiers = """\
Development Status :: 4 - Beta
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
      version          = '1.0.2',
      author           = 'Filipe Fernandes',
      author_email     = 'ocefpaf@gmail.com',
      maintainer       = 'Filipe Fernandes',
      maintainer_email = 'ocefpaf@gmail.com',
      url              = 'http://ocefpaf.tiddlyspot.com/#python-seawater',
      description      = 'Seawater Libray for Python',
      long_description = """\
Pythoni Seawater: Codiname Jörmungandr
--------------------------------------

This module is a translation of the original SEAWATER-3.2 MATLAB toolkit routines for calculating the properties of sea water. It consists of a self contained library and is extremely easy to use. This version requires NumPy.
""",
      download_url     = 'http://pypi.python.org/packages/source/s/seawater/seawater-1.0.2.tar.gz',
      packages         = ['seawater'],
      classifiers      = filter(None, classifiers.split("\n")),
      platforms        = ['any'],
      license          = 'MIT',
      keywords         = 'oceanography seawater',
      requires         = 'NumPy'
    )
