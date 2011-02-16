#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

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
      author           = 'Filipe Fernandes, Eric Firing, Ådlandsvik Bjørn',
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
      classifiers      = filter(None, classifiers.split("\n")),
      platforms        = 'any',
      packages = find_packages(),
      package_data     = {'':['data/*.npz']},
      zip_safe = False,
      license          = 'MIT',
      keywords         = 'oceanography seawater',
    )
