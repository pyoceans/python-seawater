#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

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

setup(name = 'seawater',
      version = '2.0.1',
      packages = ['seawater','seawater/gibbs','seawater/csiro', 'seawater/extras', 'seawater/extras/waves','seawater/extras/sw_extras'],
      package_data = {'':['gibbs/data/*.npz']},
      license = 'MIT',
      description = 'Seawater Libray for Python',
      long_description = open('README.txt').read(),
      author = 'Filipe Fernandes, Eric Firing, Ådlandsvik Bjørn',
      author_email = 'ocefpaf@gmail.com',
      maintainer = 'Filipe Fernandes',
      maintainer_email = 'ocefpaf@gmail.com',
      url = 'http://pypi.python.org/pypi/seawater/',
      download_url = 'http://pypi.python.org/packages/source/s/seawater/',
      classifiers = filter(None, classifiers.split("\n")),
      platforms = 'any',
      zip_safe = False,
      keywords = 'oceanography seawater',
    )