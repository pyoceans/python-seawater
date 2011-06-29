#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

from distutils.command.sdist import sdist
import os

class sdist_hg(sdist):

    user_options = sdist.user_options + [
            ('dev', None, "Add a dev marker")
            ]

    def initialize_options(self):
        sdist.initialize_options(self)
        self.dev = 0

    def run(self):
        if self.dev:
            suffix = '.dev%d' % self.get_tip_revision()
            self.distribution.metadata.version += suffix
        sdist.run(self)

    def get_tip_revision(self, path=os.getcwd()):
        from mercurial.hg import repository
        from mercurial.ui import ui
        from mercurial import node
        repo = repository(ui(), path)
        tip = repo.changelog.tip()
        return repo.changelog.rev(tip)


try: # Python 3
  from distutils.command.build_py import build_py_2to3 as build_py
except ImportError: # Python 2
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

setup(name = 'seawater',
      version = '2.0.2',
      packages = ['seawater', 'seawater/gibbs', 'seawater/csiro',
                  'seawater/extras', 'seawater/extras/waves',
                  'seawater/extras/sw_extras', 'seawater/test',
                  'seawater/gibbs3'],
      package_data = {'':['gibbs/data/*.npz']},
      license = 'LICENSE.txt',
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
      cmdclass = {'build_py': build_py},
      #cmdclass = {'sdist': sdist_hg}, # NOTE: python setup.py sdist --dev
      keywords = ['oceanography', 'seawater'],
     )
