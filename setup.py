#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

classifiers = """\
Development Status :: beta
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

setup(name='seawater',
    version='1.0-3.2',
    description='Seawater Libray for Python',
    long_description = """\
    This python 1.0-3.2 is a translation of the original SEAWATER 3.2 MATLAB toolkit routines for calculating the
    properties of sea water. They are a self contained library and are extremely easy to use.""",
    author='Filipe Fernandes',
    author_email='ocefpaf@gmail.com',
    url='http://http://ocefpaf.tiddlyspot.com/',
    packages=['seawater'],
    classifiers=filter(None, classifiers.split("\n")),
    keywords='oceanography seawater',
    download_url = 'http://cheeseshop.python.org/packages/source/seawater.tar.gz',
    license = 'MIT',
    platforms = ['any'],
     )