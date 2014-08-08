#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools.command.test import test as TestCommand

import io
import sys

import re
VERSIONFILE = "seawater/__init__.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.txt', 'CHANGES.txt')


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['--strict', '--verbose',
                          '--tb=long', 'seawater/test']
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

source = 'http://pypi.python.org/packages/source/s/seawater'
download_url = '%s/seawater-%s.tar.gz' % (source, verstr)


README = open('README.txt').read()
CHANGES = open('CHANGES.txt').read()
LICENSE = open('LICENSE.txt').read()

setup(
    name='seawater',
    version=verstr,
    url='http://pypi.python.org/pypi/seawater/',
    license='MIT',
    author='Filipe Fernandes',
    tests_require=['pytest'],
    install_requires=['numpy'],
    extras_require={'testing': ['pytest', 'scipy', 'oct2py']},
    cmdclass={'test': PyTest},
    author_email='ocefpaf@gmail.com',
    description='Seawater Library for Python',
    long_description=long_description,
    packages=['seawater'],
    include_package_data=True,
    platforms='any',
    test_suite='seawater.test.test_result_comparison',
    classifiers=['Programming Language :: Python',
                 'Development Status :: 6 - Mature',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: MIT License',
                 'Operating System :: OS Independent',
                 'Topic :: Scientific/Engineering',
                 ],
    use_2to3=True,
    maintainer='Filipe Fernandes',
    maintainer_email='ocefpaf@gmail.com',
    download_url=download_url,
    keywords=['oceanography', 'seawater'],
    )
