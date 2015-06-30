import os
import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['--verbose']
        self.test_suite = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


def extract_version(module='seawater'):
    version = None
    fdir = os.path.dirname(__file__)
    fnme = os.path.join(fdir, module, '__init__.py')
    with open(fnme) as fd:
        for line in fd:
            if (line.startswith('__version__')):
                _, version = line.split('=')
                # Remove quotation characters.
                version = version.strip()[1:-1]
                break
    return version


rootpath = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    return open(os.path.join(rootpath, *parts), 'r').read()


long_description = '{}\n{}'.format(read('README.rst'), read('CHANGES.txt'))
LICENSE = read('LICENSE.txt')

source = 'http://pypi.python.org/packages/source/s/seawater'

setup(name='seawater',
      version=extract_version(),
      license=LICENSE,
      long_description=long_description,
      classifiers=['Programming Language :: Python',
                   'Development Status :: 6 - Mature',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: MIT License',
                   'Operating System :: OS Independent',
                   'Topic :: Scientific/Engineering',
                   ],
      description='Seawater Library for Python',
      author='Filipe Fernandes',
      maintainer='Filipe Fernandes',
      maintainer_email='ocefpaf@gmail.com',
      url='https://github.com/pyoceans/python-seawater/',
      download_url='https://pypi.python.org/pypi/seawater',
      platforms='any',
      keywords=['oceanography', 'seawater'],
      install_requires=['numpy'],
      tests_require=['pytest', 'scipy', 'oct2py'],
      cmdclass=dict(test=PyTest),
      packages=['seawater'],
      include_package_data=True,
      )
