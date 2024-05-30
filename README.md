# python-seawater

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11395.svg)](https://zenodo.org/records/11395)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/seawater)
[![Tests](https://github.com/pyoceans/python-seawater/actions/workflows/tests.yml/badge.svg)](https://github.com/pyoceans/python-seawater/actions/workflows/tests.yml)
![PyPI - License](https://img.shields.io/pypi/l/seawater)

[![No Maintenance Intended](https://unmaintained.tech/badge.svg)](http://unmaintained.tech/)

This is a Python re-write of the CSIRO seawater toolbox
([SEAWATER-3.3](https://www.cmar.csiro.au/datacentre/ext_docs/seawater.html))
for calculating the properties of sea water. The package uses the
formulas from Unesco's joint panel on oceanographic tables and
standards, UNESCO 1981 and UNESCO 1983 (EOS-80).

The EOS-80 library is considered now obsolete; it is provided here for
compatibility with old scripts, and to allow a smooth transition to the
new [TEOS-10](https://www.teos-10.org/).

## Warning

The Python version default output unit for sw.dist is *km* instead of
*nm*.

Here we assume pressure as the first dimension, i.e. M pressure by N
positions (See the table below). The Matlab version does some guessing
at this that we simply ignore to avoid confusions.

  | **P** | **S**   | **T**   |
  |------:|:-------:|:-------:|
  | 10    | 34.5487 | 28.7856 |
  | 50    | 34.7275 | 28.4329 |
  | 125   | 34.8605 | 22.8103 |
  | 250   | 34.6810 | 10.2600 |
  | 600   | 34.5680 | 6.8863  |
  | 1000  | 34.5600 | 4.4036  |

The current version was tested against the Matlab seawater v3.3
reproducing all functions and results from that release.
