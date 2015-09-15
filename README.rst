python-seawater
===============


.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.11395.png
   :target: http://dx.doi.org/10.5281/zenodo.11395
   :alt: DOI
.. image:: http://img.shields.io/pypi/v/seawater.svg?style=flat
   :target: https://pypi.python.org/pypi/seawater
   :alt: Version_status
.. image:: http://img.shields.io/pypi/dm/seawater.svg?style=flat
   :target: https://pypi.python.org/pypi/seawater
   :alt: Downloads
.. image:: http://img.shields.io/travis/pyoceans/python-seawater/master.svg?style=flat
   :target: https://travis-ci.org/pyoceans/python-seawater
   :alt: Build_status
.. image:: http://img.shields.io/badge/license-MIT-blue.svg?style=flat
   :target: https://github.com/pyoceans/python-seawater/blob/master/LICENSE.txt
   :alt: license
.. image:: http://bottlepy.org/docs/dev/_static/Gittip.png
   :target: https://gratipay.com/~ocefpaf/
   :alt: Gittip


This is a Python re-write of the CSIRO seawater toolbox
(`SEAWATER-3.3 <http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm>`__)
for calculating the properties of sea water. The package uses the
formulas from Unesco's joint panel on oceanographic tables and
standards, UNESCO 1981 and UNESCO 1983 (EOS-80).

The EOS-80 library is considered now obsolete; it is provided here for
compatibility with old scripts, and to allow a smooth transition to the
new `TEOS-10 <http://www.teos-10.org/>`__.

Warning
-------

The Python version default output unit for sw.dist is 'km' instead of
'nm'.

Here we assume pressure as the first dimension, i.e. M pressure by N
positions (See the table below).  The MatlabTM version does some guessing
at this that we simply ignore to avoid confusions.

+---------+-----------+-----------+
| **P**   | **S**     | **T**     |
+=========+===========+===========+
| 10      | 34.5487   | 28.7856   |
+---------+-----------+-----------+
| 50      | 34.7275   | 28.4329   |
+---------+-----------+-----------+
| 125     | 34.8605   | 22.8103   |
+---------+-----------+-----------+
| 250     | 34.6810   | 10.2600   |
+---------+-----------+-----------+
| 600     | 34.5680   | 6.8863    |
+---------+-----------+-----------+
| 1000    | 34.5600   | 4.4036    |
+---------+-----------+-----------+

The current version was tested against the MatlabTM seawater v3.3 reproducing
all functions and results from that release.

More information at http://pythonhosted.org/seawater
