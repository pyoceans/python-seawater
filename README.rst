python-seawater
===============

.. image:: https://badge.fury.io/py/seawater.png
   :target: http://badge.fury.io/py/seawater
.. image:: https://api.travis-ci.org/pyoceans/python-seawater.png?branch=master
   :target: https://travis-ci.org/pyoceans/python-seawater
.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.11395.png
   :target: http://dx.doi.org/10.5281/zenodo.11395
.. image:: http://bottlepy.org/docs/dev/_static/Gittip.png
   :target: https://www.gittip.com/ocefpaf/

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
