Python Seawater
===============

Introduction:
-------------

This python package contains a python translation for two MatlabTM Toolboxes.

(1) The `CSIRO seawater
toolbox <http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm>`_
(SEAWATER-3.3) for calculating the properties of sea water. Uses the formulas
from Unesco's joint panel on oceanographic tables and standards, UNESCO 1981 and
UNESCO 1983 (EOS-80).

The EOS-80 library is considered now obsolete; it is provided here for
compatibility with old scripts, and to allow a smooth transition to the
new TEOS-10.

(2) The `Gibbs Sea Water <http://www.teos-10.org/software.htm>`_ (GSW v2.0).

A oceanographic toolbox of the International Thermodynamic Equation of
Seawater - 2010 or TEOS-10.

Contains the functions for evaluating the thermodynamic properties of pure water
(using IAPWS-09) and seawater (using IAPWS-08 for the saline part).

The author has no intention to do things in a "pythonic-way", it is just a
"work around" from someone that couldn't afford MatlabTM anymore. Most of the
functions are a direct translation from MatlabTM to python, with the exception
for library._delta_SA, which was completely re-written with more clarity and
efficiennce by Eric Firing.

For those coming from matlab try the (test/matlab-test.py) script, it
is replicate the testing from the original MatlabTM toolbox. However,
Bjørn Ådlandsvik has written a much better "pythonic" testing suite
(test/test_all.py).

More information:
    http://pypi.python.org/pypi/seawater/