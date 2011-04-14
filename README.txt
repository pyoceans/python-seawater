===============
Python Seawater
===============

This python package contains a python re-write of two MatlabTM Toolboxes and some
extra fun functions. Typical usage often looks like this::

    >>> import seawater.gibbs as gsw
    >>> gsw.SA_from_SP(SP, p, lon, lat)

    >>> import seawater.csiro as sw
    >>> sw.test()

Main modules
============

gibbs
-----

The `Gibbs Sea Water <http://www.teos-10.org/software.htm>`_ (GSW v2.0).

A oceanographic toolbox of the International Thermodynamic Equation of
Seawater - 2010 or TEOS-10.

Contains the functions for evaluating the thermodynamic properties of pure water
(using IAPWS-09) and seawater (using IAPWS-08 for the saline part).

csiro
-----

The `CSIRO seawater toolbox <http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm>`_
(SEAWATER-3.3) for calculating the properties of sea water. The package uses the
formulas from Unesco's joint panel on oceanographic tables and standards, UNESCO
1981 and UNESCO 1983 (EOS-80) .

The EOS-80 library is considered now obsolete; it is provided here for
compatibility with old scripts, and to allow a smooth transition to the
new TEOS-10.

extra
-----

Just some extra functions and a wave calculator.


Thanks
======

* Bjørn Ådlandsvik - Testing unit and several bug fixes
* Eric Firing - Support for masked arrays, re-write of _delta_SA
* Trevor J. McDougall (and all of SCOR/IAPSO WG127) for making available the Matlab and Fortran versions of this software

Acknowledgments
---------------

* SCOR/IAPSO WG127. Most of the gibbs module is derived from the GSW Oceanographic Toolbox of TEOS-10.

* CSIRO. The csiro.py is a re-write of their seawater matlab toolbox.


More information:
    http://pypi.python.org/pypi/seawater/

Note on versioning:

The MAJOR.MINOR.MICRO will be used to represent:

MAJOR == The matlab version from the TEOS-10 Group

MINOR == Significant changes made in the python version

MICRO == Bug fixes only
