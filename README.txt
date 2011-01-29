===============
Python Seawater
===============

Introduction:
-------------

This python package contains a python translation for two Matlab user contributed toolboxes. (1) the `sewater <http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm>`_ (EOS-80) and (2) the `gibbs <http://www.teos-10.org/software.htm>`_ seawater (TEOS-10).


(1) Based on the original CSIRO Matlab package (SEAWATER-3.3) for calculating the properties of sea water. The package uses the formulas from Unesco's joint panel on oceanographic tables and standards, UNESCO 1981 and UNESCO 1983.

The SeaWater library of EOS-80 seawater properties is considered now obsolete; it is provided here for compatibility with old scripts and to allow a smooth transition to the new TEOS-10.

(2) Based on the Gibbs SeaWater (GSW v2.0) Oceanographic Toolbox of the International Thermodynamic Equation Of Seawater - 2010, (TEOS-10); contains the TEOS-10 subroutines for evaluating the thermodynamic properties of pure water (using IAPWS-09) and seawater (using IAPWS-08 for the saline part).

The author has no intention to do things in a "pythonic-way", it is just a "work around" from someone that couldn't afford Matlab anymore.

CAVEAT: These modules do not adhere to strict basic-SI units but rather oceanographic units are adopted.

GSW toolbox was rewritten in OO approach, i.e., there is a SaTePr class that contains all (SA,t,p) functions as methods.

Ex.:
        >>> from seawater.gibbs import SaTePr
        >>> SA = [34.5075, 34.7165, 34.8083, 34.8465, 34.8636, 34.8707, 34.8702]
        >>> t  = [27.9620, 4.4726, 2.1178, 1.6031, 1.4601, 1.4753, 1.5998]
        >>> p  = [0., 1010., 2025., 3045., 4069., 5098., 6131.]
        >>> STP = SaTePr(SA, t, p)
        >>> STP.beta_const_pt()

For more information click `here. <http://packages.python.org/seawater/>`_