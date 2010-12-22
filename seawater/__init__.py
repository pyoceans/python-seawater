# -*- coding: utf-8 -*-

"""
Python Seawater
===============

This python package contains:
1) Python translation of the original CSIRO Matlab package (SEAWATER-3.2) for calculating the properties of sea water. The package uses the formulas from Unesco's joint panel on oceanographic tables and standards, UNESCO 1981 and UNESCO 1983.

2) Python translation of the original Gibbs-SeaWater Matlab package (GSW v2.0). The Gibbs-SeaWater (GSW) Oceanographic Toolbox contains the TEOS-10 subroutines for evaluating the thermodynamic properties of pure water (using IAPWS-09) and seawater (using IAPWS-08 for the saline part).

The author has no intention to do things in a "pythonic-way", it is just a "work around" from someone that couldn't afford Matlab anymore.

CAVEAT: These modules do not adhere to strict basic-SI units but rather oceanographic units are adopted.

The present version is 1.2.1, released 24 December 2010.
"""

__authors__    = ['Filipe Fernandes']
__copyright__  = ["CSIRO"]
__credits__    = ["Filipe Fernandes", "Lindsay Pender","Phil Morgan"]
__license__    = ["CSIRO"]
__version__    = ["1.2.1"]
__maintainer__ = ["Filipe Fernandes"]
__email__      = ["ocefpaf@gmail.com"]
__status__     = ["Production"]
__all__        = ["seawater","gibbs","csiro","extras"] #FIXME
__obs__        = ["csiro: Translated from matlab CSIRO seawater toolbox version 3.2"]
__obs__        = ["gibbs: Translated from matlab GSW seawater toolbox version 2.0"]
__web__        = ['http://ocefpaf.tiddlyspot.com/']
__created__    = ["14-Jan-2010"]
__modified__   = ["24-Dec-2010"]
__all__        = ["csiro"]