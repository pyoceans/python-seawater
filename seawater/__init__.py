# -*- coding: utf-8 -*-

"""
Python Seawater
===============

This module contains a translation of the original CSIRO Matlab package (SEAWATER-3.2) for calculating the properties of sea water.
It consists of a self contained library easy to use. The only requirent is NumPy.

The author has no intention to do things in a "pythonic-way", it is just a "work around" from someone that couldn't afford Matlab anymore.

The package uses the formulas from Unesco's joint panel on oceanographic tables and standards, UNESCO 1981 and UNESCO 1983.

The present version is 1.0.5, released 23 August 2010.
"""

__authors__    = ['Filipe Fernandes']
__copyright__  = ["CSIRO"]
__credits__    = ["Filipe Fernandes", "Lindsay Pender","Phil Morgan"]
__license__    = ["CSIRO"]
__version__    = ["1.0.5"]
__maintainer__ = ["Filipe Fernandes"]
__email__      = ["ocefpaf@gmail.com"]
__status__     = ["Production"]
__all__        = ["seawater","csiro","extras"]
__obs__        = ["Translated from matlab CSIRO seawater toolbox Version 3.2"]
__web__        = ['http://ocefpaf.tiddlyspot.com/']
__created__    = ["14-Jan-2010"]
__modified__   = ["17-Aug-2010"]

import csiro
import extras

#TODO:
#if __name__ == '__main__':
    #import doctest
    #doctest.testmod()

#TODO: make all numpy arrays automatic like pylab
