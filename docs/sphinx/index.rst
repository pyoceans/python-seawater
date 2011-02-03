.. seawater documentation master file, created by
   sphinx-quickstart on Tue Aug 10 16:47:25 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to seawater's documentation!
====================================

Python Seawater
===============

Introduction:
-------------

This python package contains a python translation for two MatlabTM Toolboxes.

(1) The `CSIRO seawater toolbox <http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm>`_   (SEAWATER-3.3) for calculating the properties of sea water. Uses the formulas from Unesco's joint panel on oceanographic tables and standards, UNESCO 1981 and UNESCO 1983 (EOS-80).

The EOS-80 library is considered now obsolete; it is provided here for compatibility with old scripts, and to allow a smooth transition to the new TEOS-10.

(2) The `Gibbs Sea Water <http://www.teos-10.org/software.htm>`_ (GSW v2.0).

A oceanographic toolbox of the International Thermodynamic Equation of Seawater - 2010 or TEOS-10.

Contains the functions for evaluating the thermodynamic properties of pure water (using IAPWS-09) and seawater (using IAPWS-08 for the saline part).

The author has no intention to do things in a "pythonic-way", it is just a "work around" from someone that couldn't afford MatlabTM anymore.

GSW toolbox was rewritten in OO approach, i.e., there is a Gibbs class that contains all (SA,t,p) functions as methods.

Ex.:
        >>> from seawater.gibbs import Gibbs
        >>> SA = [34.5075, 34.7165, 34.8083, 34.8465, 34.8636, 34.8707, 34.8702]
        >>> t  = [27.9620, 4.4726, 2.1178, 1.6031, 1.4601, 1.4753, 1.5998]
        >>> p  = [0., 1010., 2025., 3045., 4069., 5098., 6131.]
        >>> STP = Gibbs(SA, t, p)
        >>> STP.beta_const_pt()

For more information:
    http://pypi.python.org/pypi/seawater/


Modules:
--------

.. toctree::
   :maxdepth: 4


gibbs Seawater Documentation
============================

This page contains the Csiro Module documentation.

The :mod:`gibbs` module
-----------------------

.. automodule:: seawater.gibbs
    :members:
    :undoc-members:
    :show-inheritance:


csiro Seawater Documentation
============================

This page contains the Csiro Module documentation.

The :mod:`csiro` module
-----------------------

.. automodule:: seawater.csiro
    :members:
    :undoc-members:
    :show-inheritance:


extras Documentation
====================

This page contains the Extras Package documentation.

Waves Documentation
===================

This page contains the Waves Package documentation.

The :mod:`waves` Package
------------------------

.. automodule:: seawater.extras.waves
    :members:
    :undoc-members:
    :show-inheritance:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
