Python Seawater
===============

The CSIRO seawater toolbox ([SEAWATER-3.3](http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm))
for calculating the properties of sea water.  The package uses the formulas
from Unesco's joint panel on oceanographic tables and standards, UNESCO 1981
and UNESCO 1983 (EOS-80).

The EOS-80 library is considered now obsolete;  it is provided here for
compatibility with old scripts, and to allow a smooth transition to the new
TEOS-10.

Notes
-----
The Python version default output unit for sw.dist is 'km' instead of  'nm'.

Another difference is that the Python version takes pressure as the first
dimension, i.e. M pressure by N positions.  The MatlabTM version does some
guessing at this that we simply ignore to confusions.  That is because matlab
have "two" types of 1D arrays (row or column), while python has just "1D"
arrays.  To form a "column" array one must actually create a 2D array M by 1.

|    P      |     S      |    T       |
|:---------:|:----------:|:----------:|
|    10     |   34.5487  |   28.7856  |
|    50     |   34.7275  |   28.4329  |
|   125     |   34.8605  |   22.8103  |
|   250     |   34.6810  |   10.2600  |
|   600     |   34.5680  |    6.8863  |
|  1000     |   34.5600  |    4.4036  |
|     .     |         .  |         .  |
|     .     |         .  |         .  |
|     .     |         .  |         .  |

Check out the `test_octave.py` routine to test the Python library against an
available MatlabTM library inside Python via the oct2py package.  The current
version was tested against seawater v3.3.

More information:
    http://pythonhosted.org/seawater
