#    This file is part of coolvib
#
#        coolvib is free software: you can redistribute it and/or modify
#        it under the terms of the GNU General Public License as published by
#        the Free Software Foundation, either version 3 of the License, or
#        (at your option) any later version.
#
#        coolvib is distributed in the hope that it will be useful,
#        but WITHOUT ANY WARRANTY; without even the implied warranty of
#        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#        GNU General Public License for more details.
#
#        You should have received a copy of the GNU General Public License
#        along with coolvib.  If not, see <http://www.gnu.org/licenses/>.
"""
Example 2: Calculating the lifetime along a given mode: Co on Cu(100)
=====================================================================

**ATTENTION** This example/functionality is broken and needs to be fixed

In this example we calculate a single lifetime and relaxation rate along 
a given vibrational normal mode for CO on Cu(100), namely the internal stretch mode.
Example 3 shows how to calculate finite difference displacements along a normal mode.

In order to run this example you need to either run example 3 with FHI-aims 
or download the dataset from following URL and extract it in the example folder


Code
-----

.. literalinclude:: ../examples/example2/example2.py


Detailed Explanation
--------------------

All settings are equivalent to the tensor case in example 1. The only 
difference is that now there is only one spectral function and routines associated 
with calculating friction are called friction (rather than friction_tensor).




"""
