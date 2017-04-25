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
Example 3: Setting up finite difference calculation for CO on Pd(111)
=====================================================================

This example uses tools in coolvib/tools/ to perform finite difference displacements 
with the help of ASE for the tensor and normal mode case.

Code
-----

.. literalinclude:: ../examples/example3/example3.py



Detailed Explanation
--------------------

In the above code we first initialize an ASE atoms object and an ASE 
FHI-aims calculator object. The output string is essential to guarantee that 
all input data is being dumped to file for each displacement. The tool 
coolvib/tools/fhiaims_add_outputoptions_to_control.py helps to identify this 
complicated string.

After setting the calculator we initialize the finite_difference object and run it.
The result is a Hessian calculation that also dumps the necessary friction input files 
into folders.

mode_displacement does the same for a given normal mode


"""
