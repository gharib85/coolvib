"""
test_aims_parsing.py

tests the different routines in parser/aims.py
"""

import numpy as np
from vibcooling.parser.aims import *

filename = 'OUTPUT'

fermi_level, kpoint_weights = aims_read_fermi_and_kpoints(filename)

print fermi_level

print kpoint_weights

nkpts = len(kpoint_weights)

eigenvalues, psi, occ = aims_read_eigenvalues_and_coefficients('./',spin=True)


print eigenvalues

print psi

print eigenvalues.shape
print psi.shape

H, S = aims_read_HS('./',spin=True)

print H.shape

print S.shape

print H[0,0,0,0]
print S[0,0,0]

