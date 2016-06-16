"""
test_aims_parsing.py

tests the different routines in parser/aims.py
"""

import numpy as np
from coolvib.parser.aims import *
from ase.all import *

atoms = read('aims_parser_test/geometry.in')
cell = atoms.cell

filename = 'aims_parser_test/OUTPUT'

fermi_level, kpoint_weights = aims_read_fermi_and_kpoints(filename, cell)

print fermi_level

print kpoint_weights

nkpts = len(kpoint_weights)

eigenvalues, psi, occ, orb_pos = aims_read_eigenvalues_and_coefficients(fermi_level, './aims_parser_test/', spin=True)

print eigenvalues
print eigenvalues.shape
print psi.shape

H, S = aims_read_HS('./aims_parser_test/',spin=True)

print H.shape

print S.shape

nspin = 0
nk = 0

psi1 = psi[nk,nspin,:,:]
S1 = S[nk,:,:]
H1 = H[nk,nspin,:,:]
S_inv = np.linalg.inv(S1)
S_invH = np.dot(S_inv,H1)

import scipy.linalg as la
E, V = la.eigh(H1,S1)
print E*27.2114

for i in range(len(psi1)):
    print np.dot(psi1[i].conjugate(),np.dot(S_invH,psi1[i]))/np.dot(psi1[i].conjugate(),psi1[i])*27.2114

