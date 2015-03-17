from vibcooling.parser.siesta import *

filename = 'siesta_parser_test/MyM'

fermi_level, eigenvalues = siesta_read_eigenvalues(filename)

print 'fermi_level ', fermi_level, 
print 'eigenvalue shape ',eigenvalues.shape

kpt_array = siesta_read_kpoints(filename)

print 'kpts'
print kpt_array.shape
print kpt_array

psi = siesta_read_coefficients(filename)

print 'psi'
print psi.shape

from time import time
start = time()
#kpt_array = kpt_array[0:4]
#print kpt_array
H,S = siesta_read_HSX(kpt_array, filename,debug=1)
print time() - start

print 'H'
print H.shape
print 'S'
print S.shape
