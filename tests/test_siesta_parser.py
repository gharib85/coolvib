from vibcooling.parser.siesta import *

filename = 'siesta_parser_test/MyM'

cell = siesta_read_struct_out(filename)
print cell

fermi_level, eigenvalues = siesta_read_eigenvalues(filename)

print 'fermi_level ', fermi_level, 
print 'eigenvalue shape ',eigenvalues.shape

kpt_array = siesta_read_kpoints(filename)

print 'kpts'
print kpt_array.shape
#print kpt_array

psi = siesta_read_coefficients(filename)

print 'psi'
print psi.shape

from time import time
start = time()
#kpt_array = kpt_array[:5]
#print kpt_array
H,S = siesta_read_HSX(kpt_array, cell, filename,debug=1)
print 'time ', time() - start

print 'H'
print H.shape
print 'S'
print S.shape

nspin = 0
nk = 0 
i = 0
psi1 = psi[nk,nspin,i,:]
H1 = H[nk,nspin,:,:]
S1 = S[nk,:,:]
S1_inv = np.linalg.inv(S1)
e = eigenvalues[nk,nspin,i]
print e

print np.dot(psi1.conjugate().transpose(),np.dot(H1,psi1))
print e*np.dot(psi1.conjugate().transpose(),np.dot(S1,psi1))

