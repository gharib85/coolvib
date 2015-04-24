import numpy as np

a=np.loadtxt('first_kpoint')
coeff=a[:,0] +1j*a[:,1]
print np.dot(coeff, coeff.conj())
#print a

