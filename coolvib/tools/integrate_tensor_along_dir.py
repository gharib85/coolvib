#!/usr/bin/python
"""
integrate_tensor_along_dir.py takes a given spectral function file 
calculated for all cartesian directions and calculates the 
friction tensor and uses it to determine the lifetime along a projected 
direction

python integrate_tensor_along_dir.py <filename> <mode v1x v1y v1z v2x v2y v2z ......>

Arguments:

    <filename>: str
        Name of the file containing the spectral functions

    <mode>: floats
        displacement vector of the mode along which the lifetime shall 
        be calculated


Code Example::

    python integrate_tensor_along_dir.py spectral_function.out 0. 0. 1.0 0. 0.0 -1.0

"""

if __name__=="__main__":
    pass
else:
    exit()

from sys import argv, exit
import numpy as np
from math import factorial
from coolvib.routines import *

f = open(str(argv[1]))#'nacs-spectrum.out')

mode = []
for i in argv[2:]:
    mode.append(float(i))

mode= np.array(mode)
mode /= np.linalg.norm(mode)

n = int(f.readline().split()[-1])
dx = float(f.readline().split()[-1])
nbins = int(f.readline().split()[-1])

xaxis = np.zeros(nbins)
for i,xx in enumerate(xaxis):
    xaxis[i] = i*dx
nspectra = (n+1)*n/2

spectrum = np.zeros([nspectra,nbins],dtype=np.complex)

for i in range(nspectra):
    f.readline()
    f.readline()
    f.readline()
    for b in range(nbins):
        tmp = f.readline().split()
        spectrum[i,b] = float(tmp[-2])+1.0j*float(tmp[-1])

friction_tensor = np.zeros([n,n,5],dtype=np.complex)

x0 = 0.0

#window
ss = np.linspace(0.05,1.0,20)

names = ['square','gauss','fermi','sine','lorentz']
print 'window   square     gauss     fermi     sine     lorentz'
for s in ss:
    k = 0
    for i in range(n):
        for j in range(i,n):
            norm = 0.
            friction = 0.+1.0j*0.
            for b in range(nbins):
                delta = square(xaxis[b],x0,s) 
                norm += delta
                friction += delta*xaxis[b]*spectrum[k,b]
            friction /= norm
            friction_tensor[i,j,0] = friction
            friction_tensor[j,i,0] = np.conjugate(friction)
            norm = 0.
            friction = 0.+1.0j*0.
            for b in range(nbins):
                delta = gauss(xaxis[b],x0,s) 
                norm += delta
                friction += delta*xaxis[b]*spectrum[k,b]
            friction /= norm
            friction_tensor[i,j,1] = friction
            friction_tensor[j,i,1] = np.conjugate(friction)
            norm = 0.
            friction = 0.+1.0j*0.
            for b in range(nbins):
                delta = fermi(xaxis[b],x0,s) 
                norm += delta
                friction += delta*xaxis[b]*spectrum[k,b]
            friction /= norm
            friction_tensor[i,j,2] = friction
            friction_tensor[j,i,2] = np.conjugate(friction)
            norm = 0.
            friction = 0.+1.0j*0.
            for b in range(nbins):
                delta = sine(xaxis[b],x0,s) 
                norm += delta
                friction += delta*xaxis[b]*spectrum[k,b]
            friction /= norm
            friction_tensor[i,j,3] = friction
            friction_tensor[j,i,3] = np.conjugate(friction)
            norm = 0.
            friction = 0.+1.0j*0.
            for b in range(nbins):
                delta = lorentz(xaxis[b],x0,s) 
                norm += delta
                friction += delta*xaxis[b]*spectrum[k,b]
            friction /= norm
            friction_tensor[i,j,4] = friction
            friction_tensor[j,i,4] = np.conjugate(friction)
            k = k + 1

    #throw away imaginary parts
    friction_tensor = np.array(friction_tensor,dtype=np.float)
    for f in range(5):
         string += str(float( 1./np.dot(e.transpose(),np.dot(friction_tensor[:,:,f],e)))) +' '
    print string
    
