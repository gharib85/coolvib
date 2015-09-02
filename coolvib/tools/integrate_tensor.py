#!/usr/bin/python
"""
integrate_tensor.py takes a given spectral function file 
calculated for all cartesian directions and calculates the 
friction tensor and analyzes the window convergence

python integrate_tensor.py <filename>

Arguments:

    <filename>: str
        Name of the file containing the spectral functions

Code Example::

    python integrate_tensor.py spectral_function.out

"""

from sys import argv

import numpy as np
from math import factorial

from coolvib.routines import *
if __name__=="__main__":


    f = open(str(argv[1]))#'nacs-spectrum.out')

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
                    delta = gaussian(xaxis[b],x0,s) 
                    norm += delta
                    friction += delta*xaxis[b]*spectrum[k,b]
                friction /= norm
                friction_tensor[i,j,1] = friction
                friction_tensor[j,i,1] = np.conjugate(friction)
                norm = 0.
                friction = 0.+1.0j*0.
                for b in range(nbins):
                    delta = squashed_fermi(xaxis[b],x0,s) 
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
                    delta = lorentzian(xaxis[b],x0,s) 
                    norm += delta
                    friction += delta*xaxis[b]*spectrum[k,b]
                friction /= norm
                friction_tensor[i,j,4] = friction
                friction_tensor[j,i,4] = np.conjugate(friction)
                k = k + 1

        #throw away imaginary parts
        friction_tensor = np.array(friction_tensor,dtype=np.float)
        #string = str(s) + ' '
        #for i in friction_tensor[0,0]:
            #string += str(1./i) + ' '
        #print string 
        #continue

        string = str(s) + ' '
        for e in range(n):
            for f in range(5):
                E,V = np.linalg.eig(friction_tensor[:,:,f])
                E = 1./E
                sort = np.argsort(E)
                E = np.array(E[sort],dtype=np.float)
                string += str(E[e]) + ' '
        print string

    
else:
    pass
