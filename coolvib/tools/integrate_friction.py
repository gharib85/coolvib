ps = 1. #2.41888432650E-5
bohr =1. #0.529177210  
mass = 1. #1822.6806

from sys import argv

import numpy as np
#import matplotlib.pyplot as plt
from math import pi, factorial
#import scipy as sp

pi = 1.

def square(x,x0,s):
    if np.abs(x-x0)>s:
        return 0.0
    else:
        return 1./s

def gauss(x,x0,s):
    return (1./np.sqrt(2.))*\
            np.exp(-0.5*(x-x0)*(x-x0)/(s*s))/(s*np.sqrt(np.pi))

def fermi(x,x0,s):
    y = (x/s)*np.sqrt(2./np.pi)+np.sqrt(0.5)
    return 2.*(y/s)*np.sqrt(2./np.pi)*np.exp(0.5-(y*y))

def lorentz(x,x0,s):
    return (1.0/np.pi)*((0.5*s)/ \
        ((x-x0)*(x-x0)+(0.5*s)*(0.5*s)))

def sine(x,x0,s):
    t = 2./s
    if abs(x-x0)<=0.000001:
        return 1.0
    else:
        return (np.sin((x-x0)*t/2.)*np.sin((x-x0)*t/2.))/\
                ((x-x0)*(x-x0)*t*t/4.0)


f = open(str(argv[1]))#'nacs-spectrum.out')

n = int(f.readline().split()[-1])
dx = float(f.readline().split()[-1])
nbins = int(f.readline().split()[-1])

xaxis = np.zeros(nbins)
for i,xx in enumerate(xaxis):
    xaxis[i] = i*dx*bohr
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

ss = np.linspace(0.05,1.0,20)
#ss /= 27.2114
#window
print n
print dx
print nbins
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
                #friction += delta*spectrum[k,b]
                friction += delta*xaxis[b]*spectrum[k,b]
            friction /= norm
            friction_tensor[i,j,0] = pi*friction/ps/mass
            friction_tensor[j,i,0] = np.conjugate(pi*friction/ps/mass)
            norm = 0.
            friction = 0.+1.0j*0.
            for b in range(nbins):
                delta = gauss(xaxis[b],x0,s) 
                norm += delta
                friction += delta*xaxis[b]*spectrum[k,b]
                #friction += delta*spectrum[k,b]
            friction /= norm
            friction_tensor[i,j,1] = pi*friction/ps/mass
            friction_tensor[j,i,1] = np.conjugate(pi*friction/ps/mass)
            norm = 0.
            friction = 0.+1.0j*0.
            for b in range(nbins):
                delta = fermi(xaxis[b],x0,s) 
                norm += delta
                friction += delta*xaxis[b]*spectrum[k,b]
                #friction += delta*spectrum[k,b]
            friction /= norm
            friction_tensor[i,j,2] = pi*friction/ps/mass
            friction_tensor[j,i,2] = np.conjugate(pi*friction/ps/mass)
            norm = 0.
            friction = 0.+1.0j*0.
            for b in range(nbins):
                delta = sine(xaxis[b],x0,s) 
                norm += delta
                friction += delta*xaxis[b]*spectrum[k,b]
                #friction += delta*spectrum[k,b]
            friction /= norm
            friction_tensor[i,j,3] = pi*friction/ps/mass
            friction_tensor[j,i,3] = np.conjugate(pi*friction/ps/mass)
            norm = 0.
            friction = 0.+1.0j*0.
            for b in range(nbins):
                delta = lorentz(xaxis[b],x0,s) 
                norm += delta
                friction += delta*xaxis[b]*spectrum[k,b]
                #friction += delta*spectrum[k,b]
            friction /= norm
            friction_tensor[i,j,4] = pi*friction/ps/mass
            friction_tensor[j,i,4] = np.conjugate(pi*friction/ps/mass)
            k = k + 1

    friction_tensor = np.array(friction_tensor,dtype=np.float)
    #string = str(s) + ' '
    #for i in friction_tensor[0,0]:
        #string += str(1./i) + ' '
    #print string 
    #continue

    #string = str(s*27.2114) + ' '
    string = str(s) + ' '
    for e in range(n):
    	for f in range(5):
	    #friction_tensor = np.array(friction_tensor,dtype=np.float)
            E,V = np.linalg.eig(friction_tensor[:,:,f])
	    E = 1./E
            sort = np.argsort(E)
	    E = np.array(E[sort],dtype=np.float)
            #V = V[sort]
	    string += str(E[e]) + ' '
            #print V[e]
    print string

    
