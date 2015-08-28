ps = 2.41888432650E-5
bohr =0.529177210  
mass = 1822.6806

import numpy as np
import matplotlib.pyplot as plt
from math import pi, factorial
import scipy as sp

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
        return t/2*pi  #1.0
    else:
        return 2/pi*(np.sin((x-x0)*t/2.)*np.sin((x-x0)*t/2.))/\
                ((x-x0)*(x-x0)*t)

kpts = [2,4,6,8,10,12,16,20]

for kk in kpts:
    f = open('nacs_spectrum_{0}_{0}_1.out'.format(kk))
    n = int(f.readline().split()[-1])
    dx = float(f.readline().split()[-1])
    big_nbins = int(float(f.readline().split()[-1]))

    nbins = int(big_nbins/1)#10)



    xaxis = np.zeros(nbins)
    for i,xx in enumerate(xaxis):
        xaxis[i] = i*dx*bohr
    nspectra = (n+1)*n/2

    spectrum = np.zeros([nspectra,nbins])+0j

    for i in range(nspectra):
        f.readline()
        f.readline()
        f.readline()
        print i
        for b in range(nbins):
    #        print f.readline().split()[-1]
            spectrum[i,b] = complex(f.readline().split()[-1])
        for b in range(nbins, big_nbins):
            f.readline()

    friction_tensor = np.zeros([n,n,5])#+0j

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
                friction = 0.
                for b in range(nbins):
                    delta = square(xaxis[b],x0,s) 
                    norm += delta
                    friction += delta*spectrum[k,b]
                friction /= norm
                friction_tensor[i,j,0] = friction#/ps/mass*pi
                friction_tensor[j,i,0] = friction#/ps/mass
                norm = 0.
                friction = 0.
                for b in range(nbins):
                    delta = gauss(xaxis[b],x0,s) 
                    norm += delta
                    friction += delta*spectrum[k,b]
                friction /= norm
                friction_tensor[i,j,1] = friction#/ps/mass
                friction_tensor[j,i,1] = friction#/ps/mass
                norm = 0.
                friction = 0.
                for b in range(nbins):
                    delta = fermi(xaxis[b],x0,s) 
                    norm += delta
                    friction += delta*spectrum[k,b]
                friction /= norm
                friction_tensor[i,j,2] = friction#/ps/mass
                friction_tensor[j,i,2] = friction#/ps/mass
                norm = 0.
                friction = 0.
                for b in range(nbins):
                    delta = sine(xaxis[b],x0,s) 
                    norm += delta
                    friction += delta*spectrum[k,b]
                friction /= norm
                friction_tensor[i,j,3] = friction#/ps/mass
                friction_tensor[j,i,3] = friction#/ps/mass
                norm = 0.
                friction = 0.
                for b in range(nbins):
                    delta = lorentz(xaxis[b],x0,s) 
                    norm += delta
                    friction += delta*spectrum[k,b]
                friction /= norm
                friction_tensor[i,j,4] = friction#/ps/mass
                friction_tensor[j,i,4] = friction#/ps/mass
                k = k +1

        string = str(s) + ' '
        for f in range(5):
            e_vals,e_vecs = np.linalg.eig(friction_tensor[:,:,f])
            e_vals = 1./e_vals*1e12
            for e in e_vals:
                string += ' {:.6E} '.format(e) + ' '
        print string
        with open('integrate_friction_{0}_{0}_1.out'.format(kk), "a") as file:
            file.write(string + "\n")
