import numpy as np
from vibcooling.parser.siesta import *
import os, sys
import math
from lmfit import minimize, Parameters, Parameter, report_fit
from time import time
import matplotlib.pyplot as plt

hbar=6.58211928E-16 #eV*s
pi=3.14159265359

start=time()
#print 'time ', time() - start
path = os.getcwd()

"""
For equilibrium geometry reading: 
    -Coefficients 
    -Fermi Level 
    -Eigenvalues
    -kpoints
    -Total Energy
    -Cell parameters
"""

fermi_energy, eigenvalues = siesta_read_eigenvalues(path+"/eq/MyM")
n_spin=len(eigenvalues[0,:,0])
fermi_level= float(fermi_energy)
kpoints_weights = siesta_read_kpoints(path+"/eq/MyM")
psi = siesta_read_coefficients(path+"/eq/MyM")
total_energy_eq= float(siesta_read_total_energy(path+"/eq/input"))
print "total_energy_eq = ",total_energy_eq
cell=siesta_read_struct_out(path+"/eq/MyM")
dq=0.05
'''
"""
For displaced (forward and backward) geometry reading:
    -Hamiltonian
    -Overlap
    -Total Energy
"""
reduced_mass = 6.860525
total_energy_plus= float(siesta_read_total_energy(path+"/disp_1/input"))
total_energy_minus= float(siesta_read_total_energy(path+"/disp_-1/input"))
print 'total_energy_plus = ',total_energy_plus , 'total_energy_minus = ', total_energy_minus
#H_plus,S_plus = siesta_read_HSX(kpoints, cell, path+"/disp_1/MyM",debug=1)
#H_minus,S_minus = siesta_read_HSX(kpoints, cell, path+"/disp_-1/MyM",debug=1)
"""
Parameters for Tully Method
"""
def parabola(params,t,data):
    omega = params['omega'].value
    y_shift = params['y_shift'].value
    model = 0.5*t*t*omega*omega+y_shift
    return data - model
fitting_data= np.array([(-dq,total_energy_minus),(0.0,total_energy_eq),(dq, total_energy_plus)])
#print fitting_data[:,0]

params = Parameters()
params.add('omega',value= 100)
params.add("y_shift",value = 4000)
result= minimize(parabola,params,args=(fitting_data[:,0],fitting_data[:,1]))
conversion_factor = 3.276498E+03 #conversion from sqrt(eV/amu)/ang to cm-1
omega = params['omega'].value*np.sqrt(2)
omega_hz= params['omega'].value*98226949774380.3
print 'omega = ', omega*conversion_factor/2/pi
print omega_hz
#print 'omega by hand =', (math.sqrt(total_energy_plus/2+total_energy_minus/2-total_energy_eq)/dq)*3.276498E+03

x = np.linspace(-dq,dq,100)
plt.plot()
plt.plot(fitting_data[:,0],fitting_data[:,1],label='data points',ms="s")
plt.plot(x, 0.5*params['omega'].value*params['omega'].value*x*x+params['y_shift'].value, label='fit')
plt.legend()
plt.show()
#print 'y shift = ', params['y_shift'].value


#window_width
#step for finite differences
#constant= pi/hbar/omega/omega
"""
Calculate numerical differences
"""

delta = hbar*omega_hz
window = delta*0.2
'''
gamma=0
product=0
hit_the_window=0
energy=[]
for k in range(len(kpoints_weights)):
    for s in range(n_spin):
        for i in range(len(eigenvalues[k,s,:])):
            if eigenvalues[k,s,i] > fermi_level:
                pass
            else:
                for f in range(len(eigenvalues[k,s,:])):
                    if eigenvalues[k,s,f] < fermi_level:
                        pass
                    else:
                        energy.append(eigenvalues[k,s,f]-eigenvalues[k,s,i])
                        #excit += *kpoints_weights[k]
                        #hit_the_window+=1     
data=np.asarray(sorted(i for i in energy if i <=5))

print max(data)
print len(data)


def gaussian(x, x_mean, broadening):
    
    gaussian_val = np.exp(-0.5*(((x-x_mean) / (broadening))**2));
    return gaussian_val

def dos_binning(eigenvalues,broadening=0.01, bin_width=0.01,
        coeffs=None):
    """
    performs binning for a given set of eigenvalues and 
    optionally weight coeffs.
    """

    if coeffs is None:
        coeffs = np.ones(len(eigenvalues))

    lowest_e = float(min(eigenvalues)) - 0.05
    highest_e = float(max(eigenvalues)) + 0.05
    num_bins = int((highest_e-lowest_e)/bin_width)
    print num_bins
    x_axis = np.zeros([num_bins])
    data = np.zeros([num_bins])

    #setting up x-axis
    for i in range(num_bins):
        x_axis[i] = lowest_e + i * bin_width
    #    print x_axis[i]
    #get DOS
    for e in eigenvalues:
        for i in range(num_bins):
            data[i] += gaussian(x_axis[i],e,broadening)#*coeffs[ei]
            #print data[i]
#    data /= len(eigenvalues)

    return x_axis, data
x,y = dos_binning(data)


plt.plot(x,y)
plt.legend()
plt.show()
np.savetxt('test.txt', np.vstack([x,y]).T,fmt=['%0.8f','%0.8f']) 
