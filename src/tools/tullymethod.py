import numpy as np
from vibcooling.parser.siesta import *
import os, sys
import math
from lmfit import minimize, Parameters, Parameter, report_fit
from time import time
import matplotlib.pyplot as plt

hbar=6.58211928E-16 #eV*s
pi=3.14159265359
"""
#finding the harmonic frequency from the curvature of the PES 
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
"""
"""
Hz          cm-1    
1.000000    3.335641E-11
eV          J   
1.000000    1.602177E-19  
Na      
6.022141E+23    
Hz          sqrt(ev)/Angst  
1.000000    1.018051E-14    
sqrt(eV)/Angst  Hz                  cm-1
1.000000        98226949774380.300000   3.276498E+03
"""

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

"""
For displaced (forward and backward) geometry reading:
    -Hamiltonian
    -Overlap
    -Total Energy
"""
dq=0.05
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
'''
x = np.linspace(-dq,dq,100)
plt.plot()
plt.plot(fitting_data[:,0],fitting_data[:,1],label='data points',ms="s")
plt.plot(x, 0.5*params['omega'].value*params['omega'].value*x*x+params['y_shift'].value, label='fit')
plt.legend()
plt.show()
#print 'y shift = ', params['y_shift'].value
'''
#window_width
#step for finite differences
#constant= pi/hbar/omega/omega
"""
Calculate numerical differences
"""

delta = hbar*omega_hz
window = delta*0.2
gamma=0
product=0
hit_the_window=0
for k in range(len(kpoints_weights)):
    for s in range(n_spin):
        for i in range(len(eigenvalues[k,s,:])):
            if eigenvalues[k,s,i] > fermi_level:
                pass
            else:
                for f in range(len(eigenvalues[k,s,:])):
                    if eigenvalues[k,s,f] < fermi_level:
                        pass
                    elif (delta-window) <=  (eigenvalues[k,s,f]-eigenvalues[k,s,i]) <= (delta+window):
                        H_q=(H_plus[k,s,:,:]-H_minus[k,s,:,:])/2/dq
                        S_q=(S_plus[k,:,:]-S_minus[k,:,:])/2/dq
                        product= np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(H_q-fermi_level*S_q,psi[k,s,f,:]))
                        gamma += np.square(np.absolute(product))*kpoints_weights[k]
                        hit_the_window+=1     
gamma = gamma*pi/hbar/omega_hz/omega_hz*1E+12
print gamma

