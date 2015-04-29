import numpy as np
from vibcooling.parser.siesta import *
import os, sys
import math
from lmfit import minimize, Parameters, Parameter, report_fit
from time import time
import matplotlib.pyplot as plt
import pylab
from tully_function import *

hbar=6.58211928E-16 #eV*s
pi=3.14159265359
ryd = 13.605698066 #eV
bohr = 0.529177249 #Ang

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

def parabola(params,t,data):
    omega = params['omega'].value
    y_shift = params['y_shift'].value
    model = 0.5*t*t*omega*omega+y_shift
    return data - model

kmin=3
kmax=4
window_number = 10
path = os.getcwd()
path_1=path
out_info=np.zeros([kmax-kmin,window_number,8])
for index in range(kmin,kmax):
    start=time()
    path= path_1+"/{0}".format(index)
    print path
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
    fermi_level= float(fermi_energy)
    n_spin=len(eigenvalues[0,:,0])
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
    total_energy_plus= float(siesta_read_total_energy(path+"/disp_1/input"))
    total_energy_minus= float(siesta_read_total_energy(path+"/disp_-1/input"))
    print 'total_energy_plus = ',(total_energy_plus-total_energy_eq) , 'total_energy_minus = ', (total_energy_minus-total_energy_eq)
    H_plus,S_plus = siesta_read_HSX(kpoints_weights, cell, path+"/disp_1/MyM",debug=1)
    H_minus,S_minus = siesta_read_HSX(kpoints_weights, cell, path+"/disp_-1/MyM",debug=1)
    end=time()
    print "H and S read in", "{:8.1f}".format(end-start), "  sec"
    """
    Parameters for Tully Method
    """
    fitting_data= np.array([(-dq,total_energy_minus),(0.0,total_energy_eq),(dq, total_energy_plus)])
    params = Parameters()
    params.add('omega',value= 100)
    params.add("y_shift",value = 4000)
    result= minimize(parabola,params,args=(fitting_data[:,0],fitting_data[:,1]))
    conversion_factor = 3.276498E+03 #conversion from sqrt(eV/amu)/ang to cm-1
    omega = params['omega'].value#*np.sqrt(2)
    omega_hz= params['omega'].value*98226949774380.3 # in s-1 (Hz)
    """
    Calculate numerical differences
    """
    window_interp=([])
    gamma_interp=([])
    for window_width in range(1,window_number): #vary window size from 10 to 50% of delta
        window = hbar*omega_hz*window_width/10
        print "         window", window
        d={     "kpoints_weights":kpoints_weights,\
                "n_spin":n_spin,\
                "eigenvalues":eigenvalues,\
                "fermi_level":fermi_level,\
                "omega":omega,\
                "window":window,\
                "psi":psi,\
                "H_plus":H_plus,\
                "H_minus":H_minus,\
                "S_plus":S_plus,\
                "S_minus":S_minus,\
                "dq":dq,\
                }
        gamma, hit_the_window=tully_routine(**d)
        end1=time()
        print "         Lifetime calculated in", "{:8.1f}".format(end1-end), '   sec'
        if hit_the_window ==0:
            lifetime = 0
        else:
            window_interp.append(window)
            gamma_interp.append(gamma)
            lifetime = 1/gamma*1e12
            print "         lifetime        ", lifetime
        out_info[index-kmin,window_width-1,:]= index, (total_energy_plus-total_energy_eq),(total_energy_minus-total_energy_eq), \
                            hbar*omega_hz, window, hit_the_window, end1-start, lifetime
        np.savetxt(path_1+'/tullymethod.out', np.reshape(out_info,((kmax-kmin)*window_number,8)),fmt='%3.4f')
    params = Parameters()
#    print window_interp
#    print gamma_interp
    params.add('omega',value= 100000000)
    params.add("y_shift",value = min(np.array(gamma_interp)), max=min(np.array(gamma_interp)))
    result= minimize(parabola,params,args=(np.array(window_interp),np.array(gamma_interp)))
    gamma_interpolated = params['y_shift'].value
    lifetime_interpolated = 1/gamma_interpolated*1e12
    print "         exrapolated lifetime is          ", lifetime_interpolated, "ps"
    plt.plot(np.array(window_interp),np.array(gamma_interp),marker='o', label='{0} {0} 1'.format(index,index))
    x=np.linspace(min(np.array(window_interp)/1.1), max(np.array(window_interp)*1.1), 20)
    y=np.zeros([len(x)])
    for x_index,x_value in enumerate(x):
        y[x_index]=params['y_shift'].value + np.square(params['omega'].value)*np.square(x_value)/2
    plt.plot(x,y, label = "fit")
    plt.legend()
#    plt.ion()
#    plt.draw()
    plt.show()







