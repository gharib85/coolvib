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


def gaussian(x, x_mean, broadening):
    
    gaussian_val = np.exp(-0.5*(((x-x_mean) / (broadening))**2));
    return gaussian_val

bin_width=0.01
broadening=0.05

path_1=path
for index in range(2,8):
    path= path_1+"/{0}".format(index)
    print path
    fermi_energy, eigenvalues = siesta_read_eigenvalues(path+"/eq/MyM")
    lowest_e = float(min(eigenvalues[0,0,:])) - 0.05
    highest_e = float(max(eigenvalues[0,0,:])) + 0.05
    num_bins = int((highest_e-lowest_e)/bin_width)
    x_axis = np.zeros([num_bins])
    data_excit_all = np.zeros([num_bins])
    data_dos_all = np.zeros([num_bins])
    for i in range(num_bins):
        x_axis[i] = lowest_e + i * bin_width
    n_spin=len(eigenvalues[0,:,0])
    fermi_level= float(fermi_energy)
    kpoints_weights = siesta_read_kpoints(path+"/eq/MyM")
    psi = siesta_read_coefficients(path+"/eq/MyM")
    total_energy_eq= float(siesta_read_total_energy(path+"/eq/input"))
    print "total_energy_eq = ",total_energy_eq
    cell=siesta_read_struct_out(path+"/eq/MyM")
    dq=0.05
    gamma=0
    product=0
    hit_the_window=0
    dos=[]
    excit=[]
    for k in range(len(kpoints_weights)):
        for s in range(n_spin):
            i_homo = 0
            i_lumo = 0
            for i,e1 in enumerate(eigenvalues[k,s,:]):
                if e1< fermi_level:
                    i_homo = i
                elif e1>fermi_level:
                    i_lumo=i
                    break
#                i_lumo+i_homo+1
            for i,e1 in enumerate(eigenvalues[k,s,:]):
                e_dos=e1-fermi_level
#                print "e_dos", e_dos
                bin_min=min(range(len(x_axis)), key=lambda i: abs(x_axis[i]-(e_dos-0.2)))
                bin_max=min(range(len(x_axis)), key=lambda i: abs(x_axis[i]-(e_dos+0.2)))
                for bin in range(bin_min,bin_max):
                    data_dos_all[bin] += gaussian(x_axis[bin],e_dos,broadening)*kpoints_weights[k,3]
                if e1> fermi_level:
                    pass
                else:
                    for f,e2 in enumerate(eigenvalues[k,s,i_lumo:]):
                        e = e2-e1
                        if e <= 1.0:
                            print "e", e
                            bin_min=min(range(len(x_axis)), key=lambda i: abs(x_axis[i]-(e-0.2)))
                            bin_max=min(range(len(x_axis)), key=lambda i: abs(x_axis[i]-(e+0.2)))
                            for bin in range(bin_min,bin_max):
                                data_excit_all[bin] += gaussian(x_axis[bin],e,broadening/2)*kpoints_weights[k,3]
#                                print data_excit_all[bin]
    np.savetxt(path_1+'/excit_scan_layers_{0}.txt'.format(index), np.vstack([x_axis,data_excit_all]).T,fmt=['%0.8f','%0.8f']) 
    np.savetxt(path_1+'/dos_scan_layers_{0}.txt'.format(index), np.vstack([x_axis,data_dos_all]).T,fmt=['%0.8f','%0.8f']) 



