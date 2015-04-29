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

bin_width=0.005
broadening=0.05
truncatevalue=1.0
total_energy=[[]]

path_1=path
for index in range(1,31):
    path= path_1+"/{0}".format(index)
    print path
    fermi_energy, eigenvalues = siesta_read_eigenvalues(path+"/eq/MyM")
    fermi_level= float(fermi_energy)
    n_spin=len(eigenvalues[0,:,0])
    kpoints_weights = siesta_read_kpoints(path+"/eq/MyM")
#    psi = siesta_read_coefficients(path+"/eq/MyM")
    total_energy_eq= float(siesta_read_total_energy(path+"/eq/input"))
    print "total_energy_eq = ",total_energy_eq
    total_energy.append(total_energy_eq)
    cell=siesta_read_struct_out(path+"/eq/MyM")
    dq=0.05
    gamma=0
    product=0
    hit_the_window=0
    dos=[]
    excit=[]
    orb_min=min(range(len(eigenvalues[0,0,:])), key=lambda j: abs(eigenvalues[0,0,j]-(fermi_level-1)))
    orb_max=min(range(len(eigenvalues[0,0,:])), key=lambda j: abs(eigenvalues[0,0,j]-(fermi_level+1)))
    eigs=eigenvalues[0,0,orb_min:orb_max]-fermi_level
    lowest_e = float(min(eigs)) - 0.1
    highest_e = float(max(eigs)) + 0.1
    num_bins = int((highest_e-lowest_e)/bin_width)
    x_axis = np.zeros([num_bins])
    num_bins = int((highest_e-lowest_e)/bin_width)
    data_excit_all = np.zeros([num_bins])
    data_dos_all = np.zeros([num_bins])
#    print eigs
#   print num_bins
    for i in range(num_bins):
        x_axis[i] = lowest_e + i * bin_width

    for k in range(len(kpoints_weights)):
        for s in range(n_spin):
            eigs=eigenvalues[k,s,orb_min:orb_max]-fermi_level
            i_homo = 0
            i_lumo = 0
            for i,e1 in enumerate(eigs):
                e_dos=e1
                for bin in range(num_bins):
                    data_dos_all[bin] += gaussian(x_axis[bin],e_dos,broadening)*kpoints_weights[k,3]
                if e1< 0:
                    i_homo = i
            i_lumo=i_homo+1
            for i,ei in enumerate(eigs[:i_lumo]):
                for f,ef in enumerate(eigs[i_lumo:]):
                        e = ef-ei
                        for bin in range(num_bins):
                            data_excit_all[bin] += gaussian(x_axis[bin],e,broadening/2)*kpoints_weights[k,3]
    np.savetxt(path_1+'/excit_scan_kpoints_{0}.txt'.format(index), np.vstack([x_axis,data_excit_all]).T,fmt=['%0.8f','%0.8f']) 
    np.savetxt(path_1+'/dos_scan_kpoints_{0}.txt'.format(index), np.vstack([x_axis,data_dos_all]).T,fmt=['%0.8f','%0.8f']) 
#np.savetxt(path_1+'/scan_kpoints_energies.txt'.format(index), np.vstack(total_energy).T,fmt=['%0.8f','%0.8f'])


