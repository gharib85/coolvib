"""

TEMPORARY FUNCTION BUCKET

TODO
    Needs to be cleaned out and structured
"""
#TODO split this file into different scripts that 
#each calculate DOS etc. 

import numpy as np
import os, sys
import math
from time import time
try:
    import matplotlib.pyplot as plt
except:
    raise ImportError('Cannot import matplotlib!')

broadening_default = 0.1
bin_width_default = 0.01

def gaussian(x, x_mean, broadening):
    
    gaussian_val = np.exp(-0.5*(((x-x_mean) / (broadening))**2));
    return gaussian_val

def dos_binning(eigenvalues,broadening=broadening_default, bin_width=bin_width_default,
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
    x_axis = np.zeros([num_bins])
    data = np.zeros([num_bins])

    #setting up x-axis
    for i in range(num_bins):
        x_axis[i] = lowest_e + i * bin_width
    #get DOS
    for ei,e in enumerate(eigenvalues):
        for i in range(num_bins):
            data[i] += gaussian(x_axis[i],e,broadening)*coeffs[ei]
#    data /= len(eigenvalues)

    return x_axis, data

def calculate_dos(eigenvalues, kweights=None, broadening=broadening_default, 
        bin_width=bin_width_default, fermi_shift = None):
    """
    calculate the Density of States
    """

    nkpts, nspin, nstates = eigenvalues.shape
 
    if kweights is None:
        kweights = np.ones(nkpts)
    if fermi_shift is None:
        fermi_shift = 0.0
    
    energies = []
    weights = []

    for s in range(nspin):
        for k in range(nkpts):
            for i in range(nstates):
                excit = eigenvalues[k,s,i] - fermi_shift
                energies.append(excit)
                weights.append(kweights[k])
    
    x_axis, dos = dos_binning(energies, broadening, bin_width, weights)
    return x_axis, dos

def calculate_spectrum(eigenvalues, fermi_level, cutoff = 5.0, 
        kweights=None, broadening=broadening_default, bin_width=bin_width_default):
    """
    Calculate the excitation spectrum between excitations within k points
    """
    nkpts, nspin, nstates = eigenvalues.shape
 
    if kweights is None:
        kweights = np.ones(nkpts)
    
    energies = []
    weights = []

    for s in range(nspin):
        for k in range(nkpts):
            for i in range(nstates):
                if eigenvalues[k,s,i]<fermi_level-cutoff:
                    orb_min = i
                if eigenvalues[k,s,i]<fermi_level:
                    orb_fermi = i+1
                if eigenvalues[k,s,i]<fermi_level+cutoff:
                    orb_max = i
            for i in range(orb_min, orb_fermi):
                for f in range(orb_fermi, orb_max):
                    excit = eigenvalues[k,s,f] - eigenvalues[k,s,i]
                    if excit < cutoff:
                        energies.append(excit)
                        weights.append(kweights[k])
    
    x_axis, spectrum = dos_binning(energies, broadening, bin_width, weights)
    return x_axis, spectrum


def calculate_spectrum_scattering(eigenvalues, fermi_level, cutoff = 5.0, 
        kweights=None, broadening=0.01, bin_width=0.01):
    """
    Calculate the excitation spectrum between excitations between different k points
    """
    nkpts, nspin, nstates = eigenvalues.shape
 
    if kweights is None:
        kweights = np.ones(nkpts)
    
    energies = []
    weights = []

    for s in range(nspin):
        for k1 in range(nkpts):
            for k2 in range(nkpts):
                for i in range(nstates):
                    if eigenvalues[k1,s,i]<fermi_level-cutoff:
                        orb_min = i
                    if eigenvalues[k1,s,i]<fermi_level:
                        orb_fermi1 = i+1
                    if eigenvalues[k2,s,i]<fermi_level:
                        orb_fermi2 = i
                    if eigenvalues[k2,s,i]<fermi_level+cutoff:
                        orb_max = i
                for i in range(orb_min, orb_fermi1):
                    for f in range(orb_fermi2, orb_max):
                        excit = eigenvalues[k2,s,f] - eigenvalues[k1,s,i]
                        if excit < cutoff:
                            energies.append(excit)
                            weights.append(kweights[k1]*kweights[k2])
    
    x_axis, spectrum = dos_binning(energies, broadening, bin_width, weights)
    return x_axis, spectrum

