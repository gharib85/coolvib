#!/usr/bin/python
"""
plot_spectral_function.py plots the spectral function 
of a given component using matplotlib

python plot_spectral_function.py <filename> <component no.>

Arguments:

    <filename>: str
        Name of the file holding the spectral functions

    <component no.>: int
        Index of the spectral function to be plotted

"""

from sys import argv, exit

import numpy as np
from math import pi, factorial

if __name__=="__main__":

    try:
        import matplotlib.pyplot as plt
    except:
        raise ImportError('Cannot import matplotlib')

    f = open(str(argv[1]))#'nacs-spectrum.out')
    index = argv[2]


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


    plt.plot(xaxis, spectrum[index,:].real, xaxis, spectrum[index,:].imag)
    plt.xlabel('Excitation energy in eV')
    plt.ylabel('Spectral function eV^-1 * ps^-1')
    plt.show()

else:
    pass
