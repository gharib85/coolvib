#    This file is part of coolvib
#
#        coolvib is free software: you can redistribute it and/or modify
#        it under the terms of the GNU General Public License as published by
#        the Free Software Foundation, either version 3 of the License, or
#        (at your option) any later version.
#
#        coolvib is distributed in the hope that it will be useful,
#        but WITHOUT ANY WARRANTY; without even the implied warranty of
#        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#        GNU General Public License for more details.
#
#        You should have received a copy of the GNU General Public License
#        along with coolvib.  If not, see <http://www.gnu.org/licenses/>.

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
            if len(tmp)>2:
                spectrum[i,b] = float(tmp[-2])+1.0j*float(tmp[-1])
            else:
                spectrum[i,b] = float(tmp[-1])



    plt.plot(xaxis, spectrum[index,:].real, xaxis, spectrum[index,:].imag)
    plt.xlabel('Excitation energy in eV')
    plt.ylabel('Spectral function eV^-1 * ps^-1')
    plt.show()

else:
    pass
