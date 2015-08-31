from sys import argv

import numpy as np
import matplotlib.pyplot as plt
from math import pi, factorial



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
plt.show()
