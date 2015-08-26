import numpy as np
from vibcooling.parser.siesta import *
from vibcooling.tools.spectrum  import *

kmin=4
kmax=5

index = 4
path = os.getcwd()
path_1=path
path= path_1+"/{0}".format(index)
print path

fermi_energy, eigenvalues = siesta_read_eigenvalues(path+"/eq/MyM")
fermi_level= float(fermi_energy)
kpoints_weights = siesta_read_kpoints(path+"/eq/MyM")
kweights = kpoints_weights[:,3]

x, y = calculate_spectrum(eigenvalues, fermi_level, cutoff=5.0,  keights)
#x, y = calculate_dos(eigenvalues, keights, fermi_shift=fermi_level)

np.savetxt('spectrum.txt', np.vstack([x,y]).T,fmt=['%16.8f','%16.8f'])
