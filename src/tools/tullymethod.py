import numpy as np
from siesta import 
#from siesta import *
#from parser import *

#finding the harmonic frequency from the curvature of the PES 


e_fermi,eigenvalues = siesta_read_eigenvalues("MyM")
gamma=0
constant = 
window = 
delta = 
for k in range len(kpoints_weights):
    for i in n_spin:
        for e_i in eigenvalues[k,s,:] if e_i < e_fermi:
            for e_f in eigenvalues[k,s,:] if e_f > e_fermi:
                if e_f-e_in >= delta+window and e_f-e_in <= delta-window:
                    gamma=gamma+coeffs[k,s,



    
    gamma=gamma + gamma*kpoints_weights[k]
     



