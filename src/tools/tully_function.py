import numpy as np
from vibcooling.parser.siesta import *

hbar=6.58211928E-16 #eV*s
pi=3.14159265359
conversion_factor = 98226949774380.3 #for omega from sqrt(eV/amu)/Ang to Hz
def tully_routine(kpoints_weights,n_spin,eigenvalues,fermi_level,omega,window,psi,H_plus,H_minus,S_plus,S_minus,dq):
    delta = hbar*omega*conversion_factor # in Hz 
#    print "     delta",delta
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
                        elif (delta-window/2) <=  (eigenvalues[k,s,f]-eigenvalues[k,s,i]) <= (delta+window/2):
                            H_q=(H_plus[k,s,:,:]-H_minus[k,s,:,:])/2/dq
                            S_q=(S_plus[k,:,:]-S_minus[k,:,:])/2/dq
                            product= np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(H_q-fermi_level*S_q,psi[k,s,f,:]))
#                             gamma += np.square(np.absolute(product))*kpoints_weights[k,3]*np.square(eigenvalues[k,s,f]-eigenvalues[k,s,i])
                            gamma += np.square(np.absolute(product))*kpoints_weights[k,3]*(eigenvalues[k,s,f]-eigenvalues[k,s,i])
                            hit_the_window+=1     
    print "         hit_the_window", hit_the_window
    #gamma *=pi/hbar/omega/omega
    gamma *=pi/hbar/omega/omega#/window
    return gamma, hit_the_window





