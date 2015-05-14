import numpy as np
from vibcooling.constants import hbar_1, pi, sqrt_pi, conversion_factor

from math import exp

def generate_delta_approximation(dirac_method, fermi_centered, delta, window):
    """
    function generator that returns a function, which 
    yields true or false if its argument is inside the given 
    window defined by minimum and maximum.
    """
    def dirac_weight_window(e):
        return 1./window
    def dirac_weigh_gaussian(e):
        gaussian = exp(-(e*e)/(window*window))/(window*sqrt_pi)
        return gaussian 
    
    if dirac_method is 'window':
        win = window
        dirac_weight = dirac_weight_window
    if dirac_method is 'gaussian':
        win = window*4.
        dirac_weight = dirac_weight_gaussian
    
    if fermi_centered:
        minimum=0
        maximum=win
    else:
        minimum=delta-win/2.
        maximum=delta+win/2.
    
    def e_decider(e):
        if minimum < e < maximum:
            return True
        else:
            return False
    
    return e_decider, dirac_weight

def lifetime_routine_momentum(kpoints_weights,n_spin,eigenvalues,
        fermi_level,omega,window,psi,
        H_plus,H_minus,S_plus,S_minus,
        dq,atomic_positions,cell,basis_pos,
        fermi_centered = False,
        dirac_method='window', #'gaussian'
        ):
    """
    This routine calculates the vibrational lifetime for all excitations between different k points 
    with the deltafunction centered over the excitation hbar*omega.

    In addition to other input, this routine depends on atomic_positions with dimension (Na*3) 
    and basis_pos, an array with length N_basis that contains the indices of atoms on which the 
    basis functions sit.
    
    the flag fermi_centered decides if the delta function is centered over omega or 
    over fermi_energy
    """
    
    #TODO implement basis_pos in parser
    delta = hbar*omega*conversion_factor # in Hz 
    if fermi_centered:
        inside_window, dirac_weight = generate_delta_approximation(dirac_method, fermi_centered, delta, window)

    nbasis = psi.shape[-1]
    #prepare N vector
    kpts = kweights[:,:3]
    nmin = kpts[np.nonzero(kpts)].min()
    N = np.array(np.round(kpts/nmin)),dtype=np.int)
    for n, nvec in enumerate(N):
        N[n] = np.dot(nvec, cell)
    
    gamma=0
    product=0
    hit_the_window=0
    for s in range(n_spin):
        #first loop over k points from which to excite
        for k1,kvec1,wk1 in enumerate(zip(kpoints_weights[:,:3],kpoints_weights[:,4])):
            eigenvalues_k1 = eigenvalues[k1,s,:]
            psi_k1 = psi[k1,s,:,:]
            #identify loop boundaries
            for ei,e in enumerate(eigenvalues_k1):
                if e<fermi_energy-delta:
                    orb_min_k1 = ei
                    if e<fermi_energy:
                        orb_max_k1 = ei
            #second loop over k points into which to excite
            for k2,kvec2,wk2 in enumerate(zip(kpoints_weights[:,:3],kpoints_weights[:,4])):
                eigenvalues_k2 = eigenvalues[k2,s,:]
                psi_k2 = psi[k2,s,:,:]
                #identify loop boundaries
                for ei,e in enumerate(eigenvalues_k2):
                    if e>fermi_energy:
                        orb_min_k2 = ei
                        break
                for ei,e in enumerate(eigenvalues_k2):
                    if e<fermi_energy+delta:
                        orb_max_k2 = ei
                        break
                wk = wk1*wk2 
                gamma_tmp = 0.0
                #first loop over states from which to excite
                for i in range(orb_min_k1,orb_max_k1):
                    #second loop over states into which to excite
                    for f in range(orb_min_k2,orb_max_k2):
                        e = eigenvalues_k2[f] - eigenvalues_k1[i]
                        if inside_window(e):
                            #HERE GOES the Hamiltonian that is now defined in the real 
                            #auxiliary supercell defined by the translated images of 
                            #the original unit cell
                            #loop over all real space images of the unit cell 
                            product = 0.0
                            for Nvec1 in N: 
                                for Nvec2 in N:
                                    for a in range(nbasis):
                                        RNvec1 = atomic_positions[basis_pos[a]]+Nvec1
                                        for b in range(nbasis):
                                            RNvec2 = atomic_positions[basis_pos[b]]+Nvec2
                                            tmp = 0.0
                                            #Fourier transformation of G matrix element
                                            for kl,kvec in enumerate(kpts): 
                                                H_q=(H_plus[kl,s,a,b]-H_minus[kl,s,a,b])/2/dq
                                                S_q=(S_plus[kl,a,b]-S_minus[kl,a,b])/2/dq
                                                phase = exp(-1j*np.dot(kvec,RNvec1-RNvec2)) 
                                                tmp += phase*(H_q-fermi_level*S_q)
                                            #*scattering phase shift
                                            tmp *= exp(1j*(np.dot(kvec2,RNvec2)-np.dot(kvec1,RNvec1)))
                                            #wavefunc coeffs
                                            tmp *= psi_k1[ns,k1,i,a].conjugate()*psi_k2[ns,k2,f,b]
                                            product += tmp
                                        #end of a, b basis loops
                                #end of N1,N2 loops
                            gamma_tmp += np.square(np.absolute(product))
                            gamma_tmp /= e
                            #times dirac delta weight, currently a constant, but could be 
                            #a gaussian prefactor
                            gamma_tmp *= dirac_weight(e)
                            hit_the_window+=1     
                            #EVERYTHING ENDS
                    #end of e2 loop
                #end of e1 loop
                gamma += gamma_tmp*wk
            #end of k2 loop
        #end of k1 loop
    #end of spin loop
    print "         hit_the_window", hit_the_window
    gamma *=pi/hbar_1/window#/window
    return gamma, hit_the_window





