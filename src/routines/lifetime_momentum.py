import numpy as np
from vibcooling.constants import hbar, hbar_1, pi, sqrt_pi, conversion_factor
import vibcooling.routines.lifetime_ext as lt_ext
from math import exp

from time import time

def generate_delta_approximation(dirac_method, fermi_centered, delta, window):
    """
    function generator that returns a function, which 
    yields true or false if its argument is inside the given 
    window defined by minimum and maximum.
    """
    def dirac_weight_window(exc):
        #if fermi_centered:
            #return (exc*exc)/(window*window)
        #else:
            return exc/window
    def dirac_weight_gaussian(exc):
        gaussian = exp(-((exc-delta)*(exc-delta))/(window*window))/(window*sqrt_pi)
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
        if minimum<0.:
            maximum = maximum-minimum
            minimum=0.
    
    def e_decider(exc):
        if minimum < exc < maximum:
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
    over fermi_level
    """

    delta = hbar*omega*conversion_factor # from Hz to eV
    inside_window, dirac_weight = generate_delta_approximation(dirac_method, fermi_centered, delta, window)

    nbasis = psi.shape[-1]
    #prepare N vector
    kpts = kpoints_weights[:,:3]
    kweights = kpoints_weights[:,3]
    nmin = kpts[np.nonzero(kpts)].min()
    N = np.array(np.round(kpts/nmin),dtype=np.float)
    for n, nvec in enumerate(N):
        N[n] = np.dot(nvec, cell)

    gamma=0
    hit_the_window=0
    for s in range(n_spin):
        #first loop over states from image cells from which to excite
        for n1, Nvec1 in enumerate(N):
            eigenvalues_k1 = eigenvalues[n1,s,:]
            #identify loop boundaries
            for ei,e in enumerate(eigenvalues_k1):
                if e<fermi_level-delta:
                    orb_min_k1 = ei
                if e<fermi_level:
                    orb_max_k1 = ei
            #second loop over states from image cells into which to excite
            for n2, Nvec2 in enumerate(N):
                eigenvalues_k2 = eigenvalues[n2,s,:]
                #identify loop boundaries
                for ei,e in enumerate(eigenvalues_k2):
                    if e>fermi_level:
                        orb_min_k2 = ei
                        break
                for ei,e in enumerate(eigenvalues_k2):
                    if e>fermi_level+delta:
                        orb_max_k2 = ei
                        break
                product = 0.0
                #first loop over states from which to excite
                for i in range(orb_min_k1,orb_max_k1):
                    #second loop over states into which to excite
                    for f in range(orb_min_k2,orb_max_k2):
                        e = eigenvalues_k2[f] - eigenvalues_k1[i]
                        if inside_window(e):
                            #CYTHON FUNCTION FOR NONADIABATIC COUPLING
                            nacs = lt_ext.calculate_nonadiabatic_coupling_HS_kpts(
                                    i, f, psi[:,s,:,:], Nvec1, Nvec2, 
                                    H_plus[:,s,:,:], H_minus[:,s,:,:], S_plus, S_minus, 
                                    atomic_positions, basis_pos, kpts,kweights, 
                                    dq, fermi_level, e,
                                    )
                            #times dirac delta weight, currently a constant, but could be 
                            #a gaussian prefactor
                            nacs *= dirac_weight(e)
                            hit_the_window+=1     
                            product += nacs 
                            #EVERYTHING ENDS
                    ####COMING BACK FROM CYTHON
                    #end of e2 loop
                #end of e1 loop
                gamma += product
            #end of k2 loop
        #end of k1 loop
    #end of spin loop
    gamma *=pi*hbar_1 /len(N)
    return gamma, hit_the_window


def lifetime_routine(kpoints_weights,n_spin,eigenvalues,
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
    over fermi_level
    """

    delta = hbar*omega*conversion_factor # from Hz to eV
    inside_window, dirac_weight = generate_delta_approximation(dirac_method, fermi_centered, delta, window)

    nbasis = psi.shape[-1]
    #prepare N vector
    kpts = kpoints_weights[:,:3]

    gamma=0
    hit_the_window=0
    for s in range(n_spin):
        #first loop over states from image cells from which to excite
        for k in range(len(kpts)):
            wk = kpoints_weights[k,3]
            eigenvalues_k1 = eigenvalues[k,s,:]
            for ei,e in enumerate(eigenvalues_k1):
                if e<fermi_level-delta:
                    orb_min = ei
                if e<fermi_level:
                    orb_max = ei
            #second loop over states from image cells into which to excite
            product = 0.0
            #first loop over states from which to excite
            for i in range(orb_min,orb_max):
                for f in range(orb_min,orb_max):
                    e = eigenvalues_k1[f] - eigenvalues_k1[i]
                    if inside_window(e):
                        hit_the_window +=1
                        #CYTHON FUNCTION FOR NONADIABATIC COUPLING
                        H_q=(H_plus[k,s,:,:]-H_minus[k,s,:,:])/2/dq
                        S_q=(S_plus[k,:,:]-S_minus[k,:,:])/2/dq
                        nacs= np.dot(psi[k,s,i,:].conjugate().transpose(),
                                np.dot(H_q-fermi_level*S_q,psi[k,s,f,:]))
                        nacs = np.square(np.absolute(nacs))
                        nacs /= (e*e)
                        #times dirac delta weight, currently a constant, but could be 
                        #a gaussian prefactor
                        nacs *= dirac_weight(e)
                        product += nacs 
                        #EVERYTHING ENDS
                ####COMING BACK FROM CYTHON
                #end of e2 loop
            #end of e1 loop
            gamma += wk*product
        #end of k1 loop
    #end of spin loop
    gamma *=pi*hbar_1
    return gamma, hit_the_window

