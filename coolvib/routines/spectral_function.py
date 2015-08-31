"""
spectral_function.py
"""

import numpy as np
from math import sqrt, pi

from coolvib.routines import fermi_occ
from coolvib.routines import delta_function 
from coolvib.routines import discretize_peak
from coolvib.constants import hplanck

#def calculate_spectral_function_mode()


def calculate_spectral_function_tensor(
        fermi_energy,
        eigenvalues,
        kpoints,
        psi,
        first_order_H,
        first_order_S,
        masses,
        pop,
        **kwargs):
    """
    Calculates spectral functions for all cartesian directions and 
    couplings between directions.

    Parameters:

    fermi_energy: float
        Fermi Level in eV

    eigenvalues: np.array
        numpy array of eigenvalues with dimensions given by 
        n_kpoints, n_spin, n_states

    kpoints: np.array
        array of kpoints with shape (n_kpoints, 4) where the 
        first three columns give the x, y and z components of 
        the kvector in reciprocal space and the fourth column 
        gives the kpoint weight.

    psi: np.array
        array of wavefunction coefficients with its shape given 
        by n_kpoints, n_spin, n_states and n_basisfunctions

    first_order_H: np.array
        cartesian derivative of the hamiltonian matrix with its 
        6 dimension given by n_atoms, 3 cartesian directions, 
        n_kpoints, n_spin, n_basis, n_basis

    first_order_S: np.array
        cartesian derivative of the overlap matrix with its 
        6 dimension given by n_atoms, 3 cartesian directions, 
        n_kpoints, n_spin, n_basis, n_basis
        
    masses: np.array or list
        list of masses for the n_atoms active atoms


    Returns:

    x_axis: np.array
        array with the discretized x axis values for the spectral function

    spectral_function:: np.array
        2-dimensional array 

    """

    #ANALYSE KWARGS

    default_keywords = {
            'discretization_type' : 'gaussian',
            'discretization_broadening' : 0.05,
            'discretization_length' : 0.01,
            'max_energy' : 3.0,
            'temperature': 300,
            }

    keys = {}
    for key in default_keywords.keys():
        if  key in kwargs:
            keys[key] = kwargs[key] 
        else:
            keys[key] = default_keywords[key]

    #build coupling matrix G
    
    n_atoms, n_cart, n_kpts, n_spin, n_states, n_basis = first_order_H.shape
    n_dim = n_atoms*n_cart

    kweights = kpoints[:,-1]

    G = np.zeros([n_dim,n_kpts,n_spin,n_states, n_basis],dtype=np.complex)
    counter = 0
    for a in range(n_atoms):
        for c in range(n_cart):
            for s in range(n_spin):
                G[counter,:,s,:,:] = first_order_H[a,c,:,s,:,:] - fermi_energy*first_order_S[a,c,:,:,:]
            counter += 1

    # calculate spectral components

    n_axis = int(keys['max_energy']/keys['discretization_length'])
    x_axis = np.array( [keys['discretization_length']*i for i in range(n_axis)] )
    ef = fermi_energy
    delta_method = keys['discretization_type']
    sigma = keys['discretization_broadening']
    T = keys['temperature']

    spectral_function = np.zeros([n_dim*(n_dim+1)/2,n_axis],dtype=np.complex)

    counter = 0
    for d in range(n_dim):
        for d2 in range(d,n_dim):
            print 'Calculating spectral function for components {0} and {1}'.format(d,d2)
            for s in range(n_spin):
                for k in range(n_kpts):
                    print 's ', s, 'k ', k
                    wk = kweights[k]
                    orb_min = 0
                    orb_lumo = 0
                    orb_homo = 0
                    orb_max = 0
                    for ei,e in enumerate(eigenvalues[k,s,:]):
                        occ = fermi_occ(e,ef,T)
                        # occ = pop[k,s,ei]
                        if e<ef-2.00*keys['max_energy']:
                            orb_min = ei
                        if occ>=0.999:
                            orb_lumo = ei
                        if occ>= 0.001:
                            orb_homo = ei
                        if e<ef+2.00*keys['max_energy']:
                            orb_max = ei
                    for i in range(orb_min,orb_homo):
                        for f in range(orb_lumo, orb_max):
                            e = eigenvalues[k,s,f] - eigenvalues[k,s,i]
                            occ =fermi_occ(eigenvalues[k,s,i],ef,T) - fermi_occ(eigenvalues[k,s,f],ef,T)
                            # occ =pop[k,s,i] - pop[k,s,f]
                            if e>0.0 and e<=1.0*keys['max_energy'] and occ>=1.E-5:
                                #calculate transition strength

                                nacs1 = np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(G[d,k,s,:,:],\
                                        psi[k,s,f,:]))
                                nacs2 = np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(G[d2,k,s,:,:],\
                                        psi[k,s,f,:]))
                                nacs = np.dot(nacs1.conjugate().transpose(),nacs2)
                                nacs /= (e*e)
                                nacs *= wk
                                nacs *= occ*(3.-n_spin)
                                # print e, ' ' , fermi_occ(eigenvalues[k,s,i],ef,T)-fermi_occ(eigenvalues[k,s,f],ef,T), ' ', nacs
                                # print e, ' ' , occ, ' ', nacs
                                spectral_function[counter,:] += discretize_peak(e,nacs, x_axis, sigma, delta_method)
            
            spectral_function[counter,:] *= (pi*hplanck) #/ sqrt(masses[d/3]*masses[d2/3])
            counter += 1


    return x_axis, spectral_function


def calculate_spectral_function_tensor_q(
        fermi_energy,
        eigenvalues,
        kpoints,
        psi,
        first_order_H,
        first_order_S,
        masses,
        basis_pos,
        **kwargs):
    """
    calculates spectral functions for all directions
    """

    #ANALYSE KWARGS

    default_keywords = {
            'discretization_type' : 'gaussian',
            'discretization_broadening' : 0.05,
            'discretization_length' : 0.01,
            'max_energy' : 3.0,
            'temperature': 300,
            }

    keys = {}
    for key in default_keywords.keys():
        if key in kwargs:
            keys[key] = kwargs[key] 
        else:
            keys[key] = default_keywords[key]

    #build coupling matrix G
    
    n_atoms, n_cart, n_kpts, n_spin, n_states, n_basis = first_order_H.shape
    n_dim = n_atoms*n_cart

    kweights = kpoints[:,-1]

    G = np.zeros([n_dim,n_kpts,n_spin,n_states, n_basis],dtype=np.complex)
    counter = 0
    for a in range(n_atoms):
        for c in range(n_cart):
            for s in range(n_spin):
                G[counter,:,s,:,:] = first_order_H[a,c,:,s,:,:] - fermi_energy*first_order_S[a,c,:,:,:]
            counter += 1

    # calculate spectral components

    n_axis = int(keys['max_energy']/keys['discretization_length'])
    x_axis = np.array( [keys['discretization_length']*i for i in range(n_axis)] )

    ef = fermi_energy
    delta_method = keys['discretization_type']
    sigma = keys['discretization_broadening']
    T = keys['temperature']

    spectral_function = np.zeros([n_dim*(n_dim+1)/2,n_axis],dtype=np.complex)

    counter = 0
    for d in range(n_dim):
        for d2 in range(d,n_dim):
            print 'Calculating spectral function for components {0} and {1}'.format(d,d2)
            for s in range(n_spin):
                for k in range(n_kpts):
                    wk = kweights[k]
                    for ei,e in enumerate(eigenvalues[k,s,:]):
                        occ = fermi_occ(e,ef,T)
                        #TODO ####WHAT ABOUT NO ORBITALS BELOW 3*keys
                        if e<ef-3.00*keys['max_energy']:
                            orb_min = ei
                        if occ>=0.999:
                            orb_lumo = ei
                        if occ>= 0.001:
                            orb_homo = ei
                        if e<ef+3.00*keys['max_energy']:
                            orb_max = ei

                    for i in range(orb_min,orb_homo):
                        for f in range(orb_lumo, orb_max):
                            e = eigenvalues[k,s,f] - eigenvalues[k,s,i]
                            if e>0.0 and e<1.25*keys['max_energy']:
                                #calculate transition strength

                                nacs1 = np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(G[d,k,s,:,:],psi[k,s,f,:]))
                                nacs2 = np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(G[d2,k,s,:,:],psi[k,s,f,:]))
                                nacs = np.dot(nacs1.conjugate().transpose(),nacs2)
                                nacs /= (e*e)
                                nacs *= wk
                                nacs *= fermi_occ(eigenvalues[k,s,i],ef,T) - fermi_occ(eigenvalues[k,s,f],ef,T)
                                spectral_function[counter,:] += discretize_peak(e,nacs, x_axis, sigma, delta_method)

            spectral_function[counter,:] *= (hplanck*pi) / sqrt(masses[d/3]*masses[d2/3])
            counter += 1


    return x_axis, spectral_function


