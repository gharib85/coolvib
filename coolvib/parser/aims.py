"""
Parsing routines for FHI-Aims 

Copyright Reinhard J. Maurer, Yale University, 04/08/2015
"""

import numpy as np
import os

def aims_read_fermi_and_kpoints(filename,cell=None):
    """
    Reads the fermi_level and the kpoint information from the 
    main OUTPUT file of FHI-Aims. In order for this to work one needs 
    to set the keyword 'output ki_point_list
    
    written by Reinhard J. Maurer, Yale University, 04/08/2015
    """

    if cell is None:
        raise ValueError('Please supply unit cell for aims_read_fermi_kpoints')
    rcell = 2.*np.pi*np.linalg.inv(cell.T)

    with open(filename,'r') as f:
        
        kweights_exist = False
        fermi_level_exists = False

        content = f.readlines()
        for l, line in enumerate(content):
            if 'k_point_list' in line:
                print 'Found k_point_list keyword, extracting kpoint_weights'
                kweights_exist = True

            if '| Chemical potential (Fermi level) in eV  ' in line \
                    and not fermi_level_exists:
                print 'Found Fermi Level, extracting fermi level'
                fermi_level_exists = True

        if fermi_level_exists and kweights_exist:

            kpoints = []
            fermi_level = 0.0
            k_points_are_done = False
            for l, line in enumerate(content):
                if '| K-points in task ' in line and not k_points_are_done:
                    nkpts = int(line.split()[-1])
                    n_line = l+1
                    while 'K-points in task  ' in content[n_line]:
                        nkpts += int(content[n_line].split()[-1])
                        n_line += 1
                    for nk in range(nkpts):
                        kpoint = []
                        kpoint.append(float(content[n_line].split()[4]))
                        kpoint.append(float(content[n_line].split()[5]))
                        kpoint.append(float(content[n_line].split()[6]))
                        kpoint.append(float(content[n_line].split()[9]))
                        kpoints.append(kpoint)
                        n_line += 1
                    k_points_are_done = True

                if '| Chemical potential (Fermi level) in eV      ' in line:
                    fermi_level = float(line.split()[-1])
            
            kpoint_weights = np.array(kpoints)
            #transform k points fromr reciprocal to absolute values
            for i in range(len(kpoint_weights)):
                kpoint_weights[i,:3] = np.dot(kpoint_weights[i,:3],rcell)
            return (fermi_level, kpoint_weights)

        else:
            raise RuntimeError('Cannot extract fermi level and/or kpoint_weights \n \
                    Either k_point_list keyword not set, or run is not converged')


def aims_read_eigenvalues_and_coefficients(fermi_level, directory='./', spin=False):
    """
    This routine reads eigenvalues and eigenvectors from aims output
    generated via the keywords 'output band x1 y1 z1 x2 y2 z2 2' and 
    'output eigenvectors'.

    arguments are the directory where the 
    KS_eigenvectors_XX.band_X.kpt_X.out files lie, the value of the fermi_level, 
    and a logical specifying if the calculation was spin-collinear or without spin.

    the routine returns the eigenvalues, eigenvectors, and occupations
    in this order as np.arrays

    written by Reinhard J. Maurer, Yale University, 04/08.2015
    """

    #name_base = directory+'/KS_eigenvectors'
    name_base = 'KS_eigenvectors'
    #spin 0 is dn spin 1 is up
    spin_postfix = ['_dn', '_up']
    band_basis = 2

    eigenvector_files = []
    #how many files exist
    for file in os.listdir(directory):
        if file.startswith(name_base):
            eigenvector_files.append(directory+'/'+file)

    if spin:
        n_kpts = len(eigenvector_files)/2
        n_spin = 2
    else:
        n_kpts = len(eigenvector_files)
        n_spin = 1

    name_base = directory + '/' + name_base

    #now we need to find out how many basis functions are used
    if spin:
        filename = name_base + '_dn.band_1.kpt_1.out'
    else:
        filename = name_base + '.band_1.kpt_1.out'
    with open(filename, 'r') as f:
        lines = f.readlines()
        n_basis = len(lines) - 7
        n_states = len(lines[4].split()[3:])
    #print n_basis, n_kpts, n_spin
    
    eigenvalues = np.zeros([n_kpts, n_spin, n_states])
    occ = np.zeros([n_kpts, n_spin, n_states])
    psi = np.zeros([n_kpts, n_spin, n_states, n_basis],dtype=complex)
    orbital_pos = np.zeros(n_basis,dtype=np.int)
    #loop through all files
    for s in range(n_spin):
        prefix = name_base
        if spin:
            prefix += spin_postfix[s]
        #this now assumes 2 kpts per 'band'
        for k in range(n_kpts):
            band = int(k/band_basis)+1
            kp = k%band_basis+1
            filename = prefix + '.band_{0}.kpt_{1}.out'.format(band,kp)
            #opening file
            with open(filename,'r') as f:
                print 'Reading eigenvalues and psi from {0} '.format(filename)
                lines = f.readlines()
                #line 4 contains the eigenvalues
                eigenvalues[k,s,:] = np.array(lines[4].split()[3:]).astype(np.float)
                #line 5 contains the occupations
                occ[k,s,:] = np.array(lines[5].split()[3:]).astype(np.float)
                nline = 7
                for i in range(n_basis):
                    read_psi = np.array(lines[nline+i].split()[6:]).astype(np.float).reshape(-1,2)
                    psi[k,s,:,i] = read_psi[:,0] + 1j*read_psi[:,1]
                    orbital_pos[i] = int(lines[nline+i].split()[1]) -1

    return eigenvalues+fermi_level, psi, occ, orbital_pos

def aims_read_HS(directory='./', spin=False):
    """
    This routine extracts hamiltonian and overlap matrix from FHI-Aims 
    outputfiles generated with the keyword 'output band ....' and 
    output 'output hamiltonian_matrix' and 'output overlap_matrix'

    the function takes as arguments the directory where these files can be found 
    and a logical specifying the spin-mode.

    the function returns the Hamiltonian and overlap matrix as np.arrays

    written by Reinhard J. Maurer, Yale University, 04/08.2015
    """

    name_base_H = 'KS_hamiltonian_matrix'
    name_base_S = 'KS_overlap_matrix'
    band_basis = 2

    H_files = []
    S_files = []
    #how many files exist
    for file in os.listdir(directory):
        if file.startswith(name_base_H):
            H_files.append(file)
        if file.startswith(name_base_S):
            S_files.append(file)

    assert len(H_files) == len(S_files)

    n_kpts = len(H_files)
    if spin:
        n_spin = 2
    else:
        n_spin = 1

    #now we need to find out how many basis functions are used
    filename = directory+'/'+name_base_S + '.band_1.kpt_1.out'
    with open(filename, 'r') as f:
        lines = f.readlines()
        n_basis = len(lines) -2 
    
    H = np.zeros([n_kpts, n_spin, n_basis, n_basis],dtype=complex)
    S = np.zeros([n_kpts, n_basis, n_basis],dtype=complex)

    for k in range(n_kpts):
        band = int(k/band_basis)+1
        kp = k%band_basis+1
        filename_H = directory + '/'+name_base_H + '.band_{0}.kpt_{1}.out'.format(band,kp)
        filename_S = directory + '/'+name_base_S + '.band_{0}.kpt_{1}.out'.format(band,kp)
        #opening files
        with open(filename_S,'r') as f:
            print 'Reading overlap_matrix from {0} '.format(filename_S)
            s = np.loadtxt(filename_S).reshape([n_basis,n_basis,2])
            S[k, :, :] = s[:,:,0] + 1j*s[:,:,1]
        with open(filename_H,'r') as f:
            print 'Reading hamiltonian matrix from {0} '.format(filename_H)
            h = np.loadtxt(filename_H)
            if spin:
                h_dn = h[:n_basis].reshape([n_basis,n_basis,2])
                h_up = h[n_basis:].reshape([n_basis,n_basis,2])
                H[k, 0, :, :] = h_dn[:,:,0] + 1j*h_dn[:,:,1]
                H[k, 1, :, :] = h_up[:,:,0] + 1j*h_up[:,:,1]
            else:
                h = h.reshape([n_basis,n_basis,2])
                H[k, 0, :, :] = h[:,:,0] + 1j*h[:,:,1]

    return H*27.211384500, S

