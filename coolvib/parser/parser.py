"""
Parent module and routine for parsing input from 
different QM codes

    all routines assume the following filestructure:
    all data is contained in a main directory 'path' containing subdirectories
    with following syntax:
    a{0}_c{1}_plus/
    a{0}_c{1}_minus/
    eq/

written by Reinhard J. Maurer, Yale University, 03/17/2015
"""
from coolvib.parser.siesta import *
from coolvib.parser.aims import *
import numpy as np
from os.path import join as pathjoin

def parse_siesta(model, path='./',seed='MyM',atoms=[1], incr=0.01,debug=0):
    """
    This subroutine returns all necessary information stored in the siesta output files
    and returns the fermi level, the eigenvalues, the kpoints and the hamiltonian and the 
    overlap matrices

    Input
        path : string
        seed : string
        atoms : list
        incr : float
        debug : 0 / 1

    Output
        eigenvalues
        fermi_level
        kpoints_weights
        psi
        basis_pos
        cell
        atomic_positions
        H_q
        S_q
    """

    fermi_energy, eigenvalues = siesta_read_eigenvalues(path+'/eq/'+seed)
    n_spin = eigenvalues.shape[1]
    kpoints_weights = siesta_read_kpoints(path+'/eq'+seed)
    psi, basis_pos = siesta_read_coefficients(path+'/eq/'+seed)
    n_states = psi.shape[2]
    n_basis = psi.shape[3]
    cell, atomic_positions = siesta_read_struct_out(path+'/eq/'+seed)

    H_q = np.zeros([len(atoms),3,len(kpoints_weights), n_spin, n_states, n_basis],dtype=np.complex)
    S_q = np.zeros([len(atoms),3,len(kpoints_weights), n_spin, n_states, n_basis],dtype=np.complex)

    for a in atoms:
        for c in range(3):
            H_plus, S_plus = siesta_read_HSX(kpoints_weights, path+'/a{0}c{1}+/'.format(int(a),int(c))+seed,debug=debug)
            H_minus, S_minus = siesta_read_HSX(kpoints_weights, path+'/a{0}c{1}-/'.format(int(a),int(c))+seed,debug=debug)
            H_q[a,c,:,:,:,:] = (H_plus - H_minus)/2.0/incr
            S_q[a,c,:,:,:] = (S_plus - S_minus)/2.0/incr

    model.eigenvalues = eigenvalues
    model.fermi_energy = fermi_energy
    model.kpoints = kpoints_weights
    model.psi = psi
    model.basis_pos = basis_pos
    model.first_order_H = H_q
    model.first_order_S = S_q
    

def parse_aims(model, cell, spin=True, path='./', filename='aims.out', atoms=[1], incr=0.01, debug=0):
    """
    This subroutine returns all necessary information stored in the siesta output files
    and returns the fermi level, the eigenvalues, the kpoints and the hamiltonian and the 
    overlap matrices

    Input
        cell : np.array
        spin : boolean
        path : string
        filename : string
        atoms : list
        incr : float
        debug : 0 / 1

    Output
        eigenvalues
        fermi_level
        kpoints_weights
        psi
        basis_pos
        occ
        H_q
        S_q

    """

    fermi_level, kpoints_weights = aims_read_fermi_and_kpoints(pathjoin(path+'/eq', filename), cell)
    eigenvalues, psi, occ, basis_pos = \
            aims_read_eigenvalues_and_coefficients(fermi_level, directory=pathjoin(path,'eq'),spin=spin)
    if spin:
        n_spin=2
    else:
        n_spin=1
    n_states = psi.shape[2]
    n_basis = psi.shape[3]

    H_q = np.zeros([len(atoms),3,len(kpoints_weights), n_spin, n_states, n_basis],dtype=np.complex)
    S_q = np.zeros([len(atoms),3,len(kpoints_weights), n_spin, n_states, n_basis],dtype=np.complex)

    for a in atoms:
        for c in range(3):
            H_plus, S_plus = aims_read_HS(path+'/a{0}c{1}+/'.format(int(a),int(c)),spin=spin)
            H_minus, S_minus = aims_read_HS(path+'/a{0}c{1}-/'.format(int(a),int(c)),spin=spin)
            H_q[a,c,:,:,:,:] = (H_plus - H_minus)/2.0/incr
            S_q[a,c,:,:,:] = (S_plus - S_minus)/2.0/incr

    model.eigenvalues = eigenvalues
    model.fermi_energy = fermi_level
    model.kpoints = kpoints_weights
    model.psi = psi
    model.basis_pos = basis_pos
    model.occ = occ
    model.first_order_H = H_q
    model.first_order_S = S_q



