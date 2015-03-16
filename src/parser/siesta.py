"""
Parsing routine for Siesta
"""

import numpy as np
import struct
from scipy.io import FortranFile


def siesta_read_eigenvalues(filename):
    """
    This routine reads the number of K points, the eigenvalues and the \
    and the number of states

    Input: the name of the eigenvalue file
    Output: Fermi_level, Eigenvalues i

    The Eigenvalues is a 3D numpy array:
    1) K point
    2) Spin
    3) State
    """

    with open("%s.EIG" % filename, "r") as f:
        content=f.readlines()
        fermi_level=content[0].split()[0]
        n_bands= int(content[1].split()[0])
        n_spin = int(content[1].split()[1])
        n_kpoints = int(content[1].split()[2])
        eigenvalues=np.zeros((n_kpoints,n_spin,n_bands))
    
        energy=[]
        counter=1.0
        line_number=2
        for k in range(0,n_kpoints):
            energy+=content[line_number].split()[1:]
            for n in range(1, int(np.ceil(n_bands*n_spin/10.0))):
                line_number+=1
                energy+=content[line_number].split()
            for s in range (n_spin):
                for i in range(0, n_bands):
                    eigenvalues[k,s,i] = energy[i+n_bands*s]
            line_number+=1
            energy=[]
    
        return(fermi_level,eigenvalues)    


def siesta_read_kpoints(filename):
    """
    This routine reads the weights of K points and returns a 2 d array with dimensions:
    1) index of the k point 
    2) x,y,z coordinate and weight
    """

    with open("%s.KP" % filename, "r") as f:
        content=f.readlines()
        n_kpoints=int(content[0].split()[0])
        kpoints_weights=np.zeros((n_kpoints,4))
        line_number=1
        while line_number <=len(content)-1:
            kpoints_weights[line_number-1]= content[line_number].split()[1:] 
            line_number+=1
    return kpoints_weights


def siesta_read_coefficients(filename):
    """
    This routine reads siesta eigenvectors from MyM.WFSX binary
    using the FortranFile magical binary reader
    and returns an array with the following dimensions:
    1) the k-point index
    2) the spin index
    3) the state index (nwflist)
    4) the coefficient for each basis set index (nuotot)
    """

    f= FortranFile("%s.WFSX" % filename)
    nk,gamma= f.readInts()
    print "nk, gamma  = ", nk,gamma
    nspin = int(f.readInts())
    print "nspin =",nspin
    nuotot = int(f.readInts())
    print "nuotot =",nuotot
    f.readRecord()
    psi=np.zeros((nk,nspin,nuotot,nuotot))
    psi = psi +0j
    for iik in range (1,nk+1):  #for each k-point
        for iispin in range (1,nspin+1):  #for each spin
            f.readRecord()
            ispin = int(f.readInts())
            nwflist =int(f.readInts())
            for iw in range(1,nwflist+1):  # for each state (nwflist = total number of states)
                indwf=f.readInts()
                energy=f.readReals('d')    # we first read the energy of the state 
                read_psi = f.readReals()   # and all the orbital coefficients (real value, followed by the imaginary value
                read_psi=np.reshape(read_psi, (nwflist,2))  # reshape it 
                psi[iik-1,iispin-1,iw-1,:]=read_psi[:,0]+1j*read_psi[:,1]  # and make a row of complex numbers
    return psi



"""
To call the routine from this file just type 
from siesta py import function
function(a,b)
"""
