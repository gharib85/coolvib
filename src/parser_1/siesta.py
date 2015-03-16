"""
Parsing routine for Siesta
"""

import numpy as np


def siesta_read_eigenvalues(filename):
    """
    This routine reads the number of K points, the eigenvalues and the \
    and the number of states

    Input: the name of the eigenvalue file
    Output: Fermi_level, Eigenvalues i

    The Eigenvalues is a 3D numpy array:
    1) Bands
    2) K points
    3) Spin
    """

    with open("MyM.EIG", "r") as f:
        content=f.readlines()
        fermi_level=content[0].split()[0]
        n_bands= int(content[1].split()[0])
        n_spin = int(content[1].split()[1])
        n_kpoints = int(content[1].split()[2])
        eigenvalues=np.zeros((n_bands,n_kpoints,n_spin))

        energy=[]
        counter=1.0
        line_number=2
        for k in range(0,n_kpoints):
            energy+=content[line_number].split()[1:]
            for n in range(1, int(np.ceil(n_bands/10.0))):
                line_number+=1
                energy+=content[line_number].split()
            for i in range(0, n_bands):
                eigenvalues[i,k,0] = energy[i]
            line_number+=1
            energy=[]

    return(fermi_level,eigenvalues)
    


def siesta_read_kpoints(filename):
    """
    This routine reads the weights of K points
    """

    with open("MyM.KP", "r") as f:
        content=f.readlines()
        n_kpoints=int(content[0].split()[0])
        kpoints_weights=np.zeros((n_kpoints))
        line_number=1
        while line_number <=len(content)-1:
            kpoints_weights[line_number-1]=content[line_number].split()[4]
            line_number+=1

    return kpoints_weights


def siesta_read_eigenvectors(filename):
    """
    This routine reads siesta eigenvectors
    """
    return



"""
To call the routine from this file just type 
from siesta py import function
function(a,b)
"""
