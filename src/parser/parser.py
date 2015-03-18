"""
Parent module and routine for parsing input from 
different QM codes

written by Mikhail Askerka and Reinhard J. Maurer, Yale University, 03/17/2015
"""

from __future__ import print_function
import os
import glob


default_keywords = {
        'filename' : 'MyM',
        'qm_code' : 'Siesta',
 #       'atoms_disp': [0],

        }


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--qm_code","-qm",  help = "The electronic structure code you are using", \
        default = default_keywords['qm_code'],type =str)
parser.add_argument("--filename","-f", help = "The name of the calculaton file", \
        default = default_keywords['filename'], type = str)
args = parser.parse_args()


#vibcool -f "filename" -q "Siesta" 

def parser_read_input(keywords_dict):
    """
    This subroutine returns all necessary information stored in the siesta output files
    and returns the fermi level, the eigenvalues, the kpoints and the hamiltonian and the 
    overlap matrices
    """
    qm_code = args.qm_code
    filename = args.filename
    if qm_code is 'SIESTA' or "Siesta" or "siesta":
        from siesta import siesta_read_eigenvalues, \
                siesta_read_kpoints, siesta_read_coefficients
        fermi_level, eigenvalues = siesta_read_eigenvalues(filename)
        kpts_array = siesta_read_kpoints(filename)
        coeffs = siesta_read_coefficients(filename)
        #read H and S of eq.
        #H, S = siesta_read_hs(kpts_array, filename)
        #Hdim1,Hdim2,Hdim3,Hdim4 = H.shape
        #Hd = np.zeros[Hdim1,Hdim2,Hdim3,Hdim4,len(atoms_disp)*3]
    else:
        raise NotImplementedError('Other QM codes not yet implemented!')

    return eigenvalues, fermi_level, kpts_array, coeffs #, Hd, Sd  

stuff = parser_read_input(args)



