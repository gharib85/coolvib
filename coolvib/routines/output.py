#    This file is part of coolvib
#
#        coolvib is free software: you can redistribute it and/or modify
#        it under the terms of the GNU General Public License as published by
#        the Free Software Foundation, either version 3 of the License, or
#        (at your option) any later version.
#
#        coolvib is distributed in the hope that it will be useful,
#        but WITHOUT ANY WARRANTY; without even the implied warranty of
#        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#        GNU General Public License for more details.
#
#        You should have received a copy of the GNU General Public License
#        along with coolvib.  If not, see <http://www.gnu.org/licenses/>.
"""
output.py

contains functions for output
"""
import numpy as np
from coolvib.constants import time_to_ps

def print_spectral_function(x_axis, spectral_function, n_dim, filename='nacs-spectrum.out'):
    """
    This function outputs the spectral function(s) into a file.
    """

    n_axis = len(x_axis)
    de = x_axis[1]-x_axis[0]
    file = open(filename, mode='w')

    string = ''
    string += "No. of components {0} \n".format(n_dim)
    string += "Discretization length {0} \n".format(de)
    string += "Number of bins   {0} \n".format(n_axis)

    if n_dim==1:
        spectral_function = [spectral_function]

    file.write(string)

    counter = 0
    for i in range(n_dim):
        for j in range(i,n_dim):
            string = "Friction component for   {0}  {1} \n".format(i,j)
            string +="Excitation energy in eV   Coupling element in 1/eV*s\n ========================================== \n"
            file.write(string)
            for ei, e in enumerate(x_axis):
                file.write("{0:9.4f} {1:16.6E} {2:16.6E} \n".format(e, \
                        (spectral_function[counter, ei]).real/time_to_ps,
                        (spectral_function[counter, ei]).imag/time_to_ps))
            counter += 1

    file.close()


def print_jmol_friction_eigenvectors(atoms,active_atoms, 
        eigenvectors, eigenvalues, filename='friction_eigenvectors.jmol'):
    """
    prints friction eigenvectors into a jmol format file
    """

    file = open(filename, mode='w')

    n_atoms = len(atoms)

    for n in range(len(eigenvalues)):
        file.write("{0:6d}\n".format(n_atoms))
        file.write("Mode # {0:3d} f = {1:16.5F} ps \n".format(n,1./(eigenvalues[n]/time_to_ps)) )
        
        c = 0
        for a in range(n_atoms):
            symbol = atoms[a].symbol
            pos = atoms[a].position
            if a in active_atoms:
                file.write("{0} {1:12.6f} {2:12.6f} {3:12.6f} {4:12.6f} {5:12.6f} {6:12.6f} \n"\
                        .format(symbol,pos[0],pos[1],pos[2],\
                        # eigenvectors[n,c*3].real,eigenvectors[n,c*3+1].real,eigenvectors[n,c*3+2].real ))
                        np.sign(eigenvectors[n,c*3])*np.abs(eigenvectors[n,c*3]).real,\
                        np.sign(eigenvectors[n,c*3+1])*np.abs(eigenvectors[n,c*3+1]).real,\
                        np.sign(eigenvectors[n,c*3+2])*np.abs(eigenvectors[n,c*3+2]).real))
                c += 1
            else:
                file.write("{0} {1:12.6f} {2:12.6f} {3:12.6f} {4:12.6f} {5:12.6f} {6:12.6f} \n"\
                        .format(symbol,pos[0],pos[1],pos[2],\
                        0.,0.,0. ))


def print_matrix(matrix):
    """
    fancy printing of matrix
    """
    print '      '
    print '      '.join(["{0:14d}  ".format(i) for i in range(len(matrix[1])) ])
    for i, element in enumerate(matrix):
        print "{0:6d}".format(i),''.join(['{0:6.4e} {1:6.4e} '.format(y.real,y.imag) for y in element] ) 
    print '      '

def plot_spectral_function(x_axis, spectral_function):
    """
    make plot of spectral function
    """
    try:
        import matplotlib.pyplot as plt
    except:
        raise ImportError('Could not find matplotlib')

    plt.plot(x_axis, spectral_function.real/time_to_ps, x_axis, spectral_function.imag/time_to_ps)

    plt.show()
    
