"""
friction_tensor.py
"""

from coolvib.routines import evaluate_delta_function
from coolvib.routines.output import print_matrix
from coolvib.constants import ps

import numpy as np
import numpy.linalg as la

def calculate_tensor(n_dim, x_axis, spectral_function, **kwargs):
    """
    calculates friction tensor from given spectral function
    """
    
    ##ANALYSE KWARGS

    default_keywords = {
            'delta_function_type' : 'gaussian',
            'delta_function_width' : 0.60,
            'max_energy' : 3.0,
            'perturbing_energy' : 0.0,  #centered around fermi
            }

    keys = {}
    for key in default_keywords.keys():
        if hasattr(kwargs,key):
            keys[key] = kwargs[key] 
        else:
            keys[key] = default_keywords[key]

    sigma = keys['delta_function_width']
    delta_method = keys['delta_function_type']
    perturbing_energy = keys['perturbing_energy']
    friction_tensor = np.zeros([n_dim,n_dim],dtype=np.complex)
    
    c = 0
    for d in range(n_dim):
        for d2 in range(d,n_dim):
            
            friction_tensor[d,d2] = \
                    evaluate_delta_function(x_axis, spectral_function[c], perturbing_energy ,sigma, delta_method)
            friction_tensor[d2,d] = friction_tensor[d,d2].conjugate()

            c += 1


    return friction_tensor 


def analyse_tensor(friction_tensor):
    """
    performs a spectral analysis of friction tensor
    """

    eigenvalues, eigenvectors = la.eigh(friction_tensor)

    #printing
    print 'Friction Tensor' 
    print_matrix(friction_tensor/ps)
    print ' '
    print 'Friction Eigenvalues in 1/ps'
    print ' '.join( ['{0:10d}'.format(y) for y in range(len(eigenvalues)) ])
    print ' '.join( ['{0:10.4f}'.format(y) for y in eigenvalues/ps ]  )
    print 'Friction Eigenvectors'
    print_matrix(eigenvectors)
    print 'Principal lifetimes in ps'
    print ' '.join( ['{0:10d}'.format(y) for y in range(len(eigenvalues)) ])
    print ' '.join( ['{0:10.4f}'.format(y) for y in ps/eigenvalues ]  )

    return eigenvectors, eigenvalues
    
    #do stuff
