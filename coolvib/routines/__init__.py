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
Routines __init__.py
"""

# import coolvib.routines.spectral_function as spectral_functions
# import coolvib.routines.friction_tensor as friction_tensor
from coolvib.constants import k_b

# calculate_spectral_function_tensor = spectral_functions.calculate_spectral_function_tensor
# calculate_spectral_function_tensor_q = spectral_functions.calculate_spectral_function_tensor_q
# calculate_tensor = friction_tensor.calculate_tensor

from math import exp, sin, sqrt, pi
import numpy as np

# from sys import float_info
# min_float = float_info.min

def fermi_occ(e, e0, T=300):
    """
    calculates occupation following Fermi-Dirac statistics
    """
    beta = 1.0/(T*k_b)
    try:
        return 1./(exp((e-e0)*beta)+1.)
    except OverflowError:
        return 0.0

#DIRAC delta approximations
def square(x,x0,s):
    if abs(x-x0)>s:
        return 0.0
    else:
        return 1./s

def gaussian(x,x0,s):
    try:
        return 1./sqrt(2*pi*s*s)*exp(-0.5*(((x-x0) / (s))**2))
    except OverflowError:
        return 0.0

def squashed_fermi(x,x0,s):
    y = (x/s)*sqrt(2./pi)+sqrt(0.5)
    try:
        return 2.*(y/s)*sqrt(2./pi)*exp(0.5-(y*y))
    except OverflowError:
        return 0.0

def lorentzian(x,x0,s):
    return (1.0/pi)*((0.5*s)/ \
        ((x-x0)*(x-x0)+(0.5*s)*(0.5*s)))

def sine(x,x0,s):
    t = 1./s
    if abs(x-x0)<0.00001:
        sine_val = t/pi
    else:
        sine_val = (sin(t*(x-x0)))**2/(x-x0)**2/t/pi
    return sine_val
    
def delta_function(x,x0,s, method):
    """
    function generator that returns a function, which 
    yields true or false if its argument is inside the given 
    window defined by minimum and maximum.
    """
    epsilon = 1E-100
    from math import isnan

    if method is 'gaussian':
        dirac_weight = gaussian(x,x0,s)
    elif method is 'square':
        dirac_weight = square(x,x0,s)
    elif method is 'sine':
        dirac_weight = sine(x,x0,s)
    elif method is 'lorentzian':
        dirac_weight = lorentzian(x,x0,s)
    elif method is 'sine_VSB':
        dirac_weight = sine_VSB(x,x0,s)
    elif method is 'squashed_fermi':
        dirac_weight = squashed_fermi(x,x0,s)
    else:
        raise NotImplementedError('delta method {0} is unknown'.format(method))

    if dirac_weight<epsilon:
        dirac_weight = 0.0
    return dirac_weight


def discretize_peak(e, nacs, x_axis, sigma, delta_method):
    """
    takes a delta peak with intensity and width and discretizes it 
    on a grid
    """
    de = x_axis[1] - x_axis[0]
    spectrum = np.zeros(len(x_axis),dtype=np.complex)

    norm = 0.0
    for i,x in enumerate(x_axis):
        delta = delta_function(x,e,sigma,delta_method)
        norm += delta
        spectrum[i] += delta*nacs
    norm *= de
     # if norm <= 1E-30:
        # return np.zeros(len(x_axis),dtype=np.complex)
     # else:
    return spectrum/norm


def evaluate_delta_function(x_axis, f, x0, sigma, delta_method):
    """
    evaluates the first moment of a function times a delta function
    """

    norm = 0.0
    result = 0.0
    for i, x in enumerate(x_axis):
        delta = delta_function(x,x0,sigma, delta_method)
        norm += delta
        result += f[i]*delta

    result/= norm

    return result

