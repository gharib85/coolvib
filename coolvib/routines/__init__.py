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
        return 1./(exp((e-e0)*beta)+1)
    except OverflowError:
        return 0.0

#DIRAC delta approximations
def square(x,x0,s):
    if abs(x-x0)>s:
        return 0.0
    else:
        return 1./s

def gaussian(x,x0,s):
    return 1./sqrt(2*pi*s*s)*exp(-0.5*(((x-x0) / (s))**2))

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
    if abs(x-x_mean)<0.00001:
        sine_val = 1./s/pi
    else:
        sine_val = (sin(s*(x-x0)))**2/(x-x0)**2/s/pi
    return sine_val
    
def sine_VSB(x, x0, s): 
    if abs(x-x0)<0.00001:
        sine_VSB_val = 1./2./pi*s
    else: 
        x_1 = (x-x0)*(2/s)    #t/2h_bar i.e. t=4h_bar/s
        sine_VSB_val=2./pi/s*(sin(x_1))**2./x_1**2.
    return sine_VSB_val


def delta_function(x,x0,s, method):
    """
    function generator that returns a function, which 
    yields true or false if its argument is inside the given 
    window defined by minimum and maximum.
    """
    
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
        result += f[i]*x*delta

    result /= norm

    return result

