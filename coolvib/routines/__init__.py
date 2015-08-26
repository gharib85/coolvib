"""
Routines __init__.py
"""


def gaussian(x, x_mean, broadening): 
    gaussian_val = 1/np.sqrt(2*pi*broadening*broadening)*np.exp(-0.5*(((x-x_mean) / (broadening))**2));
    return gaussian_val

def square(broadening): 
    square_val = 1/broadening;
    return square_val

def sine(x, x_mean, broadening): 
    if np.abs(x-x_mean)<0.00001:
        sine_val = 1/broadening/pi
    else:
        sine_val = (np.sin(broadening*(x-x_mean)))**2/(x-x_mean)**2/broadening/pi;
    return sine_val

def sine_VSB(x, x_mean, broadening): 
    if np.abs(x-x_mean)<0.00001:
        sine_VSB_val = 1/2/pi*broadening
    else: 
        x_1 = (x-x_mean)*(2/broadening)    #t/2h_bar i.e. t=4h_bar/broadening
        sine_VSB_val=2/pi/broadening*(np.sin(x_1))**2/x_1**2
    return sine_VSB_val

def lorentzian(x, x_mean, broadening): 
    lorentzian_val = 1/(1+((x-x_mean)/broadening)**2)/broadening/pi;
    return lorentzian_val

def generate_delta_approximation(dirac_method, fermi_centered, delta, window):
    """
    function generator that returns a function, which 
    yields true or false if its argument is inside the given 
    window defined by minimum and maximum.
    """
    def dirac_weight_window(exc):
        #if fermi_centered:
            #return (exc*exc)/(window*window)
        #else:
            return 1./window
    def dirac_weight_gaussian(exc):
        gaussian = exp(-((exc-delta)*(exc-delta))/(2*window*window))/(window*sqrt_2pi)
        return gaussian 
    
    if dirac_method is 'window':
        win = window
        dirac_weight = dirac_weight_window
    if dirac_method is 'gaussian':
        win = window*6.
        dirac_weight = dirac_weight_gaussian
    
    if fermi_centered:
        minimum=0
        maximum=win
    else:
        minimum=delta-win/2.
        maximum=delta+win/2.
        if minimum<0.:
            maximum = maximum-minimum
            minimum=0.
    
    def e_decider(exc):
        if minimum < exc < maximum:
            return True
        else:
            return False
    
    return e_decider, dirac_weight
