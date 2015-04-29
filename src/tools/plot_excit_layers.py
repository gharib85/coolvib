from __future__ import division


import numpy as np
import scipy
#import scipy.optimize as optimize
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from lmfit import minimize, Parameters, Parameter, report_fit
import matplotlib.pyplot as plt

#reading the datapoints 

for index in range(2,8):
    x = np.loadtxt('excit_scan_layers_{0}.txt'.format(index),usecols=range(0,1))
    y = np.loadtxt('excit_scan_layers_{0}.txt'.format(index),usecols=range(1,2))
    plt.plot(x,y, label='{0} layers'.format(index))
plt.xlim(0,1)
plt.ylim(0,4)
plt.legend(loc='upper left')
plt.show()

for index in range(2,8):
    x = np.loadtxt('dos_scan_layers_{0}.txt'.format(index),usecols=range(0,1))
    y = np.loadtxt('dos_scan_layers_{0}.txt'.format(index),usecols=range(1,2))
    plt.plot(x,y, label='{0} layers'.format(index))
plt.xlim(-1,1)
plt.ylim(0, 4)
plt.legend(loc='upper left')
plt.show()
