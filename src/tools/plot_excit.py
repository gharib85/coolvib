from __future__ import division


import numpy as np
import scipy
#import scipy.optimize as optimize
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from lmfit import minimize, Parameters, Parameter, report_fit
import matplotlib.pyplot as plt

#reading the datapoints 

for index in range(1,13):
    x = np.loadtxt('test_scan_kpoints_{0}.txt'.format(index),usecols=range(0,1))
    y = np.loadtxt('test_scan_kpoints_{0}.txt'.format(index),usecols=range(1,2))
    plt.plot(x,y, label='{0} {1} 1'.format(index,index))
plt.xlim(0,0.6)
plt.ylim(0, 40)
plt.legend(loc='upper left')
plt.show()

