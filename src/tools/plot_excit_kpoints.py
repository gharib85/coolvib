from __future__ import division
import numpy as np
import scipy
#import scipy.optimize as optimize
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from lmfit import minimize, Parameters, Parameter, report_fit
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as colors


num_plots = 30
#reading the datapoints 

#print colors.cnames
colormap = plt.cm.gist_ncar
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])

for index in range(1,31):
    x = np.loadtxt('excit_scan_kpoints_{0}.txt'.format(index),usecols=range(0,1))
    y = np.loadtxt('excit_scan_kpoints_{0}.txt'.format(index),usecols=range(1,2))
    plt.plot(x,y, label='{0} {0} 1'.format(index,index))
plt.xlim(0,0.5)
plt.ylim(0,0.5)
plt.legend(loc='upper left')
plt.show()

colormap = plt.cm.gist_ncar
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])

for index in range(1,31):
    x = np.loadtxt('dos_scan_kpoints_{0}.txt'.format(index),usecols=range(0,1))
    y = np.loadtxt('dos_scan_kpoints_{0}.txt'.format(index),usecols=range(1,2))
    plt.plot(x,y, label='{0} {0} 1'.format(index,index))
plt.xlim(-0.5,0.5)
plt.ylim(0, 2)
plt.legend(loc='upper left')
plt.show()
