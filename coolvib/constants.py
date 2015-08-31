"""
This file contains all conversion factors and constants

coolvib works with the same units as the 
Atomic Simulation Environment

energy : eV
length : Angstrom
mass : a.m.u.
force : eV/Angstorm
time : Angstrom*sqrt(amu/e)
action, angular momentum : sqrt(e*amu)*Angstrom

"""

"""
Hz          cm-1    
1.000000    3.335641E-11
eV          J   
1.000000    1.602177E-19  
Na      
6.022141E+23    
Hz          sqrt(ev)/Angst  
1.000000    1.018051E-14    
sqrt(eV)/Angst  Hz                  cm-1
1.000000        98226949774380.300000   3.276498E+03
"""

from math import pi

conversion_factor = 98226935315503.17 #sqrt(amu/e)*Ang to seconds value taken from ASE, CODATA02
#conversion_factor = 98226949774380.3 #for omega from sqrt(eV/amu)/Ang to Hz
time_to_ps = 1./(conversion_factor*1E-12)
hplanck = 0.40623411488 # in units of sqrt(e*amu)*Angstrom
hbar = 0.064654167 #in unit of sqrt(e*amu)*Angstrom
#hbar=6.58211928E-16 #eV*s
hbar_1=3.990313E+13  #amu*ang2/s

#hbar_1=3.990313E+13/2./pi  #amu*ang2/s

sqrt_pi = 1.772453851

sqrt_2pi = 2.506628275 

#ryd = 13.605698066 #eV old value
ryd = 13.6056978277587 #eV from ASE, CODATA02

#k_b = 8.6173324E-5 #eV / K old
k_b = 8.61738569226E-5 #eV / K from ASE, CODATA02

#bohr = 0.529177249 #Ang old
bohr = 0.5291772575 #Ang from ASE, CODATA02

# me = 1822.88839 # mass of the atom relative to electron old
me = 1822.88853006 # mass of the atom relative to electron from ASE, CODATA02

