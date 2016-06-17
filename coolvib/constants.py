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
This file contains all conversion factors and constants

coolvib works with the same units as the 
`Atomic Simulation Environment <https://wiki.fysik.dtu.dk/ase/>`_

energy : eV
length : Angstrom
mass : a.m.u.
force : eV/Angstorm
time : Angstrom*sqrt(amu/e)
action, angular momentum : sqrt(e*amu)*Angstrom


+--------------+------------------------+---------------------+
|Hz            |   cm-1                 |                     |
+--------------+------------------------+---------------------+
|1.000000      |   3.335641E-11         |                     |
+--------------+------------------------+---------------------+
|eV            |   J                    |                     |
+--------------+------------------------+---------------------+
|1.000000      |       1.602177E-19     |                     |
+--------------+------------------------+---------------------+
|Hz            |   sqrt(ev)/Angst       |                     |
+--------------+------------------------+---------------------+
|1.000000      |   1.018051E-14         |                     |
+--------------+------------------------+---------------------+
|sqrt(eV)/Ang  |  Hz                    |   cm-1              |
+--------------+------------------------+---------------------+
|1.000000      |  98226949774380.300000 |  3.276498E+03       |
+--------------+------------------------+---------------------+
"""

from math import pi

conversion_factor = 98226935315503.17 #sqrt(amu/e)*Ang to seconds value taken from ASE, CODATA02
#conversion_factor = 98226949774380.3 #for omega from sqrt(eV/amu)/Ang to Hz
time_to_ps = 1./(conversion_factor*1E-12)
hplanck = 0.40623411488 # in units of sqrt(e*amu)*Angstrom
hbar = 0.064654167 #in unit of sqrt(e*amu)*Angstrom

hbar_1=3.990313E+13/2./pi  #amu*ang2/s

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

