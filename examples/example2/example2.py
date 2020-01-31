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
example 2: 
Calculating lifetime along a given mode displacement for 
CO on Pd(111) internal stretch calculated with FHI-Aims
"""
import numpy as np
from ase.io import read
import os, sys

import coolvib
from scipy import linalg as LA

#####DEFINE SYSTEM######

system = read('mode/mode_eq/geometry.in')
cell = system

model = coolvib.workflow_mode(system, code='aims')

finite_difference_incr = 0.0025

keywords = {
    'discretization_type' : 'gaussian',
    'discretization_broadening' : 0.01,
    'discretization_length' : 0.01,
    'max_energy' : 7,
    'temperature' : 300,
    'delta_function_type': 'gaussian',
    'delta_function_width': 0.60,
    'perturbing_energy' : 0.0,
        }

print('workflow initialized and keywords set')

######READ QM INPUT DATA###

model.read_input_data(spin=False, 
        path='mode/', 
        prefix='mode', 
        filename='aims.out', 
        incr=finite_difference_incr)
print('successfully read QM input data')

# model.calculate_spectral_function(mode='default', **keywords)
# print('successfully calculated spectral_function')
# model.print_spectral_function('nacs-spectrum.out')
# model.calculate_friction_from_spectrum(**keywords)

print('successfully calculated friction tensor')
model.calculate_friction(**keywords)

model.analyse_friction()




