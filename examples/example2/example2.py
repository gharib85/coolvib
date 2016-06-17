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
example 1: 
Calculating 6x6 friction tensor for 
CO on Cu(100) calculated with FHI-Aims
"""
import numpy as np
from ase.all import *
import os, sys

import coolvib

from scipy import linalg as LA

#####DEFINE SYSTEM######


system = read('CO-on-Cu100/eq/geometry.in')
cell = system

active_atoms = [3,4] # meaning we have two atoms - C - 0 and O - 1

model = coolvib.workflow_tensor(system, code='aims', active_atoms=active_atoms)
print model.atoms.get_masses()[active_atoms]

finite_difference_incr = 0.001

keywords = {
    'discretization_type' : 'gaussian',
    'discretization_broadening' : 0.05,
    'discretization_length' : 0.01,
    'max_energy' : 6.00,
    'temperature' : 300,
    'delta_function_type': 'gaussian',
    'delta_function_width': 0.60,
    'perturbing_energy' : 0.0,
        }

print 'workflow initialized and keywords set'

######READ QM INPUT DATA###

model.read_input_data(spin=True, path='./CO-on-Cu100', filename='OUTPUT', active_atoms=active_atoms, incr=finite_difference_incr)
print 'successfully read QM input data'


model.calculate_spectral_function(mode='default', **keywords)
print 'successfully calculated spectral_function'

model.print_spectral_function('nacs-spectrum.out')

model.calculate_friction_tensor(**keywords)
print 'successfully calculated friction tensor'

model.analyse_friction_tensor()




