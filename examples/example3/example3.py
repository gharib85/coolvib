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
We perform a finite-difference calculation using ASE and FHI-aims
to generate matrix elements for coolvib"""

import numpy as np
import os, sys
from ase.io import read
from ase.calculators.aims import Aims
import coolvib
from coolvib.tools.finite_difference import finite_difference
from coolvib.tools.mode_displacement import mode_displacement

input1='start.in'
ab=read(input1)

#CALCULATOR
calc=Aims(xc='PBE',
    command = 'mpirun -np 4 <insert_path>/bin/aims.160210.mpi.x > aims.out',
    species_dir='<insert path>/aimsfiles/species_defaults/light',
    occupation_type = ['gaussian',0.1],
    sc_iter_limit = 100,
    #spin = 'collinear',
    relativistic = ['atomic_zora','scalar'],
    #default_initial_moment = 0,
    sc_accuracy_etot=1e-6,
    sc_accuracy_eev=1e-3,
    sc_accuracy_rho=1e-6,
    sc_accuracy_forces=1e-4,
    load_balancing = True,
    k_grid = [16,16,1],
    restart_aims='wvfn.dat',
    output = ["eigenvectors", "k_point_list", "hamiltonian_matrix",  "overlap_matrix",  "band     0.0000    0.0625    0.0000    0.0000    0.0625    0.0000  2 ",  "band     0.0000    0.1875    0.0000    0.0000    0.1875    0.0000  2 ",  "band     0.0000    0.3750    0.0000    0.0000    0.3750    0.0000  2 ",  "band     0.0000    0.5000    0.0000    0.0000    0.5000    0.0000  2 ",  "band     0.0625    0.0000    0.0000    0.0625    0.0000    0.0000  2 ",  "band     0.0625    0.1250    0.0000    0.0625    0.1250    0.0000  2 ",  "band     0.0625    0.3125    0.0000    0.0625    0.3125    0.0000  2 ",  "band     0.0625    0.5000    0.0000    0.0625    0.5000    0.0000  2 ",  "band     0.0625    0.6250    0.0000    0.0625    0.6250    0.0000  2 ",  "band     0.0625    0.8125    0.0000    0.0625    0.8125    0.0000  2 ",  "band     0.1250    0.0000    0.0000    0.1250    0.0000    0.0000  2 ",  "band     0.1250    0.1250    0.0000    0.1250    0.1250    0.0000  2 ",  "band     0.1250    0.3125    0.0000    0.1250    0.3125    0.0000  2 ",  "band     0.1250    0.5000    0.0000    0.1250    0.5000    0.0000  2 ",  "band     0.1250    0.6250    0.0000    0.1250    0.6250    0.0000  2 ",  "band     0.1250    0.8125    0.0000    0.1250    0.8125    0.0000  2 ",  "band     0.1875    0.0000    0.0000    0.1875    0.0000    0.0000  2 ",  "band     0.1875    0.1250    0.0000    0.1875    0.1250    0.0000  2 ",  "band     0.1875    0.3125    0.0000    0.1875    0.3125    0.0000  2 ",  "band     0.1875    0.5000    0.0000    0.1875    0.5000    0.0000  2 ",  "band     0.1875    0.6250    0.0000    0.1875    0.6250    0.0000  2 ",  "band     0.1875    0.8125    0.0000    0.1875    0.8125    0.0000  2 ",  "band     0.2500    0.0000    0.0000    0.2500    0.0000    0.0000  2 ",  "band     0.2500    0.1250    0.0000    0.2500    0.1250    0.0000  2 ",  "band     0.2500    0.3125    0.0000    0.2500    0.3125    0.0000  2 ",  "band     0.2500    0.5000    0.0000    0.2500    0.5000    0.0000  2 ",  "band     0.2500    0.6250    0.0000    0.2500    0.6250    0.0000  2 ",  "band     0.2500    0.8125    0.0000    0.2500    0.8125    0.0000  2 ",  "band     0.3125    0.0000    0.0000    0.3125    0.0000    0.0000  2 ",  "band     0.3125    0.1250    0.0000    0.3125    0.1250    0.0000  2 ",  "band     0.3125    0.3125    0.0000    0.3125    0.3125    0.0000  2 ",  "band     0.3125    0.5000    0.0000    0.3125    0.5000    0.0000  2 ",  "band     0.3125    0.6250    0.0000    0.3125    0.6250    0.0000  2 ",  "band     0.3125    0.8125    0.0000    0.3125    0.8125    0.0000  2 ",  "band     0.3750    0.0000    0.0000    0.3750    0.0000    0.0000  2 ",  "band     0.3750    0.1250    0.0000    0.3750    0.1250    0.0000  2 ",  "band     0.3750    0.3125    0.0000    0.3750    0.3125    0.0000  2 ",  "band     0.3750    0.5000    0.0000    0.3750    0.5000    0.0000  2 ",  "band     0.3750    0.6250    0.0000    0.3750    0.6250    0.0000  2 ",  "band     0.3750    0.8125    0.0000    0.3750    0.8125    0.0000  2 ",  "band     0.4375    0.0000    0.0000    0.4375    0.0000    0.0000  2 ",  "band     0.4375    0.1250    0.0000    0.4375    0.1250    0.0000  2 ",  "band     0.4375    0.3125    0.0000    0.4375    0.3125    0.0000  2 ",  "band     0.4375    0.5000    0.0000    0.4375    0.5000    0.0000  2 ",  "band     0.4375    0.6250    0.0000    0.4375    0.6250    0.0000  2 ",  "band     0.4375    0.8125    0.0000    0.4375    0.8125    0.0000  2 ",  "band     0.5000    0.0000    0.0000    0.5000    0.0000    0.0000  2 ",  "band     0.5000    0.1250    0.0000    0.5000    0.1250    0.0000  2 ",  "band     0.5000    0.3125    0.0000    0.5000    0.3125    0.0000  2 ",  "band     0.5000    0.5000    0.0000    0.5000    0.5000    0.0000  2 ",  "band     0.5000    0.8125    0.0000    0.5000    0.8125    0.0000  2 ",  "band     0.5625    0.3125    0.0000    0.5625    0.3125    0.0000  2 ",  "band     0.5625    0.8125    0.0000    0.5625    0.8125    0.0000  2 ",  "band     0.6250    0.3125    0.0000    0.6250    0.3125    0.0000  2 ",  "band     0.6250    0.8125    0.0000    0.6250    0.8125    0.0000  2 ",  "band     0.6875    0.3125    0.0000    0.6875    0.3125    0.0000  2 ",  "band     0.6875    0.8125    0.0000    0.6875    0.8125    0.0000  2 ",  "band     0.7500    0.3125    0.0000    0.7500    0.3125    0.0000  2 ",  "band     0.7500    0.8125    0.0000    0.7500    0.8125    0.0000  2 ",  "band     0.8125    0.3125    0.0000    0.8125    0.3125    0.0000  2 ",  "band     0.8125    0.8125    0.0000    0.8125    0.8125    0.0000  2 ",  "band     0.8750    0.3125    0.0000    0.8750    0.3125    0.0000  2 ",  "band     0.8750    0.8125    0.0000    0.8750    0.8125    0.0000  2 ",  "band     0.9375    0.3125    0.0000    0.9375    0.3125    0.0000  2 ",  "band     0.9375    0.8125    0.0000    0.9375    0.8125    0.0000  2 ",  ], 
    )

ab.set_calculator(calc)

indices1=[1,2]
f = finite_difference(ab,indices=indices1)

f.run()
f.summary()
f.write_jmol()

#Internal stretch mode
mode = f.get_mode(-1)

print mode 

#FOLLOW A STRETCH MODE

mode = np.array([
    [0.,0.,0.],
    [0.,0.,-1.0],
    [0.,0.,1.0],
     ])
mode /= np.linalg.norm(mode)

mode_displacement(ab, mode, name='mode', disp=0.01)


