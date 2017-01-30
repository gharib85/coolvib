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
from coolvib.routines.output import print_matrix

input1='start.in'
ab=read(input1)

#CALCULATOR
basisset = 'light'

#setting up Aims calculator in ASE
calc=Aims(xc='PBE',
    # run_command = 'mpirun -np 4 <insert_path>/bin/aims.160210.mpi.x > aims.out',
    run_command = 'mpirun -np 4 $HOME/Work/codes/QM-codes/FHI-aims/FHI-AIMS_Current/bin/aims.170122.mpi.x > aims.out',
    # species_dir='<insert path>/aimsfiles/species_defaults/'+basisset',
    species_dir='/home/reini/Work/codes/QM-codes/FHI-aims/FHI-AIMS_Current/aimsfiles/species_defaults/'+basisset,
    occupation_type = ['gaussian',0.1],
    sc_iter_limit = 100,
    #spin = 'collinear',
    relativistic = ['atomic_zora','scalar'],
    #default_initial_moment = 0,
    sc_accuracy_etot=1e-6,
    sc_accuracy_eev=0.001,
    sc_accuracy_rho=1e-5,
    sc_accuracy_forces=1e-3,
    load_balancing = True,
    empty_states=10,
    k_grid = [8,8,1],
    restart_aims='wvfn.dat',
   output = ["eigenvectors",  "k_point_list",  "hamiltonian_matrix",  "overlap_matrix",  "band     0.0000    0.0000    0.0000    0.0000    0.1250    0.0000  2 ",  "band     0.0000    0.2500    0.0000    0.0000    0.3750    0.0000  2 ",  "band     0.0000    0.5000    0.0000    0.1250    0.0000    0.0000  2 ",  "band     0.1250    0.1250    0.0000    0.1250    0.2500    0.0000  2 ",  "band     0.1250    0.5000    0.0000    0.1250    0.6250    0.0000  2 ",  "band     0.1250    0.7500    0.0000    0.2500    0.0000    0.0000  2 ",  "band     0.2500    0.1250    0.0000    0.2500    0.2500    0.0000  2 ",  "band     0.2500    0.5000    0.0000    0.2500    0.6250    0.0000  2 ",  "band     0.2500    0.7500    0.0000    0.3750    0.0000    0.0000  2 ",  "band     0.3750    0.1250    0.0000    0.3750    0.2500    0.0000  2 ",  "band     0.3750    0.5000    0.0000    0.3750    0.6250    0.0000  2 ",  "band     0.3750    0.7500    0.0000    0.5000    0.0000    0.0000  2 ",  "band     0.5000    0.1250    0.0000    0.5000    0.2500    0.0000  2 ",  "band     0.5000    0.5000    0.0000    0.5000    0.6250    0.0000  2 ",  "band     0.6250    0.1250    0.0000    0.6250    0.6250    0.0000  2 ",  "band     0.7500    0.1250    0.0000    0.7500    0.6250    0.0000  2 ",  "band     0.8750    0.1250    0.0000    0.8750    0.6250    0.0000  2 ",  ], 
    )

ab.set_calculator(calc)

indices1=[3,4]
f = finite_difference(ab,indices=indices1, delta=0.0025)

f.run()
f.summary()
f.write_jmol()

#printing the Hessian
print 'Vibrational hessian just for reference'
print_matrix(f.modes)

# for i in range(len(f.hnu)):
    # string = ''
    # mode = f.get_mode(i)[indices1]
    # for j in range(len(vib.hnu)):
        # string += ' {0:14.8f} '.format(modes[i,j])
    # print string

######
# Now we have generated all input files necessary to run example1
#####


#Calculate finite difference along a given mode

print ' Now calculating finite difference along a given normal mode'
#follow the internal stretch mode

#You can define your own mode displacement or use one calculated from 
#a Hessian calculation. It is important to use displacements in 
#mass-weighted cartesian coordinates, 
#that means, e*sqrt(mass)

# mode = np.array([
    # [0.,0.,0.],
    # [0.,0.,0.],
    # [0.,0.,0.],
    # [0.,0.,-1.0],
    # [0.,0.,1.0],
     # ])

#We pick the right mode from the hessian
#in this case the internal stretch mode
mode = f.get_mode(-1)
#we normalize
mode /= np.linalg.norm(mode)
print mode
#we mass-weight
# print np.diag(ab.get_masses())
# w = np.diag(np.sqrt(ab.get_masses()))
# print w
# mode = np.dot(mode,w)
# print mode

mode_displacement(ab, mode, name='mode', disp=0.0025)


