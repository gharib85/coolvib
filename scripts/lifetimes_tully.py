#!/usr/bin/python
from ase import *
from ase.lattice.surface import *
from ase.visualize import view
from string import Template
import os, sys
import numpy as np
from numpy import *
from ase.calculators.siesta import Siesta
from ase.data import chemical_symbols
from ase.units import Rydberg, fs
from ase.io.siesta import read_rho, read_fdf, read_struct
from ase.io.cube import read_cube_data
import os
from os.path import join, isfile, islink, getmtime
from cmath import exp
import array
from  shutil import *
#Script is written for a slab containing only one type of metal
#Things you can change:
#1) coordinates and type of slab and adsorbate: lines 23-32
#2) Monkhorst-Pack grid: lines 55-57

#Create Cu slab
slab = fcc111('Cu',
              a=3.6149,       # Cu lattice constant
              size=(1,1,3),   #3-layer slab in 1x1 configuration
              vacuum=10.0)

#create a CO molecule
adsorbate= Atoms([Atom('C',[0., 0., 0.]),
           Atom('O',[0., 0., 1.15])])

add_adsorbate(slab,adsorbate,height=1.9,position='ontop')

adsorbate_indices = slab.get_tags().tolist().count(0)
#print adsorbate_indices

file_name = 'displacement'

#definition of calculator = SIESTA

calc = Siesta(
        xc = "PBE",
        width = 0.1,
        kpts = [2,2,1],
        label = "siesta",
        write_fdf= True)
calc.set_fdf("WriteWaveFunctions", True)
calc.set_fdf("WriteDenchar", True)
calc.set_fdf("WriteWaveFunctions", True)
calc.set_fdf("WriteEigenvalues", True)

#setting up the grid"
grid = np.int8(np.zeros((3,3)))
grid[0,0]= 1
grid[1,1]= 1
grid[2,2]= 1



slab.set_calculator(calc)

#Getting the matrix of chemical species
sym_ind = np.unique(slab.get_atomic_numbers(),return_index=True)[1]
sym = [slab.get_chemical_symbols()[i] for i in sort(sym_ind)]
num_ind = np.unique(slab.get_atomic_numbers(),return_index=True)[1]
num = [slab.get_atomic_numbers()[i] for i in sort(num_ind)]
ind = np.asarray(range(1,len(np.unique(slab.get_chemical_symbols()))+1))
block_ChemicalSpeciesLabel = vstack((ind,num,sym)).T


'''
#Getting the matrix of coordinates with metal atoms frozen

slab_xyz=slab.get_positions()
j=1
k= [[0,0,0]]
cartesians = [[]]
for i in range(1, len(slab_xyz)+1):
    if i <= len(slab)-len(adsorbate): 
        np.append(np.append(slab_xyz[i-1].tolist(),j), [[0,0,0]])  
        #print np.append(np.append(slab_xyz[i-1].tolist(),j), [[0,0,0]])  #np.append(slab_xyz[i-1],np.array(j),np.array(0 ,0, 0))
    else:
        j+=1
        print np.append(slab_xyz[i-1].tolist(),j)
        #print np.append(slab_xyz[i-1].tolist(),j)
print cartesians

#print "\n".join("   ".join(str(element) for element in row) for row in slab.get_positions())

'''



###### vibdispbuilder
path = os.getcwd()
#print path

dr= 0.0001 	# value of displacement for numerical differences
a=0  		# counter for the loop


#Create a dictionary to be written into a template. This can be easily spread over to other parameters

d={ "number_of_atoms": len(slab.get_positions()),\
    "number_of_species": len(np.unique(slab.get_chemical_symbols())),\
    "block_lattice_vectors":\
    "\n".join("   ".join(str(element) for element in row) for row in slab.get_cell()),\
    "block_kgrid_Monkhorst_Pack":\
    "  0.0\n".join("   ".join(str(element) for element in row) for row in grid)+"  0.0",\
    "number_constrained": len(slab)-len(adsorbate),\
    "block_ChemicalSpeciesLabel":\
    "\n".join("   ".join(str(element) for element in row) for row in block_ChemicalSpeciesLabel),\
    "block_AtomicCoordinatesAndAtomicSpecies":\
    "\n".join("   ".join(str(element) for element in row) for row in slab.get_positions())}


#Reading the template file
#equilibrium geometry calculaton

directory = path + "/"+file_name+"_eq"
if os.path.exists(directory):
    pass
else:
    os.mkdir(directory)
file_path = directory + "/input.fdf"
temp=open("template.fdf", "r+")
src = Template(temp.read())
with open(file_path, "w+") as f:
    f.write(src.substitute(d))
with open(file_path, "r+") as f:
    content=f.readlines()
    k=content.index('%block AtomicCoordinatesAndAtomicSpecies\n')
    for i in range(k+1, k+len(slab)-len(adsorbate)+1):
        content[i] = content[i][0:-1]+ '  1 0 0 0  \n'
    m=2
    for i in range(k+len(slab)-len(adsorbate)+1,k+len(slab)+1):
        content[i] = content[i][0:-1]+ ' '+ str(m)+'\n'
        m+=1
with open(file_path, "w+") as f:
    for line in content:
        f.write(line)
copyfile("submit.pbs", directory +"/submit.pbs")

for i in np.unique(slab.get_chemical_symbols()):  #getting the pseudopotentials
    copyfile(i+ ".psf", directory +"/"+i+".psf") 

cartesians = slab.get_positions()
#displacements




for index in range(len(adsorbate)):
    a=0
    for xyz in "xyz":
        cartesians[index+len(slab)-len(adsorbate),a] += dr
        directory = path + "/displacement_{0}_{1}".format(adsorbate.get_chemical_symbols()[index], xyz)
        if os.path.exists(directory):
            pass
        else:
            os.mkdir(directory)
        file_path = directory + "/input.fdf"
        copyfile("submit.pbs", directory +"/submit.pbs") #getting the submit script
        for i in np.unique(slab.get_chemical_symbols()):  #getting the pseudopotentials
            copyfile(i+ ".psf", directory +"/"+i+".psf") 
        temp=open("template.fdf", "r+")
        src = Template(temp.read())
        with open(file_path, "w+") as f:
            d={ "number_of_atoms": len(slab.get_positions()),\
                "number_of_species": len(np.unique(slab.get_chemical_symbols())),\
                "block_lattice_vectors":\
                "\n".join("   ".join(str(element) for element in row) for row in slab.get_cell()),\
                "block_kgrid_Monkhorst_Pack":\
                "  0.0\n".join("   ".join(str(element) for element in row) for row in grid) +"  0.0",\
                "number_constrained": len(slab)-len(adsorbate),\
                "block_ChemicalSpeciesLabel":\
                "\n".join("   ".join(str(element) for element in row) for row in block_ChemicalSpeciesLabel),\
                "block_AtomicCoordinatesAndAtomicSpecies":\
                "\n".join("   ".join(str(element) for element in row) for row in cartesians)}
            f.write(src.substitute(d))
        with open(file_path, "r+") as f:
            content=f.readlines()
            k=content.index('%block AtomicCoordinatesAndAtomicSpecies\n')
            for i in range(k+1, k+len(slab)-len(adsorbate)+1):
                content[i] = content[i][0:-1]+ '  1 0 0 0  \n'
            m=2
            for i in range(k+len(slab)-len(adsorbate)+1,k+len(slab)+1):
                content[i] = content[i][0:-1]+ ' '+ str(m)+'\n'
                m+=1
        with open(file_path, "w+") as f:
            for line in content:
                f.write(line)
        cartesians[index+len(slab)-len(adsorbate),a] -= dr
        a += 1




