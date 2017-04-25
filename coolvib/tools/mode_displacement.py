import os
from ase.utils import opencew
from os import remove
from os.path import isfile, getsize
import numpy as np

def mode_displacement(atoms, mode, name='mode', disp=0.001):
    """
    This class is based on the ASE Vibrations class and only 
    slightly modified to print Hamiltonian and overlap matrix output into 
    specific folders. With this data coolvib can calculate the finite difference 
    nonadiabatic coupling elements.
    """

    filename = name + '.eq.pckl'
    fd = opencew(filename)
    if fd is not None:
        e0 = atoms.get_potential_energy()
        # print e0
        fd.write('calculated')
        try:
            os.mkdir(name+'_eq')
        except:
            pass
        os.system('mv *.out '+name+'_eq/')

    p = atoms.positions.copy()
    m = atoms.get_masses() 
    for i,a in enumerate(mode):
        mode[i] = a/np.sqrt(m[i])
    #print mode

    for sign in [-1,1]:
        filename = ('%s.%s.pckl' % (name, ' +-'[sign]))
        filename2 = ('%s_%s' %
                    (name, ' +-'[sign]))
        if (isfile(filename) and getsize(filename)==0):
            remove(filename)
        fd = opencew(filename)
        if fd is not None:
            displ = sign * mode * disp #*np.sqrt(m)
            atoms.positions = p + displ
            #print atoms.positions
            atoms.calc.calculate(atoms)
            fd.write('calculated')
            try:
                os.mkdir(filename2)
            except:
                pass
            os.system('mv *.out '+filename2+'/')
            atoms.positions = p

