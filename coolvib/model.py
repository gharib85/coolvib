"""
model.py

defines data class that holds all input
"""

import coolvib.parser.parser as parser
from coolvib.parser import codes
from ase.atoms import Atoms

class model():
    def __init__(self, atoms, code=None, active_atoms=None):
        """
        Base data container for all calculations
        """
        if isinstance(atoms, Atoms):
            self.atoms = atoms
        else:
            raise ValueError('model needs an instance of an ase Atoms class')
        self.code = code
        self.active_atoms = atoms_list

        self.masses = atoms.get_masses()
        self.

    def set_code(self, code):

        if code in codes:
            self.code = code
        else:
            print 'This code is not supported, choose one of ', codes

    def set_active_atoms(self, active_atoms=None):

        self.active_atoms = active_atoms

    def read_input_data(self, *args, **kwargs):

        if self.code is None:
            print 'Please set code first using .set_code'

        method_to_call = getattr(parser,codes[self.code])
        method_to_call(self, *args, **kwargs)

    def calculate_friction_tensor(self,):
        """
        calculates the friction tensor for a given window and delta function
        """

        calculate_spec = True
        if hasattr(self, 'spectral_function'):
            if self.spectral_function.emax >= emax:
                pass
            else:
                calcualte_spec = True
                print 'recalculating spectral function'
            calculate_spec = False
        
        #first calculate spectral function
        if calculate_spec:
            self.calculate_spectral_function(*args, **kwargs)

        ###calculate the tensor by applying the delta function
