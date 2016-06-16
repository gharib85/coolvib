"""
workflow_mode.py

Defines a class **workflow_mode** that acts as a data container and workflow 
guideline to calculate the non-adiabatic friction induced lifetime 
of a given vibrational mode.
"""
import coolvib
import coolvib.parser as parser
import coolvib.routines.spectral_function as spectral
import coolvib.routines.output as output
from ase.atoms import Atoms
import numpy as np

class workflow_mode():
    def __init__(self, atoms, code='aims', mode=None):
        """
        workflow_mode object.
        
        The workflow_tensor objet is a basic data container and workflow guideline 
        for all calculations necessary to generate a friction tensor and 
        calculate vibrational lifetimes due to non-adiabatic couplings. 

        Parameters:

        atoms: ase.atoms class 
            defines the atoms, molecules and or surfaces for 
            which the calculation is performed.

        code: str 
            A string that defines from which quantum chemistry package the input 
            is received.

        mode: list
            A list or np.array of displacements that defines the mode along which the 
            lifetime is to be calculated.
           
        Examples:

            This is how this class is initialized::

                from ase.all import Atoms
                a = Atoms("CO", [(0.,0.,0.), (0.,0.,1.0)])
                model = workflow_mode(a, code='siesta', mode=[0.,0.,1.0, 0.,0.,-1.0])

        """
        if isinstance(atoms, Atoms):
            self.atoms = atoms
        else:
            raise ValueError('model needs an instance of an ase Atoms class')
        self.code = code
        self.mode = np.array(mode)
        self.mode = self.mode/np.linalg.norm(self.mode)

        self.input_read=False

    def set_code(self, code):
        """
        Sets the code keyword
        """

        if code in codes:
            self.code = code
        else:
            print 'This code is not supported, choose one of ', codes

    def set_mode(self, mode=None):
        """
        Sets the active_atoms list
        """

        self.mode = mode

    def read_input_data(self, **kwargs):
        """
        Reads all necessary input data including first order matrices
        wavefunctions, eigenvalues, fermi_level, etc.

        This is a wrapper function for specific parser routines, depending 
        on which code is used. For detailed input explanations please see 
        the documentation in the :py:coolvib.parser: package.

        """

        if self.code is None:
            print 'Please set code first using .set_code'

        method_to_call =getattr(parser,parser.codes[self.code]+'_mode')
        self = method_to_call(self, **kwargs)

        self.input_read = True

    def calculate_spectral_function(self,mode='default', **kwargs):
        """
        initiates calculation of the spectral functions and 
        takes the same arguments as coolvib.routines.calculate_spectral_function
     
        Parameters:

        mode : str
            'default' or 'momentum' for q=0 or q not 0 versions
        
        """

        masses = self.atoms.get_masses()[self.active_atoms]
        #CASE CODE
        # if parser.code_type[self.code] is 'local':
            # #CASE mode
            # if mode is 'default':
                # self.x_axis, self.spectral_function = spectral.calculate_spectral_function_mode(
                        # self.fermi_energy,
                        # self.eigenvalues,
                        # self.kpoints,
                        # self.psi,
                        # self.first_order_H,
                        # self.first_order_S,
                        # masses,
                        # self.occ,
                        # **kwargs)
            # elif mode is 'momentum':
                # self.x_axis, self.spectral_function = spectral.calculate_spectral_function_mode_q(
                        # self.fermi_energy,
                        # self.eigenvalues,
                        # self.kpoints,
                        # self.psi,
                        # self.first_order_H,
                        # self.first_order_S,
                        # masses,
                        # self.basis_pos,
                        # **kwargs)

        # else:
            #plane wave stuff
        raise NotImplementedError('calculate_spectral_function: code is not supported!')


    def print_spectral_function(self, filename='nacs-spectrum.out'):
        """
        prints spectral function to a file

        Parameters:

        filename: str
            name of the output file, default is nacs-spectrum.out
        """

        output.print_spectral_function(
                self.x_axis, self.spectral_function,1,
                filename)
    
